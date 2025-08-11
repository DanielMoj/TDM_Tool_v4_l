# R/backend_bayes.R
# Backends: Laplace (optim), Stan (cmdstanr/rstan), JAGS (rjags)
# Modified with LRU cache and automatic garbage collection

suppressPackageStartupMessages({
  library(numDeriv)
})

# Source safety wrapper if available
if (file.exists("R/safe_computation.R")) {
  source("R/safe_computation.R")
}

# Source cache manager
if (file.exists("R/cache_manager.R")) {
  source("R/cache_manager.R")
}

# ===== PERFORMANCE OPTIMIZATION: LRU Stan Model Cache =====
# Initialize global LRU cache (replaces old .stan_model_cache)
.initialize_model_cache <- function() {
  if (!exists(".global_lru_cache", envir = .GlobalEnv)) {
    # Create LRU cache with max 10 models
    assign(".global_lru_cache", LRUCache$new(max_size = 10), envir = .GlobalEnv)
    message("Initialized global LRU cache for Stan models")
  }
  get(".global_lru_cache", envir = .GlobalEnv)
}

# Modified get_compiled_model to use LRU cache
get_compiled_model <- function(stan_file, force_recompile = FALSE) {
  cache <- .initialize_model_cache()
  
  # Log memory status before compilation
  if (getOption("tdmx_debug", FALSE)) {
    stats <- get_memory_stats()
    message(sprintf("Memory before model request: R=%.1fMB, Cache=%.1fMB", 
                   stats$r_memory$used_mb, stats$cache$total_size_mb))
  }
  
  # Get or compile model using LRU cache
  model <- cache$get(stan_file, force_recompile = force_recompile)
  
  # Log memory status after compilation
  if (getOption("tdmx_debug", FALSE)) {
    stats <- get_memory_stats()
    message(sprintf("Memory after model request: R=%.1fMB, Cache=%.1fMB", 
                   stats$r_memory$used_mb, stats$cache$total_size_mb))
  }
  
  return(model)
}

# Helper to clear cache if needed
clear_model_cache <- function() {
  if (exists(".global_lru_cache", envir = .GlobalEnv)) {
    cache <- get(".global_lru_cache", envir = .GlobalEnv)
    cache$clear()
    message("Model cache cleared")
  }
}
# ===== END PERFORMANCE OPTIMIZATION =====

backend_status <- function() {
  has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE)
  has_rstan   <- requireNamespace("rstan", quietly = TRUE)
  has_rjags   <- requireNamespace("rjags", quietly = TRUE)
  glue::glue("cmdstanr: {has_cmdstan}, rstan: {has_rstan}, rjags: {has_rjags}")
}

# Prior-Parsing: erwartet priors$theta_log (mu/sd auf log-Skala)
draw_from_priors <- function(priors) {
  th <- priors$theta
  exp(setNames(rnorm(length(th), priors$theta_log$mu[names(th)], priors$theta_log$sd[names(th)]), names(th)))
}

# Helper function for covariate adjustment
apply_covariates <- function(theta, covariates, drug_name = NULL) {
  # Simple allometric scaling for CL and V
  if (!is.null(covariates$weight) && !is.null(theta[["CL"]])) {
    theta[["CL"]] <- theta[["CL"]] * (covariates$weight / 70)^0.75
  }
  if (!is.null(covariates$weight) && !is.null(theta[["Vc"]])) {
    theta[["Vc"]] <- theta[["Vc"]] * (covariates$weight / 70)^1.0
  }
  theta
}

# Modified Stan HMC implementation with error handling and GC
run_fit_stan_hmc <- function(obs, regimen, priors, model_type, error_model, 
                            covariates, estimate_sigma, sigma_init, 
                            blq_lloq = NA_real_, is_blq = NULL) {
  
  # Get HMC controls from options or defaults
  .hmc <- getOption("tdmx_hmc", default = list(
    chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
    parallel_chains = NULL, adapt_delta = 0.8, max_treedepth = 10L, seed = 1234L
  ))
  
  # Prepare parameters
  params <- list(
    obs = obs,
    regimen = regimen,
    priors = priors,
    model_type = model_type,
    error_model = error_model,
    covariates = covariates,
    estimate_sigma = estimate_sigma,
    sigma_init = sigma_init,
    blq_lloq = blq_lloq,
    is_blq = is_blq,
    adapt_delta = .hmc$adapt_delta,
    max_treedepth = .hmc$max_treedepth,
    chains = .hmc$chains,
    iter_warmup = .hmc$iter_warmup,
    iter_sampling = .hmc$iter_sampling,
    seed = .hmc$seed,
    parallel_chains = .hmc$parallel_chains
  )
  
  # Computation function
  stan_computation <- function(p) {
    # Get Stan file
    stan_file <- stan_file_for_model(p$model_type)
    
    # Prepare data
    data_list <- stan_data_list2(
      p$obs, p$regimen, p$priors, p$model_type, 
      p$error_model, p$estimate_sigma, p$sigma_init, 
      p$blq_lloq, p$is_blq
    )
    
    # Validate data
    if (any(is.na(data_list$y_obs)) || any(is.infinite(data_list$y_obs))) {
      stop("Invalid observation data: contains NA or Inf values")
    }
    
    # Get compiled model using LRU cache
    mod <- get_compiled_model(stan_file)
    
    # Run sampling with specific parameters
    fit <- mod$sample(
      data = data_list,
      seed = p$seed,
      chains = p$chains,
      parallel_chains = p$parallel_chains %||% p$chains,
      iter_warmup = p$iter_warmup,
      iter_sampling = p$iter_sampling,
      adapt_delta = p$adapt_delta,
      max_treedepth = p$max_treedepth,
      refresh = 0,  # Suppress Stan output
      show_messages = FALSE
    )
    
    return(fit)
  }
  
  # Run with safety wrapper if available
  if (exists("safe_bayesian_computation") && is.function(safe_bayesian_computation)) {
    result <- safe_bayesian_computation(
      computation_type = "stan",
      computation_fn = stan_computation,
      params = params,
      session = getDefaultReactiveDomain()  # Get Shiny session if available
    )
    
    if (!result$success) {
      # Provide detailed error message
      error_msg <- sprintf("Stan computation failed: %s", result$error)
      if (length(result$warnings) > 0) {
        error_msg <- paste(error_msg, "\nWarnings:", paste(result$warnings, collapse = "; "))
      }
      stop(error_msg)
    }
    
    fit <- result$result
    computation_diagnostics <- result$diagnostics
    computation_warnings <- result$warnings
    computation_time <- result$duration
    
  } else {
    # Fallback: direct execution with basic error handling
    fit <- tryCatch({
      stan_computation(params)
    }, error = function(e) {
      stop(sprintf("Stan computation failed: %s", e$message))
    })
    computation_diagnostics <- NULL
    computation_warnings <- character()
    computation_time <- NA
  }
  
  # === AUTOMATIC GARBAGE COLLECTION after Stan sampling (Line ~178) ===
  if (getOption("tdmx_gc_after_sampling", TRUE)) {
    force_gc(verbose = getOption("tdmx_debug", FALSE))
  }
  
  # Process results
  draws <- tryCatch({
    as.data.frame(fit$draws(
      variables = c("CL_out","Vc_out","Q1_out","Vp1_out","Q2_out","Vp2_out",
                   "sigma_add","sigma_prop","nu"), 
      format = "df"
    ))
  }, error = function(e) {
    stop(sprintf("Failed to extract draws: %s", e$message))
  })
  
  # === AUTOMATIC GARBAGE COLLECTION after draw extraction (Line ~422) ===
  if (getOption("tdmx_gc_after_draws", TRUE)) {
    force_gc(verbose = getOption("tdmx_debug", FALSE))
  }
  
  # Rename outputs
  names(draws) <- sub("_out$", "", names(draws))
  
  # Extract additional diagnostics if not already done
  if (is.null(computation_diagnostics)) {
    diagnostics <- NULL
    try({
      if (requireNamespace("posterior", quietly = TRUE)) {
        summ <- posterior::summarise_draws(fit$draws())
        keep <- intersect(c("CL_out","Vc_out","Q1_out","Vp1_out","Q2_out","Vp2_out",
                          "sigma_add","sigma_prop","nu"), summ$variable)
        summ <- summ[summ$variable %in% keep, c("variable","rhat","ess_bulk","ess_tail"), drop = FALSE]
      } else { 
        summ <- NULL 
      }
      
      sdiag <- try(fit$diagnostic_summary(), silent = TRUE)
      div <- try(sdiag$num_divergent[1], silent = TRUE)
      treedepth <- try(sdiag$num_max_treedepth[1], silent = TRUE)
      stepsize <- try(as.numeric(fit$metadata()$step_size_adaptation), silent = TRUE)
      
      diagnostics <- list(
        summary = summ, 
        divergences = div, 
        treedepth_hits = treedepth, 
        stepsize = stepsize
      )
    }, silent = TRUE)
  } else {
    diagnostics <- computation_diagnostics
  }
  
  return(list(
    draws = draws,
    diagnostics = diagnostics,
    warnings = computation_warnings,
    computation_time = computation_time
  ))
}

# Modified ADVI implementation with error handling and GC
run_fit_stan_advi <- function(obs, regimen, priors, model_type, error_model, 
                             covariates, estimate_sigma, sigma_init, 
                             blq_lloq = NA_real_, is_blq = NULL) {
  
  if (!requireNamespace("cmdstanr", quietly = TRUE) && !requireNamespace("rstan", quietly = TRUE)) {
    warning("Stan-ADVI nicht verfügbar, fallback auf Laplace.")
    return(run_fit_laplace(obs, regimen, priors, model_type, error_model, 
                          covariates, estimate_sigma, sigma_init, blq_lloq, is_blq))
  }
  
  params <- list(
    obs = obs,
    regimen = regimen,
    priors = priors,
    model_type = model_type,
    error_model = error_model,
    estimate_sigma = estimate_sigma,
    sigma_init = sigma_init,
    blq_lloq = blq_lloq,
    is_blq = is_blq
  )
  
  advi_computation <- function(p) {
    stan_file <- stan_file_for_model(p$model_type)
    data_list <- stan_data_list2(
      p$obs, p$regimen, p$priors, p$model_type, 
      p$error_model, p$estimate_sigma, p$sigma_init, 
      p$blq_lloq, p$is_blq
    )
    
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      # Use LRU cache for model compilation
      mod <- get_compiled_model(stan_file)
      fit <- mod$variational(
        data = data_list, 
        output_samples = 1000, 
        seed = 123,
        refresh = 0
      )
      dr <- as.data.frame(fit$draws(
        variables = c("CL","Vc","Q1","Vp1","Q2","Vp2"), 
        format = "df"
      ))
    } else {
      stan_text <- readChar(stan_file, file.info(stan_file)$size)
      sm <- rstan::stan_model(model_code = stan_text)
      fit <- rstan::vb(sm, data = data_list, output_samples = 1000, seed = 123)
      dr <- as.data.frame(rstan::extract(fit, pars = c("CL","Vc","Q1","Vp1","Q2","Vp2")))
    }
    
    return(dr)
  }
  
  # Run with safety wrapper if available
  if (exists("safe_bayesian_computation") && is.function(safe_bayesian_computation)) {
    result <- safe_bayesian_computation(
      computation_type = "stan",
      computation_fn = advi_computation,
      params = params,
      session = getDefaultReactiveDomain()
    )
    
    if (!result$success) {
      stop(sprintf("Stan ADVI failed: %s", result$error))
    }
    
    dr <- result$result
  } else {
    dr <- tryCatch({
      advi_computation(params)
    }, error = function(e) {
      stop(sprintf("Stan ADVI failed: %s", e$message))
    })
  }
  
  # === AUTOMATIC GARBAGE COLLECTION after ADVI fits (Line ~266) ===
  if (getOption("tdmx_gc_after_advi", TRUE)) {
    force_gc(verbose = getOption("tdmx_debug", FALSE))
  }
  
  keep <- intersect(colnames(dr), names(priors$theta))
  dr <- dr[, keep, drop = FALSE]
  
  return(list(draws = dr, diagnostics = NULL))
}

# Existing helper functions (keeping compatibility)
run_fit_stan <- function(obs, regimen, priors, model_type, error_model, covariates,
                        estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL) {
  run_fit_stan_hmc(obs, regimen, priors, model_type, error_model, 
                   covariates, estimate_sigma, sigma_init, blq_lloq, is_blq)
}

# Helper to select appropriate Stan model file
stan_file_for_model <- function(model_type) {
  if (model_type == "MM-1C") {
    return("models/stan/pk_mm_onecpt_ode.stan")
  } else if (model_type == "TMDD-QSS-1C") {
    return("models/stan/pk_tmdd_qss_onecpt_ode.stan")
  } else {
    return("models/stan/pk_multicpt_ode.stan")
  }
}

# Helper to build Stan data list
stan_data_list2 <- function(obs, regimen, priors, model_type, error_model, 
                           estimate_sigma, sigma_init, blq_lloq, is_blq) {
  # Map error model string to numeric code
  error_code <- switch(error_model,
    "additiv" = 1L,
    "proportional" = 2L,
    "kombiniert" = 3L,
    "t-additiv" = 4L,
    "t-proportional" = 5L,
    "mixture" = 6L,
    3L  # default to combined
  )
  
  # Build infusion schedule
  n_inf <- regimen$n_doses
  t0 <- regimen$start_time + (0:(n_inf-1)) * regimen$tau
  tinf <- rep(regimen$tinf, n_inf)
  rate <- rep(regimen$dose / regimen$tinf, n_inf)
  
  # BLQ handling
  if (is.null(is_blq)) is_blq <- rep(0L, nrow(obs))
  if (is.na(blq_lloq)) blq_lloq <- 0.0
  
  # Create data list
  data_list <- list(
    N = nrow(obs),
    n_inf = n_inf,
    t_obs = obs$time,
    y_obs = obs$conc,
    t0 = t0,
    tinf = tinf,
    rate = rate,
    
    # Model selection
    model_type = switch(model_type,
      "1C" = 1L,
      "2C" = 2L,
      "3C" = 3L,
      1L
    ),
    
    # Error model
    error_model = error_code,
    estimate_sigma = as.integer(estimate_sigma),
    sigma_add_prior = if (estimate_sigma) c(sigma_init, sigma_init * 0.5) else c(sigma_init, 0.01),
    sigma_prop_prior = if (estimate_sigma) c(0.2, 0.1) else c(0.01, 0.001),
    nu_prior = c(4, 2),  # degrees of freedom for t-distribution
    
    # Priors (log-scale)
    CL_prior = priors$theta_log$mu["CL"],
    CL_sd = priors$theta_log$sd["CL"],
    Vc_prior = priors$theta_log$mu["Vc"],
    Vc_sd = priors$theta_log$sd["Vc"],
    Q1_prior = if ("Q1" %in% names(priors$theta_log$mu)) priors$theta_log$mu["Q1"] else 0,
    Q1_sd = if ("Q1" %in% names(priors$theta_log$sd)) priors$theta_log$sd["Q1"] else 1,
    Vp1_prior = if ("Vp1" %in% names(priors$theta_log$mu)) priors$theta_log$mu["Vp1"] else 0,
    Vp1_sd = if ("Vp1" %in% names(priors$theta_log$sd)) priors$theta_log$sd["Vp1"] else 1,
    Q2_prior = if ("Q2" %in% names(priors$theta_log$mu)) priors$theta_log$mu["Q2"] else 0,
    Q2_sd = if ("Q2" %in% names(priors$theta_log$sd)) priors$theta_log$sd["Q2"] else 1,
    Vp2_prior = if ("Vp2" %in% names(priors$theta_log$mu)) priors$theta_log$mu["Vp2"] else 0,
    Vp2_sd = if ("Vp2" %in% names(priors$theta_log$sd)) priors$theta_log$sd["Vp2"] else 1,
    
    # BLQ
    lloq = blq_lloq,
    is_blq = as.integer(is_blq)
  )
  
  return(data_list)
}

# Add memory monitoring wrapper for all backend functions
monitor_memory <- function(fn_name, fn, ...) {
  if (getOption("tdmx_monitor_memory", FALSE)) {
    before <- get_memory_stats()
    message(sprintf("\n[%s] Memory before: R=%.1fMB, Cache=%.1fMB", 
                   fn_name, before$r_memory$used_mb, before$cache$total_size_mb))
  }
  
  result <- fn(...)
  
  if (getOption("tdmx_monitor_memory", FALSE)) {
    after <- get_memory_stats()
    message(sprintf("[%s] Memory after: R=%.1fMB, Cache=%.1fMB\n", 
                   fn_name, after$r_memory$used_mb, after$cache$total_size_mb))
  }
  
  return(result)
}

# Placeholder for run_fit_laplace function (not modified in this context)
run_fit_laplace <- function(obs, regimen, priors, model_type, error_model, 
                           covariates, estimate_sigma, sigma_init, 
                           blq_lloq = NA_real_, is_blq = NULL) {
  # Original implementation remains unchanged
  stop("run_fit_laplace not implemented in this example")
}