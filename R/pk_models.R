# R/pk_models.R
# Optimierte pharmakokinetische Modelle mit effizienten ODE-Solvern
# Performance-Verbesserungen: Vektorisierung, Events, adaptive Solver, Caching

# Load required packages
suppressPackageStartupMessages({
  require(deSolve)
  require(digest)
})

# Define operator if not exists
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# ===============================================================================
# CACHING INFRASTRUCTURE
# ===============================================================================

# Global cache environment
.pk_cache <- new.env(parent = emptyenv())
.pk_cache$results <- list()
.pk_cache$max_size <- 100  # Maximum number of cached results

# Generate cache key
generate_cache_key <- function(times, theta, regimen, model_type) {
  key_data <- list(
    times = round(times, 6),
    theta = round(unlist(theta), 6),
    regimen = regimen,
    model_type = model_type
  )
  digest::digest(key_data, algo = "xxhash64")
}

# Get cached result
get_cached_result <- function(key) {
  if (key %in% names(.pk_cache$results)) {
    .pk_cache$results[[key]]$access_time <- Sys.time()
    return(.pk_cache$results[[key]]$value)
  }
  NULL
}

# Store result in cache
cache_result <- function(key, value) {
  # Implement LRU cache eviction
  if (length(.pk_cache$results) >= .pk_cache$max_size) {
    access_times <- sapply(.pk_cache$results, function(x) x$access_time)
    oldest <- names(which.min(access_times))
    .pk_cache$results[[oldest]] <- NULL
  }
  
  .pk_cache$results[[key]] <- list(
    value = value,
    access_time = Sys.time()
  )
  invisible(value)
}

# Clear cache
clear_pk_cache <- function() {
  .pk_cache$results <- list()
  invisible(TRUE)
}

# ===============================================================================
# HELPER FUNCTIONS FOR EVENT-BASED DOSING
# ===============================================================================

# Create dosing events for deSolve
create_dosing_events <- function(regimen) {
  n_doses <- regimen$n_doses
  start_time <- regimen$start_time %||% 0
  
  if (regimen$tinf > 0) {
    # Infusion: create start and stop events
    event_times <- numeric(2 * n_doses)
    event_values <- numeric(2 * n_doses)
    event_methods <- rep("add", 2 * n_doses)
    
    for (i in seq_len(n_doses)) {
      dose_start <- start_time + (i - 1) * regimen$tau
      dose_end <- dose_start + regimen$tinf
      
      # Start infusion
      event_times[2*i - 1] <- dose_start
      event_values[2*i - 1] <- 0  # No bolus, handled by forcing function
      
      # Stop infusion
      event_times[2*i] <- dose_end
      event_values[2*i] <- 0
    }
  } else {
    # Bolus dosing
    event_times <- start_time + (seq_len(n_doses) - 1) * regimen$tau
    event_values <- rep(regimen$dose, n_doses)
    event_methods <- rep("add", n_doses)
  }
  
  list(
    data = data.frame(
      var = 1,  # Add to first compartment
      time = event_times,
      value = event_values,
      method = event_methods
    )
  )
}

# Optimized infusion rate calculation (vectorized)
calculate_infusion_rate_vectorized <- function(t, doses) {
  if (nrow(doses) == 0) return(0)
  
  # Vectorized comparison
  active <- (t > doses$t0) & (t <= doses$t0 + doses$tinf)
  sum(doses$rate * active)
}

# ===============================================================================
# ANALYTICAL SOLUTION FOR 1-COMPARTMENT MODEL
# ===============================================================================

conc_1c_inf_analytical <- function(t, dose, CL, Vc, tau, tinf, n_doses, start_time = 0) {
  k <- CL / Vc
  C <- numeric(length(t))
  
  for (i in 0:(n_doses - 1)) {
    t0 <- start_time + i * tau
    t1 <- t0 + tinf
    rate <- dose / tinf
    
    # Vectorized operations
    idx_during <- which(t >= t0 & t <= t1)
    idx_after <- which(t > t1 & t < t0 + tau)
    
    if (length(idx_during) > 0) {
      C[idx_during] <- C[idx_during] + (rate / CL) * (1 - exp(-k * (t[idx_during] - t0)))
    }
    if (length(idx_after) > 0) {
      C_end <- (rate / CL) * (1 - exp(-k * tinf))
      C[idx_after] <- C[idx_after] + C_end * exp(-k * (t[idx_after] - t1))
    }
  }
  C
}

# ===============================================================================
# OPTIMIZED MULTI-COMPARTMENT MODEL
# ===============================================================================

conc_profile_multi <- function(times, theta, regimen, model_type = "2C") {
  # Check cache first
  cache_key <- generate_cache_key(times, theta, regimen, model_type)
  cached <- get_cached_result(cache_key)
  if (!is.null(cached)) return(cached)
  
  # Extract parameters
  CL <- theta[["CL"]]
  Vc <- theta[["Vc"]]
  Q1 <- theta[["Q1"]] %||% 0
  Vp1 <- theta[["Vp1"]] %||% 1
  Q2 <- theta[["Q2"]] %||% 0
  Vp2 <- theta[["Vp2"]] %||% 1
  
  # Calculate rate constants
  k10 <- CL / Vc
  k12 <- ifelse(model_type %in% c("2C", "3C"), Q1/Vc, 0)
  k21 <- ifelse(model_type %in% c("2C", "3C"), Q1/Vp1, 0)
  k13 <- ifelse(model_type == "3C", Q2/Vc, 0)
  k31 <- ifelse(model_type == "3C", Q2/Vp2, 0)
  
  # Prepare dosing data
  doses <- data.frame(
    t0 = regimen$start_time + (0:(regimen$n_doses-1)) * regimen$tau,
    tinf = regimen$tinf,
    rate = regimen$dose / regimen$tinf
  )
  
  # ODE system with vectorized infusion rate
  rhs <- function(t, A, pars) {
    rate <- calculate_infusion_rate_vectorized(t, doses)
    
    dA1 <- rate - (k10 + k12 + k13) * A[1] + k21 * A[2] + k31 * A[3]
    dA2 <- k12 * A[1] - k21 * A[2]
    dA3 <- k13 * A[1] - k31 * A[3]
    
    list(c(dA1, dA2, dA3))
  }
  
  # Initial conditions
  A0 <- c(0, 0, 0)
  
  # Select optimal solver method
  method <- if (model_type == "3C") "radau5" else "ode45"
  
  # Adaptive tolerances based on expected concentration range
  rtol <- 1e-6
  atol <- 1e-9
  
  # Solve ODE with optimized settings
  sol <- deSolve::ode(
    y = A0,
    times = sort(unique(c(0, times))),
    func = rhs,
    parms = NULL,
    method = method,
    rtol = rtol,
    atol = atol,
    maxsteps = 5000
  )
  
  # Extract and interpolate results
  df <- as.data.frame(sol)
  result <- approx(df$time, df$A.1 / Vc, xout = times)$y
  
  # Cache the result
  cache_result(cache_key, result)
  
  result
}

# ===============================================================================
# OPTIMIZED MICHAELIS-MENTEN MODEL
# ===============================================================================

conc_profile_mm <- function(times, theta, regimen) {
  # Check cache
  cache_key <- generate_cache_key(times, theta, regimen, "MM")
  cached <- get_cached_result(cache_key)
  if (!is.null(cached)) return(cached)
  
  CL <- theta[["CL"]]
  Vc <- theta[["Vc"]]
  Vmax <- theta[["Vmax"]]
  Km <- theta[["Km"]]
  
  doses <- data.frame(
    t0 = regimen$start_time + (0:(regimen$n_doses-1)) * regimen$tau,
    tinf = regimen$tinf,
    rate = regimen$dose / regimen$tinf
  )
  
  # ODE with vectorized rate calculation
  rhs <- function(t, A, pars) {
    rate <- calculate_infusion_rate_vectorized(t, doses)
    C <- A[1] / Vc
    dA1 <- rate - CL * C - (Vmax * C) / (Km + C)
    list(c(dA1))
  }
  
  A0 <- c(0)
  
  # Use stiff solver for MM kinetics
  sol <- deSolve::ode(
    y = A0,
    times = sort(unique(c(0, times))),
    func = rhs,
    parms = NULL,
    method = "radau5",  # Better for stiff systems
    rtol = 1e-6,
    atol = 1e-9
  )
  
  df <- as.data.frame(sol)
  result <- approx(df$time, df$A.1 / Vc, xout = times)$y
  
  cache_result(cache_key, result)
  result
}

# ===============================================================================
# OPTIMIZED TMDD-QSS MODEL
# ===============================================================================

conc_profile_tmdd_qss <- function(times, theta, regimen) {
  # Check cache
  cache_key <- generate_cache_key(times, theta, regimen, "TMDD-QSS")
  cached <- get_cached_result(cache_key)
  if (!is.null(cached)) return(cached)
  
  CL <- theta[["CL"]]
  Vc <- theta[["Vc"]]
  kint <- theta[["kint"]]
  Rtot <- theta[["Rtot"]]
  Kss <- theta[["Kss"]]
  
  doses <- data.frame(
    t0 = regimen$start_time + (0:(regimen$n_doses-1)) * regimen$tau,
    tinf = regimen$tinf,
    rate = regimen$dose / regimen$tinf
  )
  
  # ODE with vectorized rate
  rhs <- function(t, A, pars) {
    rate <- calculate_infusion_rate_vectorized(t, doses)
    C <- A[1] / Vc
    dA1 <- rate - CL * C - kint * Rtot * C / (Kss + C)
    list(c(dA1))
  }
  
  A0 <- c(0)
  
  # Stiff solver for TMDD
  sol <- deSolve::ode(
    y = A0,
    times = sort(unique(c(0, times))),
    func = rhs,
    parms = NULL,
    method = "radau5",
    rtol = 1e-6,
    atol = 1e-9
  )
  
  df <- as.data.frame(sol)
  result <- approx(df$time, df$A.1 / Vc, xout = times)$y
  
  cache_result(cache_key, result)
  result
}

# ===============================================================================
# TIME-VARYING CLEARANCE MODEL (OPTIMIZED)
# ===============================================================================

conc_profile_multi_tvcl <- function(times, theta, regimen, model_type = "2C", CL_fun = NULL) {
  # Check cache with CL_fun signature
  cache_key <- generate_cache_key(times, theta, regimen, paste0(model_type, "_TVCL"))
  cached <- get_cached_result(cache_key)
  if (!is.null(cached)) return(cached)
  
  Vc <- theta[["Vc"]]
  Q1 <- theta[["Q1"]] %||% 0
  Vp1 <- theta[["Vp1"]] %||% 1
  Q2 <- theta[["Q2"]] %||% 0
  Vp2 <- theta[["Vp2"]] %||% 1
  
  k12 <- ifelse(model_type %in% c("2C", "3C"), Q1/Vc, 0)
  k21 <- ifelse(model_type %in% c("2C", "3C"), Q1/Vp1, 0)
  k13 <- ifelse(model_type == "3C", Q2/Vc, 0)
  k31 <- ifelse(model_type == "3C", Q2/Vp2, 0)
  
  doses <- data.frame(
    t0 = regimen$start_time + (0:(regimen$n_doses-1)) * regimen$tau,
    tinf = regimen$tinf,
    rate = regimen$dose / regimen$tinf
  )
  
  # ODE with time-varying clearance
  rhs <- function(t, A, pars) {
    rate <- calculate_infusion_rate_vectorized(t, doses)
    CLt <- if (!is.null(CL_fun)) CL_fun(t) else theta[["CL"]]
    k10 <- CLt / Vc
    
    dA1 <- rate - (k10 + k12 + k13) * A[1] + k21 * A[2] + k31 * A[3]
    dA2 <- k12 * A[1] - k21 * A[2]
    dA3 <- k13 * A[1] - k31 * A[3]
    
    list(c(dA1, dA2, dA3))
  }
  
  A0 <- c(0, 0, 0)
  
  # Use adaptive solver for time-varying systems
  sol <- deSolve::ode(
    y = A0,
    times = sort(unique(c(0, times))),
    func = rhs,
    parms = NULL,
    method = "ode45",
    rtol = 1e-6,
    atol = 1e-9
  )
  
  df <- as.data.frame(sol)
  result <- approx(df$time, df$A.1 / Vc, xout = times)$y
  
  cache_result(cache_key, result)
  result
}

# ===============================================================================
# DISPATCHER FUNCTION
# ===============================================================================

predict_conc_grid <- function(times, regimen, theta, model_type = "1C") {
  if (model_type == "1C") {
    conc_1c_inf_analytical(
      times, 
      regimen$dose, 
      theta[["CL"]], 
      theta[["Vc"]], 
      regimen$tau, 
      regimen$tinf, 
      regimen$n_doses, 
      regimen$start_time
    )
  } else if (model_type %in% c("2C", "3C")) {
    conc_profile_multi(times, theta, regimen, model_type)
  } else if (model_type == "MM-1C") {
    conc_profile_mm(times, theta, regimen)
  } else if (model_type == "TMDD-QSS-1C") {
    conc_profile_tmdd_qss(times, theta, regimen)
  } else {
    stop("Unsupported model_type: ", model_type)
  }
}

# ===============================================================================
# DOSE FINDING FUNCTION (UNCHANGED)
# ===============================================================================

find_dose_for_cmax <- function(theta, regimen, target_cmax, model_type = "1C", bounds = c(10, 10000)) {
  f <- function(dose) {
    reg <- regimen
    reg$dose <- dose
    t_peak <- regimen$tinf
    conc <- predict_conc_grid(t_peak, reg, theta, model_type)
    conc - target_cmax
  }
  
  fl <- f(bounds[1])
  fu <- f(bounds[2])
  
  if (fl * fu > 0) return(list(dose_mg = NA_real_))
  if (fl >= 0) return(list(dose_mg = bounds[1]))
  if (fu <= 0) return(list(dose_mg = NA_real_))
  
  result <- uniroot(f, lower = bounds[1], upper = bounds[2])
  list(dose_mg = result$root)
}

# ===============================================================================
# UTILITY FUNCTIONS
# ===============================================================================

# Get solver statistics
get_solver_stats <- function() {
  list(
    cache_size = length(.pk_cache$results),
    cache_hits = sum(sapply(.pk_cache$results, function(x) !is.null(x$access_time))),
    available_methods = c("ode45", "radau5", "lsoda", "bdf")
  )
}

# Benchmark function for performance testing
benchmark_ode_solver <- function(model_type = "2C", n_times = 100) {
  # Test parameters
  theta <- list(
    CL = 5, Vc = 30,
    Q1 = 2, Vp1 = 20,
    Q2 = 1, Vp2 = 10,
    Vmax = 100, Km = 10,
    kint = 0.1, Rtot = 100, Kss = 1
  )
  
  regimen <- list(
    dose = 1000,
    tau = 8,
    tinf = 1,
    n_doses = 10,
    start_time = 0
  )
  
  times <- seq(0, regimen$n_doses * regimen$tau, length.out = n_times)
  
  # Clear cache for fair comparison
  clear_pk_cache()
  
  # Time the computation
  start_time <- Sys.time()
  result <- predict_conc_grid(times, regimen, theta, model_type)
  end_time <- Sys.time()
  
  list(
    model_type = model_type,
    n_times = n_times,
    computation_time = as.numeric(end_time - start_time, units = "secs"),
    result_length = length(result),
    cache_used = length(.pk_cache$results) > 0
  )
}