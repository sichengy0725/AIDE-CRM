# ============================================================
# Random generation of monotone dose-toxicity scenarios
# User specifies:
#   target      = target toxicity probability
#   target_diff = desired average adjacent probability difference
# ============================================================


# -----------------------------
# Helper functions
# -----------------------------

logit <- function(p) {
  qlogis(p)
}

inv_logit <- function(x) {
  plogis(x)
}


# Convert desired probability difference around target
# into a mean spacing on the logit scale.
#
# We solve for mu such that:
#
#   [logit^{-1}(logit(theta) + mu) -
#    logit^{-1}(logit(theta) - mu)] / 2 = target_diff
#
get_logit_spacing <- function(target, target_diff) {
  
  if (target <= 0 || target >= 1) {
    stop("target must be between 0 and 1.")
  }
  
  if (target_diff <= 0 || target_diff >= 0.5) {
    stop("target_diff must be between 0 and 0.5.")
  }
  
  f <- function(mu) {
    upper <- inv_logit(logit(target) + mu)
    lower <- inv_logit(logit(target) - mu)
    (upper - lower) / 2 - target_diff
  }
  
  uniroot(f, lower = 1e-8, upper = 50)$root
}


# Draw from a positive normal distribution by rejection sampling
rposnorm <- function(n, mean, sd) {
  
  out <- numeric(n)
  i <- 1
  
  while (i <= n) {
    x <- rnorm(1, mean = mean, sd = sd)
    
    if (x > 0) {
      out[i] <- x
      i <- i + 1
    }
  }
  
  out
}


# -----------------------------
# Generate one scenario
# -----------------------------

generate_one_dose_toxicity_curve <- function(
    ndose,
    target = 0.30,
    target_diff = 0.15,
    mtd = NULL,
    spacing_sd = NULL,
    mtd_jitter = NULL,
    digits = 2,
    max_try = 10000
) {
  
  if (ndose < 2) {
    stop("ndose must be at least 2.")
  }
  
  if (target <= 0 || target >= 1) {
    stop("target must be between 0 and 1.")
  }
  
  if (target_diff <= 0 || target_diff >= 0.5) {
    stop("target_diff must be between 0 and 0.5.")
  }
  
  # Randomly choose the true MTD if not specified
  if (is.null(mtd)) {
    mtd <- sample(seq_len(ndose), size = 1)
  }
  
  if (mtd < 1 || mtd > ndose) {
    stop("mtd must be between 1 and ndose.")
  }
  
  # Mean adjacent spacing on the logit scale
  spacing_mean <- get_logit_spacing(
    target = target,
    target_diff = target_diff
  )
  
  # Default spacing variability
  if (is.null(spacing_sd)) {
    spacing_sd <- spacing_mean / 3
  }
  
  # Toxicity probability at the MTD is generated near target
  if (is.null(mtd_jitter)) {
    mtd_jitter <- target_diff / 2
  }
  
  lower_mtd <- max(1e-6, target - mtd_jitter)
  upper_mtd <- min(1 - 1e-6, target + mtd_jitter)
  
  for (attempt in seq_len(max_try)) {
    
    # Generate probability at the MTD
    p_mtd <- runif(1, lower_mtd, upper_mtd)
    
    eta <- rep(NA_real_, ndose)
    eta[mtd] <- logit(p_mtd)
    
    # Generate lower doses
    if (mtd > 1) {
      for (j in (mtd - 1):1) {
        delta <- rposnorm(
          n = 1,
          mean = spacing_mean,
          sd = spacing_sd
        )
        eta[j] <- eta[j + 1] - delta
      }
    }
    
    # Generate higher doses
    if (mtd < ndose) {
      for (j in (mtd + 1):ndose) {
        delta <- rposnorm(
          n = 1,
          mean = spacing_mean,
          sd = spacing_sd
        )
        eta[j] <- eta[j - 1] + delta
      }
    }
    
    ptox <- inv_logit(eta)
    ptox_round <- round(ptox, digits)
    
    # Check that the rounded curve is still monotone
    rounded_monotone <- all(diff(ptox_round) > 0)
    
    # Check that the rounded curve still has the intended MTD
    rounded_mtd <- which.min(abs(ptox_round - target))
    rounded_mtd_ok <- rounded_mtd == mtd
    
    if (rounded_monotone && rounded_mtd_ok) {
      return(list(
        ptox = ptox_round,
        mtd = mtd,
        target = target,
        target_diff = target_diff,
        spacing_mean_logit = spacing_mean,
        spacing_sd_logit = spacing_sd,
        attempt = attempt
      ))
    }
  }
  
  stop("Failed to generate a valid rounded curve. Try reducing spacing_sd or mtd_jitter.")
}


# -----------------------------
# Generate many scenarios
# -----------------------------

generate_many_dose_toxicity_curves <- function(
    nscenario,
    ndose,
    target = 0.30,
    target_diff = 0.15,
    spacing_sd = NULL,
    mtd_jitter = NULL,
    digits = 2,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  curves <- matrix(NA_real_, nrow = nscenario, ncol = ndose)
  mtd_vec <- integer(nscenario)
  attempt_vec <- integer(nscenario)
  
  for (s in seq_len(nscenario)) {
    
    res <- generate_one_dose_toxicity_curve(
      ndose = ndose,
      target = target,
      target_diff = target_diff,
      spacing_sd = spacing_sd,
      mtd_jitter = mtd_jitter,
      digits = digits
    )
    
    curves[s, ] <- res$ptox
    mtd_vec[s] <- res$mtd
    attempt_vec[s] <- res$attempt
  }
  
  colnames(curves) <- paste0("Dose", seq_len(ndose))
  
  out <- data.frame(
    Scenario = seq_len(nscenario),
    MTD = mtd_vec,
    curves,
    Attempt = attempt_vec,
    row.names = NULL
  )
  
  return(out)
}


# ============================================================
# Example
# ============================================================
scenarios <- generate_many_dose_toxicity_curves(
  nscenario = 5,
  ndose = 8,
  target = 0.30,
  target_diff = 0.2
)

print(scenarios)


