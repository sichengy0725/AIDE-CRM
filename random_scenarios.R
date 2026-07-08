# ============================================================
# Random generation of monotone dose-toxicity scenarios
# User specifies:
#   target      = target toxicity probability
#   target_diff = legacy desired average adjacent probability difference
#   target_diff_below / target_diff_above =
#     optional desired adjacent probability differences below/above the MTD
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


validate_target_diff <- function(target, target_diff, side = c("symmetric", "below", "above")) {
  
  side <- match.arg(side)
  
  if (!is.numeric(target_diff) || length(target_diff) != 1L || !is.finite(target_diff)) {
    stop("target_diff values must be single finite numbers.")
  }
  
  if (side == "symmetric") {
    if (target_diff <= 0 || target_diff >= 0.5) {
      stop("target_diff must be between 0 and 0.5.")
    }
    return(invisible(TRUE))
  }
  
  max_diff <- if (side == "below") target else 1 - target
  
  if (target_diff <= 0 || target_diff >= max_diff) {
    stop(
      sprintf(
        "target_diff_%s must be between 0 and %.4f for target = %.4f.",
        side,
        max_diff,
        target
      )
    )
  }
  
  invisible(TRUE)
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
  
  validate_target_diff(target, target_diff, side = "symmetric")
  
  f <- function(mu) {
    upper <- inv_logit(logit(target) + mu)
    lower <- inv_logit(logit(target) - mu)
    (upper - lower) / 2 - target_diff
  }
  
  uniroot(f, lower = 1e-8, upper = 50)$root
}


get_logit_spacing_one_sided <- function(target, target_diff, side = c("below", "above")) {
  
  side <- match.arg(side)
  
  if (target <= 0 || target >= 1) {
    stop("target must be between 0 and 1.")
  }
  
  validate_target_diff(target, target_diff, side = side)
  
  target_eta <- logit(target)
  
  f <- if (side == "below") {
    function(mu) {
      target - inv_logit(target_eta - mu) - target_diff
    }
  } else {
    function(mu) {
      inv_logit(target_eta + mu) - target - target_diff
    }
  }
  
  uniroot(f, lower = 1e-8, upper = 50)$root
}


resolve_target_diff_inputs <- function(
    target,
    target_diff,
    target_diff_below = NULL,
    target_diff_above = NULL
) {
  
  legacy_spacing <- is.null(target_diff_below) && is.null(target_diff_above)
  
  if (legacy_spacing) {
    if (is.null(target_diff)) {
      stop("target_diff must be specified when target_diff_below and target_diff_above are NULL.")
    }
    
    spacing_mean <- get_logit_spacing(
      target = target,
      target_diff = target_diff
    )
    
    return(list(
      legacy_spacing = TRUE,
      target_diff = target_diff,
      target_diff_below = target_diff,
      target_diff_above = target_diff,
      spacing_mean_below = spacing_mean,
      spacing_mean_above = spacing_mean
    ))
  }
  
  if (is.null(target_diff_below)) {
    if (is.null(target_diff)) {
      stop("target_diff_below must be specified when target_diff is NULL.")
    }
    target_diff_below <- target_diff
  }
  
  if (is.null(target_diff_above)) {
    if (is.null(target_diff)) {
      stop("target_diff_above must be specified when target_diff is NULL.")
    }
    target_diff_above <- target_diff
  }
  
  spacing_mean_below <- get_logit_spacing_one_sided(
    target = target,
    target_diff = target_diff_below,
    side = "below"
  )
  
  spacing_mean_above <- get_logit_spacing_one_sided(
    target = target,
    target_diff = target_diff_above,
    side = "above"
  )
  
  list(
    legacy_spacing = FALSE,
    target_diff = target_diff,
    target_diff_below = target_diff_below,
    target_diff_above = target_diff_above,
    spacing_mean_below = spacing_mean_below,
    spacing_mean_above = spacing_mean_above
  )
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
    target_diff_below = NULL,
    target_diff_above = NULL,
    mtd = NULL,
    spacing_sd = NULL,
    spacing_sd_below = NULL,
    spacing_sd_above = NULL,
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
  
  # Randomly choose the true MTD if not specified
  if (is.null(mtd)) {
    mtd <- sample(seq_len(ndose), size = 1)
  }
  
  if (mtd < 1 || mtd > ndose) {
    stop("mtd must be between 1 and ndose.")
  }
  
  # Mean adjacent spacing on the logit scale. If target_diff_below/above
  # are omitted, this preserves the legacy symmetric spacing behavior.
  target_diff_inputs <- resolve_target_diff_inputs(
    target = target,
    target_diff = target_diff,
    target_diff_below = target_diff_below,
    target_diff_above = target_diff_above
  )
  
  spacing_mean_below <- target_diff_inputs$spacing_mean_below
  spacing_mean_above <- target_diff_inputs$spacing_mean_above
  
  # Default spacing variability
  if (!is.null(spacing_sd) &&
      (!is.numeric(spacing_sd) || length(spacing_sd) != 1L || !is.finite(spacing_sd) || spacing_sd <= 0)) {
    stop("spacing_sd must be a single positive finite number.")
  }
  
  if (is.null(spacing_sd_below)) {
    spacing_sd_below <- if (is.null(spacing_sd)) spacing_mean_below / 3 else spacing_sd
  }
  
  if (is.null(spacing_sd_above)) {
    spacing_sd_above <- if (is.null(spacing_sd)) spacing_mean_above / 3 else spacing_sd
  }
  
  if (!is.numeric(spacing_sd_below) ||
      length(spacing_sd_below) != 1L ||
      !is.finite(spacing_sd_below) ||
      spacing_sd_below <= 0) {
    stop("spacing_sd_below must be a single positive finite number.")
  }
  
  if (!is.numeric(spacing_sd_above) ||
      length(spacing_sd_above) != 1L ||
      !is.finite(spacing_sd_above) ||
      spacing_sd_above <= 0) {
    stop("spacing_sd_above must be a single positive finite number.")
  }
  
  # Toxicity probability at the MTD is generated near target
  if (is.null(mtd_jitter)) {
    nearest_target_diff <- if (mtd == 1L) {
      target_diff_inputs$target_diff_above
    } else if (mtd == ndose) {
      target_diff_inputs$target_diff_below
    } else {
      min(target_diff_inputs$target_diff_below, target_diff_inputs$target_diff_above)
    }
    
    mtd_jitter <- nearest_target_diff / 2
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
          mean = spacing_mean_below,
          sd = spacing_sd_below
        )
        eta[j] <- eta[j + 1] - delta
      }
    }
    
    # Generate higher doses
    if (mtd < ndose) {
      for (j in (mtd + 1):ndose) {
        delta <- rposnorm(
          n = 1,
          mean = spacing_mean_above,
          sd = spacing_sd_above
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
        target_diff_below = target_diff_inputs$target_diff_below,
        target_diff_above = target_diff_inputs$target_diff_above,
        spacing_mean_logit = if (target_diff_inputs$legacy_spacing) spacing_mean_below else NA_real_,
        spacing_mean_logit_below = spacing_mean_below,
        spacing_mean_logit_above = spacing_mean_above,
        spacing_sd_logit = if (target_diff_inputs$legacy_spacing &&
                                identical(spacing_sd_below, spacing_sd_above)) spacing_sd_below else NA_real_,
        spacing_sd_logit_below = spacing_sd_below,
        spacing_sd_logit_above = spacing_sd_above,
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
    target_diff_below = NULL,
    target_diff_above = NULL,
    spacing_sd = NULL,
    spacing_sd_below = NULL,
    spacing_sd_above = NULL,
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
      target_diff_below = target_diff_below,
      target_diff_above = target_diff_above,
      spacing_sd = spacing_sd,
      spacing_sd_below = spacing_sd_below,
      spacing_sd_above = spacing_sd_above,
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
  
  attr(out, "target") <- target
  attr(out, "target_diff") <- target_diff
  attr(out, "target_diff_below") <- if (is.null(target_diff_below)) target_diff else target_diff_below
  attr(out, "target_diff_above") <- if (is.null(target_diff_above)) target_diff else target_diff_above
  
  return(out)
}


# -----------------------------
# Save generated scenarios
# -----------------------------

fmt_scenario_param <- function(x, digits = 4) {
  x <- round(as.numeric(x), digits)
  out <- format(x, scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  out[out == "-0"] <- "0"
  out <- gsub("-", "m", out)
  out <- gsub("\\.", "p", out)
  out
}


make_random_scenario_filename <- function(
    ndose,
    nscenario,
    target,
    target_diff_below,
    target_diff_above,
    out_dir = "scenario_sets"
) {
  file.path(
    out_dir,
    paste0(
      "random_scenarios",
      "_ndose", as.integer(ndose),
      "_nscenario", as.integer(nscenario),
      "_target", fmt_scenario_param(target),
      "_tdiffbelow", fmt_scenario_param(target_diff_below),
      "_tdiffabove", fmt_scenario_param(target_diff_above),
      ".csv"
    )
  )
}


save_generated_scenarios_csv <- function(scenarios, file) {
  if (!is.data.frame(scenarios)) {
    stop("scenarios must be a data.frame.")
  }

  out_dir <- dirname(file)
  if (nzchar(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  utils::write.csv(scenarios, file = file, row.names = FALSE)
  invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
}


# ============================================================
# User settings / run
# ============================================================
ndose <- 5L
nscenario <- 10000L
target <- 0.30
target_diff_below <- 0.15
target_diff_above <- 0.15

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L) {
  if (length(args) != 5L) {
    stop(
      "Usage: Rscript random_scenarios.R ",
      "<ndose> <nscenario> <target> <target_diff_below> <target_diff_above>"
    )
  }

  ndose <- as.integer(args[1])
  nscenario <- as.integer(args[2])
  target <- as.numeric(args[3])
  target_diff_below <- as.numeric(args[4])
  target_diff_above <- as.numeric(args[5])
}

scenarios <- generate_many_dose_toxicity_curves(
  nscenario = nscenario,
  ndose = ndose,
  target = target,
  target_diff_below = target_diff_below,
  target_diff_above = target_diff_above
)

scenario_file <- make_random_scenario_filename(
  ndose = ndose,
  nscenario = nscenario,
  target = target,
  target_diff_below = target_diff_below,
  target_diff_above = target_diff_above
)

saved_file <- save_generated_scenarios_csv(scenarios, scenario_file)

cat("Saved generated scenarios to:", saved_file, "\n")
print(scenarios)


