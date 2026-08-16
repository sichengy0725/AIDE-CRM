aide_beta_prior <- function(x, name) {
  x <- as.numeric(x); if (length(x) != 2L || any(!is.finite(x)) || any(x <= 0)) stop(name, " must be a positive beta prior.")
  x
}

aide_default_toxicity_skeleton <- function(ndose, target) {
  # A weak, target-calibrated CRM-style skeleton anchors the initial dose at a
  # genuinely low toxicity level instead of the 0.5 mean of Beta(1,1).
  exp(seq(log(.05), log(min(.80, max(.40, 2 * target))), length.out = ndose))
}

aide_phase12_require_jags <- function() {
  if (!requireNamespace("rjags", quietly = TRUE) || !requireNamespace("coda", quietly = TRUE))
    stop("JAGS fitting requires the R packages 'rjags' and 'coda'. No fallback model is used.")
}

aide_phase12_jags_file <- function(name) {
  path <- file.path(getOption("tite_aide_root"), "inst", "jags", name)
  if (!file.exists(path)) stop("Required JAGS model file is unavailable: ", path)
  path
}

aide_expand_beta_prior <- function(prior, ndose, name) {
  prior <- as.numeric(prior)
  if (length(prior) == 2L) prior <- rep(prior, ndose)
  if (length(prior) != 2L * ndose || any(!is.finite(prior)) || any(prior <= 0))
    stop(name, " must be a positive length-2 prior or an ndose-by-2 prior.")
  matrix(prior, ncol = 2L, byrow = TRUE)
}

aide_beta_binomial_futility <- function(interim, efficacy, ndose) {
  ## Futility is deliberately separate from the delayed-outcome efficacy
  ## model used for utility.  It has the specified Beta(1, 1) posterior and
  ## therefore uses only outcomes that have already been ascertained.
  if (length(ndose) != 1L || !is.finite(ndose) || ndose < 1L ||
      ndose != as.integer(ndose)) {
    stop("ndose must be a positive integer.")
  }
  ndose <- as.integer(ndose)
  if (!is.data.frame(interim) || !all(c("dose", "y", "ascertained") %in% names(interim))) {
    stop("interim efficacy data must contain dose, y, and ascertained columns.")
  }
  if (length(efficacy$threshold) != 1L || !is.finite(efficacy$threshold) ||
      efficacy$threshold <= 0 || efficacy$threshold >= 1) {
    stop("efficacy$threshold must be a scalar in (0, 1).")
  }
  if (length(efficacy$futility_cutoff) != 1L || !is.finite(efficacy$futility_cutoff) ||
      efficacy$futility_cutoff <= 0 || efficacy$futility_cutoff >= 1) {
    stop("efficacy$futility_cutoff must be a scalar in (0, 1).")
  }
  if (length(efficacy$min_eff_n_for_futility) != 1L ||
      !is.finite(efficacy$min_eff_n_for_futility) ||
      efficacy$min_eff_n_for_futility < 0 ||
      efficacy$min_eff_n_for_futility != as.integer(efficacy$min_eff_n_for_futility)) {
    stop("efficacy$min_eff_n_for_futility must be a non-negative integer.")
  }
  if (!nrow(interim)) {
    n <- y <- integer(ndose)
  } else {
    dose <- as.integer(interim$dose)
    if (any(!is.finite(interim$dose)) || any(is.na(dose)) ||
        any(dose < 1L | dose > ndose)) {
      stop("interim efficacy data contain an invalid dose.")
    }
    ascertained <- !is.na(interim$ascertained) & as.logical(interim$ascertained)
    observed_y <- interim$y[ascertained]
    if (any(is.na(observed_y)) || any(!observed_y %in% c(0L, 1L))) {
      stop("Ascertainable efficacy outcomes must be coded 0 or 1.")
    }
    n <- tabulate(dose[ascertained], nbins = ndose)
    y <- tabulate(dose[ascertained & interim$y == 1L], nbins = ndose)
  }

  ## p_E,j | D ~ Beta(1 + y_j, 1 + n_j - y_j).  Pending efficacy records
  ## contribute to neither count and hence cannot affect this probability.
  prob_below_threshold <- stats::pbeta(
    efficacy$threshold,
    shape1 = 1 + y,
    shape2 = 1 + n - y
  )
  futile <- n >= as.integer(efficacy$min_eff_n_for_futility) &
    prob_below_threshold > efficacy$futility_cutoff
  list(
    n_fully_ascertained_by_dose = as.integer(n),
    y_fully_ascertained_by_dose = as.integer(y),
    prob_below_threshold = prob_below_threshold,
    futile = futile,
    prior = c(1, 1)
  )
}

aide_fit_toxicity <- function(interim, config, ndose) {
  # Skeleton TITE-CRM fit.  All MCMC draws are transient and immediately
  # reduced to posterior means/probabilities returned to the event engine.
  aide_phase12_require_jags()
  skeleton <- config$toxicity$skeleton %||% aide_default_toxicity_skeleton(ndose, config$toxicity$target)
  if (length(skeleton) != ndose) stop("toxicity$skeleton must have one entry per dose.")
  tx <- config$toxicity
  carry_prior <- aide_beta_prior(tx$carryover_prior, "toxicity$carryover_prior")
  dat <- interim
  if (!nrow(dat)) dat <- data.frame(dose = 1L, previous_dose = 1L, ipde = 0L, y = 0L, weight = 0)
  dat$y <- ifelse(is.na(dat$y), 0L, as.integer(dat$y)); dat$weight[dat$y == 1L] <- 1
  jags_data <- list(N = nrow(dat), J = ndose, y = as.integer(dat$y), dose = as.integer(dat$dose),
    ipde = as.integer(dat$ipde), w = as.numeric(dat$weight), q = as.numeric(skeleton), beta_mean = tx$beta_prior_mean,
    a_r = carry_prior[1], b_r = carry_prior[2])
  additive <- identical(tx$model, "additive_alpha")
  if (additive) {
    jags_data$previous_dose <- ifelse(is.na(dat$previous_dose), 1L, as.integer(dat$previous_dose))
    jags_data$tau_beta <- 1 / tx$beta_prior_sd^2
    model_file <- aide_phase12_jags_file("previous_dose_additive_CRM_TITE.bug")
    monitors <- c("p", "theta_ipde", "alpha")
  } else {
    jags_data$tau_alpha <- 1 / tx$beta_prior_sd^2
    model_file <- aide_phase12_jags_file("random_CRM_TITE.bug")
    monitors <- c("p", "theta_ipde", "r")
  }
  jm <- rjags::jags.model(model_file, data = jags_data, n.chains = tx$n_chains, n.adapt = tx$n_adapt, quiet = TRUE)
  if (tx$n_burnin > 0L) stats::update(jm, n.iter = tx$n_burnin, progress.bar = "none")
  draws <- as.matrix(rjags::coda.samples(jm, variable.names = monitors, n.iter = tx$n_iter, thin = tx$thin, progress.bar = "none"))
  p_cols <- paste0("p[", seq_len(ndose), "]"); ipde_cols <- paste0("theta_ipde[", seq_len(ndose), "]")
  if (!all(p_cols %in% colnames(draws)) || !all(ipde_cols %in% colnames(draws))) stop("The toxicity JAGS model did not return all dose summaries.")
  carry_draws <- if (additive) draws[, "alpha"] else draws[, "r"]
  prob_overtox <- colMeans(draws[, p_cols, drop = FALSE] > tx$target)
  prob_ipde_overtox <- colMeans(draws[, ipde_cols, drop = FALSE] > tx$target)
  list(model = tx$model, p_regular_mean = colMeans(draws[, p_cols, drop = FALSE]),
       carryover_mean = rep(mean(carry_draws), ndose), p_ipde_mean = colMeans(draws[, ipde_cols, drop = FALSE]),
       prob_overtox_by_dose = prob_overtox, prob_ipde_overtox_by_dose = prob_ipde_overtox,
       eliminated = prob_overtox > tx$cutoff,
       skeleton = skeleton, posterior_carryover_mean = mean(carry_draws))
}

aide_fit_efficacy <- function(interim, config, ndose) {
  aide_phase12_require_jags(); ef <- config$efficacy
  futility <- aide_beta_binomial_futility(interim, ef, ndose)
  regular_prior <- aide_expand_beta_prior(ef$prior, ndose, "efficacy$prior")
  carry_prior <- aide_expand_beta_prior(ef$carryover_prior, ndose, "efficacy$carryover_prior")
  dat <- interim
  if (!nrow(dat)) dat <- data.frame(dose = 1L, previous_dose = 1L, ipde = 0L, y = NA_integer_, ascertained = FALSE, weight = 0)
  dat$y <- ifelse(is.na(dat$y), 0L, as.integer(dat$y)); dat$ascertained <- as.integer(dat$ascertained)
  model <- ef$model; additive <- model %in% c("previous_dose_additive", "dose_specific_previous_dose_additive")
  model_file <- switch(model,
    dose_specific_carryover = "beta_binomial_tite_dose_specific_carryover.jags",
    shared_carryover = "beta_binomial_tite_shared_carryover.jags",
    previous_dose_additive = "previous_dose_additive_tite_beta_binomial_efficacy.jags",
    dose_specific_previous_dose_additive = "previous_dose_additive_dose_specific_tite_beta_binomial_efficacy.jags")
  jags_data <- list(N = nrow(dat), ndose = ndose, eff = as.integer(dat$y), dose = as.integer(dat$dose),
    ipde = as.integer(dat$ipde),
    eff_ascertained = as.integer(dat$ascertained), eff_weight = as.numeric(dat$weight), zeros = rep.int(0L, nrow(dat)),
    eps = 1e-12, zero_trick_constant = 1, a_regular = regular_prior[, 1L], b_regular = regular_prior[, 2L])
  if (additive) { jags_data$previous_dose <- ifelse(is.na(dat$previous_dose), 1L, as.integer(dat$previous_dose)); jags_data$a_alpha <- if (model == "previous_dose_additive") carry_prior[1L, 1L] else carry_prior[, 1L]; jags_data$b_alpha <- if (model == "previous_dose_additive") carry_prior[1L, 2L] else carry_prior[, 2L]; monitors <- c("p_regular", "alpha")
  } else if (model == "shared_carryover") { jags_data$a_carry <- carry_prior[1L, 1L]; jags_data$b_carry <- carry_prior[1L, 2L]; monitors <- c("p_regular", "r_e")
  } else { jags_data$a_carry <- carry_prior[, 1L]; jags_data$b_carry <- carry_prior[, 2L]; monitors <- c("p_regular", "r_e") }
  jm <- rjags::jags.model(aide_phase12_jags_file(model_file), data = jags_data, n.chains = ef$n_chains, n.adapt = ef$n_adapt, quiet = TRUE)
  if (ef$n_burnin > 0L) stats::update(jm, n.iter = ef$n_burnin, progress.bar = "none")
  draws <- as.matrix(rjags::coda.samples(jm, variable.names = monitors, n.iter = ef$n_iter, thin = ef$thin, progress.bar = "none"))
  p_cols <- paste0("p_regular[", seq_len(ndose), "]"); if (!all(p_cols %in% colnames(draws))) stop("The efficacy JAGS model did not return all regular-dose draws.")
  p_draws <- draws[, p_cols, drop = FALSE]
  if (additive) { alpha_cols <- if (model == "previous_dose_additive") "alpha" else paste0("alpha[", seq_len(ndose), "]"); carry_draws <- if (length(alpha_cols) == 1L) matrix(draws[, alpha_cols], ncol = ndose, nrow = nrow(draws)) else draws[, alpha_cols, drop = FALSE]; ipde_draws <- matrix(pmin(1, p_draws + carry_draws * p_draws), nrow = nrow(p_draws), ncol = ncol(p_draws))
  } else { r_cols <- if (model == "shared_carryover") "r_e" else paste0("r_e[", seq_len(ndose), "]"); carry_draws <- if (length(r_cols) == 1L) matrix(draws[, r_cols], ncol = ndose, nrow = nrow(draws)) else draws[, r_cols, drop = FALSE]; ipde_draws <- carry_draws + (1 - carry_draws) * p_draws }
  list(model = model, p_regular_mean = colMeans(p_draws), carryover_mean = colMeans(carry_draws),
       p_ipde_mean = colMeans(ipde_draws),
       prob_below_threshold = futility$prob_below_threshold,
       futile = futility$futile,
       n_fully_ascertained_by_dose = futility$n_fully_ascertained_by_dose,
       y_fully_ascertained_by_dose = futility$y_fully_ascertained_by_dose,
       futility_prior = futility$prior)
}

aide_update_futility <- function(efficacy_fit, previous) list(
  futile = previous | efficacy_fit$futile,
  n_fully_ascertained_by_dose = efficacy_fit$n_fully_ascertained_by_dose,
  y_fully_ascertained_by_dose = efficacy_fit$y_fully_ascertained_by_dose,
  posterior_futility_probabilities = efficacy_fit$prob_below_threshold
)

aide_compute_utility <- function(p_e, p_t, settings) {
  utility_type <- as.integer(settings$type)
  scores <- as.numeric(settings$scores)
  if (length(scores) != 4L || any(!is.finite(scores)))
    stop("utility$scores must contain four finite values: u00, u01, u10, and u11.")
  names(scores) <- c("u00", "u01", "u10", "u11")
  switch(
    as.character(utility_type),
    "1" = p_e,
    "2" = p_e - settings$lambda_T * p_t,
    "3" = scores[["u00"]] * (1 - p_t) * (1 - p_e) +
      scores[["u01"]] * (1 - p_t) * p_e +
      scores[["u10"]] * p_t * (1 - p_e) +
      scores[["u11"]] * p_t * p_e,
    stop("utility$type must be 1, 2, or 3.")
  )
}

aide_toxicity_recommendation <- function(current_dose, toxicity_fit, config, eliminated) {
  ndose <- length(toxicity_fit$p_regular_mean)
  if (eliminated[1L]) return(list(action = "stop", recommended_dose = NA_integer_, stop = TRUE))
  if (eliminated[current_dose] || toxicity_fit$prob_overtox_by_dose[current_dose] > config$toxicity$cutoff)
    return(list(action = "de_escalate", recommended_dose = max(1L, current_dose - 1L), stop = FALSE))
  local <- seq.int(max(1L, current_dose - 1L), min(ndose, current_dose + 1L))
  local <- local[!eliminated[local]]
  if (!length(local)) return(list(action = "stop", recommended_dose = NA_integer_, stop = TRUE))
  recommended <- aide_select_lowest_tied(-abs(toxicity_fit$p_regular_mean - config$toxicity$target), local)
  list(action = aide_action_from_next(recommended, current_dose), recommended_dose = recommended, stop = FALSE)
}
