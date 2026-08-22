validation_script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg)) return(dirname(normalizePath(file_arg[1L], winslash = "/")))
  source_files <- vapply(sys.frames(), function(frame) {
    value <- frame$ofile %||% ""
    if (length(value)) as.character(value)[1L] else ""
  }, character(1))
  source_file <- tail(source_files[nzchar(source_files)], 1L)
  if (length(source_file)) dirname(normalizePath(source_file, winslash = "/")) else normalizePath(getwd(), winslash = "/")
}

validation_source_tite_aide <- function() {
  validation_dir <- validation_script_directory()
  source(file.path(dirname(validation_dir), "TITE-AIDE-Rebuild", "TITE-AIDE.R"))
  invisible(validation_dir)
}

validation_complete_toxicity_model <- function(model) {
  if (model == "discount_r") {
    return("model {
      alpha ~ dnorm(beta_mean, tau_alpha)
      g_r_1 ~ dgamma(a_r, 1)
      g_r_2 ~ dgamma(b_r, 1)
      r <- g_r_1 / (g_r_1 + g_r_2)
      for (j in 1:J) {
        p[j] <- pow(q[j], exp(alpha))
        theta_ipde[j] <- r + (1 - r) * p[j]
      }
      for (i in 1:N) {
        theta_base[i] <- (1 - ipde[i]) * p[dose[i]] + ipde[i] * theta_ipde[dose[i]]
        y[i] ~ dbern(theta_base[i])
      }
    }")
  }
  if (model == "multicycle_additive") {
    return("model {
      beta ~ dnorm(beta_mean, tau_beta)
      g_alpha_1 ~ dgamma(a_r, 1)
      g_alpha_2 ~ dgamma(b_r, 1)
      alpha <- g_alpha_1 / (g_alpha_1 + g_alpha_2)
      for (j in 1:J) {
        p[j] <- pow(q[j], exp(beta))
        theta_ipde[j] <- min(1, p[j] + alpha * p[j])
      }
      state_tox[1] <- 0
      for (i in 1:N) {
        state_tox[i + 1] <- p[dose[i]] + alpha * state_tox[previous_state_index[i]]
        theta[i] <- min(1, state_tox[i + 1])
        y[i] ~ dbern(theta[i])
      }
    }")
  }
  "model {
    beta ~ dnorm(beta_mean, tau_beta)
    g_alpha_1 ~ dgamma(a_r, 1)
    g_alpha_2 ~ dgamma(b_r, 1)
    alpha <- g_alpha_1 / (g_alpha_1 + g_alpha_2)
    for (j in 1:J) {
      p[j] <- pow(q[j], exp(beta))
      theta_ipde[j] <- min(1, p[j] + alpha * p[max(1, j - 1)])
    }
    for (i in 1:N) {
      theta_regular[i] <- p[dose[i]]
      theta_retreated[i] <- min(1, p[dose[i]] + alpha * p[previous_dose[i]])
      theta_base[i] <- (1 - ipde[i]) * theta_regular[i] + ipde[i] * theta_retreated[i]
      y[i] ~ dbern(theta_base[i])
    }
  }"
}

validation_complete_efficacy_model <- function(model) {
  if (model == "dose_specific_carryover") {
    return("model {
      for (j in 1:ndose) {
        g_regular_1[j] ~ dgamma(a_regular[j], 1)
        g_regular_2[j] ~ dgamma(b_regular[j], 1)
        p_regular[j] <- g_regular_1[j] / (g_regular_1[j] + g_regular_2[j])
        g_carry_1[j] ~ dgamma(a_carry[j], 1)
        g_carry_2[j] ~ dgamma(b_carry[j], 1)
        r_e[j] <- g_carry_1[j] / (g_carry_1[j] + g_carry_2[j])
        p_ipde[j] <- r_e[j] + (1 - r_e[j]) * p_regular[j]
      }
      for (i in 1:N) {
        theta_eff[i] <- (1 - ipde[i]) * p_regular[dose[i]] + ipde[i] * p_ipde[dose[i]]
        eff[i] ~ dbern(theta_eff[i])
      }
    }")
  }
  if (model == "shared_carryover") {
    return("model {
      g_carry_1 ~ dgamma(a_carry, 1)
      g_carry_2 ~ dgamma(b_carry, 1)
      r_e <- g_carry_1 / (g_carry_1 + g_carry_2)
      for (j in 1:ndose) {
        g_regular_1[j] ~ dgamma(a_regular[j], 1)
        g_regular_2[j] ~ dgamma(b_regular[j], 1)
        p_regular[j] <- g_regular_1[j] / (g_regular_1[j] + g_regular_2[j])
        p_ipde[j] <- r_e + (1 - r_e) * p_regular[j]
      }
      for (i in 1:N) {
        theta_eff[i] <- (1 - ipde[i]) * p_regular[dose[i]] + ipde[i] * p_ipde[dose[i]]
        eff[i] ~ dbern(theta_eff[i])
      }
    }")
  }
  if (model == "previous_dose_additive") {
    return("model {
      g_alpha_1 ~ dgamma(a_alpha, 1)
      g_alpha_2 ~ dgamma(b_alpha, 1)
      alpha <- g_alpha_1 / (g_alpha_1 + g_alpha_2)
      for (j in 1:ndose) {
        g_regular_1[j] ~ dgamma(a_regular[j], 1)
        g_regular_2[j] ~ dgamma(b_regular[j], 1)
        p_regular[j] <- g_regular_1[j] / (g_regular_1[j] + g_regular_2[j])
      }
      for (i in 1:N) {
        theta_ipde[i] <- min(1, p_regular[dose[i]] + alpha * p_regular[previous_dose[i]])
        theta_eff[i] <- (1 - ipde[i]) * p_regular[dose[i]] + ipde[i] * theta_ipde[i]
        eff[i] ~ dbern(theta_eff[i])
      }
    }")
  }
  if (model == "multicycle_additive") {
    return("model {
      g_alpha_1 ~ dgamma(a_alpha, 1)
      g_alpha_2 ~ dgamma(b_alpha, 1)
      alpha <- g_alpha_1 / (g_alpha_1 + g_alpha_2)
      for (j in 1:ndose) {
        g_regular_1[j] ~ dgamma(a_regular[j], 1)
        g_regular_2[j] ~ dgamma(b_regular[j], 1)
        p_regular[j] <- g_regular_1[j] / (g_regular_1[j] + g_regular_2[j])
      }
      state_eff[1] <- 0
      for (i in 1:N) {
        state_eff[i + 1] <- p_regular[dose[i]] + alpha * state_eff[previous_state_index[i]]
        theta_eff[i] <- min(1, state_eff[i + 1])
        eff[i] ~ dbern(theta_eff[i])
      }
    }")
  }
  "model {
    for (j in 1:ndose) {
      g_alpha_1[j] ~ dgamma(a_alpha[j], 1)
      g_alpha_2[j] ~ dgamma(b_alpha[j], 1)
      alpha[j] <- g_alpha_1[j] / (g_alpha_1[j] + g_alpha_2[j])
      g_regular_1[j] ~ dgamma(a_regular[j], 1)
      g_regular_2[j] ~ dgamma(b_regular[j], 1)
      p_regular[j] <- g_regular_1[j] / (g_regular_1[j] + g_regular_2[j])
    }
    for (i in 1:N) {
      theta_ipde[i] <- min(1, p_regular[dose[i]] + alpha[dose[i]] * p_regular[previous_dose[i]])
      theta_eff[i] <- (1 - ipde[i]) * p_regular[dose[i]] + ipde[i] * theta_ipde[i]
      eff[i] ~ dbern(theta_eff[i])
    }
  }"
}

validation_fit_complete_toxicity <- function(interim, config, ndose) {
  tx <- config$toxicity
  skeleton <- tx$skeleton %||% aide_default_toxicity_skeleton(ndose, tx$target)
  carry_prior <- as.numeric(tx$carryover_prior)
  dat <- interim
  multicycle <- identical(tx$model, "multicycle_additive")
  if (multicycle) dat <- aide_multicycle_history_data(dat)
  dat$y <- as.integer(dat$y)
  jags_data <- list(
    N = nrow(dat), J = ndose, y = dat$y, dose = as.integer(dat$dose),
    ipde = as.integer(dat$ipde), q = as.numeric(skeleton),
    beta_mean = tx$beta_prior_mean, a_r = carry_prior[1L], b_r = carry_prior[2L]
  )
  if (multicycle) {
    jags_data$previous_state_index <- as.integer(dat$previous_state_index)
    jags_data$tau_beta <- 1 / tx$beta_prior_sd^2
    monitors <- c("p", "theta_ipde", "alpha")
  } else if (tx$model == "additive_alpha") {
    jags_data$previous_dose <- as.integer(dat$previous_dose)
    jags_data$tau_beta <- 1 / tx$beta_prior_sd^2
    monitors <- c("p", "theta_ipde", "alpha")
  } else {
    jags_data$tau_alpha <- 1 / tx$beta_prior_sd^2
    monitors <- c("p", "theta_ipde", "r")
  }
  model_connection <- textConnection(validation_complete_toxicity_model(tx$model))
  jm <- rjags::jags.model(model_connection, data = jags_data, n.chains = tx$n_chains, n.adapt = tx$n_adapt, quiet = TRUE)
  close(model_connection)
  if (tx$n_burnin > 0L) stats::update(jm, n.iter = tx$n_burnin, progress.bar = "none")
  draws <- as.matrix(rjags::coda.samples(jm, variable.names = monitors, n.iter = tx$n_iter, thin = tx$thin, progress.bar = "none"))
  p_cols <- paste0("p[", seq_len(ndose), "]")
  ipde_cols <- paste0("theta_ipde[", seq_len(ndose), "]")
  carry_draws <- if (tx$model %in% c("additive_alpha", "multicycle_additive")) draws[, "alpha"] else draws[, "r"]
  list(
    model = tx$model,
    p_regular_mean = colMeans(draws[, p_cols, drop = FALSE]),
    carryover_mean = rep(mean(carry_draws), ndose),
    p_ipde_mean = colMeans(draws[, ipde_cols, drop = FALSE]),
    prob_overtox_by_dose = colMeans(draws[, p_cols, drop = FALSE] > tx$target),
    prob_ipde_overtox_by_dose = colMeans(draws[, ipde_cols, drop = FALSE] > tx$target),
    eliminated = colMeans(draws[, p_cols, drop = FALSE] > tx$target) > tx$cutoff,
    skeleton = skeleton,
    posterior_carryover_mean = mean(carry_draws)
  )
}

validation_fit_complete_efficacy <- function(interim, config, ndose) {
  ef <- config$efficacy
  regular_prior <- aide_expand_beta_prior(ef$prior, ndose, "efficacy$prior")
  carry_prior <- aide_expand_beta_prior(ef$carryover_prior, ndose, "efficacy$carryover_prior")
  dat <- interim
  multicycle <- identical(ef$model, "multicycle_additive")
  if (multicycle) dat <- aide_multicycle_history_data(dat)
  jags_data <- list(
    N = nrow(dat), ndose = ndose, eff = as.integer(dat$y), dose = as.integer(dat$dose),
    ipde = as.integer(dat$ipde), a_regular = regular_prior[, 1L], b_regular = regular_prior[, 2L]
  )
  additive <- ef$model %in% c("previous_dose_additive", "dose_specific_previous_dose_additive", "multicycle_additive")
  if (multicycle) {
    jags_data$previous_state_index <- as.integer(dat$previous_state_index)
    jags_data$a_alpha <- carry_prior[1L, 1L]
    jags_data$b_alpha <- carry_prior[1L, 2L]
    monitors <- c("p_regular", "alpha")
  } else if (additive) {
    jags_data$previous_dose <- as.integer(dat$previous_dose)
    jags_data$a_alpha <- if (ef$model == "previous_dose_additive") carry_prior[1L, 1L] else carry_prior[, 1L]
    jags_data$b_alpha <- if (ef$model == "previous_dose_additive") carry_prior[1L, 2L] else carry_prior[, 2L]
    monitors <- c("p_regular", "alpha")
  } else if (ef$model == "shared_carryover") {
    jags_data$a_carry <- carry_prior[1L, 1L]
    jags_data$b_carry <- carry_prior[1L, 2L]
    monitors <- c("p_regular", "r_e")
  } else {
    jags_data$a_carry <- carry_prior[, 1L]
    jags_data$b_carry <- carry_prior[, 2L]
    monitors <- c("p_regular", "r_e")
  }
  model_connection <- textConnection(validation_complete_efficacy_model(ef$model))
  jm <- rjags::jags.model(model_connection, data = jags_data, n.chains = ef$n_chains, n.adapt = ef$n_adapt, quiet = TRUE)
  close(model_connection)
  if (ef$n_burnin > 0L) stats::update(jm, n.iter = ef$n_burnin, progress.bar = "none")
  draws <- as.matrix(rjags::coda.samples(jm, variable.names = monitors, n.iter = ef$n_iter, thin = ef$thin, progress.bar = "none"))
  p_cols <- paste0("p_regular[", seq_len(ndose), "]")
  p_draws <- draws[, p_cols, drop = FALSE]
  if (additive) {
    alpha_cols <- if (ef$model %in% c("previous_dose_additive", "multicycle_additive")) "alpha" else paste0("alpha[", seq_len(ndose), "]")
    carry_draws <- if (length(alpha_cols) == 1L) matrix(draws[, alpha_cols], ncol = ndose, nrow = nrow(draws)) else draws[, alpha_cols, drop = FALSE]
    ipde_draws <- matrix(pmin(1, p_draws + carry_draws * p_draws), nrow = nrow(p_draws), ncol = ncol(p_draws))
  } else {
    r_cols <- if (ef$model == "shared_carryover") "r_e" else paste0("r_e[", seq_len(ndose), "]")
    carry_draws <- if (length(r_cols) == 1L) matrix(draws[, r_cols], ncol = ndose, nrow = nrow(draws)) else draws[, r_cols, drop = FALSE]
    ipde_draws <- carry_draws + (1 - carry_draws) * p_draws
  }
  full_n <- full_y <- integer(ndose)
  for (j in seq_len(ndose)) {
    x <- interim[interim$dose == j, , drop = FALSE]
    full_n[j] <- nrow(x)
    full_y[j] <- sum(x$y)
  }
  below <- colMeans(p_draws < ef$threshold)
  list(
    model = ef$model,
    p_regular_mean = colMeans(p_draws),
    carryover_mean = colMeans(carry_draws),
    p_ipde_mean = colMeans(ipde_draws),
    prob_below_threshold = below,
    futile = full_n >= ef$min_eff_n_for_futility & below > ef$futility_cutoff,
    n_fully_ascertained_by_dose = full_n
  )
}

validation_decision_state <- function(admin, ndose, current_dose) {
  list(
    admin = admin,
    eliminated = rep(FALSE, ndose),
    futile = rep(FALSE, ndose),
    current_dose = current_dose,
    stage = list(active = "one_stage")
  )
}

validation_final_selection <- function(admin, toxicity_fit, efficacy_fit, config) {
  tried <- sort(unique(admin$dose))
  safe_tried <- tried[!toxicity_fit$eliminated[tried]]
  MTD <- safe_tried[which.min(abs(toxicity_fit$p_regular_mean[safe_tried] - config$toxicity$target))]
  utility <- aide_compute_utility(efficacy_fit$p_regular_mean, toxicity_fit$p_regular_mean, config$utility)
  obd_candidates <- safe_tried[!efficacy_fit$futile[safe_tried] & safe_tried <= MTD]
  OBD <- obd_candidates[which.max(efficacy_fit$p_regular_mean[obd_candidates])]
  list(MTD = MTD, OBD = OBD, utilities = utility)
}

validation_fit_summary <- function(fit, endpoint, method, toxicity_model, efficacy_model) {
  data.frame(
    toxicity_model = toxicity_model,
    efficacy_model = efficacy_model,
    endpoint = endpoint,
    method = method,
    dose = seq_along(fit$p_regular_mean),
    p_regular_mean = fit$p_regular_mean,
    p_ipde_mean = fit$p_ipde_mean,
    carryover_mean = fit$carryover_mean,
    stringsAsFactors = FALSE
  )
}
