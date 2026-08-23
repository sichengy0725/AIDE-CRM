aide_integral_controls <- function(config) {
  x <- config$integration
  list(
    r_grid_n = x$r_grid_n,
    beta_lower = x$beta_lower,
    beta_upper = x$beta_upper,
    rel_tol = x$rel_tol,
    eps = 1e-12
  )
}

aide_integral_beta_grid <- function(a, b, n) {
  stats::qbeta((seq_len(n) - 0.5) / n, a, b)
}

aide_integral_value <- function(f, lower, upper, rel_tol) {
  stats::integrate(f, lower = lower, upper = upper, rel.tol = rel_tol)$value
}

aide_integral_toxicity_density <- function(beta, r, dat, skeleton, prior_mean, prior_sd, eps) {
  vapply(beta, function(x) {
    p <- skeleton^exp(x)
    theta_ipde <- r + (1 - r) * p
    theta <- ifelse(dat$ipde == 1L, theta_ipde[dat$dose], p[dat$dose])
    theta <- pmin(1 - eps, pmax(eps, dat$weight * theta))
    exp(stats::dnorm(x, prior_mean, prior_sd, log = TRUE) +
          sum(dat$y * log(theta) + (1 - dat$y) * log1p(-theta)))
  }, numeric(1))
}

aide_fit_toxicity_integral <- function(interim, config, ndose) {
  tx <- config$toxicity
  ctl <- aide_integral_controls(config)
  skeleton <- tx$skeleton %||% aide_default_toxicity_skeleton(ndose, tx$target)
  dat <- interim
  if (!nrow(dat)) dat <- data.frame(dose = 1L, ipde = 0L, y = 0L, weight = 0)
  dat$y <- ifelse(is.na(dat$y), 0L, as.integer(dat$y))
  dat$weight[dat$y == 1L] <- 1
  r_grid <- aide_integral_beta_grid(tx$carryover_prior[1L], tx$carryover_prior[2L], ctl$r_grid_n)
  marginal <- numeric(ctl$r_grid_n)
  regular <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)
  ipde <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)
  overtox <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)
  ipde_overtox <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)

  for (g in seq_len(ctl$r_grid_n)) {
    r <- r_grid[g]
    density <- function(beta) aide_integral_toxicity_density(
      beta, r, dat, skeleton, tx$beta_prior_mean, tx$beta_prior_sd, ctl$eps
    )
    marginal[g] <- aide_integral_value(density, ctl$beta_lower, ctl$beta_upper, ctl$rel_tol)
    for (j in seq_len(ndose)) {
      regular[g, j] <- aide_integral_value(
        function(beta) density(beta) * skeleton[j]^exp(beta),
        ctl$beta_lower, ctl$beta_upper, ctl$rel_tol
      ) / marginal[g]
      ipde[g, j] <- aide_integral_value(
        function(beta) density(beta) * (r + (1 - r) * skeleton[j]^exp(beta)),
        ctl$beta_lower, ctl$beta_upper, ctl$rel_tol
      ) / marginal[g]
      overtox[g, j] <- aide_integral_value(
        function(beta) density(beta) * as.numeric(skeleton[j]^exp(beta) > tx$target),
        ctl$beta_lower, ctl$beta_upper, ctl$rel_tol
      ) / marginal[g]
      ipde_overtox[g, j] <- aide_integral_value(
        function(beta) density(beta) * as.numeric(r + (1 - r) * skeleton[j]^exp(beta) > tx$target),
        ctl$beta_lower, ctl$beta_upper, ctl$rel_tol
      ) / marginal[g]
    }
  }

  posterior_r <- marginal / sum(marginal)
  prob_overtox <- colSums(regular * 0 + overtox * posterior_r)
  prob_ipde_overtox <- colSums(ipde * 0 + ipde_overtox * posterior_r)
  list(
    model = "integral_discount_r",
    p_regular_mean = colSums(regular * posterior_r),
    carryover_mean = rep(sum(r_grid * posterior_r), ndose),
    p_ipde_mean = colSums(ipde * posterior_r),
    prob_overtox_by_dose = prob_overtox,
    prob_ipde_overtox_by_dose = prob_ipde_overtox,
    eliminated = aide_toxicity_elimination(prob_overtox, tx$cutoff),
    skeleton = skeleton,
    posterior_carryover_mean = sum(r_grid * posterior_r)
  )
}

aide_integral_efficacy_density <- function(p, r, rows, eps) {
  vapply(p, function(x) {
    theta <- ifelse(rows$ipde == 1L, r + (1 - r) * x, x)
    theta <- pmin(1 - eps, pmax(eps, theta))
    nonresponse_power <- rows$ascertained * (1 - rows$y) + (1 - rows$ascertained) * rows$weight
    exp(sum(rows$ascertained * rows$y * log(theta) + nonresponse_power * log1p(-theta)))
  }, numeric(1))
}

aide_fit_efficacy_integral <- function(interim, config, ndose) {
  ef <- config$efficacy
  ctl <- aide_integral_controls(config)
  regular_prior <- aide_expand_beta_prior(ef$prior, ndose, "efficacy$prior")
  carry_prior <- aide_expand_beta_prior(ef$carryover_prior, ndose, "efficacy$carryover_prior")
  dat <- interim
  if (!nrow(dat)) dat <- data.frame(dose = 1L, ipde = 0L, y = 0L, ascertained = 0L, weight = 0)
  dat$y <- ifelse(is.na(dat$y), 0L, as.integer(dat$y))
  dat$ascertained <- as.integer(dat$ascertained)
  r_grid <- aide_integral_beta_grid(carry_prior[1L, 1L], carry_prior[1L, 2L], ctl$r_grid_n)
  marginal_by_r <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)
  regular_by_r <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)
  ipde_by_r <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)
  below_by_r <- matrix(0, nrow = ctl$r_grid_n, ncol = ndose)

  for (g in seq_len(ctl$r_grid_n)) {
    for (j in seq_len(ndose)) {
      r <- r_grid[g]
      rows <- dat[dat$dose == j, , drop = FALSE]
      density <- function(p) {
        aide_integral_efficacy_density(p, r, rows, ctl$eps) *
          stats::dbeta(p, regular_prior[j, 1L], regular_prior[j, 2L])
      }
      marginal_by_r[g, j] <- aide_integral_value(density, ctl$eps, 1 - ctl$eps, ctl$rel_tol)
      regular_by_r[g, j] <- aide_integral_value(
        function(p) density(p) * p, ctl$eps, 1 - ctl$eps, ctl$rel_tol
      ) / marginal_by_r[g, j]
      ipde_by_r[g, j] <- aide_integral_value(
        function(p) density(p) * (r + (1 - r) * p), ctl$eps, 1 - ctl$eps, ctl$rel_tol
      ) / marginal_by_r[g, j]
      below_by_r[g, j] <- aide_integral_value(
        density, ctl$eps, min(ef$threshold, 1 - ctl$eps), ctl$rel_tol
      ) / marginal_by_r[g, j]
    }
  }

  log_marginal <- rowSums(log(marginal_by_r))
  posterior_r <- exp(log_marginal - max(log_marginal))
  posterior_r <- posterior_r / sum(posterior_r)
  full_n <- full_y <- integer(ndose)
  for (j in seq_len(ndose)) {
    x <- interim[interim$dose == j & interim$ascertained, , drop = FALSE]
    full_n[j] <- nrow(x)
    full_y[j] <- sum(x$y)
  }
  below <- colSums(below_by_r * posterior_r)
  list(
    model = "integral_shared_carryover",
    p_regular_mean = colSums(regular_by_r * posterior_r),
    carryover_mean = rep(sum(r_grid * posterior_r), ndose),
    p_ipde_mean = colSums(ipde_by_r * posterior_r),
    prob_below_threshold = below,
    futile = full_n >= ef$min_eff_n_for_futility & below > ef$futility_cutoff,
    n_fully_ascertained_by_dose = full_n
  )
}
