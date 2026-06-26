## ============================================================
## AIDE-CRM helper with five CRM backends
##   1) fixed/r_fixed : discount CRM with fixed r through JAGS
##   2) random        : discount CRM with one random r through JAGS
##   3) level         : discount CRM with level-specific random r through JAGS
##   4) alpha_crm     : deterministic alpha-CRM using alpha grid + integrate()
##   5) cumu_crm      : logistic cumulative-dose CRM through JAGS
##
## User-prespecified inputs:
##   - skeleton      : prior CRM toxicity skeleton, length ndose
##   - dose_values   : actual/scaled dose values for alpha cumulative dose
##   - dose_scores   : dose scores for logistic cumulative CRM
##   - all priors    : exposed in crm_fit(), crm_move(), select.mtd.crm()
##
## Data requirements:
##   Required columns: dose and y, or dose and y_obs
##   Optional columns : type, id, cycle, arrival_time/t_start/t_arrival/t_eval
##   For alpha_crm, arrival time is used if available; otherwise row order is used.
## ============================================================

## ------------------------------------------------------------
## Shared small utilities
## ------------------------------------------------------------
crm_inv_logit <- function(x) 1 / (1 + exp(-x))

crm_r_model_choices <- function() c("fixed", "r_fixed", "random", "level", "alpha_crm", "cumu_crm")

crm_normalize_r_model <- function(r_model) {
  r_model <- match.arg(r_model, choices = crm_r_model_choices())
  if (r_model == "r_fixed") "fixed" else r_model
}

crm_validate_skeleton <- function(skeleton, ndose) {
  if (length(skeleton) != ndose) {
    stop("skeleton must have length ndose.")
  }
  if (any(!is.finite(skeleton)) || any(skeleton <= 0 | skeleton >= 1)) {
    stop("skeleton must contain probabilities in (0, 1).")
  }
  invisible(TRUE)
}

crm_default_dose_values <- function(ndose) seq_len(ndose)

crm_default_dose_scores <- function(skeleton) {
  ## Matches the user's original model.option 0/1 convention:
  ## doses <- logit(prior) - 3
  stats::qlogis(skeleton) - 3
}

crm_prepare_dat <- function(dat,
                            ndose,
                            time_col = NULL,
                            require_time = FALSE) {
  if (is.null(dat) || nrow(dat) == 0L) return(dat)

  if (!"dose" %in% names(dat)) {
    stop("dat must contain column: dose.")
  }

  y_col <- if ("y" %in% names(dat)) {
    "y"
  } else if ("y_obs" %in% names(dat)) {
    "y_obs"
  } else {
    stop("dat must contain either column y or y_obs.")
  }

  out <- data.frame(
    dose = as.integer(dat$dose),
    y = as.integer(dat[[y_col]]),
    stringsAsFactors = FALSE
  )

  if (any(is.na(out$y)) || any(!out$y %in% c(0L, 1L))) {
    stop("dat$y or dat$y_obs must be 0/1 with no missing values.")
  }
  if (any(is.na(out$dose)) || any(out$dose < 1L | out$dose > ndose)) {
    stop("dat$dose must be integer dose levels from 1 to ndose.")
  }

  out$type <- if ("type" %in% names(dat)) as.character(dat$type) else "new"
  out$id <- if ("id" %in% names(dat)) dat$id else seq_len(nrow(dat))

  if ("cycle" %in% names(dat)) {
    out$cycle <- dat$cycle
  } else if ("ncycle" %in% names(dat)) {
    out$cycle <- dat$ncycle
  } else {
    out$cycle <- ave(seq_len(nrow(dat)), out$id, FUN = seq_along)
  }

  if (!is.null(time_col) && time_col %in% names(dat)) {
    out$arrival_time <- as.numeric(dat[[time_col]])
  } else {
    candidate_time_cols <- c("arrival_time", "t_start", "t_arrival", "t_eval")
    found_time_cols <- candidate_time_cols[candidate_time_cols %in% names(dat)]
    if (length(found_time_cols) > 0L) {
      out$arrival_time <- as.numeric(dat[[found_time_cols[1L]]])
    } else {
      if (require_time) {
        warning(
          "No arrival-time column found. alpha_crm will use within-dataset row order as arrival_time. ",
          "Pass time_col=... if you want calendar-time carryover."
        )
      }
      out$arrival_time <- as.numeric(seq_len(nrow(dat)) - 1L)
    }
  }

  if (any(!is.finite(out$arrival_time))) {
    stop("arrival_time/time_col contains non-finite values.")
  }

  out <- out[order(out$id, out$cycle, out$arrival_time), , drop = FALSE]
  rownames(out) <- NULL
  out
}

crm_make_result <- function(p_hat,
                            target,
                            skeleton,
                            r_model,
                            post = NULL,
                            theta_ipde_hat = NULL,
                            r_hat = NA_real_,
                            prob_overtox = NA_real_,
                            stop = NA_integer_,
                            extra = list()) {
  ndose <- length(skeleton)
  if (is.null(theta_ipde_hat)) theta_ipde_hat <- rep(NA_real_, ndose)
  out <- list(
    p_hat = as.numeric(p_hat),
    theta_ipde_hat = as.numeric(theta_ipde_hat),
    r_hat = as.numeric(r_hat),
    prob_overtox = as.numeric(prob_overtox),
    stop = as.integer(stop),
    post = post,
    r_model = r_model,
    target = target,
    skeleton = skeleton
  )
  c(out, extra)
}

## ------------------------------------------------------------
## Dispatch fitter: call this when the CRM backend is user-selected
## ------------------------------------------------------------
crm_fit <- function(dat,
                    ndose,
                    skeleton,
                    target = 0.30,
                    cutoff = 0.95,
                    r_model = crm_r_model_choices(),
                    r_carry = 0.10,
                    a_r = 1,
                    b_r = 9,
                    alpha_sd = 2,
                    model_file = NULL,
                    fixed_model_file = "fix_CRM.bug",
                    random_model_file = "random_CRM.bug",
                    level_model_file = "random_CRM_level.bug",
                    n_chains = 2,
                    n_adapt = 500,
                    n_burnin = 500,
                    n_iter = 2000,
                    thin = 1,
                    seed = NULL,
                    dose_values = NULL,
                    dose_scores = NULL,
                    time_col = NULL,
                    alpha_grid = seq(0.01, 0.99, length.out = 61),
                    alpha_T = 28,
                    theta_prior_mean = 0,
                    theta_prior_sd = sqrt(2),
                    alpha_L = 8,
                    alpha_rel_tol = 1e-8,
                    alpha_eps = 1e-12,
                    alpha_n_draw_prior = 5000,
                    cumu_model_file = NULL,
                    cumu_beta0_mean = stats::qlogis(target),
                    cumu_beta0_prec = 4,
                    cumu_beta0_df = 1,
                    cumu_beta1_shape = 5.83,
                    cumu_beta1_rate = 1.21,
                    cumu_beta2_rate = 1,
                    cumu_include_current = FALSE) {

  r_model <- crm_normalize_r_model(r_model)
  crm_validate_skeleton(skeleton, ndose)

  if (r_model %in% c("fixed", "random", "level")) {
    return(crm_fit_discount(
      dat = dat,
      ndose = ndose,
      skeleton = skeleton,
      target = target,
      cutoff = cutoff,
      r_model = r_model,
      r_carry = r_carry,
      a_r = a_r,
      b_r = b_r,
      alpha_sd = alpha_sd,
      model_file = model_file,
      fixed_model_file = fixed_model_file,
      random_model_file = random_model_file,
      level_model_file = level_model_file,
      n_chains = n_chains,
      n_adapt = n_adapt,
      n_burnin = n_burnin,
      n_iter = n_iter,
      thin = thin,
      seed = seed
    ))
  }

  if (r_model == "alpha_crm") {
    if (is.null(dose_values)) dose_values <- crm_default_dose_values(ndose)
    return(crm_fit_alpha_integrate(
      dat = dat,
      ndose = ndose,
      skeleton = skeleton,
      dose_values = dose_values,
      target = target,
      cutoff = cutoff,
      T = alpha_T,
      eps = alpha_eps,
      alpha_grid = alpha_grid,
      theta_prior_mean = theta_prior_mean,
      theta_prior_sd = theta_prior_sd,
      L = alpha_L,
      rel.tol = alpha_rel_tol,
      n_draw_prior = alpha_n_draw_prior,
      time_col = time_col,
      seed = seed
    ))
  }

  if (r_model == "cumu_crm") {
    if (is.null(dose_scores)) dose_scores <- crm_default_dose_scores(skeleton)
    return(crm_fit_cumu_jags(
      dat = dat,
      ndose = ndose,
      skeleton = skeleton,
      dose_scores = dose_scores,
      target = target,
      cutoff = cutoff,
      model_file = cumu_model_file,
      beta0_mean = cumu_beta0_mean,
      beta0_prec = cumu_beta0_prec,
      beta0_df = cumu_beta0_df,
      beta1_shape = cumu_beta1_shape,
      beta1_rate = cumu_beta1_rate,
      beta2_rate = cumu_beta2_rate,
      include_current = cumu_include_current,
      n_chains = n_chains,
      n_adapt = n_adapt,
      n_burnin = n_burnin,
      n_iter = n_iter,
      thin = thin,
      seed = seed
    ))
  }
}

## ============================================================
## 1/2) Discount CRM through JAGS: fixed r and random r
## ============================================================
crm_fit_discount <- function(dat,
                             ndose,
                             skeleton,
                             target = 0.30,
                             cutoff = 0.95,
                             cutoff.eli = NULL,
                             r_model = c("fixed", "r_fixed", "random", "level"),
                             r_carry = 0.10,
                             a_r = 1,
                             b_r = 9,
                             alpha_sd = 2,
                             model_file = NULL,
                             fixed_model_file = "fix_CRM.bug",
                             random_model_file = "random_CRM.bug",
                             level_model_file = "random_CRM_level.bug",
                             n_chains = 2,
                             n_adapt = 500,
                             n_burnin = 500,
                             n_iter = 2000,
                             thin = 1,
                             seed = NULL) {

  r_model <- crm_normalize_r_model(r_model)
  if (!is.null(cutoff.eli)) cutoff <- cutoff.eli
  crm_validate_skeleton(skeleton, ndose)

  if (is.null(model_file)) {
    model_file <- switch(
      r_model,
      fixed  = fixed_model_file,
      random = random_model_file,
      level  = level_model_file
    )
  }

  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required for CRM. Install/load JAGS and rjags first.")
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for CRM posterior samples.")
  }
  if (alpha_sd <= 0 || !is.finite(alpha_sd)) {
    stop("alpha_sd must be positive.")
  }
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff >= 1) {
    stop("cutoff/cutoff.eli must be a scalar in (0, 1).")
  }
  if (r_model == "level" && ndose != 5L) {
    stop("The level-specific random r model random_CRM_level.bug is written for exactly 5 dose levels.")
  }
  if (r_model == "fixed" && (length(r_carry) != 1L || r_carry < 0 || r_carry >= 1)) {
    stop("For fixed r, r_carry must be a scalar in [0, 1).")
  }

  dat2 <- crm_prepare_dat(dat, ndose)
  if (is.null(dat2) || nrow(dat2) == 0L) {
    stop("CRM requires at least one observed assignment in dat.")
  }

  y_vec <- as.integer(dat2$y)
  dose_vec <- as.integer(dat2$dose)
  ipde_vec <- as.integer(dat2$type == "retreat")

  jags_data <- list(
    N = length(y_vec),
    J = ndose,
    y = y_vec,
    dose = dose_vec,
    ipde = ipde_vec,
    q = as.numeric(skeleton),
    tau_alpha = 1 / alpha_sd^2
  )

  monitors <- c("alpha", "p", "theta_ipde")

  if (r_model == "fixed") {
    jags_data$r <- r_carry
  } else if (r_model == "random") {
    jags_data$a_r <- a_r
    jags_data$b_r <- b_r
    monitors <- c(monitors, "r")
  } else if (r_model == "level") {
    jags_data$a_r <- a_r
    jags_data$b_r <- b_r
    monitors <- c(monitors, paste0("r", 2:ndose))
  }

  if (!file.exists(model_file)) {
    stop("Cannot find JAGS model file: ", model_file)
  }

  jm <- rjags::jags.model(
    file = model_file,
    data = jags_data,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )

  stats::update(jm, n.iter = n_burnin, progress.bar = "none")

  samp <- rjags::coda.samples(
    model = jm,
    variable.names = monitors,
    n.iter = n_iter,
    thin = thin,
    progress.bar = "none"
  )

  post <- as.matrix(samp)

  p_cols <- paste0("p[", seq_len(ndose), "]")
  if (!all(p_cols %in% colnames(post))) {
    stop("Posterior samples do not contain p[1],...,p[J]. Check the JAGS model file.")
  }
  p_hat <- colMeans(post[, p_cols, drop = FALSE])

  theta_cols <- paste0("theta_ipde[", seq_len(ndose), "]")
  theta_hat <- rep(NA_real_, ndose)
  if (all(theta_cols %in% colnames(post))) {
    theta_hat <- colMeans(post[, theta_cols, drop = FALSE])
  }

  r_hat <- if (r_model == "fixed") {
    as.numeric(r_carry)
  } else if (r_model == "random") {
    as.numeric(mean(post[, "r"]))
  } else {
    r_cols <- paste0("r", 2:ndose)
    if (!all(r_cols %in% colnames(post))) {
      stop("Posterior samples do not contain r2,...,rJ. Check the level-specific CRM model file.")
    }
    out <- rep(NA_real_, ndose)
    names(out) <- paste0("D", seq_len(ndose))
    out[2:ndose] <- colMeans(post[, r_cols, drop = FALSE])
    out
  }
  prob_overtox <- mean(post[, p_cols[1L]] > target)
  stop_flag <- as.integer(prob_overtox > cutoff)

  crm_make_result(
    p_hat = p_hat,
    target = target,
    skeleton = skeleton,
    r_model = r_model,
    post = post,
    theta_ipde_hat = theta_hat,
    r_hat = r_hat,
    prob_overtox = prob_overtox,
    stop = stop_flag,
    extra = list(
      prob_p1_over_target = prob_overtox,
      earlystop = stop_flag,
      eliminated = if (stop_flag == 1L) rep(1L, ndose) else rep(0L, ndose),
      model_file = model_file
    )
  )
}

## ============================================================
## 3) Deterministic alpha-CRM via alpha-grid + integrate()
## ============================================================
interp_S_vec <- function(D, d_grid, s_grid) {
  K <- length(d_grid)
  out <- numeric(length(D))

  idx_lo <- D <= d_grid[1]
  out[idx_lo] <- s_grid[1]

  idx_hi <- D >= d_grid[K]
  if (any(idx_hi)) {
    t <- (D[idx_hi] - d_grid[K - 1L]) / (d_grid[K] - d_grid[K - 1L])
    out[idx_hi] <- pmin(1, s_grid[K - 1L] + t * (s_grid[K] - s_grid[K - 1L]))
  }

  idx_mid <- !(idx_lo | idx_hi)
  if (any(idx_mid)) {
    Dm <- D[idx_mid]
    k <- findInterval(Dm, d_grid)
    k <- pmax(1L, pmin(k, K - 1L))
    t <- (Dm - d_grid[k]) / (d_grid[k + 1L] - d_grid[k])
    out[idx_mid] <- s_grid[k] + t * (s_grid[k + 1L] - s_grid[k])
  }
  out
}

compute_Dij_alpha <- function(tmp, alpha, T, d_grid) {
  tmp <- tmp[order(tmp$id, tmp$cycle, tmp$arrival_time), ]
  Dij <- numeric(nrow(tmp))
  idx_by_id <- split(seq_len(nrow(tmp)), tmp$id)

  for (idx in idx_by_id) {
    sub <- tmp[idx, , drop = FALSE]
    for (r in seq_len(nrow(sub))) {
      t_j <- sub$arrival_time[r]
      s <- 0
      for (rp in seq_len(r)) {
        s <- s + d_grid[sub$dose[rp]] * alpha^((t_j - sub$arrival_time[rp]) / T)
      }
      Dij[idx[r]] <- s
    }
  }
  Dij
}

logpost_theta_given_alpha_vec <- function(theta,
                                          y,
                                          S_at_D,
                                          eps = 1e-12,
                                          theta_prior_mean = 0,
                                          theta_prior_sd = sqrt(2)) {
  lp <- stats::dnorm(theta, theta_prior_mean, theta_prior_sd, log = TRUE)
  et <- exp(theta)
  ll <- vapply(et, function(et1) {
    p <- S_at_D^et1
    p <- pmin(1 - eps, pmax(eps, p))
    sum(y * log(p) + (1 - y) * log1p(-p))
  }, numeric(1))
  lp + ll
}

crm_fit_alpha_integrate <- function(dat,
                                    ndose,
                                    skeleton,
                                    dose_values = NULL,
                                    target = 0.30,
                                    cutoff = 0.95,
                                    T = 28,
                                    eps = 1e-12,
                                    alpha_grid = seq(0.01, 0.99, length.out = 61),
                                    theta_prior_mean = 0,
                                    theta_prior_sd = sqrt(2),
                                    L = 8,
                                    rel.tol = 1e-8,
                                    n_draw_prior = 5000,
                                    time_col = NULL,
                                    seed = NULL) {

  crm_validate_skeleton(skeleton, ndose)
  if (is.null(dose_values)) dose_values <- crm_default_dose_values(ndose)
  if (length(dose_values) != ndose) stop("dose_values must have length ndose.")
  if (any(!is.finite(dose_values))) stop("dose_values must be finite.")
  if (is.unsorted(dose_values, strictly = TRUE)) {
    stop("dose_values must be strictly increasing for alpha_crm interpolation.")
  }
  if (length(alpha_grid) < 1L || any(alpha_grid <= 0 | alpha_grid >= 1)) {
    stop("alpha_grid must contain values in (0, 1).")
  }
  if (theta_prior_sd <= 0 || !is.finite(theta_prior_sd)) {
    stop("theta_prior_sd must be positive.")
  }
  if (is.null(dat) || nrow(dat) == 0L) {
    theta <- stats::rnorm(n_draw_prior, theta_prior_mean, theta_prior_sd)
    posttox <- vapply(seq_along(skeleton), function(k) mean(skeleton[k]^exp(theta)), 0.0)
    prob_overtox <- mean(skeleton[1L]^exp(theta) > target)
    stop_flag <- as.integer(prob_overtox > cutoff)
    return(crm_make_result(
      p_hat = posttox,
      target = target,
      skeleton = skeleton,
      r_model = "alpha_crm",
      prob_overtox = prob_overtox,
      stop = stop_flag,
      extra = list(
        alpha_grid = alpha_grid,
        alpha_weight = rep(1 / length(alpha_grid), length(alpha_grid)),
        dose_values = dose_values,
        MTD = which.min(abs(posttox - target)),
        prob_p1_over_target = prob_overtox,
        earlystop = stop_flag,
        eliminated = if (stop_flag == 1L) rep(1L, ndose) else rep(0L, ndose)
      )
    ))
  }

  tmp <- crm_prepare_dat(dat, ndose, time_col = time_col, require_time = TRUE)
  y <- as.integer(tmp$y)

  lower <- theta_prior_mean - L * theta_prior_sd
  upper <- theta_prior_mean + L * theta_prior_sd

  theta_star <- log(log(target) / log(skeleton[1L]))

  K <- length(skeleton)
  A <- length(alpha_grid)

  m_alpha <- numeric(A)
  m_lt_alpha <- numeric(A)
  cache_S <- vector("list", A)

  for (a in seq_len(A)) {
    alpha <- alpha_grid[a]
    Dij <- compute_Dij_alpha(tmp, alpha, T, dose_values)
    S_at_D <- interp_S_vec(Dij, dose_values, skeleton)
    S_at_D <- pmin(1 - eps, pmax(eps, S_at_D))
    cache_S[[a]] <- S_at_D

    f_den <- function(theta) {
      exp(logpost_theta_given_alpha_vec(
        theta = theta,
        y = y,
        S_at_D = S_at_D,
        eps = eps,
        theta_prior_mean = theta_prior_mean,
        theta_prior_sd = theta_prior_sd
      ))
    }

    m_alpha[a] <- stats::integrate(f_den, lower, upper, rel.tol = rel.tol)$value

    up2 <- min(theta_star, upper)
    m_lt_alpha[a] <- if (up2 <= lower) {
      0
    } else {
      stats::integrate(f_den, lower, up2, rel.tol = rel.tol)$value
    }
  }

  if (any(!is.finite(m_alpha)) || sum(m_alpha) <= 0) {
    stop("alpha_crm integration failed: marginal likelihood is zero/non-finite. Try smaller L or looser rel.tol.")
  }

  w <- m_alpha / sum(m_alpha)
  prob_overtox <- sum(w * (m_lt_alpha / m_alpha))
  stop_flag <- as.integer(prob_overtox > cutoff)

  posttox <- numeric(K)
  for (k in seq_len(K)) {
    Ek_alpha <- numeric(A)
    for (a in seq_len(A)) {
      S_at_D <- cache_S[[a]]
      f_num <- function(theta) {
        base <- exp(logpost_theta_given_alpha_vec(
          theta = theta,
          y = y,
          S_at_D = S_at_D,
          eps = eps,
          theta_prior_mean = theta_prior_mean,
          theta_prior_sd = theta_prior_sd
        ))
        base * (skeleton[k]^exp(theta))
      }
      num <- stats::integrate(f_num, lower, upper, rel.tol = rel.tol)$value
      Ek_alpha[a] <- num / m_alpha[a]
    }
    posttox[k] <- sum(w * Ek_alpha)
  }

  crm_make_result(
    p_hat = posttox,
    target = target,
    skeleton = skeleton,
    r_model = "alpha_crm",
    prob_overtox = prob_overtox,
    stop = stop_flag,
    extra = list(
      alpha_grid = alpha_grid,
      alpha_weight = w,
      dose_values = dose_values,
      MTD = which.min(abs(posttox - target)),
      prob_p1_over_target = prob_overtox,
      earlystop = stop_flag,
      eliminated = if (stop_flag == 1L) rep(1L, ndose) else rep(0L, ndose)
    )
  )
}

## ============================================================
## 4) Logistic cumulative-dose CRM through JAGS
## ============================================================
crm_write_cumu_model_file <- function(model_file = NULL) {
  model_string <- "
model {
  for (i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- beta0 + beta1 * assigned_d[i] + beta2 * cumu_d[i]
  }

  beta0 ~ dt(beta0_mean, beta0_prec, beta0_df)
  beta1 ~ dgamma(beta1_shape, beta1_rate)
  beta2 ~ dexp(beta2_rate)

  for (j in 1:J) {
    logit(p_base[j]) <- beta0 + beta1 * dose_scores[j]
  }
}
"
  if (is.null(model_file)) {
    model_file <- tempfile(pattern = "modelCRM_logistic_cumu_", fileext = ".bug")
  }
  writeLines(model_string, con = model_file)
  model_file
}

crm_compute_cumu_d <- function(dat2,
                               dose_scores,
                               include_current = FALSE) {
  dat2 <- dat2[order(dat2$id, dat2$cycle, dat2$arrival_time), , drop = FALSE]
  out <- numeric(nrow(dat2))
  idx_by_id <- split(seq_len(nrow(dat2)), dat2$id)

  for (idx in idx_by_id) {
    sub <- dat2[idx, , drop = FALSE]
    cum <- 0
    for (r in seq_len(nrow(sub))) {
      current_score <- dose_scores[sub$dose[r]]
      if (include_current) {
        out[idx[r]] <- cum + current_score
      } else {
        out[idx[r]] <- cum
      }
      cum <- cum + current_score
    }
  }
  out
}

crm_fit_cumu_jags <- function(dat,
                              ndose,
                              skeleton,
                              dose_scores = NULL,
                              target = 0.30,
                              cutoff = 0.95,
                              model_file = NULL,
                              beta0_mean = stats::qlogis(target),
                              beta0_prec = 4,
                              beta0_df = 1,
                              beta1_shape = 5.83,
                              beta1_rate = 1.21,
                              beta2_rate = 1,
                              include_current = FALSE,
                              n_chains = 2,
                              n_adapt = 500,
                              n_burnin = 500,
                              n_iter = 2000,
                              thin = 1,
                              seed = NULL) {

  crm_validate_skeleton(skeleton, ndose)
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required for cumu_crm. Install/load JAGS and rjags first.")
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for cumu_crm posterior samples.")
  }
  if (is.null(dose_scores)) dose_scores <- crm_default_dose_scores(skeleton)
  if (length(dose_scores) != ndose) stop("dose_scores must have length ndose.")
  if (any(!is.finite(dose_scores))) stop("dose_scores must be finite.")

  dat2 <- crm_prepare_dat(dat, ndose)
  if (is.null(dat2) || nrow(dat2) == 0L) {
    stop("cumu_crm requires at least one observed assignment in dat.")
  }

  assigned_d <- dose_scores[dat2$dose]
  cumu_d <- crm_compute_cumu_d(dat2, dose_scores, include_current = include_current)

  if (is.null(model_file)) {
    model_file <- crm_write_cumu_model_file(NULL)
  } else if (!file.exists(model_file)) {
    model_file <- crm_write_cumu_model_file(model_file)
  }

  jags_data <- list(
    N = length(dat2$y),
    J = ndose,
    y = as.integer(dat2$y),
    assigned_d = as.numeric(assigned_d),
    cumu_d = as.numeric(cumu_d),
    dose_scores = as.numeric(dose_scores),
    beta0_mean = beta0_mean,
    beta0_prec = beta0_prec,
    beta0_df = beta0_df,
    beta1_shape = beta1_shape,
    beta1_rate = beta1_rate,
    beta2_rate = beta2_rate
  )

  jm <- rjags::jags.model(
    file = model_file,
    data = jags_data,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )

  stats::update(jm, n.iter = n_burnin, progress.bar = "none")

  samp <- rjags::coda.samples(
    model = jm,
    variable.names = c("beta0", "beta1", "beta2", "p_base"),
    n.iter = n_iter,
    thin = thin,
    progress.bar = "none"
  )

  post <- as.matrix(samp)
  p_cols <- paste0("p_base[", seq_len(ndose), "]")
  if (!all(p_cols %in% colnames(post))) {
    stop("Posterior samples do not contain p_base[1],...,p_base[J]. Check cumu CRM model.")
  }

  p_hat <- colMeans(post[, p_cols, drop = FALSE])
  prob_overtox <- mean(post[, p_cols[1L]] > target)
  stop_flag <- as.integer(prob_overtox > cutoff)

  crm_make_result(
    p_hat = p_hat,
    target = target,
    skeleton = skeleton,
    r_model = "cumu_crm",
    post = post,
    prob_overtox = prob_overtox,
    stop = stop_flag,
    extra = list(
      dose_scores = dose_scores,
      cumu_include_current = include_current,
      prob_p1_over_target = prob_overtox,
      earlystop = stop_flag,
      eliminated = if (stop_flag == 1L) rep(1L, ndose) else rep(0L, ndose),
      model_file = model_file,
      beta_prior = list(
        beta0_mean = beta0_mean,
        beta0_prec = beta0_prec,
        beta0_df = beta0_df,
        beta1_shape = beta1_shape,
        beta1_rate = beta1_rate,
        beta2_rate = beta2_rate
      )
    )
  )
}

## ------------------------------------------------------------
## CRM dose-move function, analogous to boin_move()
## ------------------------------------------------------------
crm_move <- function(current_dose,
                     ndose,
                     dat,
                     target = 0.30,
                     skeleton,
                     r_model = crm_r_model_choices(),
                     r_carry = 0.10,
                     a_r = 1,
                     b_r = 9,
                     alpha_sd = 2,
                     model_file = NULL,
                     fixed_model_file = "fix_CRM.bug",
                     random_model_file = "random_CRM.bug",
                     level_model_file = "random_CRM_level.bug",
                     n_chains = 2,
                     n_adapt = 500,
                     n_burnin = 500,
                     n_iter = 2000,
                     thin = 1,
                     seed = NULL,
                     elimi = rep(0L, ndose),
                     n_trt_curr = NA_real_,
                     dose_cap = 3L,
                     no_skip = TRUE,
                     cutoff = 0.95,
                     cutoff.eli = NULL,
                     dose_values = NULL,
                     dose_scores = NULL,
                     time_col = NULL,
                     alpha_grid = seq(0.01, 0.99, length.out = 61),
                     alpha_T = 28,
                     theta_prior_mean = 0,
                     theta_prior_sd = sqrt(2),
                     alpha_L = 8,
                     alpha_rel_tol = 1e-8,
                     alpha_eps = 1e-12,
                     alpha_n_draw_prior = 5000,
                     cumu_model_file = NULL,
                     cumu_beta0_mean = stats::qlogis(target),
                     cumu_beta0_prec = 4,
                     cumu_beta0_df = 1,
                     cumu_beta1_shape = 5.83,
                     cumu_beta1_rate = 1.21,
                     cumu_beta2_rate = 1,
                     cumu_include_current = FALSE) {

  r_model <- crm_normalize_r_model(r_model)
  if (!is.null(cutoff.eli)) cutoff <- cutoff.eli
  crm_validate_skeleton(skeleton, ndose)

  next_dose <- current_dose
  action <- "stay"

  can_escalate <- function() {
    current_dose < ndose &&
      elimi[current_dose + 1L] == 0L &&
      is.finite(n_trt_curr) &&
      n_trt_curr >= dose_cap
  }

  if (is.null(dat) || nrow(dat) == 0L) {
    return(list(
      next_dose = current_dose,
      action = "stay_no_data",
      mu_hat = NA_real_,
      n_eff = 0,
      method = paste0("crm_", r_model),
      p_hat = rep(NA_real_, ndose),
      theta_ipde_hat = rep(NA_real_, ndose),
      r_hat = if (r_model == "fixed") {
        r_carry
      } else if (r_model == "level") {
        out <- rep(NA_real_, ndose)
        names(out) <- paste0("D", seq_len(ndose))
        out
      } else {
        NA_real_
      },
      prob_overtox = NA_real_,
      prob_p1_over_target = NA_real_,
      stop = 0L,
      earlystop = 0L,
      stop_trial = FALSE,
      eliminated = elimi,
      crm_selected_dose = NA_integer_
    ))
  }

  fit <- crm_fit(
    dat = dat,
    ndose = ndose,
    skeleton = skeleton,
    target = target,
    cutoff = cutoff,
    r_model = r_model,
    r_carry = r_carry,
    a_r = a_r,
    b_r = b_r,
    alpha_sd = alpha_sd,
    model_file = model_file,
    fixed_model_file = fixed_model_file,
    random_model_file = random_model_file,
    level_model_file = level_model_file,
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_burnin = n_burnin,
    n_iter = n_iter,
    thin = thin,
    seed = seed,
    dose_values = dose_values,
    dose_scores = dose_scores,
    time_col = time_col,
    alpha_grid = alpha_grid,
    alpha_T = alpha_T,
    theta_prior_mean = theta_prior_mean,
    theta_prior_sd = theta_prior_sd,
    alpha_L = alpha_L,
    alpha_rel_tol = alpha_rel_tol,
    alpha_eps = alpha_eps,
    alpha_n_draw_prior = alpha_n_draw_prior,
    cumu_model_file = cumu_model_file,
    cumu_beta0_mean = cumu_beta0_mean,
    cumu_beta0_prec = cumu_beta0_prec,
    cumu_beta0_df = cumu_beta0_df,
    cumu_beta1_shape = cumu_beta1_shape,
    cumu_beta1_rate = cumu_beta1_rate,
    cumu_beta2_rate = cumu_beta2_rate,
    cumu_include_current = cumu_include_current
  )

  if (!is.null(fit$stop) && isTRUE(fit$stop == 1L)) {
    return(list(
      next_dose = 99L,
      action = "earlystop_p1_over_toxic",
      mu_hat = fit$p_hat[current_dose],
      n_eff = nrow(dat),
      method = paste0("crm_", r_model),
      p_hat = fit$p_hat,
      theta_ipde_hat = fit$theta_ipde_hat,
      r_hat = fit$r_hat,
      prob_overtox = fit$prob_overtox,
      prob_p1_over_target = fit$prob_overtox,
      stop = 1L,
      earlystop = 1L,
      stop_trial = TRUE,
      eliminated = if (!is.null(fit$eliminated)) fit$eliminated else rep(1L, ndose),
      crm_selected_dose = NA_integer_,
      r_carry = r_carry,
      r_model = r_model,
      model_file = if (!is.null(fit$model_file)) fit$model_file else model_file
    ))
  }

  dist <- abs(fit$p_hat - target)
  dist[elimi == 1L] <- Inf

  if (all(!is.finite(dist))) {
    return(list(
      next_dose = current_dose,
      action = "stay_all_eliminated",
      mu_hat = fit$p_hat[current_dose],
      n_eff = nrow(dat),
      method = paste0("crm_", r_model),
      p_hat = fit$p_hat,
      theta_ipde_hat = fit$theta_ipde_hat,
      r_hat = fit$r_hat,
      prob_overtox = fit$prob_overtox,
      prob_p1_over_target = fit$prob_overtox,
      stop = fit$stop,
      earlystop = 0L,
      stop_trial = FALSE,
      eliminated = elimi,
      crm_selected_dose = NA_integer_,
      model_file = if (!is.null(fit$model_file)) fit$model_file else model_file
    ))
  }

  crm_selected_dose <- which.min(dist)

  if (crm_selected_dose > current_dose) {
    if (can_escalate()) {
      next_dose <- if (no_skip) current_dose + 1L else crm_selected_dose
      action <- "escalate"
    } else {
      next_dose <- current_dose
      action <- "stay_escalation_blocked"
    }
  } else if (crm_selected_dose < current_dose) {
    next_dose <- if (no_skip) current_dose - 1L else crm_selected_dose
    next_dose <- max(1L, next_dose)
    action <- "de-escalate"
  } else {
    next_dose <- current_dose
    action <- "stay"
  }

  next_dose <- max(1L, min(ndose, as.integer(next_dose)))

  list(
    next_dose = next_dose,
    action = action,
    mu_hat = fit$p_hat[current_dose],
    n_eff = nrow(dat),
    method = paste0("crm_", r_model),
    p_hat = fit$p_hat,
    theta_ipde_hat = fit$theta_ipde_hat,
    r_hat = fit$r_hat,
    prob_overtox = fit$prob_overtox,
    prob_p1_over_target = fit$prob_overtox,
    stop = fit$stop,
    earlystop = 0L,
    stop_trial = FALSE,
    eliminated = elimi,
    crm_selected_dose = as.integer(crm_selected_dose),
    r_carry = r_carry,
    r_model = r_model,
    model_file = if (!is.null(fit$model_file)) fit$model_file else model_file
  )
}

## ------------------------------------------------------------
## Optional final MTD selection using CRM posterior p_j
## ------------------------------------------------------------
select.mtd.crm <- function(target,
                           dat,
                           ndose,
                           skeleton,
                           cutoff.eli = 0.95,
                           r_model = crm_r_model_choices(),
                           r_carry = 0.10,
                           a_r = 1,
                           b_r = 9,
                           alpha_sd = 2,
                           model_file = NULL,
                           fixed_model_file = "fix_CRM.bug",
                           random_model_file = "random_CRM.bug",
                           level_model_file = "random_CRM_level.bug",
                           n_chains = 2,
                           n_adapt = 500,
                           n_burnin = 500,
                           n_iter = 5000,
                           thin = 1,
                           seed = NULL,
                           restrict_to_tried = TRUE,
                           dose_values = NULL,
                           dose_scores = NULL,
                           time_col = NULL,
                           alpha_grid = seq(0.01, 0.99, length.out = 61),
                           alpha_T = 28,
                           theta_prior_mean = 0,
                           theta_prior_sd = sqrt(2),
                           alpha_L = 8,
                           alpha_rel_tol = 1e-8,
                           alpha_eps = 1e-12,
                           alpha_n_draw_prior = 5000,
                           cumu_model_file = NULL,
                           cumu_beta0_mean = stats::qlogis(target),
                           cumu_beta0_prec = 4,
                           cumu_beta0_df = 1,
                           cumu_beta1_shape = 5.83,
                           cumu_beta1_rate = 1.21,
                           cumu_beta2_rate = 1,
                           cumu_include_current = FALSE) {

  r_model <- crm_normalize_r_model(r_model)
  crm_validate_skeleton(skeleton, ndose)

  phat_out <- rep(NA_real_, ndose)
  elimi <- rep(0L, ndose)

  if (is.null(dat) || nrow(dat) == 0L) {
    return(list(
      MTD = 99L,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = paste0("crm_", r_model),
      r_hat = if (r_model == "fixed") {
        r_carry
      } else if (r_model == "level") {
        out <- rep(NA_real_, ndose)
        names(out) <- paste0("D", seq_len(ndose))
        out
      } else {
        NA_real_
      },
      prob_overtox = NA_real_,
      prob_p1_over_target = NA_real_,
      model_file = model_file,
      earlystop = 0L,
      stop = 0L
    ))
  }

  dat2 <- crm_prepare_dat(dat, ndose, time_col = time_col)
  n <- tabulate(as.integer(dat2$dose), nbins = ndose)
  y <- tabulate(as.integer(dat2$dose[dat2$y == 1L]), nbins = ndose)

  ## Keep the original empirical overdose elimination for doses with enough data.
  for (j in seq_len(ndose)) {
    if (n[j] > 2L) {
      post_over <- 1 - stats::pbeta(target, y[j] + 1, n[j] - y[j] + 1)
      if (post_over > cutoff.eli) {
        elimi[j:ndose] <- 1L
        break
      }
    }
  }

  fit <- crm_fit(
    dat = dat2,
    ndose = ndose,
    skeleton = skeleton,
    target = target,
    cutoff = cutoff.eli,
    r_model = r_model,
    r_carry = r_carry,
    a_r = a_r,
    b_r = b_r,
    alpha_sd = alpha_sd,
    model_file = model_file,
    fixed_model_file = fixed_model_file,
    random_model_file = random_model_file,
    level_model_file = level_model_file,
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_burnin = n_burnin,
    n_iter = n_iter,
    thin = thin,
    seed = seed,
    dose_values = dose_values,
    dose_scores = dose_scores,
    time_col = time_col,
    alpha_grid = alpha_grid,
    alpha_T = alpha_T,
    theta_prior_mean = theta_prior_mean,
    theta_prior_sd = theta_prior_sd,
    alpha_L = alpha_L,
    alpha_rel_tol = alpha_rel_tol,
    alpha_eps = alpha_eps,
    alpha_n_draw_prior = alpha_n_draw_prior,
    cumu_model_file = cumu_model_file,
    cumu_beta0_mean = cumu_beta0_mean,
    cumu_beta0_prec = cumu_beta0_prec,
    cumu_beta0_df = cumu_beta0_df,
    cumu_beta1_shape = cumu_beta1_shape,
    cumu_beta1_rate = cumu_beta1_rate,
    cumu_beta2_rate = cumu_beta2_rate,
    cumu_include_current = cumu_include_current
  )

  phat_out <- fit$p_hat

  if (elimi[1L] == 1L || (!is.null(fit$stop) && isTRUE(fit$stop == 1L))) {
    return(list(
      MTD = 99L,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = paste0("crm_", r_model),
      p_hat = fit$p_hat,
      theta_ipde_hat = fit$theta_ipde_hat,
      r_hat = fit$r_hat,
      prob_overtox = fit$prob_overtox,
      prob_p1_over_target = fit$prob_overtox,
      model_file = if (!is.null(fit$model_file)) fit$model_file else model_file,
      earlystop = 1L,
      stop = 1L
    ))
  }

  admissible <- elimi == 0L
  if (restrict_to_tried) admissible <- admissible & (n > 0L)

  dist <- abs(phat_out - target)
  dist[!admissible] <- Inf

  if (all(!is.finite(dist))) {
    mtd <- NA_integer_
  } else {
    mtd <- which.min(dist)
  }

  list(
    MTD = as.integer(mtd),
    phat = phat_out,
    pj_iso = phat_out,
    eliminated = elimi,
    approx = paste0("crm_", r_model),
    p_hat = fit$p_hat,
    theta_ipde_hat = fit$theta_ipde_hat,
    r_hat = fit$r_hat,
    prob_overtox = fit$prob_overtox,
    prob_p1_over_target = fit$prob_overtox,
    model_file = if (!is.null(fit$model_file)) fit$model_file else model_file,
    earlystop = if (!is.null(fit$earlystop)) fit$earlystop else fit$stop,
    stop = fit$stop,
    fit = fit
  )
}
