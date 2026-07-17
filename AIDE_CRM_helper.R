## ============================================================
## AIDE-CRM discount-model add-on
## Requires JAGS model files:
##   fixed r:  fix_CRM.bug
##   random r: random_CRM.bug
## Assumed JAGS data names:
##   N, J, y, dose, ipde, q, tau_alpha
##   fixed model additionally uses data r
##   random model additionally uses data a_r, b_r and parameter r
## ============================================================

## ------------------------------------------------------------
## Fit CRM discount model using JAGS
## ------------------------------------------------------------
crm_fit_discount <- function(dat,
                             ndose,
                             skeleton,
                             target = 0.30,
                             cutoff = 0.95,
                             r_model = c("fixed", "random"),
                             r_carry = 0.10,
                             a_r = 1,
                             b_r = 9,
                             alpha_sd = 2,
                             fixed_model_file = "fix_CRM.bug",
                             random_model_file = "random_CRM.bug",
                             n_chains = 2,
                             n_adapt = 500,
                             n_burnin = 500,
                             n_iter = 2000,
                             thin = 1,
                             seed = NULL) {

  r_model <- match.arg(r_model)

  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required for CRM. Install/load JAGS and rjags first.")
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for CRM posterior samples.")
  }

  if (length(skeleton) != ndose) {
    stop("skeleton must have length ndose.")
  }
  if (any(!is.finite(skeleton)) || any(skeleton <= 0 | skeleton >= 1)) {
    stop("skeleton must contain probabilities in (0, 1).")
  }
  if (alpha_sd <= 0 || !is.finite(alpha_sd)) {
    stop("alpha_sd must be positive.")
  }
  if (r_model == "fixed" && (length(r_carry) != 1L || r_carry < 0 || r_carry >= 1)) {
    stop("For fixed r, r_carry must be a scalar in [0, 1).")
  }

  ## dat should contain dose, y, and type, where type is "new" or "retreat".
  if (is.null(dat) || nrow(dat) == 0L) {
    stop("CRM requires at least one observed assignment in dat.")
  }
  if (!all(c("dose", "y", "type") %in% names(dat))) {
    stop("dat must contain columns: dose, y, type.")
  }

  y_vec <- as.integer(dat$y)
  dose_vec <- as.integer(dat$dose)
  ipde_vec <- as.integer(dat$type == "retreat")

  if (any(is.na(y_vec)) || any(!y_vec %in% c(0L, 1L))) {
    stop("dat$y must be 0/1 with no missing values.")
  }
  if (any(is.na(dose_vec)) || any(dose_vec < 1L | dose_vec > ndose)) {
    stop("dat$dose must be integer dose levels from 1 to ndose.")
  }

  jags_data <- list(
    N = length(y_vec),
    J = ndose,
    y = y_vec,
    dose = dose_vec,
    ipde = ipde_vec,
    q = as.numeric(skeleton),
    tau_alpha = 1 / alpha_sd^2
  )

  model_file <- fixed_model_file
  monitors <- c("alpha", "p", "theta_ipde")

  if (r_model == "fixed") {
    jags_data$r <- r_carry
  } else {
    model_file <- random_model_file
    jags_data$a_r <- a_r
    jags_data$b_r <- b_r
    monitors <- c(monitors, "r")
  }

  if (!file.exists(model_file)) {
    stop("Cannot find JAGS model file: ", model_file)
  }

  if (!is.null(seed)) set.seed(seed)

  rng_names <- c(
    "base::Wichmann-Hill",
    "base::Marsaglia-Multicarry",
    "base::Super-Duper",
    "base::Mersenne-Twister"
  )

  # inits <- vector("list", n_chains)
  # for (ch in seq_len(n_chains)) {
  #   init_ch <- list(
  #     alpha = rnorm(1, 0, 0.25),
  #     .RNG.name = rng_names[(ch - 1L) %% length(rng_names) + 1L],
  #     .RNG.seed = if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed + ch
  #   )
  #   if (r_model == "random") {
  #     init_ch$r <- rbeta(1, a_r, b_r)
  #   }
  #   inits[[ch]] <- init_ch
  # }

  jm <- rjags::jags.model(
    file = model_file,
    data = jags_data,
    # inits = inits,
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

  r_hat <- if (r_model == "fixed") r_carry else mean(post[, "r"])
  prob_overtox <- mean(post[, p_cols[1L]] > target)
  stop_flag <- as.integer(prob_overtox > cutoff)

  list(
    p_hat = as.numeric(p_hat),
    theta_ipde_hat = as.numeric(theta_hat),
    r_hat = as.numeric(r_hat),
    post = post,
    r_model = r_model,
    target = target,
    skeleton = skeleton,
    prob_overtox = prob_overtox,
    prob_p1_over_target = prob_overtox,
    stop = stop_flag,
    earlystop = stop_flag,
    eliminated = if (stop_flag == 1L) rep(1L, ndose) else rep(0L, ndose)
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
                     r_model = c("fixed", "random"),
                     r_carry = 0.10,
                     a_r = 1,
                     b_r = 9,
                     alpha_sd = 2,
                     fixed_model_file = "fix_CRM.bug",
                     random_model_file = "random_CRM.bug",
                     n_chains = 2,
                     n_adapt = 500,
                     n_burnin = 500,
                     n_iter = 2000,
                     thin = 1,
                     seed = NULL,
                     elimi = rep(0L, ndose),
                     n_trt_curr = NA_real_,
                     dose_cap = 3L,
                     no_skip = TRUE) {

  r_model <- match.arg(r_model)

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
      r_hat = if (r_model == "fixed") r_carry else NA_real_,
      crm_selected_dose = NA_integer_
    ))
  }

  fit <- crm_fit_discount(
    dat = dat,
    ndose = ndose,
    skeleton = skeleton,
    target = target,
    r_model = r_model,
    r_carry = r_carry,
    a_r = a_r,
    b_r = b_r,
    alpha_sd = alpha_sd,
    fixed_model_file = fixed_model_file,
    random_model_file = random_model_file,
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_burnin = n_burnin,
    n_iter = n_iter,
    thin = thin,
    seed = seed
  )

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
      crm_selected_dose = NA_integer_
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
    crm_selected_dose = as.integer(crm_selected_dose),
    r_carry = r_carry,
    r_model = r_model
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
                           r_model = c("fixed", "random"),
                           r_carry = 0.10,
                           a_r = 1,
                           b_r = 9,
                           alpha_sd = 2,
                           fixed_model_file = "fix_CRM.bug",
                           random_model_file = "random_CRM.bug",
                           n_chains = 2,
                           n_adapt = 500,
                           n_burnin = 500,
                           n_iter = 5000,
                           thin = 1,
                           seed = NULL,
                           restrict_to_tried = TRUE) {

  r_model <- match.arg(r_model)

  phat_out <- rep(NA_real_, ndose)
  elimi <- rep(0L, ndose)

  if (is.null(dat) || nrow(dat) == 0L) {
    return(list(
      MTD = 99L,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = paste0("crm_", r_model)
    ))
  }

  n <- tabulate(as.integer(dat$dose), nbins = ndose)
  fit <- crm_fit_discount(
    dat = dat,
    ndose = ndose,
    skeleton = skeleton,
    target = target,
    cutoff = cutoff.eli,
    r_model = r_model,
    r_carry = r_carry,
    a_r = a_r,
    b_r = b_r,
    alpha_sd = alpha_sd,
    fixed_model_file = fixed_model_file,
    random_model_file = random_model_file,
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_burnin = n_burnin,
    n_iter = n_iter,
    thin = thin,
    seed = seed
  )

  phat_out <- fit$p_hat

  ## Over-toxicity is determined by the fitted CRM posterior, rather than
  ## dose-specific beta-binomial calculations.
  if (!is.null(fit$eliminated)) {
    elimi <- as.integer(fit$eliminated)
  } else if (!is.null(fit$stop) && isTRUE(fit$stop == 1L)) {
    elimi <- rep(1L, ndose)
  }

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
    earlystop = if (!is.null(fit$earlystop)) fit$earlystop else fit$stop,
    stop = fit$stop
  )
}
