## ============================================================
## AIDE-BOIN functions
## Supports decision / final selection methods:
##   1) original BOIN boundary: method = "boin"
##   2) approximate IPDE MLE: method = "approx1"
##   3) moment / closed-form approximation: method = "approx2"
##
## For approx1/approx2, r can be handled in two ways:
##   r_estimator = "r_fixed": use the user-specified fixed r_carry
##   r_estimator = "r_mle"  : estimate r by the smoothed MLE and plug it in
##   r_estimator = "r_adaptive": estimate one global r by posterior integration
##
## Smoothed r_MLE uses posterior means under Beta(phi/2, 1 - phi/2):
##   pR_tilde = (YR + phi/2) / (NR + 1)
##   pI_tilde = (YI + phi/2) / (NI + 1)
##   r_hat = (pI_tilde - pR_tilde) / (1 - pR_tilde), truncated to [0, 1).
## ============================================================

clip01 <- function(x) {
  pmax(0, pmin(1, x))
}

xlogy <- function(x, y) {
  ifelse(x == 0, 0, x * log(y))
}

log_pow <- function(x, a) {
  if (a == 0) return(0)
  if (x <= 0) return(-Inf)
  a * log(x)
}

logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

## ------------------------------------------------------------
## Exact mixed MLE for model 2
## Regular/new patients:     Y_R ~ Bin(N_R, p)
## IPDE/recycled patients:   Y_I ~ Bin(N_I, r + (1-r)p)
##
## The old calls passed N = N_R + N_I.  For identifiability and correct
## likelihood evaluation, this function now also accepts NR and NI.
## ------------------------------------------------------------
exact_mixed_mle_model2 <- function(yR,
                                   yI,
                                   N = NULL,
                                   r_carry = 0.1,
                                   NR = NULL,
                                   NI = NULL,
                                   tol = 1e-10) {
  if (length(r_carry) != 1L || is.na(r_carry) || r_carry < 0 || r_carry >= 1) {
    stop("r_carry must be a single value in [0, 1).")
  }
  if (is.null(NR) || is.null(NI)) {
    stop("exact_mixed_mle_model2() requires NR and NI for the mixed likelihood.")
  }
  if (any(lengths(list(yR, yI, NR, NI)) != 1L)) {
    stop("yR, yI, NR, and NI must be scalar values.")
  }

  yR <- as.numeric(yR)
  yI <- as.numeric(yI)
  NR <- as.numeric(NR)
  NI <- as.numeric(NI)

  if (!is.null(N)) {
    N <- as.numeric(N)
    if (length(N) != 1L || is.na(N)) stop("N must be a scalar if provided.")
    if (abs(N - (NR + NI)) > 1e-8) {
      stop("N must equal NR + NI when provided.")
    }
  }

  vals <- c(yR, yI, NR, NI)
  if (any(!is.finite(vals)) || any(vals < 0)) {
    stop("yR, yI, NR, and NI must be finite nonnegative values.")
  }
  if (yR > NR + 1e-8 || yI > NI + 1e-8) {
    stop("Observed toxicities cannot exceed corresponding sample sizes.")
  }

  Nstar <- NR + NI
  if (Nstar <= 0) return(NA_real_)

  ## Numerically stable binomial log-likelihood without constants.
  loglik <- function(p) {
    if (!is.finite(p) || p < 0 || p > 1) return(-Inf)

    q <- r_carry + (1 - r_carry) * p

    out <- 0

    ## Regular contribution.
    if (NR > 0) {
      if (p == 0 && yR > 0) return(-Inf)
      if (p == 1 && (NR - yR) > 0) return(-Inf)

      out <- out + if (yR == 0) 0 else yR * log(p)
      out <- out + if ((NR - yR) == 0) 0 else (NR - yR) * log1p(-p)
    }

    ## IPDE contribution.
    if (NI > 0) {
      if (q == 0 && yI > 0) return(-Inf)
      if (q == 1 && (NI - yI) > 0) return(-Inf)

      out <- out + if (yI == 0) 0 else yI * log(q)
      out <- out + if ((NI - yI) == 0) 0 else (NI - yI) * log1p(-q)
    }

    out
  }

  eps <- tol
  candidates <- c(0, 1)

  opt <- tryCatch(
    optimize(
      f = function(p) -loglik(p),
      interval = c(eps, 1 - eps),
      tol = tol
    ),
    error = function(e) NULL
  )

  if (!is.null(opt) && is.finite(opt$minimum)) {
    candidates <- c(candidates, opt$minimum)
  }

  ll <- vapply(candidates, loglik, numeric(1))
  best <- candidates[which.max(ll)]

  clip01(best)
}

## ------------------------------------------------------------
## Smoothed MLE for r under model 2
##
## Raw MLE:
##   r_hat = (YI/NI - YR/NR) / (1 - YR/NR)
##
## Here, YI/NI and YR/NR are replaced by posterior means under
## Beta(phi/2, 1 - phi/2).  Since a0 + b0 = 1,
##   p_post = (y + phi/2) / (n + 1).
## ------------------------------------------------------------
beta_binom_post_rate_phi <- function(y, n, phi) {
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0 || phi >= 1) {
    stop("phi must be a scalar in (0, 1).")
  }
  if (length(y) != 1L || length(n) != 1L ||
      !is.finite(y) || !is.finite(n) || y < 0 || n < 0) {
    stop("y and n must be finite nonnegative scalars.")
  }
  if (y > n + 1e-8) {
    stop("Observed toxicities cannot exceed sample size.")
  }

  a0 <- phi / 2
  b0 <- 1 - phi / 2

  (y + a0) / (n + a0 + b0)
}

## ------------------------------------------------------------
## Pooled regular-patient toxicity cap for r
##
## For current dose j:
##   j = 1: use regular data from dose 1 only
##   j > 1: use regular data pooled from doses j-1 and j
##
## Smoothed by Beta(phi/2, 1 - phi/2):
##   p_cap = (Y_pool + phi/2) / (N_pool + 1)
## ------------------------------------------------------------
estimate_r_cap_pooled_regular <- function(current_dose,
                                          y_regular,
                                          n_regular,
                                          phi) {
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0 || phi >= 1) {
    stop("phi must be a scalar in (0, 1).")
  }

  if (is.null(y_regular) || is.null(n_regular)) {
    return(NA_real_)
  }

  y_regular <- as.numeric(y_regular)
  n_regular <- as.numeric(n_regular)

  if (length(y_regular) != length(n_regular)) {
    stop("y_regular and n_regular must have the same length.")
  }
  if (length(y_regular) < 1L) {
    return(NA_real_)
  }
  if (current_dose < 1L || current_dose > length(y_regular)) {
    stop("current_dose must be between 1 and length(y_regular).")
  }

  if (current_dose == 1L) {
    idx <- 1L
  } else {
    idx <- c(current_dose - 1L, current_dose)
  }

  y_pool <- sum(y_regular[idx], na.rm = TRUE)
  n_pool <- sum(n_regular[idx], na.rm = TRUE)

  if (!is.finite(y_pool) || !is.finite(n_pool) ||
      y_pool < 0 || n_pool < 0 || y_pool > n_pool + 1e-8) {
    stop("Invalid pooled regular-patient counts for r truncation.")
  }

  if (n_pool <= 0) {
    return(NA_real_)
  }

  beta_binom_post_rate_phi(
    y = y_pool,
    n = n_pool,
    phi = phi
  )
}


estimate_r_mle_beta_binom <- function(yR,
                                      yI,
                                      NR,
                                      NI,
                                      phi,
                                      r_cap = NA_real_,
                                      eps = 1e-12) {
  pR_post <- beta_binom_post_rate_phi(y = yR, n = NR, phi = phi)
  pI_post <- beta_binom_post_rate_phi(y = yI, n = NI, phi = phi)

  denom <- 1 - pR_post

  if (!is.finite(denom) || denom <= eps) {
    r_hat_raw <- 1 - eps
  } else {
    r_hat_raw <- (pI_post - pR_post) / denom
    r_hat_raw <- max(0, min(1 - eps, r_hat_raw))
  }

  ## Truncate raw r_MLE by pooled adjacent regular-patient toxicity.
  ## If r_cap is unavailable, keep the raw r_MLE.
  if (length(r_cap) != 1L || !is.finite(r_cap)) {
    r_use <- r_hat_raw
  } else {
    r_cap <- max(0, min(1 - eps, r_cap))
    r_use <- min(r_hat_raw, r_cap)
  }

  list(
    r_hat = as.numeric(r_hat_raw),
    r_use = as.numeric(r_use),
    r_cap = if (is.finite(r_cap)) as.numeric(r_cap) else NA_real_,
    pR_post = as.numeric(pR_post),
    pI_post = as.numeric(pI_post)
  )
}

## ------------------------------------------------------------
## Adaptive/global posterior estimator of r
##
## Model:
##   Y_Rk ~ Bin(N_Rk, p_k)
##   Y_Ik ~ Bin(N_Ik, r + (1-r)p_k)
##
## For a candidate r, p_k is integrated out analytically under
## p_k ~ Beta(p_prior[1], p_prior[2]).  The global posterior for r is
## then integrated numerically on [0, r_max], with a transformation that
## stabilizes Beta priors whose density is singular at r = 0.
## ------------------------------------------------------------
check_adaptive_counts <- function(yR, nR, yI, nI) {
  vals <- list(yR = yR, nR = nR, yI = yI, nI = nI)
  lens <- unique(vapply(vals, length, integer(1)))
  if (length(lens) != 1L) {
    stop("yR, nR, yI, and nI must have the same length.")
  }

  flat <- unlist(vals)
  if (any(!is.finite(flat)) || any(flat < 0)) {
    stop("Counts must be finite and nonnegative.")
  }
  if (any(abs(flat - round(flat)) > 1e-8)) {
    stop("Counts must be integers.")
  }
  if (any(yR > nR + 1e-8) || any(yI > nI + 1e-8)) {
    stop("Toxicity counts cannot exceed sample sizes.")
  }

  invisible(TRUE)
}

log_marginal_dose_adaptive <- function(r,
                                       yR,
                                       nR,
                                       yI,
                                       nI,
                                       p_prior = c(1, 1)) {
  if (r < 0 || r >= 1) return(-Inf)
  if (nR + nI == 0) return(0)

  ap <- p_prior[1]
  bp <- p_prior[2]
  if (ap <= 0 || bp <= 0) stop("p_prior must be positive.")

  yI <- as.integer(round(yI))
  failures <- nR - yR + nI - yI
  j <- 0:yI

  terms <- vapply(
    j,
    function(jj) {
      lchoose(yI, jj) +
        log_pow(r, yI - jj) +
        log_pow(1 - r, nI - yI + jj) +
        lbeta(yR + jj + ap, failures + bp) -
        lbeta(ap, bp)
    },
    numeric(1)
  )

  logsumexp(terms)
}

log_posterior_r_adaptive <- function(r,
                                     yR,
                                     nR,
                                     yI,
                                     nI,
                                     r_prior = c(1, 9),
                                     p_prior = c(1, 1),
                                     r_max = 0.6) {
  if (r < 0 || r > r_max || r >= 1) return(-Inf)

  ar <- r_prior[1]
  br <- r_prior[2]
  if (ar <= 0 || br <= 0) stop("r_prior must be positive.")
  if (r_max <= 0 || r_max >= 1) stop("r_max must be in (0, 1).")

  log_prior <- stats::dbeta(r, ar, br, log = TRUE) -
    log(stats::pbeta(r_max, ar, br))

  log_m <- vapply(
    seq_along(yR),
    function(k) {
      log_marginal_dose_adaptive(
        r = r,
        yR = yR[k],
        nR = nR[k],
        yI = yI[k],
        nI = nI[k],
        p_prior = p_prior
      )
    },
    numeric(1)
  )

  log_prior + sum(log_m)
}

estimate_global_r_adaptive <- function(yR,
                                       nR,
                                       yI,
                                       nI,
                                       r_prior = c(1, 9),
                                       p_prior = c(1, 1),
                                       r_max = 0.6,
                                       plug_in = c("mean", "map"),
                                       rel.tol = 1e-6,
                                       eps = 1e-8) {
  plug_in <- match.arg(plug_in)

  yR <- as.numeric(yR)
  nR <- as.numeric(nR)
  yI <- as.numeric(yI)
  nI <- as.numeric(nI)
  check_adaptive_counts(yR, nR, yI, nI)

  ar <- r_prior[1]
  br <- r_prior[2]
  if (ar <= 0 || br <= 0) stop("r_prior must be positive.")
  if (r_max <= 0 || r_max >= 1) stop("r_max must be in (0, 1).")

  logpost_scalar <- function(r) {
    log_posterior_r_adaptive(
      r = r,
      yR = yR,
      nR = nR,
      yI = yI,
      nI = nI,
      r_prior = r_prior,
      p_prior = p_prior,
      r_max = r_max
    )
  }

  logpost_vec <- function(r) vapply(r, logpost_scalar, numeric(1))

  opt_r <- stats::optimize(
    f = logpost_scalar,
    interval = c(eps, r_max - eps),
    maximum = TRUE
  )

  candidate_r <- c(eps, opt_r$maximum, r_max - eps)
  candidate_lp <- logpost_vec(candidate_r)
  r_map <- if (ar < 1) 0 else candidate_r[which.max(candidate_lp)]

  if (ar < 1) {
    r_of_u <- function(u) r_max * u^(1 / ar)
    log_jacobian <- function(u) log(r_max / ar) + (1 / ar - 1) * log(u)
  } else {
    r_of_u <- function(u) r_max * u
    log_jacobian <- function(u) rep(log(r_max), length(u))
  }

  u_lower <- eps
  u_upper <- 1 - eps

  log_denom_u_scalar <- function(u) {
    if (u <= 0 || u >= 1) return(-Inf)
    r <- r_of_u(u)
    logpost_scalar(r) + log_jacobian(u)
  }

  log_denom_u_vec <- function(u) vapply(u, log_denom_u_scalar, numeric(1))

  opt_u <- stats::optimize(
    f = log_denom_u_scalar,
    interval = c(u_lower, u_upper),
    maximum = TRUE
  )

  candidate_u <- c(u_lower, opt_u$maximum, u_upper)
  candidate_lu <- log_denom_u_vec(candidate_u)
  finite_lu <- is.finite(candidate_lu)

  if (!any(finite_lu)) {
    r_hat <- if (plug_in == "map") r_map else opt_r$maximum
    r_hat <- max(0, min(r_max, r_hat))
    return(list(
      r_hat = r_hat,
      r_mean = NA_real_,
      r_map = r_map,
      r_sd = NA_real_,
      r_prior = r_prior,
      p_prior = p_prior,
      r_max = r_max,
      plug_in = plug_in
    ))
  }

  lu_scale <- max(candidate_lu[finite_lu])

  denom_integrand_u <- function(u) exp(log_denom_u_vec(u) - lu_scale)
  numer_integrand_u <- function(u) {
    r <- r_of_u(u)
    r * exp(log_denom_u_vec(u) - lu_scale)
  }

  denom <- tryCatch(
    stats::integrate(
      denom_integrand_u,
      lower = u_lower,
      upper = u_upper,
      rel.tol = rel.tol,
      subdivisions = 250L,
      stop.on.error = FALSE
    )$value,
    error = function(e) NA_real_
  )

  numer <- tryCatch(
    stats::integrate(
      numer_integrand_u,
      lower = u_lower,
      upper = u_upper,
      rel.tol = rel.tol,
      subdivisions = 250L,
      stop.on.error = FALSE
    )$value,
    error = function(e) NA_real_
  )

  if (!is.finite(denom) || denom <= 0 || !is.finite(numer)) {
    r_mean <- NA_real_
    r_hat <- r_map
    r_sd <- NA_real_
  } else {
    r_mean <- numer / denom
    r_hat <- if (plug_in == "mean") r_mean else r_map

    var_integrand_u <- function(u) {
      r <- r_of_u(u)
      (r - r_mean)^2 * exp(log_denom_u_vec(u) - lu_scale)
    }

    var_num <- tryCatch(
      stats::integrate(
        var_integrand_u,
        lower = u_lower,
        upper = u_upper,
        rel.tol = rel.tol,
        subdivisions = 250L,
        stop.on.error = FALSE
      )$value,
      error = function(e) NA_real_
    )

    r_sd <- if (is.finite(var_num)) sqrt(var_num / denom) else NA_real_
  }

  r_hat <- max(0, min(r_max, r_hat))

  list(
    r_hat = as.numeric(r_hat),
    r_mean = as.numeric(r_mean),
    r_map = as.numeric(r_map),
    r_sd = as.numeric(r_sd),
    r_prior = r_prior,
    p_prior = p_prior,
    r_max = r_max,
    plug_in = plug_in
  )
}

## ------------------------------------------------------------
## DLT-time generator
## ------------------------------------------------------------
gen.tite <- function(dist = 2, n, pi, alpha = 0.5, Tobs) {

  weib <- function(n, pi, pihalft) {
    shape <- log(log(1 - pi) / log(1 - pihalft)) / log(2)
    lambda <- -log(1 - pi) / (Tobs ^ shape)
    (-log(runif(n)) / lambda) ^ (1 / shape)
  }

  llogit <- function(n, pi, pihalft) {
    shape <- log((1 / (1 - pi) - 1) / (1 / (1 - pihalft) - 1)) / log(2)
    lambda <- (1 / (1 - pi) - 1) / (Tobs ^ shape)
    ((1 / runif(n) - 1) / lambda) ^ (1 / shape)
  }

  if (length(pi) != 1L || is.na(pi) || pi < 0 || pi > 1) {
    stop("pi must be a single value in [0,1].")
  }
  if (length(Tobs) != 1L || is.na(Tobs) || Tobs < 0) {
    stop("Tobs must be a single nonnegative value.")
  }
  if (length(alpha) != 1L || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be a single value in [0,1].")
  }
  if (!dist %in% c(1, 2, 3, 4)) {
    stop("dist must be 1, 2, 3, or 4.")
  }

  tox <- rep(0L, n)
  t.tox <- rep(Tobs, n)

  if (Tobs == 0) {
    tox <- rbinom(n, 1L, pi)
    t.tox <- rep(0, n)
    return(list(tox = tox, t.tox = t.tox, ntox.st = sum(tox)))
  }

  if (pi == 0) {
    return(list(tox = tox, t.tox = t.tox, ntox.st = 0L))
  }

  ## dist = 1: continuous uniform DLT time conditional on DLT
  if (dist == 1) {
    tox <- rbinom(n, 1L, pi)
    ntox.st <- sum(tox)

    if (ntox.st > 0L) {
      t.tox[tox == 1L] <- runif(ntox.st, min = 0, max = Tobs)
    }

    return(list(tox = tox, t.tox = t.tox, ntox.st = ntox.st))
  }

  ## dist = 2: Weibull time to DLT
  if (dist == 2) {
    if (alpha <= 0 || alpha >= 1) {
      stop("For dist = 2, alpha should be strictly between 0 and 1.")
    }

    eps <- 1e-12
    pi_eff <- min(max(pi, eps), 1 - eps)
    pihalft <- alpha * pi_eff

    t.raw <- weib(n, pi_eff, pihalft)

    if (pi >= 1) {
      tox <- rep(1L, n)
      t.tox <- pmin(t.raw, Tobs)
    } else {
      tox[t.raw <= Tobs] <- 1L
      t.tox[tox == 1L] <- t.raw[tox == 1L]
    }

    return(list(tox = tox, t.tox = t.tox, ntox.st = sum(tox)))
  }

  ## dist = 3: log-logistic time to DLT
  if (dist == 3) {
    if (alpha <= 0 || alpha >= 1) {
      stop("For dist = 3, alpha should be strictly between 0 and 1.")
    }

    eps <- 1e-12
    pi_eff <- min(max(pi, eps), 1 - eps)
    pihalft <- alpha * pi_eff

    t.raw <- llogit(n, pi_eff, pihalft)

    if (pi >= 1) {
      tox <- rep(1L, n)
      t.tox <- pmin(t.raw, Tobs)
    } else {
      tox[t.raw <= Tobs] <- 1L
      t.tox[tox == 1L] <- t.raw[tox == 1L]
    }

    return(list(tox = tox, t.tox = t.tox, ntox.st = sum(tox)))
  }

  ## dist = 4: discrete uniform DLT time
  if (dist == 4) {
    tox <- rbinom(n, 1L, pi)
    ntox.st <- sum(tox)

    if (ntox.st > 0L) {
      if (abs(Tobs - round(Tobs)) > .Machine$double.eps ^ 0.5) {
        stop("dist = 4 requires integer Tobs because it uses sample.int(Tobs).")
      }

      t.tox[tox == 1L] <- sample.int(
        n = as.integer(round(Tobs)),
        size = ntox.st,
        replace = TRUE
      )
    }

    return(list(tox = tox, t.tox = t.tox, ntox.st = ntox.st))
  }
}

## ------------------------------------------------------------
## BOIN continuous boundaries for approx1 / approx2
## ------------------------------------------------------------
boin_boundary <- function(phi,
                          phi1 = 0.6 * phi,
                          phi2 = 1.4 * phi) {
  if (phi <= 0 || phi >= 1) stop("phi must be in (0,1)")
  if (phi1 <= 0 || phi1 >= phi) stop("phi1 must satisfy 0 < phi1 < phi")
  if (phi2 <= phi || phi2 >= 1) stop("phi2 must satisfy phi < phi2 < 1")

  lambda_e <- log((1 - phi1) / (1 - phi)) /
    log(phi * (1 - phi1) / (phi1 * (1 - phi)))

  lambda_d <- log((1 - phi) / (1 - phi2)) /
    log(phi2 * (1 - phi) / (phi * (1 - phi2)))

  list(lambda_e = lambda_e, lambda_d = lambda_d,
       phi = phi, phi1 = phi1, phi2 = phi2)
}

## ------------------------------------------------------------
## Original BOIN tabular boundaries
## ------------------------------------------------------------
get.boundary <- function(target, ncohort, cohortsize = 3,
                         design = 3, cutoff.eli = 0.95) {
  density1 <- function(p, n, m1, m2) {
    pbinom(m1, n, p) + 1 - pbinom(m2 - 1, n, p)
  }
  density2 <- function(p, n, m1) {
    1 - pbinom(m1, n, p)
  }
  density3 <- function(p, n, m2) {
    pbinom(m2 - 1, n, p)
  }
  df <- function(p, gamma, n, target) {
    l <- (n / (1 - p)) /
      (log(p / (1 - p)) - log(target / (1 - target)))
    l <- l - (log(gamma) - n * log((1 - p) / (1 - target))) /
      p / (1 - p) /
      (log(p / (1 - p)) - log(target / (1 - target)))^2
  }

  if (!design %in% 1:5) {
    stop("design must be 1, 2, 3, 4, or 5.")
  }

  if (design == 1) {
    if (target < 0.3) {
      lambda1 <- target - 0.09
      lambda2 <- target + 0.09
    } else if (target < 0.4) {
      lambda1 <- target - 0.1
      lambda2 <- target + 0.1
    } else if (target < 0.45) {
      lambda1 <- target - 0.12
      lambda2 <- target + 0.12
    } else {
      lambda1 <- target - 0.13
      lambda2 <- target + 0.13
    }
  }

  if (design == 2 || design == 4) {
    p.saf <- target - 0.05
    p.tox <- target + 0.05
  }

  if (design == 3) {
    p.saf <- target * 0.6
    p.tox <- target * 1.4

    lambda1 <- log((1 - p.saf) / (1 - target)) /
      log(target * (1 - p.saf) / (p.saf * (1 - target)))

    lambda2 <- log((1 - target) / (1 - p.tox)) /
      log(p.tox * (1 - target) / (target * (1 - p.tox)))
  }

  if (design == 5) {
    c0 <- log(1.1) / 3
  }

  ntrt <- NULL
  b.e <- NULL
  b.d <- NULL
  elim <- NULL

  for (n in seq_len(ncohort) * cohortsize) {
    ntrt <- c(ntrt, n)

    if (design == 1 || design == 3) {
      cutoff1 <- floor(n * lambda1)
      cutoff2 <- floor(n * lambda2) + 1
    }

    if (design == 5) {
      gamma <- exp(c0 * sqrt(n))

      p.tox <- uniroot(
        df,
        c(target + 0.0001, target * 2),
        gamma = gamma,
        n = n,
        target = target
      )$root

      p.saf <- uniroot(
        df,
        c(0.0001, target - 0.0001),
        gamma = gamma,
        n = n,
        target = target
      )$root

      lambda1 <- (log((1 - p.saf) / (1 - target)) - log(gamma) / n) /
        log(target * (1 - p.saf) / (p.saf * (1 - target)))

      lambda2 <- (log((1 - target) / (1 - p.tox)) + log(gamma) / n) /
        log(p.tox * (1 - target) / (target * (1 - p.tox)))

      cutoff1 <- floor(n * lambda1)
      cutoff2 <- floor(n * lambda2) + 1
    }

    if (design == 2 || design == 4) {
      error.min <- 3

      for (m1 in 0:floor(target * n)) {
        for (m2 in ceiling(target * n):n) {
          if (design == 2) {
            error1 <- integrate(density1, lower = p.saf, upper = p.tox,
                                n, m1, m2)$value / (p.tox - p.saf)
            error2 <- integrate(density2, lower = 0, upper = p.saf,
                                n, m1)$value / p.saf
            error3 <- integrate(density3, lower = p.tox, upper = 1,
                                n, m2)$value / (1 - p.tox)
          }

          if (design == 4) {
            epsilon <- p.tox - p.saf
            error1 <- integrate(density1, lower = p.saf, upper = p.tox,
                                n, m1, m2)$value / epsilon
            error2 <- integrate(density2, lower = p.saf - epsilon,
                                upper = p.saf, n, m1)$value / epsilon
            error3 <- integrate(density3, lower = p.tox,
                                upper = p.tox + epsilon, n, m2)$value / epsilon
          }

          error <- error1 + error2 + error3

          if (error < error.min) {
            error.min <- error
            cutoff1 <- m1
            cutoff2 <- m2
          }
        }
      }
    }

    b.e <- c(b.e, cutoff1)
    b.d <- c(b.d, cutoff2)

    elimineed <- 0

    if (n < 3) {
      elim <- c(elim, NA)
    } else {
      for (ntox in 3:n) {
        if (1 - pbeta(target, ntox + 1, n - ntox + 1) > cutoff.eli) {
          elimineed <- 1
          break
        }
      }

      if (elimineed == 1) {
        elim <- c(elim, ntox)
      } else {
        elim <- c(elim, NA)
      }
    }
  }

  for (i in seq_along(b.d)) {
    if (!is.na(elim[i]) && b.d[i] > elim[i]) {
      b.d[i] <- elim[i]
    }
  }

  boundaries <- rbind(ntrt, elim, b.d, b.e)
  rownames(boundaries) <- c(
    "Number of patients treated",
    "Eliminate if # of DLT >=",
    "Deescalate if # of DLT >=",
    "Escalate if # of DLT <="
  )
  colnames(boundaries) <- rep("", ncohort)

  boundaries
}

## ------------------------------------------------------------
## Unified dose-move function
## ------------------------------------------------------------
boin_move <- function(current_dose, ndose,
                      method = c("boin", "approx1", "approx2"),
                      r_estimator = c("r_fixed", "r_mle", "r_adaptive"),
                      y_curr = NULL,
                      n_curr = NULL,
                      b.e = NULL,
                      b.d = NULL,
                      C = 3L,
                      YR = 0,
                      YI = 0,
                      NR_star = 0,
                      NI_star = 0,
                      lambda_e = NULL,
                      lambda_d = NULL,
                      phi = NULL,
                      r_carry = 0.1,
                      y_regular_all = NULL,
                      n_regular_all = NULL,
                      y_recycle_all = NULL,
                      n_recycle_all = NULL,
                      r_adaptive_prior = c(1, 9),
                      p_adaptive_prior = c(1, 1),
                      r_adaptive_max = 0.6,
                      r_adaptive_plug_in = c("mean", "map"),
                      r_adaptive_rel_tol = 1e-6,
                      elimi = rep(0L, ndose),
                      n_trt_curr = n_curr,
                      dose_cap = 3L) {

  method <- match.arg(method)
  r_estimator <- match.arg(r_estimator)
  r_adaptive_plug_in <- match.arg(r_adaptive_plug_in)

  next_dose <- current_dose
  action <- "stay"
  mu_hat <- NA_real_
  n_eff <- NA_real_
  nc <- NA_integer_
  r_hat <- NA_real_
  r_use <- NA_real_
  r_cap <- NA_real_

  can_escalate <- function() {
    current_dose < ndose &&
      elimi[current_dose + 1L] == 0L &&
      is.finite(n_trt_curr) &&
      n_trt_curr >= dose_cap
  }

  if (method == "boin") {

    if (is.null(y_curr) || is.null(n_curr) || is.null(b.e) || is.null(b.d)) {
      stop("For method = 'boin', provide y_curr, n_curr, b.e, and b.d.")
    }

    y_curr <- as.numeric(y_curr)
    n_curr <- as.numeric(n_curr)
    n_eff <- n_curr

    if (n_curr > 0) {
      mu_hat <- y_curr / n_curr
    }

    nc <- ceiling(n_curr / C)
    nc <- max(1L, min(nc, length(b.e)))

    if (n_curr <= 0) {
      next_dose <- current_dose
      action <- "stay"
    } else if (y_curr <= b.e[nc]) {
      if (can_escalate()) {
        next_dose <- current_dose + 1L
        action <- "escalate"
      } else {
        next_dose <- current_dose
        action <- "stay"
      }
    } else if (y_curr >= b.d[nc]) {
      if (current_dose > 1L) {
        next_dose <- current_dose - 1L
        action <- "de-escalate"
      } else {
        next_dose <- current_dose
        action <- "stay"
      }
    } else {
      next_dose <- current_dose
      action <- "stay"
    }

  } else {

    if (is.null(lambda_e) || is.null(lambda_d)) {
      stop("For approx1/approx2, provide lambda_e and lambda_d.")
    }

    YR <- as.numeric(YR)
    YI <- as.numeric(YI)
    NR_star <- as.numeric(NR_star)
    NI_star <- as.numeric(NI_star)

    N_star <- NR_star + NI_star
    n_eff <- N_star

    if (N_star <= 0) {

      mu_hat <- NA_real_
      r_use <- NA_real_
      r_hat <- NA_real_
      r_cap <- NA_real_

    } else {

      ##########################################################
      # For approx1/approx2, use either:
      #   r_fixed: fixed r_carry
      #   r_mle  : smoothed MLE, truncated by pooled adjacent
      #            regular-patient toxicity.
      #
      # Important:
      # When only IPDE patients are observed at the current dose
      # and no regular patients are observed, we still use r_mle.
      # We do NOT fall back to the empirical IPDE rate.
      ##########################################################

      if (r_estimator == "r_fixed") {

        if (r_carry < 0 || r_carry >= 1) {
          stop("r_carry must be in [0, 1).")
        }

        r_use <- r_carry
        r_hat <- NA_real_
        r_cap <- NA_real_

      } else if (r_estimator == "r_mle") {

        if (is.null(phi) || length(phi) != 1L ||
            !is.finite(phi) || phi <= 0 || phi >= 1) {
          stop("For r_estimator = 'r_mle', provide phi in (0, 1), typically phi = target.")
        }

        r_cap <- estimate_r_cap_pooled_regular(
          current_dose = current_dose,
          y_regular = y_regular_all,
          n_regular = n_regular_all,
          phi = phi
        )

        r_out <- estimate_r_mle_beta_binom(
          yR = YR,
          yI = YI,
          NR = NR_star,
          NI = NI_star,
          phi = phi,
          r_cap = r_cap
        )

        r_hat <- r_out$r_hat
        r_use <- r_out$r_use
        r_cap <- r_out$r_cap

      } else if (r_estimator == "r_adaptive") {

        if (is.null(y_regular_all) || is.null(n_regular_all) ||
            is.null(y_recycle_all) || is.null(n_recycle_all)) {
          stop("For r_estimator = 'r_adaptive', provide regular and IPDE count vectors.")
        }

        r_out <- estimate_global_r_adaptive(
          yR = y_regular_all,
          nR = n_regular_all,
          yI = y_recycle_all,
          nI = n_recycle_all,
          r_prior = r_adaptive_prior,
          p_prior = p_adaptive_prior,
          r_max = r_adaptive_max,
          plug_in = r_adaptive_plug_in,
          rel.tol = r_adaptive_rel_tol
        )

        r_hat <- r_out$r_hat
        r_use <- r_out$r_hat
        r_cap <- NA_real_
      }

      if (method == "approx1") {

        mu_hat <- exact_mixed_mle_model2(
          yR = YR,
          yI = YI,
          N = N_star,
          NR = NR_star,
          NI = NI_star,
          r_carry = r_use
        )

      } else if (method == "approx2") {

        mu_hat <- clip01(
          (YR + (YI - r_use * NI_star) / (1 - r_use)) / N_star
        )
      }
    }

    if (!is.finite(mu_hat)) {
      next_dose <- current_dose
      action <- "stay"
    } else if (mu_hat <= lambda_e) {
      if (can_escalate()) {
        next_dose <- current_dose + 1L
        action <- "escalate"
      } else {
        next_dose <- current_dose
        action <- "stay"
      }
    } else if (mu_hat >= lambda_d) {
      if (current_dose > 1L) {
        next_dose <- current_dose - 1L
        action <- "de-escalate"
      } else {
        next_dose <- current_dose
        action <- "stay"
      }
    } else {
      next_dose <- current_dose
      action <- "stay"
    }
  }

  next_dose <- max(1L, min(ndose, as.integer(next_dose)))

  list(
    next_dose = next_dose,
    action = action,
    mu_hat = mu_hat,
    n_eff = n_eff,
    nc = nc,
    method = method,
    YR = YR,
    YI = YI,
    NR_star = NR_star,
    NI_star = NI_star,
    STFT = NA_real_,
    pi_E = NA_real_,
    pi_D = NA_real_,
    r_carry = r_carry,
    r_estimator = r_estimator,
    r_hat = r_hat,
    r_use = r_use,
    r_cap = r_cap,
    phi = if (r_estimator == "r_mle") phi else NA_real_,
    r_adaptive_prior = if (r_estimator == "r_adaptive") r_adaptive_prior else c(NA_real_, NA_real_),
    p_adaptive_prior = if (r_estimator == "r_adaptive") p_adaptive_prior else c(NA_real_, NA_real_),
    r_adaptive_max = if (r_estimator == "r_adaptive") r_adaptive_max else NA_real_,
    r_adaptive_plug_in = if (r_estimator == "r_adaptive") r_adaptive_plug_in else NA_character_
  )
}

## ------------------------------------------------------------
## Final MTD selection supporting boin / approx1 / approx2
## ------------------------------------------------------------
select.mtd <- function(target,
                       y, n,
                       cutoff.eli = 0.95,
                       approx = c("boin", "approx1", "approx2"),
                       r_carry = 0.1,
                       r_estimator = c("r_fixed", "r_mle", "r_adaptive"),
                       phi = target,
                       y_new = NULL,
                       n_new = NULL,
                       y_recycle = NULL,
                       n_recycle = NULL,
                       r_adaptive_prior = c(1, 9),
                       p_adaptive_prior = c(1, 1),
                       r_adaptive_max = 0.6,
                       r_adaptive_plug_in = c("mean", "map"),
                       r_adaptive_rel_tol = 1e-6,
                       restrict_to_tried = TRUE) {

  approx <- match.arg(approx)
  r_estimator <- match.arg(r_estimator)
  r_adaptive_plug_in <- match.arg(r_adaptive_plug_in)
  ndose <- length(n)

  pava <- function(x, wt = rep(1, length(x))) {
    n0 <- length(x)
    if (n0 <= 1) return(x)

    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed.")
    }

    lvlsets <- seq_len(n0)

    repeat {
      viol <- as.vector(diff(x)) < 0
      if (!any(viol)) break

      i <- min(which(viol))
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1L]
      idx <- lvlsets == lvl1 | lvlsets == lvl2

      x[idx] <- sum(x[idx] * wt[idx]) / sum(wt[idx])
      lvlsets[idx] <- lvl1
    }

    x
  }

  y <- as.numeric(y)
  n <- as.numeric(n)

  if (length(y) != ndose) {
    stop("y and n must have the same length.")
  }

  if (is.null(y_new)) y_new <- y
  if (is.null(n_new)) n_new <- n
  if (is.null(y_recycle)) y_recycle <- rep(0, ndose)
  if (is.null(n_recycle)) n_recycle <- rep(0, ndose)

  y_new <- as.numeric(y_new)
  n_new <- as.numeric(n_new)
  y_recycle <- as.numeric(y_recycle)
  n_recycle <- as.numeric(n_recycle)

  if (!all(lengths(list(y_new, n_new, y_recycle, n_recycle)) == ndose)) {
    stop("y_new, n_new, y_recycle, and n_recycle must all have length length(n).")
  }

  ## BOIN-style elimination still uses raw total administrations.
  elimi <- rep(0L, ndose)

  for (j in seq_len(ndose)) {
    if (n[j] > 2L) {
      post_over <- 1 - pbeta(
        target,
        y[j] + 1,
        n[j] - y[j] + 1
      )

      if (post_over > cutoff.eli) {
        elimi[j:ndose] <- 1L
        break
      }
    }
  }

  phat_out <- rep(NA_real_, ndose)
  r_hat <- rep(NA_real_, ndose)
  r_use <- rep(NA_real_, ndose)
  r_cap <- rep(NA_real_, ndose)

  if (elimi[1L] == 1L) {
    return(list(
      MTD = 99L,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = approx,
      r_estimator = r_estimator,
      r_hat = r_hat,
      r_use = r_use,
      r_cap = r_cap
    ))
  }

  restrict_to_tried <- isTRUE(restrict_to_tried)

  if (!any(n > 0L)) {
    return(list(
      MTD = 99L,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = approx,
      r_estimator = r_estimator,
      r_hat = r_hat,
      r_use = r_use,
      r_cap = r_cap
    ))
  }

  mu_hat <- rep(NA_real_, ndose)
  n_eff <- rep(0, ndose)
  adaptive_r_out <- NULL

  if (approx != "boin" && r_estimator == "r_adaptive" &&
      any((n_new + n_recycle) > 0L)) {
    adaptive_r_out <- estimate_global_r_adaptive(
      yR = y_new,
      nR = n_new,
      yI = y_recycle,
      nI = n_recycle,
      r_prior = r_adaptive_prior,
      p_prior = p_adaptive_prior,
      r_max = r_adaptive_max,
      plug_in = r_adaptive_plug_in,
      rel.tol = r_adaptive_rel_tol
    )
  }

  for (j in seq_len(ndose)) {

    if (approx == "boin") {

      if (n[j] > 0) {
        mu_hat[j] <- y[j] / n[j]
        n_eff[j] <- n[j]
      } else if (!restrict_to_tried) {
        mu_hat[j] <- 0.5
        n_eff[j] <- 0
      }

    } else {

      YR <- y_new[j]
      YI <- y_recycle[j]
      NR <- n_new[j]
      NI <- n_recycle[j]

      N_star <- NR + NI
      n_eff[j] <- N_star

      if (N_star > 0) {

        if (r_estimator == "r_fixed") {

          if (r_carry < 0 || r_carry >= 1) {
            stop("r_carry must be in [0, 1).")
          }

          r_use[j] <- r_carry
          r_hat[j] <- NA_real_
          r_cap[j] <- NA_real_

        } else if (r_estimator == "r_mle") {

          if (length(phi) != 1L || !is.finite(phi) || phi <= 0 || phi >= 1) {
            stop("For r_estimator = 'r_mle', provide phi in (0, 1), typically phi = target.")
          }

          r_cap[j] <- estimate_r_cap_pooled_regular(
            current_dose = j,
            y_regular = y_new,
            n_regular = n_new,
            phi = phi
          )

          r_out <- estimate_r_mle_beta_binom(
            yR = YR,
            yI = YI,
            NR = NR,
            NI = NI,
            phi = phi,
            r_cap = r_cap[j]
          )

          r_hat[j] <- r_out$r_hat
          r_use[j] <- r_out$r_use
          r_cap[j] <- r_out$r_cap

        } else if (r_estimator == "r_adaptive") {

          if (is.null(adaptive_r_out)) {
            adaptive_r_out <- estimate_global_r_adaptive(
              yR = y_new,
              nR = n_new,
              yI = y_recycle,
              nI = n_recycle,
              r_prior = r_adaptive_prior,
              p_prior = p_adaptive_prior,
              r_max = r_adaptive_max,
              plug_in = r_adaptive_plug_in,
              rel.tol = r_adaptive_rel_tol
            )
          }

          r_hat[j] <- adaptive_r_out$r_hat
          r_use[j] <- adaptive_r_out$r_hat
          r_cap[j] <- NA_real_
        }

        if (approx == "approx1") {

          mu_hat[j] <- exact_mixed_mle_model2(
            yR = YR,
            yI = YI,
            N = N_star,
            NR = NR,
            NI = NI,
            r_carry = r_use[j]
          )

        } else if (approx == "approx2") {

          mu_hat[j] <- clip01(
            (YR + (YI - r_use[j] * NI) / (1 - r_use[j])) / N_star
          )
        }
      } else if (!restrict_to_tried) {
        mu_hat[j] <- 0.5
        n_eff[j] <- 0
      }
    }
  }

  highest_uneliminated <- max(which(elimi == 0L))
  highest_tried <- max(which(n != 0L))

  nadmis <- if (restrict_to_tried) {
    min(highest_uneliminated, highest_tried)
  } else {
    highest_uneliminated
  }

  mu_use <- mu_hat[seq_len(nadmis)]
  neff_use <- n_eff[seq_len(nadmis)]

  valid <- is.finite(mu_use) & is.finite(neff_use)
  if (restrict_to_tried) valid <- valid & neff_use > 0

  if (!any(valid)) {
    return(list(
      MTD = NA_integer_,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = approx,
      r_estimator = r_estimator,
      r_hat = r_hat,
      r_use = r_use,
      r_cap = r_cap
    ))
  }

  pseudo_y <- mu_use * neff_use
  pseudo_n <- neff_use
  pseudo_y <- pmax(0, pmin(pseudo_y, pseudo_n))

  a_post <- pseudo_y + 0.005
  b_post <- pseudo_n - pseudo_y + 0.005

  phat <- a_post / (a_post + b_post)

  phat.var <- a_post * b_post /
    ((a_post + b_post)^2 * (a_post + b_post + 1))

  phat_for_iso <- phat
  phat_for_iso[!valid] <- target
  phat.var[!valid] <- Inf

  wt <- 1 / phat.var
  wt[!is.finite(wt)] <- 1e-8

  phat.iso <- pava(phat_for_iso, wt = wt)
  phat.iso <- phat.iso + seq_len(nadmis) * 1e-10

  dist_to_target <- abs(phat.iso - target)
  dist_to_target[!valid] <- Inf

  selectdose <- which.min(dist_to_target)

  phat_out[seq_len(nadmis)] <- phat.iso

  list(
    MTD = as.integer(selectdose),
    phat = phat_out,
    pj_iso = phat_out,
    eliminated = elimi,
    approx = approx,
    mu_hat = mu_hat,
    n_eff = n_eff,
    r_estimator = r_estimator,
    r_hat = r_hat,
    r_use = r_use,
    r_cap = r_cap
  )
}

## ------------------------------------------------------------
## CFO utilities
## ------------------------------------------------------------
if (file.exists("CFO_tox_utils.R")) {
  source("CFO_tox_utils.R")
}
