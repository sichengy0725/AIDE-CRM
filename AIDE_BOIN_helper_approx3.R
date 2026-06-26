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

estimate_r_mle_beta_binom <- function(yR,
                                      yI,
                                      NR,
                                      NI,
                                      phi,
                                      eps = 1e-12) {
  pR_post <- beta_binom_post_rate_phi(y = yR, n = NR, phi = phi)
  pI_post <- beta_binom_post_rate_phi(y = yI, n = NI, phi = phi)
  
  denom <- 1 - pR_post
  
  if (!is.finite(denom) || denom <= eps) {
    r_hat <- 1 - eps
  } else {
    r_hat <- (pI_post - pR_post) / denom
    r_hat <- max(0, min(1 - eps, r_hat))
  }
  
  list(
    r_hat = as.numeric(r_hat),
    pR_post = as.numeric(pR_post),
    pI_post = as.numeric(pI_post)
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
                      r_estimator = c("r_fixed", "r_mle"),
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
                      elimi = rep(0L, ndose),
                      n_trt_curr = n_curr,
                      dose_cap = 3L) {
  
  method <- match.arg(method)
  r_estimator <- match.arg(r_estimator)
  
  next_dose <- current_dose
  action <- "stay"
  mu_hat <- NA_real_
  n_eff <- NA_real_
  nc <- NA_integer_
  r_hat <- NA_real_
  r_use <- NA_real_
  
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
    
    N_star <- NR_star + NI_star
    n_eff <- N_star
    
    if (N_star <= 0) {
      mu_hat <- NA_real_
    } else {
      if (r_estimator == "r_fixed") {
        if (r_carry < 0 || r_carry >= 1) {
          stop("r_carry must be in [0, 1).")
        }
        r_use <- r_carry
      } else if (r_estimator == "r_mle") {
        if (is.null(phi) || length(phi) != 1L ||
            !is.finite(phi) || phi <= 0 || phi >= 1) {
          stop("For r_estimator = 'r_mle', provide phi in (0, 1), typically phi = target.")
        }
        r_out <- estimate_r_mle_beta_binom(
          yR = YR,
          yI = YI,
          NR = NR_star,
          NI = NI_star,
          phi = phi
        )
        r_hat <- r_out$r_hat
        r_use <- r_hat
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
    r_use = r_use,
    r_hat = r_hat,
    phi = if (r_estimator == "r_mle") phi else NA_real_
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
                       r_estimator = c("r_fixed", "r_mle"),
                       phi = target,
                       y_new = NULL,
                       n_new = NULL,
                       y_recycle = NULL,
                       n_recycle = NULL) {
  
  approx <- match.arg(approx)
  r_estimator <- match.arg(r_estimator)
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
  
  if (elimi[1L] == 1L) {
    return(list(
      MTD = 99L,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = approx,
      r_estimator = r_estimator,
      r_hat = r_hat,
      r_use = r_use
    ))
  }
  
  if (!any(n > 0L)) {
    return(list(
      MTD = 99L,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = approx,
      r_estimator = r_estimator,
      r_hat = r_hat,
      r_use = r_use
    ))
  }
  
  mu_hat <- rep(NA_real_, ndose)
  n_eff <- rep(0, ndose)
  
  for (j in seq_len(ndose)) {
    
    if (approx == "boin") {
      
      if (n[j] > 0) {
        mu_hat[j] <- y[j] / n[j]
        n_eff[j] <- n[j]
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
        } else if (r_estimator == "r_mle") {
          r_out <- estimate_r_mle_beta_binom(
            yR = YR,
            yI = YI,
            NR = NR,
            NI = NI,
            phi = phi
          )
          r_hat[j] <- r_out$r_hat
          r_use[j] <- r_hat[j]
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
      }
    }
  }
  
  nadmis <- min(
    max(which(elimi == 0L)),
    max(which(n != 0L))
  )
  
  mu_use <- mu_hat[seq_len(nadmis)]
  neff_use <- n_eff[seq_len(nadmis)]
  
  valid <- is.finite(mu_use) & is.finite(neff_use) & neff_use > 0
  
  if (!any(valid)) {
    return(list(
      MTD = NA_integer_,
      phat = phat_out,
      pj_iso = phat_out,
      eliminated = elimi,
      approx = approx,
      r_estimator = r_estimator,
      r_hat = r_hat,
      r_use = r_use
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
    r_use = r_use
  )
}
