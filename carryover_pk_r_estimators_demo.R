## ============================================================
## Carryover r and dose-specific p_k estimators
##
## Implements:
##   - global posterior estimator of r with p_k integrated out
##   - plug-in pooled moment estimator of p_k
##   - plug-in exact conditional MLE estimator of p_k
##
## Main update:
##   - add tryCatch wrappers around integrate()
##   - enter browser() only when an integration/replicate error occurs
##   - remove unconditional browser() calls
## ============================================================

## -------------------------------
## Small utilities
## -------------------------------

clamp <- function(x, lower = 0, upper = 1) {
  pmin(upper, pmax(lower, x))
}

check_counts <- function(yR, nR, yI, nI) {
  x <- list(yR = yR, nR = nR, yI = yI, nI = nI)
  len <- unique(vapply(x, length, integer(1)))
  
  if (length(len) != 1L) {
    stop("yR, nR, yI, and nI must have the same length.")
  }
  
  if (any(!is.finite(unlist(x)))) {
    stop("Counts must be finite.")
  }
  
  if (any(unlist(x) < 0)) {
    stop("Counts must be nonnegative.")
  }
  
  if (any(abs(unlist(x) - round(unlist(x))) > 1e-8)) {
    stop("Counts must be integers.")
  }
  
  if (any(yR > nR) || any(yI > nI)) {
    stop("Toxicity counts cannot exceed sample sizes.")
  }
  
  invisible(TRUE)
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

make_counts <- function(dat) {
  required <- c("dose", "nR", "yR", "nI", "yI")
  missing <- setdiff(required, names(dat))
  
  if (length(missing) > 0L) {
    stop("dat is missing columns: ", paste(missing, collapse = ", "))
  }
  
  optional <- intersect(c("p_true", "q_true"), names(dat))
  dat <- dat[order(dat$dose), c(required, optional)]
  
  check_counts(dat$yR, dat$nR, dat$yI, dat$nI)
  
  dat
}

## -------------------------------
## Two conditional estimators of p_k
## -------------------------------

pk_pool_moment_scalar <- function(yR, nR, yI, nI, r) {
  if (r < 0 || r >= 1) stop("r must be in [0, 1).")
  if (nR + nI == 0) return(NA_real_)
  
  regular_part <- if (nR > 0) yR else 0
  ipde_part <- if (nI > 0) {
    nI * ((yI / nI - r) / (1 - r))
  } else {
    0
  }
  
  clamp((regular_part + ipde_part) / (nR + nI), 0, 1)
}

pk_pool_moment <- function(yR, nR, yI, nI, r) {
  check_counts(yR, nR, yI, nI)
  
  vapply(
    seq_along(yR),
    function(k) {
      pk_pool_moment_scalar(
        yR = yR[k],
        nR = nR[k],
        yI = yI[k],
        nI = nI[k],
        r = r
      )
    },
    numeric(1)
  )
}

pk_conditional_mle_scalar <- function(yR, nR, yI, nI, r) {
  if (r < 0 || r >= 1) stop("r must be in [0, 1).")
  if (nR + nI == 0) return(NA_real_)
  
  ## Special cases.
  if (nI == 0) {
    return(clamp(yR / nR, 0, 1))
  }
  
  if (nR == 0) {
    return(clamp((yI / nI - r) / (1 - r), 0, 1))
  }
  
  N <- nR + nI
  B <- nR - yR + nI - yI
  
  b <- yR * (1 - 2 * r) -
    B * r +
    yI * (1 - r)
  
  disc <- b^2 + 4 * (1 - r) * N * yR * r
  
  p_hat <- (b + sqrt(max(disc, 0))) / (2 * (1 - r) * N)
  
  clamp(p_hat, 0, 1)
}

pk_conditional_mle <- function(yR, nR, yI, nI, r) {
  check_counts(yR, nR, yI, nI)
  
  vapply(
    seq_along(yR),
    function(k) {
      pk_conditional_mle_scalar(
        yR = yR[k],
        nR = nR[k],
        yI = yI[k],
        nI = nI[k],
        r = r
      )
    },
    numeric(1)
  )
}

observed_pooled_curve <- function(yR, nR, yI, nI) {
  check_counts(yR, nR, yI, nI)
  
  n <- nR + nI
  out <- rep(NA_real_, length(n))
  
  ok <- n > 0
  out[ok] <- (yR[ok] + yI[ok]) / n[ok]
  
  out
}

## -------------------------------
## Global posterior estimator of r
## -------------------------------

## Analytic p_k-integrated marginal likelihood contribution.
## This evaluates:
## int L_k(p_k, r) pi(p_k) dp_k
## by expanding {r + (1-r)p}^{yI}.
log_marginal_dose <- function(r, yR, nR, yI, nI, p_prior = c(1, 1)) {
  if (r < 0 || r >= 1) return(-Inf)
  if (nR + nI == 0) return(0)
  
  ap <- p_prior[1]
  bp <- p_prior[2]
  
  if (ap <= 0 || bp <= 0) {
    stop("p_prior must be positive.")
  }
  
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

log_posterior_r <- function(r,
                            counts,
                            r_prior = c(1, 9),
                            p_prior = c(1, 1),
                            r_max = 0.4) {
  if (r < 0 || r > r_max || r >= 1) return(-Inf)
  
  ar <- r_prior[1]
  br <- r_prior[2]
  
  if (ar <= 0 || br <= 0) {
    stop("r_prior must be positive.")
  }
  
  if (r_max <= 0 || r_max >= 1) {
    stop("r_max must be in (0, 1).")
  }
  
  ## Truncated beta prior on [0, r_max].
  log_prior <- stats::dbeta(r, ar, br, log = TRUE) -
    log(stats::pbeta(r_max, ar, br))
  
  log_m <- vapply(
    seq_len(nrow(counts)),
    function(k) {
      log_marginal_dose(
        r = r,
        yR = counts$yR[k],
        nR = counts$nR[k],
        yI = counts$yI[k],
        nI = counts$nI[k],
        p_prior = p_prior
      )
    },
    numeric(1)
  )
  
  log_prior + sum(log_m)
}

estimate_global_r <- function(dat,
                              r_prior = c(1, 9),
                              p_prior = c(1, 1),
                              r_max = 0.4,
                              rel.tol = 1e-8,
                              eps = 1e-8,
                              debug_on_error = TRUE) {
  counts <- make_counts(dat)
  
  ar <- r_prior[1]
  br <- r_prior[2]
  
  if (ar <= 0 || br <= 0) stop("r_prior must be positive.")
  if (r_max <= 0 || r_max >= 1) stop("r_max must be in (0, 1).")
  
  logpost_scalar <- function(r) {
    log_posterior_r(
      r = r,
      counts = counts,
      r_prior = r_prior,
      p_prior = p_prior,
      r_max = r_max
    )
  }
  
  logpost_vec <- function(r) {
    vapply(r, logpost_scalar, numeric(1))
  }
  
  ## Posterior MAP for r.
  ## If ar < 1, beta prior is singular at 0, so posterior MAP is often r = 0.
  opt_r <- stats::optimize(
    f = logpost_scalar,
    interval = c(eps, r_max - eps),
    maximum = TRUE
  )
  
  candidate_r <- c(eps, opt_r$maximum, r_max - eps)
  candidate_lp <- logpost_vec(candidate_r)
  
  if (ar < 1) {
    r_map <- 0
  } else {
    r_map <- candidate_r[which.max(candidate_lp)]
  }
  
  ## ==========================================================
  ## Stable transformed integration
  ##
  ## Problem:
  ##   With Beta(ar, br), ar < 1 gives density ~ r^(ar - 1),
  ##   which is singular at r = 0.
  ##
  ## Transformation:
  ##   r = r_max * u^(1/ar), u in [0, 1].
  ##
  ## Near u = 0:
  ##   r^(ar - 1) dr/du becomes approximately constant.
  ##
  ## This removes the endpoint singularity numerically.
  ## ==========================================================
  
  if (ar < 1) {
    r_of_u <- function(u) {
      r_max * u^(1 / ar)
    }
    
    log_jacobian <- function(u) {
      log(r_max / ar) + (1 / ar - 1) * log(u)
    }
  } else {
    r_of_u <- function(u) {
      r_max * u
    }
    
    log_jacobian <- function(u) {
      rep(log(r_max), length(u))
    }
  }
  
  u_lower <- eps
  u_upper <- 1 - eps
  
  log_denom_u_scalar <- function(u) {
    if (u <= 0 || u >= 1) return(-Inf)
    
    r <- r_of_u(u)
    
    lp <- logpost_scalar(r)
    lj <- log_jacobian(u)
    
    lp + lj
  }
  
  log_denom_u_vec <- function(u) {
    vapply(u, log_denom_u_scalar, numeric(1))
  }
  
  opt_u <- stats::optimize(
    f = log_denom_u_scalar,
    interval = c(u_lower, u_upper),
    maximum = TRUE
  )
  
  candidate_u <- c(u_lower, opt_u$maximum, u_upper)
  candidate_lu <- log_denom_u_vec(candidate_u)
  
  finite_lu <- is.finite(candidate_lu)
  
  if (!any(finite_lu)) {
    message("\nAll transformed log-integrand candidate values are non-finite.")
    message("candidate_u = ", paste(candidate_u, collapse = ", "))
    message("candidate_lu = ", paste(candidate_lu, collapse = ", "))
    
    if (debug_on_error) browser()
    
    stop("Cannot scale transformed posterior integrand.")
  }
  
  lu_scale <- max(candidate_lu[finite_lu])
  
  denom_integrand_u <- function(u) {
    exp(log_denom_u_vec(u) - lu_scale)
  }
  
  numer_integrand_u <- function(u) {
    r <- r_of_u(u)
    r * exp(log_denom_u_vec(u) - lu_scale)
  }
  
  safe_integrate <- function(label, f, lower, upper) {
    tryCatch(
      {
        stats::integrate(
          f,
          lower = lower,
          upper = upper,
          rel.tol = rel.tol,
          subdivisions = 500L,
          stop.on.error = TRUE
        )
      },
      error = function(e) {
        message("\n================ TRANSFORMED INTEGRATION ERROR ================")
        message("Failed integral: ", label)
        message("Error message: ", conditionMessage(e))
        message("u lower = ", lower)
        message("u upper = ", upper)
        message("r_max = ", r_max)
        message("eps = ", eps)
        message("r_prior = ", paste(r_prior, collapse = ", "))
        message("p_prior = ", paste(p_prior, collapse = ", "))
        message("opt_r$maximum = ", opt_r$maximum)
        message("opt_r$objective = ", opt_r$objective)
        message("candidate_r = ", paste(candidate_r, collapse = ", "))
        message("candidate_lp = ", paste(candidate_lp, collapse = ", "))
        message("opt_u$maximum = ", opt_u$maximum)
        message("opt_u$objective = ", opt_u$objective)
        message("candidate_u = ", paste(candidate_u, collapse = ", "))
        message("candidate_lu = ", paste(candidate_lu, collapse = ", "))
        message("lu_scale = ", lu_scale)
        
        debug_u <- seq(lower, upper, length.out = 500)
        debug_r <- r_of_u(debug_u)
        debug_logpost <- try(logpost_vec(debug_r), silent = TRUE)
        debug_log_integrand <- try(log_denom_u_vec(debug_u), silent = TRUE)
        debug_kernel <- try(f(debug_u), silent = TRUE)
        
        if (!inherits(debug_logpost, "try-error")) {
          finite_logpost <- is.finite(debug_logpost)
          message(
            "logpost finite count = ",
            sum(finite_logpost),
            "/",
            length(debug_logpost)
          )
          
          if (any(finite_logpost)) {
            message(
              "logpost finite range = ",
              paste(range(debug_logpost[finite_logpost]), collapse = ", ")
            )
          }
        }
        
        if (!inherits(debug_log_integrand, "try-error")) {
          finite_li <- is.finite(debug_log_integrand)
          message(
            "transformed log-integrand finite count = ",
            sum(finite_li),
            "/",
            length(debug_log_integrand)
          )
          
          if (any(finite_li)) {
            message(
              "transformed log-integrand finite range = ",
              paste(range(debug_log_integrand[finite_li]), collapse = ", ")
            )
          }
        }
        
        if (!inherits(debug_kernel, "try-error")) {
          finite_kernel <- is.finite(debug_kernel)
          message(
            "kernel finite count = ",
            sum(finite_kernel),
            "/",
            length(debug_kernel)
          )
          
          if (any(finite_kernel)) {
            message(
              "kernel finite range = ",
              paste(range(debug_kernel[finite_kernel]), collapse = ", ")
            )
          }
        }
        
        message("Objects available in browser():")
        message("  counts, r_prior, p_prior, r_max, eps")
        message("  opt_r, candidate_r, candidate_lp")
        message("  opt_u, candidate_u, candidate_lu, lu_scale")
        message("  r_of_u, log_jacobian, log_denom_u_vec")
        message("  debug_u, debug_r, debug_logpost, debug_log_integrand, debug_kernel")
        message("===============================================================\n")
        
        if (debug_on_error) browser()
        
        stop(e)
      }
    )
  }
  
  denom <- safe_integrate(
    label = "denom",
    f = denom_integrand_u,
    lower = u_lower,
    upper = u_upper
  )$value
  
  numer <- safe_integrate(
    label = "numer",
    f = numer_integrand_u,
    lower = u_lower,
    upper = u_upper
  )$value
  
  if (!is.finite(denom) || denom <= 0) {
    message("\nInvalid denominator after transformed integration.")
    message("denom = ", denom)
    if (debug_on_error) browser()
    stop("Posterior normalization constant is non-positive or non-finite.")
  }
  
  r_mean <- numer / denom
  
  if (!is.finite(r_mean)) {
    message("\nInvalid r_mean.")
    message("numer = ", numer)
    message("denom = ", denom)
    if (debug_on_error) browser()
    stop("r_mean is non-finite.")
  }
  
  var_integrand_u <- function(u) {
    r <- r_of_u(u)
    (r - r_mean)^2 * exp(log_denom_u_vec(u) - lu_scale)
  }
  
  var_num <- safe_integrate(
    label = "variance numerator",
    f = var_integrand_u,
    lower = u_lower,
    upper = u_upper
  )$value
  
  r_sd <- sqrt(var_num / denom)
  
  list(
    r_mean = r_mean,
    r_map = r_map,
    r_sd = r_sd,
    log_post_at_map = logpost_scalar(ifelse(r_map == 0, eps, r_map)),
    r_prior = r_prior,
    p_prior = p_prior,
    r_max = r_max,
    counts = counts
  )
}

estimate_pk_with_global_r <- function(dat,
                                      r_prior = c(1, 9),
                                      p_prior = c(1, 1),
                                      r_max = 0.4,
                                      plug_in = c("mean", "map"),
                                      debug_on_error = TRUE) {
  plug_in <- match.arg(plug_in)
  counts <- make_counts(dat)
  
  r_fit <- estimate_global_r(
    dat = counts,
    r_prior = r_prior,
    p_prior = p_prior,
    r_max = r_max,
    debug_on_error = debug_on_error
  )
  
  r_hat <- if (plug_in == "mean") {
    r_fit$r_mean
  } else {
    r_fit$r_map
  }
  
  ## Avoid exactly r = 1.
  r_hat <- clamp(r_hat, 0, min(r_max, 1 - 1e-8))
  
  out <- counts
  out$obs_pooled <- observed_pooled_curve(out$yR, out$nR, out$yI, out$nI)
  out$p_pool_plugin <- pk_pool_moment(out$yR, out$nR, out$yI, out$nI, r_hat)
  out$p_mle_plugin <- pk_conditional_mle(out$yR, out$nR, out$yI, out$nI, r_hat)
  
  list(
    r = r_fit,
    r_hat = r_hat,
    plug_in = plug_in,
    estimates = out
  )
}

## -------------------------------
## Simulation and plotting helpers
## -------------------------------

expand_n_by_dose <- function(n, K, name) {
  if (length(n) == 1L) {
    n <- rep(n, K)
  }
  
  if (length(n) != K) {
    stop(name, " must be length 1 or length K.")
  }
  
  if (any(n < 0) || any(abs(n - round(n)) > 1e-8)) {
    stop(name, " must contain nonnegative integer counts.")
  }
  
  as.integer(n)
}

simulate_regular_ipde_counts <- function(p_true,
                                         n_regular = 20,
                                         n_ipde = 20,
                                         alpha_true = 0.1,
                                         seed = NULL,
                                         ipde_model = c("document", "aide_alpha")) {
  ipde_model <- match.arg(ipde_model)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  p_true <- as.numeric(p_true)
  K <- length(p_true)
  
  if (any(!is.finite(p_true)) || any(p_true < 0 | p_true > 1)) {
    stop("p_true must contain probabilities in [0, 1].")
  }
  
  if (length(alpha_true) != 1L || !is.finite(alpha_true) || alpha_true < 0) {
    stop("alpha_true must be a finite nonnegative scalar.")
  }
  
  nR <- expand_n_by_dose(n_regular, K, "n_regular")
  nI <- expand_n_by_dose(n_ipde, K, "n_ipde")
  
  ## In the document model, alpha_true is the simulation truth for r:
  ## q_k = alpha_true + (1 - alpha_true) p_k, capped at 1.
  ##
  ## The optional "aide_alpha" mode matches older AIDE-style simulation:
  ## q_k = p_k + alpha_true p_{k-1}, with dose 1 unchanged.
  if (ipde_model == "document") {
    q_true <- pmin(1, alpha_true + (1 - alpha_true) * p_true)
  } else {
    q_true <- p_true
    if (K >= 2L) {
      q_true[-1] <- pmin(1, p_true[-1] + alpha_true * p_true[-K])
    }
  }
  
  data.frame(
    dose = seq_len(K),
    p_true = p_true,
    q_true = q_true,
    nR = nR,
    yR = stats::rbinom(K, nR, p_true),
    nI = nI,
    yI = stats::rbinom(K, nI, q_true)
  )
}

plot_carryover_fit <- function(fit,
                               scenario_name = "scenario",
                               alpha_true = NA_real_,
                               file = NULL,
                               ylim = c(0, 1)) {
  est <- fit$estimates
  dose <- est$dose
  
  main <- paste0(
    scenario_name,
    ", alpha=",
    alpha_true,
    ", r_hat=",
    sprintf("%.3f", fit$r_hat)
  )
  
  if (!is.null(file)) {
    grDevices::png(file, width = 950, height = 650, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  graphics::plot(
    dose,
    est$p_true,
    type = "b",
    pch = 16,
    lwd = 2,
    ylim = ylim,
    xlab = "Dose level",
    ylab = "DLT probability",
    main = main
  )
  
  graphics::lines(
    dose,
    est$obs_pooled,
    type = "b",
    pch = 1,
    lwd = 2,
    lty = 3,
    col = "gray35"
  )
  
  graphics::lines(
    dose,
    est$p_pool_plugin,
    type = "b",
    pch = 17,
    lwd = 2,
    col = "#0072B2"
  )
  
  graphics::lines(
    dose,
    est$p_mle_plugin,
    type = "b",
    pch = 15,
    lwd = 2,
    col = "#D55E00"
  )
  
  graphics::abline(h = 0.30, col = "gray75", lty = 2)
  
  graphics::legend(
    "topleft",
    bty = "n",
    lwd = 2,
    pch = c(16, 1, 17, 15),
    lty = c(1, 3, 1, 1),
    col = c("black", "gray35", "#0072B2", "#D55E00"),
    legend = c(
      "True current toxicity p_k",
      "Observed pooled regular + IPDE",
      "Plug-in pooled moment p_k",
      "Plug-in exact MLE p_k"
    )
  )
}

run_carryover_estimator_demo <- function(scenarios,
                                         alpha_values = c(0, 0.1, 0.25, 0.4),
                                         n_regular = 24,
                                         n_ipde = 24,
                                         r_prior = c(1, 9),
                                         p_prior = c(1, 1),
                                         r_max = 0.6,
                                         plug_in = "mean",
                                         ipde_model = "document",
                                         seed = 20260628,
                                         plot_dir = "carryover_pk_r_demo_plots",
                                         debug_on_error = TRUE) {
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  all_rows <- list()
  idx <- 1L
  
  scenario_names <- names(scenarios)
  
  if (is.null(scenario_names) || any(!nzchar(scenario_names))) {
    scenario_names <- paste0("scenario_", seq_along(scenarios))
  }
  
  for (s in seq_along(scenarios)) {
    p_true <- scenarios[[s]]
    scenario_name <- scenario_names[s]
    
    for (a in alpha_values) {
      dat <- simulate_regular_ipde_counts(
        p_true = p_true,
        n_regular = n_regular,
        n_ipde = n_ipde,
        alpha_true = a,
        seed = seed + idx,
        ipde_model = ipde_model
      )
      
      fit <- tryCatch(
        {
          estimate_pk_with_global_r(
            dat = dat,
            r_prior = r_prior,
            p_prior = p_prior,
            r_max = r_max,
            plug_in = plug_in,
            debug_on_error = debug_on_error
          )
        },
        error = function(e) {
          message("\n================ DEMO ERROR ================")
          message("scenario = ", scenario_name)
          message("alpha_true = ", a)
          message("idx = ", idx)
          message("seed used = ", seed + idx)
          message("Error message: ", conditionMessage(e))
          print(dat)
          message("Objects available in browser():")
          message("  dat, scenario_name, a, idx")
          message("  r_prior, p_prior, r_max, plug_in")
          message("===========================================\n")
          
          if (debug_on_error) {
            browser()
          }
          
          stop(e)
        }
      )
      
      est <- fit$estimates
      
      one <- data.frame(
        scenario = scenario_name,
        alpha_true = a,
        r_hat_mean = fit$r$r_mean,
        r_hat_map = fit$r$r_map,
        r_hat_sd = fit$r$r_sd,
        dose = est$dose,
        nR = est$nR,
        yR = est$yR,
        nI = est$nI,
        yI = est$yI,
        p_true = est$p_true,
        q_true = est$q_true,
        obs_pooled = est$obs_pooled,
        p_pool_plugin = est$p_pool_plugin,
        p_mle_plugin = est$p_mle_plugin,
        stringsAsFactors = FALSE
      )
      
      all_rows[[idx]] <- one
      
      plot_file <- file.path(
        plot_dir,
        paste0(
          gsub("[^A-Za-z0-9]+", "_", scenario_name),
          "_alpha_",
          gsub("\\.", "p", a),
          ".png"
        )
      )
      
      plot_carryover_fit(
        fit = fit,
        scenario_name = scenario_name,
        alpha_true = a,
        file = plot_file
      )
      
      idx <- idx + 1L
    }
  }
  
  out <- do.call(rbind, all_rows)
  rownames(out) <- NULL
  
  out
}

summarize_replicate_estimates <- function(replicate_results) {
  if (is.null(replicate_results) || nrow(replicate_results) == 0L) {
    stop("replicate_results must contain at least one row.")
  }
  
  keys <- interaction(
    replicate_results$scenario,
    replicate_results$alpha_true,
    replicate_results$dose,
    drop = TRUE
  )
  
  out <- do.call(
    rbind,
    lapply(
      split(replicate_results, keys),
      function(x) {
        unique_r_by_rep <- unique(x[, c("replicate", "r_hat_used")])
        
        data.frame(
          scenario = x$scenario[1],
          alpha_true = x$alpha_true[1],
          dose = x$dose[1],
          n_reps = length(unique(x$replicate)),
          p_true = x$p_true[1],
          q_true = x$q_true[1],
          mean_yR = mean(x$yR),
          mean_yI = mean(x$yI),
          mean_obs_pooled = mean(x$obs_pooled, na.rm = TRUE),
          sd_obs_pooled = stats::sd(x$obs_pooled, na.rm = TRUE),
          mean_p_pool_plugin = mean(x$p_pool_plugin, na.rm = TRUE),
          sd_p_pool_plugin = stats::sd(x$p_pool_plugin, na.rm = TRUE),
          bias_p_pool_plugin = mean(x$p_pool_plugin - x$p_true, na.rm = TRUE),
          mean_p_mle_plugin = mean(x$p_mle_plugin, na.rm = TRUE),
          sd_p_mle_plugin = stats::sd(x$p_mle_plugin, na.rm = TRUE),
          bias_p_mle_plugin = mean(x$p_mle_plugin - x$p_true, na.rm = TRUE),
          mean_r_hat = mean(x$r_hat_used, na.rm = TRUE),
          sd_r_hat = stats::sd(unique_r_by_rep$r_hat_used, na.rm = TRUE),
          mean_r_hat_posterior_mean = mean(x$r_hat_mean, na.rm = TRUE),
          mean_r_hat_map = mean(x$r_hat_map, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  rownames(out) <- NULL
  out[order(out$scenario, out$alpha_true, out$dose), ]
}

plot_carryover_average <- function(summary_rows,
                                   scenario_name = NULL,
                                   alpha_true = NULL,
                                   file = NULL,
                                   ylim = c(0, 1)) {
  dat <- summary_rows
  
  if (!is.null(scenario_name)) {
    dat <- dat[dat$scenario == scenario_name, , drop = FALSE]
  }
  
  if (!is.null(alpha_true)) {
    dat <- dat[dat$alpha_true == alpha_true, , drop = FALSE]
  }
  
  if (nrow(dat) == 0L) {
    stop("No rows to plot after filtering.")
  }
  
  if (length(unique(dat$scenario)) != 1L ||
      length(unique(dat$alpha_true)) != 1L) {
    stop("summary_rows must contain exactly one scenario and alpha_true for this plot.")
  }
  
  scenario_name <- unique(dat$scenario)
  alpha_true <- unique(dat$alpha_true)
  dat <- dat[order(dat$dose), ]
  
  main <- paste0(
    scenario_name,
    ", alpha=",
    alpha_true,
    ", mean r_hat=",
    sprintf("%.3f", unique(dat$mean_r_hat)[1])
  )
  
  if (!is.null(file)) {
    grDevices::png(file, width = 950, height = 650, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  graphics::plot(
    dat$dose,
    dat$p_true,
    type = "b",
    pch = 16,
    lwd = 2,
    ylim = ylim,
    xlab = "Dose level",
    ylab = "Average DLT probability estimate",
    main = main
  )
  
  graphics::lines(
    dat$dose,
    dat$mean_obs_pooled,
    type = "b",
    pch = 1,
    lwd = 2,
    lty = 3,
    col = "gray35"
  )
  
  graphics::lines(
    dat$dose,
    dat$mean_p_pool_plugin,
    type = "b",
    pch = 17,
    lwd = 2,
    col = "#0072B2"
  )
  
  graphics::lines(
    dat$dose,
    dat$mean_p_mle_plugin,
    type = "b",
    pch = 15,
    lwd = 2,
    col = "#D55E00"
  )
  
  graphics::abline(h = 0.30, col = "gray75", lty = 2)
  
  graphics::legend(
    "topleft",
    bty = "n",
    lwd = 2,
    pch = c(16, 1, 17, 15),
    lty = c(1, 3, 1, 1),
    col = c("black", "gray35", "#0072B2", "#D55E00"),
    legend = c(
      "True current toxicity p_k",
      "Average observed pooled regular + IPDE",
      "Average plug-in pooled moment p_k",
      "Average plug-in exact MLE p_k"
    )
  )
}

plot_carryover_average_by_alpha <- function(summary_rows,
                                            alpha_true,
                                            scenario_names = NULL,
                                            file = NULL,
                                            ylim = c(0, 1),
                                            panel_ncol = 3L) {
  dat <- summary_rows[summary_rows$alpha_true == alpha_true, , drop = FALSE]
  
  if (!is.null(scenario_names)) {
    dat <- dat[dat$scenario %in% scenario_names, , drop = FALSE]
  } else {
    scenario_names <- unique(dat$scenario)
  }
  
  if (nrow(dat) == 0L) {
    stop("No rows to plot for alpha_true = ", alpha_true, ".")
  }
  
  scenario_names <- scenario_names[scenario_names %in% unique(dat$scenario)]
  panel_ncol <- as.integer(panel_ncol)
  panel_nrow <- ceiling(length(scenario_names) / panel_ncol)
  
  if (!is.null(file)) {
    grDevices::png(file, width = 1500, height = 900, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  graphics::par(
    mfrow = c(panel_nrow, panel_ncol),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 3, 0)
  )
  
  for (i in seq_along(scenario_names)) {
    scenario_name <- scenario_names[i]
    one <- dat[dat$scenario == scenario_name, , drop = FALSE]
    one <- one[order(one$dose), ]
    
    main <- paste0(
      scenario_name,
      ", mean r_hat=",
      sprintf("%.3f", unique(one$mean_r_hat)[1])
    )
    
    graphics::plot(
      one$dose,
      one$p_true,
      type = "b",
      pch = 16,
      lwd = 2,
      ylim = ylim,
      xlab = "Dose level",
      ylab = "DLT probability",
      main = main
    )
    
    graphics::lines(
      one$dose,
      one$mean_obs_pooled,
      type = "b",
      pch = 1,
      lwd = 2,
      lty = 3,
      col = "gray35"
    )
    
    graphics::lines(
      one$dose,
      one$mean_p_pool_plugin,
      type = "b",
      pch = 17,
      lwd = 2,
      col = "#0072B2"
    )
    
    graphics::lines(
      one$dose,
      one$mean_p_mle_plugin,
      type = "b",
      pch = 15,
      lwd = 2,
      col = "#D55E00"
    )
    
    graphics::abline(h = 0.30, col = "gray75", lty = 2)
    
    if (i == 1L) {
      graphics::legend(
        "topleft",
        bty = "n",
        cex = 0.75,
        lwd = 2,
        pch = c(16, 1, 17, 15),
        lty = c(1, 3, 1, 1),
        col = c("black", "gray35", "#0072B2", "#D55E00"),
        legend = c(
          "True p_k",
          "Average observed pooled",
          "Average pooled moment",
          "Average exact MLE"
        )
      )
    }
  }
  
  graphics::mtext(
    paste0("Carryover estimator average curves, alpha = ", alpha_true),
    outer = TRUE,
    cex = 1.2,
    font = 2
  )
  
  invisible(file)
}

plot_carryover_average_alpha_files <- function(summary_rows,
                                               alpha_values,
                                               scenario_names,
                                               plot_dir = "carryover_pk_r_demo_plots") {
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  files <- vapply(
    alpha_values,
    function(a) {
      plot_file <- file.path(
        plot_dir,
        paste0(
          "sce1_to_sce6_alpha_",
          gsub("\\.", "p", a),
          "_avg.png"
        )
      )
      
      plot_carryover_average_by_alpha(
        summary_rows = summary_rows,
        alpha_true = a,
        scenario_names = scenario_names,
        file = plot_file
      )
      
      plot_file
    },
    character(1)
  )
  
  names(files) <- paste0("alpha_", alpha_values)
  files
}

run_carryover_estimator_replicates <- function(scenarios,
                                               alpha_values = c(0, 0.1, 0.25, 0.4),
                                               n_reps = 100,
                                               n_regular = 24,
                                               n_ipde = 24,
                                               r_prior = c(1, 9),
                                               p_prior = c(1, 1),
                                               r_max = 0.6,
                                               plug_in = "mean",
                                               ipde_model = "document",
                                               seed = 20260628,
                                               plot_dir = "carryover_pk_r_replicate_plots",
                                               debug_on_error = TRUE) {
  if (length(n_reps) != 1L || !is.finite(n_reps) || n_reps < 1L) {
    stop("n_reps must be a positive integer.")
  }
  
  n_reps <- as.integer(n_reps)
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  all_rows <- list()
  idx <- 1L
  
  scenario_names <- names(scenarios)
  
  if (is.null(scenario_names) || any(!nzchar(scenario_names))) {
    scenario_names <- paste0("scenario_", seq_along(scenarios))
  }
  
  for (s in seq_along(scenarios)) {
    p_true <- scenarios[[s]]
    scenario_name <- scenario_names[s]
    
    for (a in alpha_values) {
      for (rep_id in seq_len(n_reps)) {
        dat <- simulate_regular_ipde_counts(
          p_true = p_true,
          n_regular = n_regular,
          n_ipde = n_ipde,
          alpha_true = a,
          seed = seed + idx,
          ipde_model = ipde_model
        )
        
        fit <- tryCatch(
          {
            estimate_pk_with_global_r(
              dat = dat,
              r_prior = r_prior,
              p_prior = p_prior,
              r_max = r_max,
              plug_in = plug_in,
              debug_on_error = debug_on_error
            )
          },
          error = function(e) {
            message("\n================ REPLICATE ERROR ================")
            message("scenario = ", scenario_name)
            message("alpha_true = ", a)
            message("replicate = ", rep_id)
            message("idx = ", idx)
            message("seed used = ", seed + idx)
            message("Error message: ", conditionMessage(e))
            print(dat)
            message("Objects available in browser():")
            message("  dat, scenario_name, a, rep_id, idx")
            message("  r_prior, p_prior, r_max, plug_in")
            message("=================================================\n")
            
            if (debug_on_error) {
              browser()
            }
            
            stop(e)
          }
        )
        
        est <- fit$estimates
        
        all_rows[[idx]] <- data.frame(
          scenario = scenario_name,
          alpha_true = a,
          replicate = rep_id,
          r_hat_used = fit$r_hat,
          r_hat_mean = fit$r$r_mean,
          r_hat_map = fit$r$r_map,
          r_hat_sd = fit$r$r_sd,
          dose = est$dose,
          nR = est$nR,
          yR = est$yR,
          nI = est$nI,
          yI = est$yI,
          p_true = est$p_true,
          q_true = est$q_true,
          obs_pooled = est$obs_pooled,
          p_pool_plugin = est$p_pool_plugin,
          p_mle_plugin = est$p_mle_plugin,
          stringsAsFactors = FALSE
        )
        
        idx <- idx + 1L
      }
    }
  }
  
  replicate_results <- do.call(rbind, all_rows)
  rownames(replicate_results) <- NULL
  
  summary_results <- summarize_replicate_estimates(replicate_results)
  
  for (scenario_name in unique(summary_results$scenario)) {
    alpha_set <- unique(summary_results$alpha_true[summary_results$scenario == scenario_name])
    
    for (a in alpha_set) {
      plot_file <- file.path(
        plot_dir,
        paste0(
          gsub("[^A-Za-z0-9]+", "_", scenario_name),
          "_alpha_",
          gsub("\\.", "p", a),
          "_avg.png"
        )
      )
      
      plot_carryover_average(
        summary_rows = summary_results,
        scenario_name = scenario_name,
        alpha_true = a,
        file = plot_file
      )
    }
  }
  
  plot_carryover_average_alpha_files(
    summary_rows = summary_results,
    alpha_values = alpha_values,
    scenario_names = scenario_names,
    plot_dir = plot_dir
  )
  
  list(
    replicate_results = replicate_results,
    summary_results = summary_results,
    n_reps = n_reps,
    plot_dir = plot_dir
  )
}

## -------------------------------
## Example test run
## -------------------------------

example_scenarios <- list(
  sce1 = c(0.30, 0.35, 0.40, 0.45, 0.50),
  sce2 = c(0.15, 0.30, 0.38, 0.45, 0.55),
  sce3 = c(0.15, 0.20, 0.30, 0.35, 0.45),
  sce4 = c(0.05, 0.10, 0.18, 0.30, 0.40),
  sce5 = c(0.07, 0.12, 0.17, 0.22, 0.30)
)

## Your original prior:
##   r_prior = c(0.3 / 2, 1 - 0.3 / 2)
## gives Beta(0.15, 0.85), which has infinite density at r = 0.
## This is mathematically proper but numerically unstable for integrate().
##
## For stable testing, I recommend this prior with the same mean 0.15:
##   Beta(1.5, 8.5), ESS = 10.
##
## You can switch back to Beta(0.15, 0.85) if you want to debug the boundary spike.

demo_results <- run_carryover_estimator_replicates(
  scenarios = example_scenarios,
  alpha_values = c(1.2, 1.5),
  n_reps = 10,
  n_regular = rep(100, 5),
  n_ipde = rep(100, 5),
  
  ## More stable version of prior mean 0.15:
  r_prior = c(0.15, 0.85),
  p_prior = c(0.3 / 2, 1 - (0.3 / 2)),
  
  ## You used r_max = 0.99.
  ## This is okay, but r near 1 can make p_k plug-in unstable.
  r_max = 0.99,
  
  plug_in = "mean",
  ipde_model = "document",
  seed = 20260628,
  plot_dir = "carryover_pk_r_demo_plots",
  debug_on_error = TRUE
)

print(demo_results$summary_results)

write.csv(
  demo_results$replicate_results,
  "carryover_pk_r_replicate_results.csv",
  row.names = FALSE
)

write.csv(
  demo_results$summary_results,
  "carryover_pk_r_summary_results.csv",
  row.names = FALSE
)

cat("\nSaved replicate table: carryover_pk_r_replicate_results.csv\n")
cat("Saved summary table: carryover_pk_r_summary_results.csv\n")
cat("Saved plots under: carryover_pk_r_demo_plots/\n")

