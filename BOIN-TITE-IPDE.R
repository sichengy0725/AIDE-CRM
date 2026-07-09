source('TITE-BOIN-2.R')
## ----------------------------------------------------------
## final MTD selection: binary version (original BOIN/TITE-BOIN)
## ----------------------------------------------------------
## ==========================================================
## Falke-style helper for alpha_ij
## alpha_ij(t) = [ current dose + decayed previous doses ] / current dose
## with decay fraction falke_alpha^(age / window), default falke_alpha = 0.5
## ==========================================================
compute_alpha_ij <- function(patient,
                             row_index,
                             t_ref,
                             window,
                             falke_alpha = 0.5) {
  if (length(row_index) == 0L) return(numeric(0))
  
  if (length(t_ref) == 1L) {
    t_ref <- rep(t_ref, length(row_index))
  }
  if (length(t_ref) != length(row_index)) {
    stop("t_ref must have length 1 or length(row_index)")
  }
  if (window <= 0) {
    stop("window must be > 0 for Falke-style alpha_ij")
  }
  if (falke_alpha <= 0 || falke_alpha > 1) {
    stop("falke_alpha must be in (0, 1]")
  }
  
  out <- rep(1, length(row_index))
  
  for (k in seq_along(row_index)) {
    rr <- row_index[k]
    
    cur_id    <- patient$id[rr]
    cur_cycle <- patient$cycle[rr]
    cur_dose  <- patient$dose[rr]
    
    if (is.na(cur_dose) || cur_dose <= 0) {
      out[k] <- NA_real_
      next
    }
    
    if (cur_cycle <= 1L) {
      out[k] <- 1
      next
    }
    
    prev_idx <- which(patient$id == cur_id & patient$cycle < cur_cycle)
    
    if (length(prev_idx) == 0L) {
      out[k] <- 1
      next
    }
    
    age <- pmax(t_ref[k] - patient$arrival_time[prev_idx], 0)
    
    carry_dose <- sum(patient$dose[prev_idx] * falke_alpha^(age / window))
    
    out[k] <- (cur_dose + carry_dose) / cur_dose
  }
  
  out
}

## ==========================================================
## Collect per-row contributions for TITE-gBOIN using the new MLE
##
## MLE:
##   p_hat = A / (A + B)
## where
##   A = sum observed DLTs
##   B = sum regular no-DLT contribution
##       + sum alpha_ij * IPDE no-DLT contribution
##
## For pending rows, no-DLT contribution is time_weight.
##
## NOTE:
## - discount_alpha is kept only for backward compatibility.
## - falke_alpha is the carryover decay parameter in Falke's form.
## ==========================================================
collect_gboin_dose_data <- function(patient,
                                    dose_level,
                                    t_now,
                                    window,
                                    discount_alpha = 0,
                                    falke_alpha = 0.5) {
  idx <- which(patient$dose == dose_level & patient$arrival_time <= t_now)
  
  if (length(idx) == 0L) {
    return(data.frame(
      id = integer(0),
      cycle = integer(0),
      dose = integer(0),
      arrival_time = numeric(0),
      eval_time = numeric(0),
      DLT_time = numeric(0),
      ipde_true = integer(0),
      y = integer(0),
      prev_dose = numeric(0),
      time_follow = numeric(0),
      time_weight = numeric(0),
      alpha_ij = numeric(0),
      observed = logical(0),
      dlt_observed = logical(0),
      no_dlt_complete = logical(0),
      pending = logical(0),
      dlt_contrib = numeric(0),
      nodlt_contrib = numeric(0),
      B_contrib = numeric(0),
      mle_den_contrib = numeric(0),
      ex_contrib = numeric(0),      # alias for compatibility
      n_eff_contrib = numeric(0)    # alias for compatibility
    ))
  }
  
  dat <- patient[idx, , drop = FALSE]
  
  ## previous-cycle dose for display/debug
  key_all  <- paste(patient$id, patient$cycle)
  prev_key <- paste(dat$id, dat$cycle - 1L)
  prev_idx <- match(prev_key, key_all)
  
  prev_dose <- rep(NA_real_, nrow(dat))
  has_prev <- !is.na(prev_idx) & dat$cycle > 1L
  prev_dose[has_prev] <- patient$dose[prev_idx[has_prev]]
  
  ## TITE follow-up
  if (window > 0) {
    time_follow <- pmax(0, pmin(t_now - dat$arrival_time, window))
    time_weight <- pmin(time_follow / window, 1)
  } else {
    time_follow <- rep(0, nrow(dat))
    time_weight <- rep(1, nrow(dat))
  }
  
  is_ipde <- (dat$ipde_true == 1L) & (dat$cycle > 1L)
  
  ## observed status at t_now
  dlt_observed    <- (dat$y == 1L) & is.finite(dat$DLT_time) & (dat$DLT_time <= t_now)
  no_dlt_complete <- (dat$y == 0L) & (dat$eval_time <= t_now)
  observed        <- dlt_observed | no_dlt_complete
  pending         <- !observed
  
  ## choose reference time for alpha_ij
  ## - observed DLT: DLT time
  ## - completed no DLT: end of assessment
  ## - pending: current decision time
  time_ref <- ifelse(
    dlt_observed,
    dat$DLT_time,
    ifelse(no_dlt_complete, dat$eval_time, t_now)
  )
  
  alpha_ij <- rep(1, nrow(dat))
  ipde_rows <- which(is_ipde)
  
  if (length(ipde_rows) > 0L) {
    alpha_ij[ipde_rows] <- compute_alpha_ij(
      patient    = patient,
      row_index  = idx[ipde_rows],
      t_ref      = time_ref[ipde_rows],
      window     = window,
      falke_alpha = falke_alpha
    )
  }
  
  ## A-part of MLE: observed DLT count
  dlt_contrib <- as.numeric(dlt_observed)
  
  ## no-DLT contribution:
  ## completed no DLT -> 1
  ## pending          -> time_weight
  ## observed DLT     -> 0
  nodlt_contrib <- ifelse(no_dlt_complete, 1, ifelse(pending, time_weight, 0))
  
  ## B-part of MLE:
  ## regular row -> nodlt_contrib
  ## IPDE row    -> alpha_ij * nodlt_contrib
  B_contrib <- ifelse(is_ipde, alpha_ij, 1) * nodlt_contrib
  
  ## denominator contribution for p_hat = A / (A + B)
  mle_den_contrib <- dlt_contrib + B_contrib
  
  out <- dat
  out$prev_dose       <- prev_dose
  out$time_follow     <- time_follow
  out$time_weight     <- time_weight
  out$alpha_ij        <- alpha_ij
  out$observed        <- observed
  out$dlt_observed    <- dlt_observed
  out$no_dlt_complete <- no_dlt_complete
  out$pending         <- pending
  out$dlt_contrib     <- dlt_contrib
  out$nodlt_contrib   <- nodlt_contrib
  out$B_contrib       <- B_contrib
  out$mle_den_contrib <- mle_den_contrib
  
  ## backward-compatible aliases
  out$ex_contrib    <- dlt_contrib
  out$n_eff_contrib <- mle_den_contrib
  
  out
}

## ==========================================================
## TITE-gBOIN move using the approximate MLE
##   p_hat = A / (A + B)
## ==========================================================
boin_move_gboin_tite <- function(current_dose,
                                 patient,
                                 t_now,
                                 ndose,
                                 window,
                                 lambda_e,
                                 lambda_d,
                                 elimi,
                                 discount_alpha = 0,
                                 falke_alpha = 0.5) {
  admissible <- which(elimi == 0L)
  if (length(admissible) == 0L) {
    return(list(
      next_dose = NA_integer_,
      action = "stop",
      mu_hat = NA_real_,
      A = NA_real_,
      B = NA_real_,
      denom = NA_real_,
      ex = NA_real_,
      n_eff = NA_real_,
      dose_data = NULL
    ))
  }
  
  max_admissible <- max(admissible)
  
  dose_data <- collect_gboin_dose_data(
    patient        = patient,
    dose_level     = current_dose,
    t_now          = t_now,
    window         = window,
    discount_alpha = discount_alpha,  # kept for compatibility
    falke_alpha    = falke_alpha
  )
  
  if (nrow(dose_data) == 0L) {
    return(list(
      next_dose = current_dose,
      action = "stay",
      mu_hat = NA_real_,
      A = 0,
      B = 0,
      denom = 0,
      ex = 0,
      n_eff = 0,
      dose_data = dose_data
    ))
  }
  
  A <- sum(dose_data$dlt_contrib, na.rm = TRUE)
  B <- sum(dose_data$B_contrib,   na.rm = TRUE)
  denom <- A + B
  
  if (denom <= 0) {
    return(list(
      next_dose = current_dose,
      action = "stay",
      mu_hat = NA_real_,
      A = A,
      B = B,
      denom = denom,
      ex = A,
      n_eff = denom,
      dose_data = dose_data
    ))
  }
  
  mu_hat <- A / denom
  
  if (mu_hat <= lambda_e) {
    cand <- current_dose + 1L
    if (cand <= max_admissible) {
      return(list(
        next_dose = cand,
        action = "escalate",
        mu_hat = mu_hat,
        A = A,
        B = B,
        denom = denom,
        ex = A,
        n_eff = denom,
        dose_data = dose_data
      ))
    } else {
      return(list(
        next_dose = current_dose,
        action = "stay",
        mu_hat = mu_hat,
        A = A,
        B = B,
        denom = denom,
        ex = A,
        n_eff = denom,
        dose_data = dose_data
      ))
    }
  }
  
  if (mu_hat >= lambda_d) {
    return(list(
      next_dose = max(current_dose - 1L, 1L),
      action = "de-escalate",
      mu_hat = mu_hat,
      A = A,
      B = B,
      denom = denom,
      ex = A,
      n_eff = denom,
      dose_data = dose_data
    ))
  }
  
  list(
    next_dose = current_dose,
    action = "stay",
    mu_hat = mu_hat,
    A = A,
    B = B,
    denom = denom,
    ex = A,
    n_eff = denom,
    dose_data = dose_data
  )
}
select_mtd_weighted <- function(patient,
                                target,
                                cutoff.eli = 0.95,
                                discount_alpha = 0,
                                falke_alpha = 0.5) {
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1) return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- 1:n
    repeat {
      viol <- diff(x) < 0
      if (!any(viol)) break
      i <- which(viol)[1]
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      idx <- lvlsets %in% c(lvl1, lvl2)
      x[idx] <- sum(x[idx] * wt[idx]) / sum(wt[idx])
      lvlsets[idx] <- lvl1
    }
    x
  }
  
  if (nrow(patient) == 0L) return(NA_integer_)
  
  dat <- patient[order(patient$id, patient$cycle, patient$arrival_time), , drop = FALSE]
  nd <- max(dat$dose)
  
  ## elimination still based on raw binary DLT counts
  n <- tabulate(dat$dose, nbins = nd)
  y_bin <- tabulate(dat$dose[dat$y == 1L], nbins = nd)
  
  elimi <- rep(0L, nd)
  for (i in seq_len(nd)) {
    if (n[i] > 2L) {
      if (1 - pbeta(target, y_bin[i] + 1, n[i] - y_bin[i] + 1) > cutoff.eli) {
        elimi[i:nd] <- 1L
        break
      }
    }
  }
  
  if (elimi[1] == 1L) return(99L)
  if (!any(n > 0L)) return(NA_integer_)
  
  ## full-data alpha_ij evaluated at observed outcome time
  idx_all <- seq_len(nrow(dat))
  is_ipde <- (dat$ipde_true == 1L) & (dat$cycle > 1L)
  
  time_ref <- ifelse(
    (dat$y == 1L) & is.finite(dat$DLT_time),
    dat$DLT_time,
    dat$eval_time
  )
  
  alpha_ij <- rep(1, nrow(dat))
  ipde_rows <- which(is_ipde)
  if (length(ipde_rows) > 0L) {
    alpha_ij[ipde_rows] <- compute_alpha_ij(
      patient = dat,
      row_index = idx_all[ipde_rows],
      t_ref = time_ref[ipde_rows],
      window = max(dat$eval_time - dat$arrival_time, na.rm = TRUE),
      falke_alpha = falke_alpha
    )
  }
  
  A_contrib <- as.numeric(dat$y == 1L)
  B_contrib <- ifelse(dat$y == 0L, ifelse(is_ipde, alpha_ij, 1), 0)
  den_contrib <- A_contrib + B_contrib
  
  A_sum <- numeric(nd)
  D_sum <- numeric(nd)
  
  tmpA <- tapply(A_contrib, dat$dose, sum)
  if (length(tmpA) > 0L) A_sum[as.integer(names(tmpA))] <- as.numeric(tmpA)
  
  tmpD <- tapply(den_contrib, dat$dose, sum)
  if (length(tmpD) > 0L) D_sum[as.integer(names(tmpD))] <- as.numeric(tmpD)
  
  nadmis <- min(max(which(elimi == 0L)), max(which(n > 0L)))
  mu_hat <- A_sum[1:nadmis] / pmax(D_sum[1:nadmis], 1e-8)
  mu_var <- pmax(mu_hat * (1 - mu_hat) / pmax(D_sum[1:nadmis], 1e-8), 1e-8)
  
  mu_iso <- pava(mu_hat, wt = 1 / mu_var)
  mu_iso <- mu_iso + seq_len(nadmis) * 1e-10
  
  which.min(abs(mu_iso - target))
}
# boin_move_gboin_tite <- function(current_dose,
#                                  patient,
#                                  t_now,
#                                  ndose,
#                                  window,
#                                  lambda_e,
#                                  lambda_d,
#                                  elimi,
#                                  discount_alpha = 0) {
#   if (discount_alpha < 0) {
#     stop("discount_alpha must be >= 0")
#   }
#   
#   admissible <- which(elimi == 0L)
#   if (length(admissible) == 0L) {
#     return(list(
#       next_dose = NA_integer_,
#       action = "stop",
#       mu_hat = NA_real_,
#       ex = NA_real_,
#       n_eff = NA_real_,
#       dose_data = NULL
#     ))
#   }
#   
#   max_admissible <- max(admissible)
#   
#   dose_data <- collect_gboin_dose_data(
#     patient = patient,
#     dose_level = current_dose,
#     t_now = t_now,
#     window = window,
#     discount_alpha = discount_alpha
#   )
#   
#   if (nrow(dose_data) == 0L) {
#     return(list(
#       next_dose = current_dose,
#       action = "stay",
#       mu_hat = NA_real_,
#       ex = 0,
#       n_eff = 0,
#       dose_data = dose_data
#     ))
#   }
#   
#   ## BOIN-style moment estimator:
#   ##   p_hat = (y_N + y_I) / (n_tilde_N + alpha_ipde * n_tilde_I)
#   ex    <- sum(dose_data$ex_contrib, na.rm = TRUE)
#   n_eff <- sum(dose_data$n_eff_contrib, na.rm = TRUE)
#   
#   if (n_eff <= 0) {
#     return(list(
#       next_dose = current_dose,
#       action = "stay",
#       mu_hat = NA_real_,
#       ex = ex,
#       n_eff = n_eff,
#       dose_data = dose_data
#     ))
#   }
#   
#   mu_hat <- ex / n_eff
#   
#   if (mu_hat <= lambda_e) {
#     cand <- current_dose + 1L
#     if (cand <= max_admissible) {
#       return(list(
#         next_dose = cand,
#         action = "escalate",
#         mu_hat = mu_hat,
#         ex = ex,
#         n_eff = n_eff,
#         dose_data = dose_data
#       ))
#     } else {
#       return(list(
#         next_dose = current_dose,
#         action = "stay",
#         mu_hat = mu_hat,
#         ex = ex,
#         n_eff = n_eff,
#         dose_data = dose_data
#       ))
#     }
#   }
#   
#   if (mu_hat >= lambda_d) {
#     return(list(
#       next_dose = max(current_dose - 1L, 1L),
#       action = "de-escalate",
#       mu_hat = mu_hat,
#       ex = ex,
#       n_eff = n_eff,
#       dose_data = dose_data
#     ))
#   }
#   
#   list(
#     next_dose = current_dose,
#     action = "stay",
#     mu_hat = mu_hat,
#     ex = ex,
#     n_eff = n_eff,
#     dose_data = dose_data
#   )
# }


collect_gboin_dose_data <- function(patient,
                                    dose_level,
                                    t_now,
                                    window,
                                    discount_alpha = 0) {
  if (discount_alpha < 0) {
    stop("discount_alpha must be >= 0")
  }
  
  ## keep old argument name for compatibility:
  ## alpha_ipde = 1 means no extra IPDE toxicity
  ## alpha_ipde > 1 inflates IPDE effective sample size
  alpha_ipde <- 1 + discount_alpha
  
  idx <- which(patient$dose == dose_level & patient$arrival_time <= t_now)
  if (length(idx) == 0L) {
    return(data.frame(
      id = integer(0),
      cycle = integer(0),
      dose = integer(0),
      arrival_time = numeric(0),
      eval_time = numeric(0),
      DLT_time = numeric(0),
      ipde_true = integer(0),
      y = integer(0),
      prev_dose = numeric(0),
      time_follow = numeric(0),
      time_weight = numeric(0),
      dose_weight = numeric(0),
      observed = logical(0),
      dlt_observed = logical(0),
      no_dlt_complete = logical(0),
      pending = logical(0),
      et = numeric(0),
      ex_contrib = numeric(0),
      m_contrib = numeric(0),
      u_contrib = numeric(0),
      n_tilde_contrib = numeric(0),
      n_eff_contrib = numeric(0)
    ))
  }
  
  dat <- patient[idx, , drop = FALSE]
  
  ## previous dose = previous cycle dose from same patient
  key_all  <- paste(patient$id, patient$cycle)
  prev_key <- paste(dat$id, dat$cycle - 1L)
  prev_idx <- match(prev_key, key_all)
  prev_dose <- patient$dose[prev_idx]
  prev_dose[dat$cycle <= 1L | is.na(prev_idx)] <- NA_real_
  
  ## TITE follow-up weight
  if (window > 0) {
    time_follow <- pmax(0, pmin(t_now - dat$arrival_time, window))
    time_weight <- pmin(time_follow / window, 1)
  } else {
    time_follow <- rep(0, nrow(dat))
    time_weight <- rep(1, nrow(dat))
  }
  
  ## identify IPDE administrations
  is_ipde <- (dat$ipde_true == 1L) & (dat$cycle > 1L) & !is.na(prev_dose)
  
  ## observed status:
  ## - DLT already occurred by t_now, OR
  ## - completed assessment with no DLT by t_now
  dlt_observed    <- (dat$y == 1L) & is.finite(dat$DLT_time) & (dat$DLT_time <= t_now)
  no_dlt_complete <- (dat$y == 0L) & (dat$eval_time <= t_now)
  observed        <- dlt_observed | no_dlt_complete
  pending         <- !observed
  
  ## simple BOIN approximation
  ## numerator: raw observed DLT count
  ex_contrib <- as.numeric(dlt_observed)
  
  ## denominator decomposition:
  ##   n_tilde = y + m + u
  m_contrib <- as.numeric(no_dlt_complete)
  u_contrib <- ifelse(pending, time_weight, 0)
  n_tilde_contrib <- ex_contrib + m_contrib + u_contrib
  
  ## IPDE rows get alpha_ipde multiplier in denominator only
  alpha_weight <- ifelse(is_ipde, alpha_ipde, 1)
  n_eff_contrib <- alpha_weight * n_tilde_contrib
  
  out <- dat
  out$prev_dose        <- prev_dose
  out$time_follow      <- time_follow
  out$time_weight      <- time_weight
  out$dose_weight      <- alpha_weight   # kept name for backward compatibility
  out$observed         <- observed
  out$dlt_observed     <- dlt_observed
  out$no_dlt_complete  <- no_dlt_complete
  out$pending          <- pending
  out$et               <- ex_contrib     # backward-compatible alias
  out$ex_contrib       <- ex_contrib
  out$m_contrib        <- m_contrib
  out$u_contrib        <- u_contrib
  out$n_tilde_contrib  <- n_tilde_contrib
  out$n_eff_contrib    <- n_eff_contrib
  
  out
}
global_elimination <- function(y, dv, t.enter, t.event, ndose,
                               target, cutoff.eli, t_now) {
  elimi <- rep(0L, ndose)
  
  for (dd in seq_len(ndose)) {
    keep <- (dv == dd) & (t.enter <= t_now)
    n_dd <- sum(keep)
    
    if (n_dd >= 3L) {
      enter_dd <- t.enter[keep]
      event_dd <- t.event[keep]
      y_dd <- y[keep]
      
      observed_dd <- (enter_dd + event_dd) <= t_now
      s_dd <- sum(y_dd[observed_dd] == 1)
      
      post_overtox <- 1 - pbeta(target, s_dd + 1, n_dd - s_dd + 1)
      if (post_overtox > cutoff.eli) {
        elimi[dd:ndose] <- 1L
        break
      }
    }
  }
  
  elimi
}

boin_tite_move <- function(current_dose,
                           ndose,
                           n,
                           s,
                           pending_follow,
                           window,
                           target,
                           lambda_e,
                           lambda_d,
                           elimi,
                           maxpen = 0.5,
                           alpha_prior = 0.5 * target,
                           beta_prior  = 1 - 0.5 * target) {
  admissible <- which(elimi == 0L)
  if (length(admissible) == 0L) {
    return(list(next_dose = NA_integer_, action = "stop"))
  }
  
  max_admissible <- max(admissible)
  
  safe_escalate <- function(dose_now, max_ok) {
    cand <- dose_now + 1L
    if (cand <= max_ok) {
      list(next_dose = cand, action = "escalate")
    } else {
      list(next_dose = dose_now, action = "stay")
    }
  }
  
  if (n == 0L) {
    return(list(next_dose = current_dose, action = "stay"))
  }
  
  pending_n <- length(pending_follow)
  n_eval <- n - pending_n
  
  ## gating rule for TITE branch
  if (n_eval < maxpen * n) {
    return(list(next_dose = current_dose, action = "stay"))
  }
  
  ## no pending -> ordinary BOIN
  if (pending_n == 0L) {
    phat <- s / n
    
    if (phat <= lambda_e) {
      return(safe_escalate(current_dose, max_admissible))
    } else if (phat >= lambda_d) {
      return(list(next_dose = max(current_dose - 1L, 1L), action = "de-escalate"))
    } else {
      return(list(next_dose = current_dose, action = "stay"))
    }
  }
  
  ## TITE-BOIN STFT rule
  STFT <- if (window > 0) sum(pmin(pending_follow, window) / window) else 0
  r <- n_eval
  
  pbar <- (s + alpha_prior) / (r + alpha_prior + beta_prior)
  obs_rate <- s / n
  
  pi_E <- if (obs_rate < target) {
    pending_n - ((1 - pbar) / pbar) * (n * lambda_e - s)
  } else {
    Inf
  }
  # browser()
  pi_D <- if (obs_rate > target) {
    pending_n - ((1 - pbar) / pbar) * (n * lambda_d - s)
  } else {
    -Inf
  }
  
  if (STFT >= pi_E) {
    return(safe_escalate(current_dose, max_admissible))
  } else if (STFT <= pi_D) {
    return(list(next_dose = max(current_dose - 1L, 1L), action = "de-escalate"))
  } else {
    return(list(next_dose = current_dose, action = "stay"))
  }
}
# ============================================================
# Run many trials using simulate_IPCRM_trial()
# Returns objects analogous to get_oc_tite_boin_benchmark_clock()
# ============================================================

simulate_IPCRM_trial <- function(
    PI, PI_ipde,
    J = length(PI),
    COHORTSIZE = 3,
    ncohort = 10,
    Kmax = 3,
    TARGET = 0.30,
    window = 28,
    arrival_rate = 1/14,
    seed = 1,
    verbose = FALSE,
    model = "IPCRM",
    parameters = NULL,
    escalation_rule_ipde = 1,
    escalation_rule_new = 1,
    TITE = FALSE,
    gBOIN = FALSE,
    discount_alpha = 0,
    maxpen = 0.5,
    cutoff.eli = 0.95,
    design = 1,
    dose_cap = 3L,
    dose_cap1 = 0
) {
  set.seed(seed)
  
  if (!design %in% c(1, 2, 3)) {
    stop("design must be either 1, 2, or 3")
  }
  if (discount_alpha < 0) {
    stop("discount_alpha must be >= 0")
  }
  
  bd <- boin_boundary(TARGET)
  lambda_e <- bd$lambda_e
  lambda_d <- bd$lambda_d
  
  ## ----------------------------------------------------------
  ## ordinary BOIN move for non-TITE branch
  ## ----------------------------------------------------------
  boin_move_non_tite <- function(current_dose, n_eval, s, ndose, elimi,
                                 lambda_e, lambda_d) {
    admissible <- which(elimi == 0L)
    if (length(admissible) == 0L) {
      return(list(next_dose = NA_integer_, action = "stop"))
    }
    
    if (n_eval <= 0L) {
      return(list(next_dose = current_dose, action = "stay"))
    }
    
    phat <- s / n_eval
    
    if (phat <= lambda_e) {
      cand <- current_dose + 1L
      if (cand <= ndose && elimi[cand] == 0L) {
        return(list(next_dose = cand, action = "escalate"))
      } else {
        return(list(next_dose = current_dose, action = "stay"))
      }
    } else if (phat >= lambda_d) {
      return(list(next_dose = max(current_dose - 1L, 1L), action = "de-escalate"))
    } else {
      return(list(next_dose = current_dose, action = "stay"))
    }
  }
  
  ## ----------------------------------------------------------
  ## collect weighted patient data at one dose and one time
  ##
  ## dose weight:
  ##   regular administration: 1
  ##   IPDE administration: d_j / (d_j + alpha * d_{j-1})
  ##
  ## quasi-Bernoulli observed ET:
  ##   regular observed DLT -> 1
  ##   IPDE observed DLT    -> dose_weight
  ##   otherwise            -> 0
  ##
  ## effective denominator contribution:
  ##   completed no-DLT -> 1
  ##   pending row      -> time_weight = follow-up / window
  ##   observed DLT     -> 0 in m-part, ET in ex-part
  ## ----------------------------------------------------------
  collect_gboin_dose_data <- function(patient, dose_level, t_now, window,
                                      discount_alpha = 0) {
    idx <- which(patient$dose == dose_level & patient$arrival_time <= t_now)
    
    if (length(idx) == 0L) {
      return(data.frame(
        id = integer(0),
        cycle = integer(0),
        dose = integer(0),
        arrival_time = numeric(0),
        eval_time = numeric(0),
        DLT_time = numeric(0),
        ipde_true = integer(0),
        ipde_ok = integer(0),
        y = integer(0),
        prev_dose = numeric(0),
        time_follow = numeric(0),
        time_weight = numeric(0),
        dose_weight = numeric(0),
        observed = logical(0),
        et = numeric(0),
        ex_contrib = numeric(0),
        m_contrib = numeric(0),
        n_eff_contrib = numeric(0)
      ))
    }
    
    dat <- patient[idx, , drop = FALSE]
    
    ## previous dose for each administration = previous cycle dose of same patient
    key_all  <- paste(patient$id, patient$cycle)
    prev_key <- paste(dat$id, dat$cycle - 1L)
    prev_idx <- match(prev_key, key_all)
    
    prev_dose <- rep(NA_real_, nrow(dat))
    has_prev <- !is.na(prev_idx) & dat$cycle > 1L
    prev_dose[has_prev] <- patient$dose[prev_idx[has_prev]]
    
    ## evaluated follow-up up to t_now
    if (window > 0) {
      eval_cutoff <- pmin(t_now, dat$eval_time)
      time_follow <- pmax(0, eval_cutoff - dat$arrival_time)
      time_follow <- pmin(time_follow, window)
      time_weight <- pmin(time_follow / window, 1)
    } else {
      time_follow <- rep(0, nrow(dat))
      time_weight <- rep(1, nrow(dat))
    }
    
    ## dose discount factor
    dose_weight <- rep(1, nrow(dat))
    is_ipde <- (dat$ipde_true == 1L) & !is.na(prev_dose)
    
    dose_weight[is_ipde] <- dat$dose[is_ipde] /
      (dat$dose[is_ipde] + discount_alpha * prev_dose[is_ipde])
    
    dose_weight <- pmin(pmax(dose_weight, 0), 1)
    
    ## observed means DLT status fully known by t_now
    observed <- dat$eval_time <= t_now
    
    ## weighted ET score
    dlt_observed <- observed & (dat$y == 1L) & (dat$DLT_time <= t_now)
    et <- ifelse(dlt_observed, dose_weight, 0)
    
    ## TITE-gBOIN pieces
    completed_no_dlt <- observed & (dat$y == 0L)
    pending <- !observed
    
    ex_contrib <- ifelse(dlt_observed, et, 0)
    m_contrib <- ifelse(completed_no_dlt, 1, ifelse(pending, time_weight, 0))
    n_eff_contrib <- ex_contrib + m_contrib
    
    out <- dat
    out$prev_dose <- prev_dose
    out$time_follow <- time_follow
    out$time_weight <- time_weight
    out$dose_weight <- dose_weight
    out$observed <- observed
    out$et <- et
    out$ex_contrib <- ex_contrib
    out$m_contrib <- m_contrib
    out$n_eff_contrib <- n_eff_contrib
    
    out
  }
  
  Nmax <- COHORTSIZE * ncohort
  ndose <- length(PI)
  
  patient <- data.frame(
    id = integer(0),
    cycle = integer(0),
    dose = integer(0),
    arrival_time = numeric(0),
    eval_time = numeric(0),
    DLT_time = numeric(0),
    ipde_true = integer(0),
    ipde_ok = integer(0),
    y = integer(0),
    MTD = integer(0)
  )
  
  active <- rep(FALSE, Nmax)
  current_dose <- rep(NA_integer_, Nmax)
  current_cycle <- rep(0L, Nmax)
  stop_reason <- rep(NA_character_, Nmax)
  
  t_last <- 0
  next_pid <- 0L
  j_recent <- 1L
  
  trial_stop <- FALSE
  stop_ipde <- FALSE
  
  for (coh in seq_len(ncohort)) {
    
    interarrival <- rexp(COHORTSIZE, rate = arrival_rate)
    arrival_time <- t_last + cumsum(interarrival)
    t_last <- tail(arrival_time, 1)
    
    for (k in seq_len(COHORTSIZE)) {
      
      t_now <- arrival_time[k]
      
      repeat {
        pending_idx <- which(patient$eval_time <= t_now & patient$ipde_ok == 1L)
        if (length(pending_idx) == 0L) break
        
        r <- pending_idx[which.min(patient$eval_time[pending_idx])]
        t_admin <- patient$eval_time[r]
        pid <- patient$id[r]
        d_cur <- patient$dose[r]
        cyc <- patient$cycle[r]
        
        if (!active[pid]) {
          patient$ipde_ok[r] <- 0L
          next
        }
        
        y_obs <- as.integer(patient$DLT_time[r] <= t_admin)
        
        if (verbose) {
          cat(sprintf(
            "[EVAL] t=%.2f pid=%d cycle=%d dose=%d y_obs=%d j_recent=%d design=%d TITE=%d gBOIN=%d\n",
            t_admin, pid, cyc, d_cur, y_obs, j_recent, design, TITE, gBOIN
          ))
        }
        
        patient$ipde_ok[r] <- 0L
        
        if (y_obs == 1L) {
          active[pid] <- FALSE
          if (is.na(stop_reason[pid])) stop_reason[pid] <- "DLT"
          next
        }
        
        if (cyc >= Kmax) {
          active[pid] <- FALSE
          if (is.na(stop_reason[pid])) stop_reason[pid] <- "maxK"
          next
        }
        
        ## elimination test on all doses BEFORE IPDE
        elimi <- global_elimination(
          y = patient$y,
          dv = patient$dose,
          t.enter = patient$arrival_time,
          t.event = patient$DLT_time - patient$arrival_time,
          ndose = J,
          target = TARGET,
          cutoff.eli = cutoff.eli,
          t_now = t_admin
        )
        
        if (elimi[1] == 1L) {
          trial_stop <- TRUE
          break
        }
        
        admissible <- which(elimi == 0L)
        if (length(admissible) == 0L) {
          trial_stop <- TRUE
          break
        }
        
        patient_current_dose <- current_dose[pid]
        if (is.na(patient_current_dose)) {
          patient_current_dose <- d_cur
        }
        
        if (design == 1) {
          ## design 1: IPDE uses global working dose j_recent
          if (is.na(j_recent) || elimi[j_recent] == 1L) {
            j_work <- max(admissible)
          } else {
            j_work <- j_recent
          }
          
          idx_cur <- which(patient$dose == j_work & patient$arrival_time <= t_admin)
          n_total_cur <- length(idx_cur)
          n_eval_cur <- sum(patient$eval_time[idx_cur] <= t_admin)
          
          if (TITE) {
            allow_ipde <- (n_total_cur >= dose_cap1) &&
              (n_eval_cur >= maxpen * n_total_cur) &&
              (elimi[j_work] == 0L)
          } else {
            allow_ipde <- (n_eval_cur >= dose_cap1) &&
              (elimi[j_work] == 0L)
          }
          
          if (!allow_ipde) {
            active[pid] <- FALSE
            if (is.na(stop_reason[pid])) stop_reason[pid] <- "no_ipde"
            next
          }
          
          if (nrow(patient) >= Nmax) {
            stop_ipde <- TRUE
            break
          }
          
          ## only perform IPDE if j_work is higher than patient's current dose
          if (j_work <= patient_current_dose) {
            active[pid] <- FALSE
            if (is.na(stop_reason[pid])) stop_reason[pid] <- "no_ipde"
            next
          }
          
          new_dose <- j_work
          new_cycle <- cyc + 1L
          
        } else {
          ## design 2 and 3: IPDE uses patient's own current dose
          j_ipde <- patient_current_dose
          
          idx_cur <- which(patient$dose == j_ipde & patient$arrival_time <= t_admin)
          n_total_cur <- length(idx_cur)
          n_eval_cur <- sum(patient$eval_time[idx_cur] <= t_admin)
          s_cur <- sum(patient$y[idx_cur] == 1L & patient$DLT_time[idx_cur] <= t_admin)
          
          if (TITE) {
            allow_ipde <- (n_total_cur >= dose_cap) &&
              (n_eval_cur >= maxpen * n_total_cur) &&
              (elimi[j_ipde] == 0L)
          } else {
            allow_ipde <- (n_eval_cur >= dose_cap) &&
              (elimi[j_ipde] == 0L)
          }
          
          if (!allow_ipde) {
            active[pid] <- FALSE
            if (is.na(stop_reason[pid])) stop_reason[pid] <- "no_ipde"
            next
          }
          
          if (nrow(patient) >= Nmax) {
            stop_ipde <- TRUE
            break
          }
          
          if (TITE) {
            if (gBOIN) {
              move_ipde <- boin_move_gboin_tite(
                current_dose = j_ipde,
                patient = patient,
                t_now = t_admin,
                ndose = J,
                window = window,
                lambda_e = lambda_e,
                lambda_d = lambda_d,
                elimi = elimi,
                discount_alpha = discount_alpha
              )
            } else {
              idx_pending_cur <- idx_cur[patient$eval_time[idx_cur] > t_admin]
              pending_follow <- if (length(idx_pending_cur) > 0L) {
                pmax(0, t_admin - patient$arrival_time[idx_pending_cur])
              } else {
                numeric(0)
              }
              
              move_ipde <- boin_tite_move(
                current_dose = j_ipde,
                ndose = J,
                n = n_total_cur,
                s = s_cur,
                pending_follow = pending_follow,
                window = window,
                target = TARGET,
                lambda_e = lambda_e,
                lambda_d = lambda_d,
                elimi = elimi,
                maxpen = maxpen
              )
            }
          } else {
            move_ipde <- boin_move_non_tite(
              current_dose = j_ipde,
              n_eval = n_eval_cur,
              s = s_cur,
              ndose = J,
              elimi = elimi,
              lambda_e = lambda_e,
              lambda_d = lambda_d
            )
          }
          
          if (is.na(move_ipde$next_dose)) {
            trial_stop <- TRUE
            break
          }
          
          ## design 2 and 3 only perform IPDE if local decision is escalate
          if (!(identical(move_ipde$action, "escalate") &&
                move_ipde$next_dose > patient_current_dose)) {
            active[pid] <- FALSE
            if (is.na(stop_reason[pid])) stop_reason[pid] <- "no_ipde"
            next
          }
          
          new_dose <- move_ipde$next_dose
          new_cycle <- cyc + 1L
        }
        
        y_new <- rbinom(1L, 1L, PI_ipde[new_dose])
        
        if (window == 0) {
          DLT_new <- if (y_new == 1L) t_admin else Inf
          eval_new <- t_admin
        } else {
          DLT_new <- if (y_new == 1L) t_admin + sample.int(window, 1) else Inf
          eval_new <- if (is.finite(DLT_new)) DLT_new else (t_admin + window)
        }
        
        patient <- rbind(
          patient,
          data.frame(
            id = pid,
            cycle = new_cycle,
            dose = new_dose,
            arrival_time = t_admin,
            eval_time = eval_new,
            DLT_time = DLT_new,
            ipde_true = 1L,
            ipde_ok = 1L,
            y = y_new,
            MTD = j_recent
          )
        )
        
        current_dose[pid] <- new_dose
        current_cycle[pid] <- new_cycle
        
        ## design 2 updates j_recent after IPDE
        if (design == 2) {
          j_recent <- new_dose
        }
        
        if (verbose) {
          cat(sprintf(
            "      -> RECYCLE pid=%d to dose=%d (cycle=%d), new eval=%.2f\n",
            pid, new_dose, new_cycle, eval_new
          ))
        }
      }
      
      if (trial_stop || stop_ipde) break
      
      ## ------------------------------------------------------
      ## new cohort / new patient allocation
      ## ------------------------------------------------------
      if (nrow(patient) == 0L) {
        j_S_curr <- 1L
      } else {
        ## elimination test on all doses BEFORE new patient allocation
        elimi <- global_elimination(
          y = patient$y,
          dv = patient$dose,
          t.enter = patient$arrival_time,
          t.event = patient$DLT_time - patient$arrival_time,
          ndose = J,
          target = TARGET,
          cutoff.eli = cutoff.eli,
          t_now = t_now
        )
        
        if (elimi[1] == 1L) {
          trial_stop <- TRUE
          break
        }
        
        admissible <- which(elimi == 0L)
        if (length(admissible) == 0L) {
          trial_stop <- TRUE
          break
        }
        
        ## new-patient decision is run on j_recent
        if (is.na(j_recent) || elimi[j_recent] == 1L) {
          j_work <- max(admissible)
        } else {
          j_work <- j_recent
        }
        
        idx_cur <- which(patient$dose == j_work & patient$arrival_time <= t_now)
        n_total_cur <- length(idx_cur)
        n_eval_cur <- sum(patient$eval_time[idx_cur] <= t_now)
        s_cur <- sum(patient$y[idx_cur] == 1L & patient$DLT_time[idx_cur] <= t_now)
        
        if (TITE) {
          allow_escalate <- (n_total_cur >= dose_cap) &&
            (n_eval_cur >= maxpen * n_total_cur)
          
          if (allow_escalate) {
            if (gBOIN) {
              move <- boin_move_gboin_tite(
                current_dose = j_work,
                patient = patient,
                t_now = t_now,
                ndose = J,
                window = window,
                lambda_e = lambda_e,
                lambda_d = lambda_d,
                elimi = elimi,
                discount_alpha = discount_alpha
              )
            } else {
              idx_pending_cur <- idx_cur[patient$eval_time[idx_cur] > t_now]
              pending_follow <- if (length(idx_pending_cur) > 0L) {
                pmax(0, t_now - patient$arrival_time[idx_pending_cur])
              } else {
                numeric(0)
              }
              
              move <- boin_tite_move(
                current_dose = j_work,
                ndose = J,
                n = n_total_cur,
                s = s_cur,
                pending_follow = pending_follow,
                window = window,
                target = TARGET,
                lambda_e = lambda_e,
                lambda_d = lambda_d,
                elimi = elimi,
                maxpen = maxpen
              )
            }
            
            if (is.na(move$next_dose)) {
              trial_stop <- TRUE
              break
            } else {
              j_S_curr <- move$next_dose
            }
          } else {
            j_S_curr <- j_work
          }
          
        } else {
          if (n_eval_cur >= dose_cap) {
            move <- boin_move_non_tite(
              current_dose = j_work,
              n_eval = n_eval_cur,
              s = s_cur,
              ndose = J,
              elimi = elimi,
              lambda_e = lambda_e,
              lambda_d = lambda_d
            )
            
            if (is.na(move$next_dose)) {
              trial_stop <- TRUE
              break
            } else {
              j_S_curr <- move$next_dose
            }
          } else {
            j_S_curr <- j_work
          }
        }
      }
      
      if (nrow(patient) >= Nmax) break
      
      next_pid <- next_pid + 1L
      pid_new <- next_pid
      
      active[pid_new] <- TRUE
      current_dose[pid_new] <- j_S_curr
      current_cycle[pid_new] <- 1L
      
      y0 <- rbinom(1L, 1L, PI[j_S_curr])
      
      if (window == 0) {
        DLT0 <- if (y0 == 1L) t_now else Inf
        eval0 <- t_now
      } else {
        DLT0 <- if (y0 == 1L) t_now + sample.int(window, 1) else Inf
        eval0 <- if (is.finite(DLT0)) DLT0 else (t_now + window)
      }
      
      patient <- rbind(
        patient,
        data.frame(
          id = pid_new,
          cycle = 1L,
          dose = j_S_curr,
          arrival_time = t_now,
          eval_time = eval0,
          DLT_time = DLT0,
          ipde_true = 0L,
          ipde_ok = 1L,
          y = y0,
          MTD = j_recent
        )
      )
      
      ## most recent assigned dose updated after each new-patient assignment
      j_recent <- j_S_curr
      
      if (verbose) {
        cat(sprintf(
          "[ARR] t=%.2f pid=%d startdose=%d eval=%.2f y(latent)=%d\n",
          t_now, pid_new, j_S_curr, eval0, y0
        ))
      }
      
      if (nrow(patient) >= Nmax) break
    }
    
    if (trial_stop || stop_ipde || nrow(patient) >= Nmax) break
  }
  
  ## ----------------------------------------------------------
  ## final summary
  ## ----------------------------------------------------------
  if (nrow(patient) == 0L) {
    t_final <- 0
    final_MTD <- NA_integer_
    posttox <- rep(NA_real_, ndose)
  } else {
    t_start <- min(patient$arrival_time[patient$ipde_true == 0L])
    t_end <- max(patient$eval_time)
    t_final <- t_end - t_start
    
    if (trial_stop) {
      final_MTD <- 99L
      posttox <- rep(NA_real_, ndose)
    } else {
      y_end <- tabulate(patient$dose[patient$y == 1L], nbins = ndose)
      n_end <- tabulate(patient$dose, nbins = ndose)
      
      if (gBOIN) {
        final_MTD <- select_mtd_weighted(
          patient = patient,
          target = TARGET,
          cutoff.eli = cutoff.eli,
          discount_alpha = discount_alpha
        )
      } else {
        final_MTD <- select_mtd_binary(
          target = TARGET,
          y = y_end,
          n = n_end,
          cutoff.eli = cutoff.eli
        )
      }
      
      posttox <- rep(NA_real_, ndose)
      nz <- which(n_end > 0L)
      posttox[nz] <- y_end[nz] / n_end[nz]
    }
  }
  
  list(
    patient = patient,
    stop_reason = if (next_pid > 0L) stop_reason[1:next_pid] else character(0),
    final_MTD = final_MTD,
    posttox = posttox,
    trial_time = t_final,
    trial_stop = trial_stop
  )
}
get_oc_sim_IPDE <- function(
    target,
    p.true,
    p.true_ipde = p.true,
    ntrial = 1000,
    seed = 1,
    ...
) {
  ndose <- length(p.true)
  
  sel_count <- integer(ndose)
  stop_count <- 0L
  na_count <- 0L
  
  # average number of treatment administrations at each dose
  n_by_dose <- numeric(ndose)
  
  # average number of IPDE administrations at each dose
  nipde_by_dose <- numeric(ndose)
  
  # average number of unique patients who ever received each dose
  unique_n_by_dose <- numeric(ndose)
  
  total_admin <- numeric(ntrial)
  total_unique <- numeric(ntrial)
  duration <- numeric(ntrial)
  
  raw_trials <- vector("list", ntrial)
  
  for (tt in seq_len(ntrial)) {
    fit <- simulate_IPCRM_trial(
      PI = p.true,
      PI_ipde = p.true_ipde,
      TARGET = target,
      seed = seed + tt - 1,
      ...
    )
    
    raw_trials[[tt]] <- fit
    pat <- fit$patient
    
    if (nrow(pat) > 0) {
      # total administrations at each dose
      n_by_dose <- n_by_dose + tabulate(pat$dose, nbins = ndose)
      
      # IPDE administrations at each dose
      if (any(pat$ipde_true == 1L)) {
        nipde_by_dose <- nipde_by_dose + tabulate(pat$dose[pat$ipde_true == 1L], nbins = ndose)
      }
      
      # unique patients who ever received each dose
      id_dose <- unique(pat[, c("id", "dose"), drop = FALSE])
      unique_n_by_dose <- unique_n_by_dose + tabulate(id_dose$dose, nbins = ndose)
      
      total_admin[tt] <- nrow(pat)
      total_unique[tt] <- length(unique(pat$id))
    } else {
      total_admin[tt] <- 0
      total_unique[tt] <- 0
    }
    
    duration[tt] <- fit$trial_time
    
    if (isTRUE(fit$trial_stop) || identical(fit$final_MTD, 99L)) {
      stop_count <- stop_count + 1L
    } else if (!is.na(fit$final_MTD) &&
               fit$final_MTD >= 1L &&
               fit$final_MTD <= ndose) {
      sel_count[fit$final_MTD] <- sel_count[fit$final_MTD] + 1L
    } else {
      na_count <- na_count + 1L
    }
  }
  
  mean_admin <- mean(total_admin)
  mean_unique <- mean(total_unique)
  
  out <- list(
    target = target,
    p.true = p.true,
    p.true_ipde = p.true_ipde,
    
    selpercent = 100 * sel_count / ntrial,
    npatients = n_by_dose / ntrial,
    nipde = nipde_by_dose / ntrial,
    totaln = mean_admin,
    percentstop = 100 * stop_count / ntrial,
    duration = mean(duration),
    
    nuniquepatients = unique_n_by_dose / ntrial,
    total_unique = mean_unique,
    percentna = 100 * na_count / ntrial,
    
    raw = raw_trials
  )
  
  out
}
# ============================================================
# Print table like picture 2 for get_oc_sim_IPDE()
# By default, % treated uses treatment administrations
# Set use_unique = TRUE to use unique patients instead
# ============================================================
print_oc_table_counts <- function(res,
                                  p.true,
                                  scenario_id,
                                  digits = 2,
                                  design = NULL,
                                  method_label = NULL,
                                  alpha_true = NULL,
                                  discount_alpha = NULL) {
  ndose <- length(p.true)
  dose_names <- paste0("Dose ", seq_len(ndose))
  
  tab <- data.frame(
    Metric = c("True DLT rate", "Selection %", "# Pts Treated", "# IPDE Doses"),
    check.names = FALSE
  )
  
  for (j in seq_len(ndose)) {
    tab[[dose_names[j]]] <- c(
      p.true[j],
      res$selpercent[j],
      res$npatients[j],
      res$nipde[j]
    )
  }
  
  tab[["# Patients"]] <- c(NA, NA, res$total_unique, NA)
  tab[["% Early Stopping"]] <- c(NA, res$percentstop, NA, NA)
  tab[["Duration"]] <- c(NA, res$duration, NA, NA)
  
  num_cols <- setdiff(names(tab), "Metric")
  tab[num_cols] <- lapply(tab[num_cols], function(x) {
    ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
  })
  
  cat("\n==============================\n")
  cat("Scenario", scenario_id, "\n")
  if (!is.null(design)) cat("Design:", design, "\n")
  if (!is.null(method_label)) cat("Method:", method_label, "\n")
  if (!is.null(alpha_true)) cat("True IPDE alpha:", alpha_true, "\n")
  if (!is.null(discount_alpha)) cat("Discount alpha:", discount_alpha, "\n")
  cat("==============================\n")
  print(tab, row.names = FALSE, right = TRUE)
}
