## ============================================================
## AIDE-CRM TITE helper and simulation wrapper
##
## This file layers a TITE event-clock simulation over the updated
## likelihoods in AIDE_CRM_helper_modified.R.
##
## Cohort assignment rule:
##   - eligible IPDE/retreat patients are used before new arrivals
##   - a recycled patient must have received their most recent dose at the
##     immediately lower dose level
##   - remaining events wait until a cohort slot opens
## ============================================================

if (!exists("crm_fit", mode = "function") ||
    !exists("crm_move", mode = "function") ||
    !exists("select.mtd.crm", mode = "function")) {
  source("AIDE_CRM_helper_modified.R")
}

aide_crm_tite_draw_dlt <- function(pi_val,
                                   t_start,
                                   window = 28,
                                   dist = 1,
                                   alpha = 0.5) {
  pi_val <- as.numeric(pi_val)
  if (!is.finite(pi_val) || pi_val < 0 || pi_val > 1) {
    stop("pi_val must be a probability in [0, 1].")
  }

  y <- stats::rbinom(1L, 1L, pi_val)
  if (window <= 0) {
    return(list(
      y = as.integer(y),
      t_tox = if (y == 1L) 0 else Inf,
      t_eval = as.numeric(t_start)
    ))
  }

  if (y == 0L) {
    return(list(
      y = 0L,
      t_tox = Inf,
      t_eval = as.numeric(t_start + window)
    ))
  }

  tox_delay <- switch(
    as.character(dist),
    "2" = {
      pihalf <- min(max(alpha * pi_val, 1e-8), pi_val * 0.999)
      shape <- log(log(1 - pi_val) / log(1 - pihalf)) / log(2)
      rate <- -log(1 - pi_val) / (window^shape)
      u <- stats::runif(1L, 0, pi_val)
      (-log(1 - u) / rate)^(1 / shape)
    },
    "3" = {
      pihalf <- min(max(alpha * pi_val, 1e-8), pi_val * 0.999)
      shape <- log((1 / (1 - pi_val) - 1) /
                     (1 / (1 - pihalf) - 1)) / log(2)
      rate <- (1 / (1 - pi_val) - 1) / (window^shape)
      u <- stats::runif(1L, 0, pi_val)
      (u / (rate * (1 - u)))^(1 / shape)
    },
    stats::runif(1L, 0, window)
  )

  tox_delay <- min(max(tox_delay, 0), window)
  list(
    y = 1L,
    t_tox = as.numeric(tox_delay),
    t_eval = as.numeric(t_start + tox_delay)
  )
}

aide_crm_tite_make_decision_dat <- function(admin,
                                            t_now,
                                            window = 28) {
  if (is.null(admin) || nrow(admin) == 0L) {
    return(admin)
  }

  dat <- admin[admin$t_start <= t_now, , drop = FALSE]
  if (nrow(dat) == 0L) {
    return(dat)
  }

  tox_time <- dat$t_start + dat$t_tox
  observed_dlt <- dat$y == 1L & tox_time <= t_now
  fully_evaluated <- dat$t_eval <= t_now
  observed <- observed_dlt | fully_evaluated

  followup <- if (window > 0) {
    pmin(window, pmax(0, t_now - dat$t_start))
  } else {
    rep(0, nrow(dat))
  }
  tite_weight <- if (window > 0) followup / window else rep(1, nrow(dat))
  tite_weight[observed] <- 1
  tite_weight[observed_dlt] <- 1
  tite_weight <- pmin(1, pmax(0, tite_weight))

  out <- data.frame(
    id = dat$id,
    dose = dat$dose,
    y = as.integer(observed_dlt),
    type = dat$type,
    t_arrival = dat$t_arrival,
    t_start = dat$t_start,
    t_tox = dat$t_tox,
    t_eval = dat$t_eval,
    ncycle = dat$ncycle,
    cycle = dat$ncycle,
    cohort = dat$cohort,
    tite_weight = tite_weight,
    followup_time = followup,
    observed = observed,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

simulate_AIDE_CRM_TITE_design <- function(
    p_true,
    p_ipde = p_true,
    target = 0.30,
    N_pat = 30L,
    Nmax_eff = 30L,
    cohortsize = 3L,
    T_assess = 28,
    cycle_max = 2L,
    arrival_rate = 0.2,
    t0 = 0,
    startdose = 1L,
    cutoff = 0.95,
    dose_cap = 3L,
    dlt_dist = 1,
    dlt_alpha = 0.5,
    seed = NULL,
    verbose = FALSE,

    crm_r_model = c("fixed", "r_fixed", "random", "alpha_crm", "cumu_crm", "ipcrm"),
    crm_skeleton,
    crm_model_file = NULL,
    crm_fixed_model_file = "fix_CRM_TITE.bug",
    crm_random_model_file = "random_CRM_TITE.bug",
    crm_level_model_file = "random_CRM_level.bug",
    crm_alpha_sd = 2,
    r_carry = 0.1,
    crm_a_r = 1,
    crm_b_r = 9,

    crm_dose_values = NULL,
    crm_time_col = "t_start",
    crm_alpha_grid = seq(0.01, 0.99, length.out = 61),
    crm_alpha_T = T_assess,
    crm_theta_prior_mean = 0,
    crm_theta_prior_sd = sqrt(2),
    crm_alpha_L = 8,
    crm_alpha_rel_tol = 1e-8,
    crm_alpha_eps = 1e-12,
    crm_alpha_n_draw_prior = 5000,

    crm_dose_scores = NULL,
    crm_cumu_model_file = "cumu_CRM_TITE.bug",
    crm_cumu_beta0_mean = stats::qlogis(target),
    crm_cumu_beta0_prec = 4,
    crm_cumu_beta0_df = 1,
    crm_cumu_beta1_shape = 5.83,
    crm_cumu_beta1_rate = 1.21,
    crm_cumu_beta2_rate = 1,
    crm_cumu_include_current = FALSE,

    crm_n_chains = 2,
    crm_n_adapt = 500,
    crm_n_burnin = 500,
    crm_n_iter = 2000,
    crm_thin = 1,
    restrict_to_tried = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  p_true <- as.numeric(p_true)
  p_ipde <- as.numeric(p_ipde)
  ndose <- length(p_true)
  if (length(p_ipde) != ndose) {
    stop("p_ipde must have the same length as p_true.")
  }
  if (length(crm_skeleton) != ndose) {
    stop("crm_skeleton must have the same length as p_true.")
  }
  if (N_pat < cohortsize) {
    stop("Need N_pat >= cohortsize.")
  }
  if (Nmax_eff < cohortsize) {
    stop("Need Nmax_eff >= cohortsize.")
  }
  if (!is.finite(arrival_rate) || arrival_rate <= 0) {
    stop("arrival_rate must be positive.")
  }

  crm_r_model <- match.arg(crm_r_model)
  crm_r_model <- crm_normalize_r_model(crm_r_model)

  arrival_times <- as.numeric(t0 + cumsum(stats::rexp(N_pat, rate = arrival_rate)))

  event_seq <- seq_len(N_pat)
  future_events <- data.frame(
    time = arrival_times,
    type = rep("new", N_pat),
    id = seq_len(N_pat),
    seq = event_seq,
    stringsAsFactors = FALSE
  )
  next_event_seq <- N_pat + 1L

  waiting <- data.frame(
    time = numeric(0),
    type = character(0),
    id = integer(0),
    seq = integer(0),
    stringsAsFactors = FALSE
  )

  admin <- data.frame(
    row_id = integer(0),
    id = integer(0),
    t_arrival = numeric(0),
    t_start = numeric(0),
    t_tox = numeric(0),
    t_eval = numeric(0),
    dose = integer(0),
    y = integer(0),
    ncycle = integer(0),
    cohort = integer(0),
    type = character(0),
    stringsAsFactors = FALSE
  )

  decision_log <- data.frame(
    time = numeric(0),
    cohort = integer(0),
    current_dose = integer(0),
    next_dose = integer(0),
    action = character(0),
    n_admin = integer(0),
    queue = integer(0),
    stringsAsFactors = FALSE
  )

  elimi <- rep(0L, ndose)
  earlystop <- 0L
  trial_stop <- FALSE
  current_dose <- as.integer(startdose)
  current_dose <- max(1L, min(ndose, current_dose))
  cohort_id <- 1L
  cohort_open <- TRUE
  cohort_enrolled <- 0L
  t_now <- t0

  schedule_event <- function(time, type, id) {
    future_events <<- rbind(
      future_events,
      data.frame(
        time = as.numeric(time),
        type = as.character(type),
        id = as.integer(id),
        seq = next_event_seq,
        stringsAsFactors = FALSE
      )
    )
    next_event_seq <<- next_event_seq + 1L
    invisible(NULL)
  }

  pop_due_events <- function(time) {
    if (nrow(future_events) == 0L) {
      return(future_events)
    }
    due <- future_events[future_events$time <= time + 1e-10, , drop = FALSE]
    future_events <<- future_events[future_events$time > time + 1e-10, , drop = FALSE]
    rownames(future_events) <<- NULL
    due[order(due$time, due$seq), , drop = FALSE]
  }

  add_to_waiting <- function(events) {
    if (nrow(events) == 0L) return(invisible(NULL))

    for (ii in seq_len(nrow(events))) {
      ev <- events[ii, , drop = FALSE]

      if (ev$type == "retreat") {
        hist <- admin[admin$id == ev$id, , drop = FALSE]
        if (nrow(hist) == 0L) next
        hist <- hist[order(hist$t_start, hist$row_id), , drop = FALSE]
        last <- hist[nrow(hist), , drop = FALSE]
        if (!(last$y == 0L &&
              last$t_eval <= ev$time + 1e-10 &&
              last$ncycle < cycle_max)) {
          next
        }
      }

      waiting <<- rbind(
        waiting,
        data.frame(
          time = ev$time,
          type = ev$type,
          id = ev$id,
          seq = ev$seq,
          stringsAsFactors = FALSE
        )
      )
    }

    waiting <<- waiting[order(waiting$time, waiting$seq), , drop = FALSE]
    rownames(waiting) <<- NULL
    invisible(NULL)
  }

  recycle_event_is_eligible <- function(ev, dose = current_dose) {
    if (ev$type != "retreat" || dose <= 1L) return(FALSE)

    hist <- admin[admin$id == ev$id, , drop = FALSE]
    if (nrow(hist) == 0L) return(FALSE)

    hist <- hist[order(hist$t_start, hist$row_id), , drop = FALSE]
    last <- hist[nrow(hist), , drop = FALSE]

    isTRUE(
      last$y == 0L &&
        last$t_eval <= ev$time + 1e-10 &&
        last$ncycle < cycle_max &&
        last$dose == dose - 1L
    )
  }

  next_waiting_index <- function() {
    if (nrow(waiting) == 0L) return(integer(0))

    retreat_idx <- which(waiting$type == "retreat")
    if (length(retreat_idx) > 0L) {
      eligible_retreat_idx <- retreat_idx[vapply(
        retreat_idx,
        function(i) recycle_event_is_eligible(waiting[i, , drop = FALSE]),
        logical(1)
      )]
      if (length(eligible_retreat_idx) > 0L) {
        return(eligible_retreat_idx[1L])
      }
    }

    new_idx <- which(waiting$type == "new")
    if (length(new_idx) > 0L) new_idx[1L] else integer(0)
  }

  has_assignable_waiting <- function() {
    length(next_waiting_index()) > 0L
  }

  assign_waiting_one <- function(time) {
    if (!cohort_open || cohort_enrolled >= cohortsize || nrow(waiting) == 0L) {
      return(FALSE)
    }
    if (nrow(admin) >= Nmax_eff) {
      return(FALSE)
    }

    waiting_index <- next_waiting_index()
    if (length(waiting_index) == 0L) return(FALSE)

    ev <- waiting[waiting_index, , drop = FALSE]
    waiting_new <- waiting[-waiting_index, , drop = FALSE]
    rownames(waiting_new) <- NULL
    waiting <<- waiting_new

    row_type <- ev$type
    if (row_type == "new") {
      row_cycle <- 1L
      pi_val <- p_true[current_dose]
    } else {
      if (!recycle_event_is_eligible(ev)) return(FALSE)

      hist <- admin[admin$id == ev$id, , drop = FALSE]
      hist <- hist[order(hist$t_start, hist$row_id), , drop = FALSE]
      last <- hist[nrow(hist), , drop = FALSE]
      row_cycle <- as.integer(last$ncycle + 1L)
      pi_val <- p_ipde[current_dose]
      row_type <- "retreat"
    }

    dlt <- aide_crm_tite_draw_dlt(
      pi_val = pi_val,
      t_start = ev$time,
      window = T_assess,
      dist = dlt_dist,
      alpha = dlt_alpha
    )

    admin <<- rbind(
      admin,
      data.frame(
        row_id = nrow(admin) + 1L,
        id = ev$id,
        t_arrival = ev$time,
        t_start = ev$time,
        t_tox = dlt$t_tox,
        t_eval = dlt$t_eval,
        dose = current_dose,
        y = dlt$y,
        ncycle = row_cycle,
        cohort = cohort_id,
        type = row_type,
        stringsAsFactors = FALSE
      )
    )

    if (dlt$y == 0L && row_cycle < cycle_max) {
      schedule_event(dlt$t_eval, "retreat", ev$id)
    }

    cohort_enrolled <<- cohort_enrolled + 1L
    if (cohort_enrolled >= cohortsize) {
      cohort_open <<- FALSE
    }

    if (verbose) {
      message(
        "t=", round(ev$time, 2),
        " cohort=", cohort_id,
        " dose=", current_dose,
        " type=", row_type,
        " id=", ev$id,
        " slot=", cohort_enrolled, "/", cohortsize
      )
    }

    TRUE
  }

  fill_open_cohort <- function(time) {
    repeat {
      did_assign <- assign_waiting_one(time)
      if (!did_assign) break
      if (!cohort_open || nrow(admin) >= Nmax_eff) break
    }
    invisible(NULL)
  }

  open_next_cohort <- function(time) {
    if (nrow(admin) == 0L) {
      cohort_open <<- TRUE
      cohort_enrolled <<- 0L
      return(invisible(NULL))
    }

    dat_dec <- aide_crm_tite_make_decision_dat(admin, time, T_assess)
    n_trt_curr <- sum(admin$dose == current_dose & admin$t_start <= time)

    move <- crm_move(
      current_dose = current_dose,
      ndose = ndose,
      dat = dat_dec,
      target = target,
      cutoff = cutoff,
      cutoff.eli = cutoff,
      skeleton = crm_skeleton,
      r_model = crm_r_model,
      r_carry = r_carry,
      a_r = crm_a_r,
      b_r = crm_b_r,
      alpha_sd = crm_alpha_sd,
      model_file = crm_model_file,
      fixed_model_file = crm_fixed_model_file,
      random_model_file = crm_random_model_file,
      level_model_file = crm_level_model_file,
      n_chains = crm_n_chains,
      n_adapt = crm_n_adapt,
      n_burnin = crm_n_burnin,
      n_iter = crm_n_iter,
      thin = crm_thin,
      seed = seed,
      elimi = elimi,
      n_trt_curr = n_trt_curr,
      dose_cap = dose_cap,
      no_skip = TRUE,
      dose_values = crm_dose_values,
      dose_scores = crm_dose_scores,
      time_col = crm_time_col,
      assessment_window = T_assess,
      decision_time = time,
      alpha_grid = crm_alpha_grid,
      alpha_T = crm_alpha_T,
      theta_prior_mean = crm_theta_prior_mean,
      theta_prior_sd = crm_theta_prior_sd,
      alpha_L = crm_alpha_L,
      alpha_rel_tol = crm_alpha_rel_tol,
      alpha_eps = crm_alpha_eps,
      alpha_n_draw_prior = crm_alpha_n_draw_prior,
      cumu_model_file = crm_cumu_model_file,
      cumu_beta0_mean = crm_cumu_beta0_mean,
      cumu_beta0_prec = crm_cumu_beta0_prec,
      cumu_beta0_df = crm_cumu_beta0_df,
      cumu_beta1_shape = crm_cumu_beta1_shape,
      cumu_beta1_rate = crm_cumu_beta1_rate,
      cumu_beta2_rate = crm_cumu_beta2_rate,
      cumu_include_current = crm_cumu_include_current
    )

    if (!is.null(move$eliminated)) {
      elimi <<- move$eliminated
    }

    decision_log <<- rbind(
      decision_log,
      data.frame(
        time = time,
        cohort = cohort_id,
        current_dose = current_dose,
        next_dose = as.integer(move$next_dose),
        action = move$action,
        n_admin = nrow(admin),
        queue = nrow(waiting),
        stringsAsFactors = FALSE
      )
    )

    if (isTRUE(move$stop_trial) || isTRUE(as.logical(move$earlystop))) {
      earlystop <<- 1L
      trial_stop <<- TRUE
      return(invisible(NULL))
    }

    current_dose <<- max(1L, min(ndose, as.integer(move$next_dose)))
    cohort_id <<- cohort_id + 1L
    cohort_enrolled <<- 0L
    cohort_open <<- TRUE
    invisible(NULL)
  }

  while (!trial_stop &&
         nrow(admin) < Nmax_eff &&
         (nrow(future_events) > 0L || nrow(waiting) > 0L)) {
    if (!has_assignable_waiting()) {
      if (nrow(future_events) == 0L) break
      future_events <- future_events[order(future_events$time, future_events$seq), , drop = FALSE]
      t_now <- future_events$time[1L]
      add_to_waiting(pop_due_events(t_now))
    }

    if (cohort_open) {
      fill_open_cohort(t_now)
    }

    repeat {
      if (trial_stop || nrow(admin) >= Nmax_eff) break
      if (!cohort_open && nrow(waiting) > 0L) {
        open_next_cohort(t_now)
        if (trial_stop) break
        fill_open_cohort(t_now)
      } else {
        break
      }
    }
  }

  if (nrow(admin) == 0L) {
    return(list(
      admin = admin,
      waiting = waiting,
      future_events = future_events,
      decision_log = decision_log,
      final = list(
        t_end = NA_real_,
        MTD = 99L,
        trial_time = NA_real_,
        eliminated = elimi,
        earlystop = 1L,
        phat = rep(NA_real_, ndose),
        pj_iso = rep(NA_real_, ndose),
        r_model = crm_r_model,
        model = "CRM_TITE"
      )
    ))
  }

  t_end <- max(admin$t_eval, na.rm = TRUE)
  dat_final <- aide_crm_tite_make_decision_dat(admin, t_end, T_assess)

  if (trial_stop || earlystop == 1L) {
    final_fit <- list(
      MTD = 99L,
      phat = rep(NA_real_, ndose),
      pj_iso = rep(NA_real_, ndose),
      eliminated = elimi,
      earlystop = 1L,
      stop = 1L
    )
  } else {
    final_fit <- select.mtd.crm(
      target = target,
      dat = dat_final,
      ndose = ndose,
      skeleton = crm_skeleton,
      cutoff.eli = cutoff,
      r_model = crm_r_model,
      r_carry = r_carry,
      a_r = crm_a_r,
      b_r = crm_b_r,
      alpha_sd = crm_alpha_sd,
      model_file = crm_model_file,
      fixed_model_file = crm_fixed_model_file,
      random_model_file = crm_random_model_file,
      level_model_file = crm_level_model_file,
      n_chains = crm_n_chains,
      n_adapt = crm_n_adapt,
      n_burnin = crm_n_burnin,
      n_iter = max(5000, crm_n_iter),
      thin = crm_thin,
      seed = seed,
      restrict_to_tried = restrict_to_tried,
      dose_values = crm_dose_values,
      dose_scores = crm_dose_scores,
      time_col = crm_time_col,
      assessment_window = T_assess,
      decision_time = t_end,
      alpha_grid = crm_alpha_grid,
      alpha_T = crm_alpha_T,
      theta_prior_mean = crm_theta_prior_mean,
      theta_prior_sd = crm_theta_prior_sd,
      alpha_L = crm_alpha_L,
      alpha_rel_tol = crm_alpha_rel_tol,
      alpha_eps = crm_alpha_eps,
      alpha_n_draw_prior = crm_alpha_n_draw_prior,
      cumu_model_file = crm_cumu_model_file,
      cumu_beta0_mean = crm_cumu_beta0_mean,
      cumu_beta0_prec = crm_cumu_beta0_prec,
      cumu_beta0_df = crm_cumu_beta0_df,
      cumu_beta1_shape = crm_cumu_beta1_shape,
      cumu_beta1_rate = crm_cumu_beta1_rate,
      cumu_beta2_rate = crm_cumu_beta2_rate,
      cumu_include_current = crm_cumu_include_current
    )
  }

  trial_time <- max(admin$t_eval, na.rm = TRUE) -
    min(admin$t_arrival, na.rm = TRUE)

  list(
    admin = admin,
    waiting = waiting,
    future_events = future_events,
    decision_log = decision_log,
    final = list(
      t_end = t_end,
      MTD = final_fit$MTD,
      trial_time = trial_time,
      eliminated = if (!is.null(final_fit$eliminated)) final_fit$eliminated else elimi,
      earlystop = if (!is.null(final_fit$earlystop)) final_fit$earlystop else earlystop,
      phat = final_fit$phat,
      pj_iso = final_fit$pj_iso,
      p_hat = if (!is.null(final_fit$p_hat)) final_fit$p_hat else final_fit$phat,
      prob_p1_over_target = if (!is.null(final_fit$prob_p1_over_target)) {
        final_fit$prob_p1_over_target
      } else {
        NA_real_
      },
      r_model = crm_r_model,
      model_file = if (!is.null(final_fit$model_file)) final_fit$model_file else crm_model_file,
      model = "CRM_TITE"
    )
  )
}

get_oc_sim_AIDE_CRM_TITE <- function(
    target,
    p.true,
    p.true_ipde = p.true,
    ntrial = 1000,
    seed = 1,
    store_raw = FALSE,
    ...
) {
  ndose <- length(p.true)
  sel_count <- integer(ndose)
  stop_count <- 0L
  na_count <- 0L

  n_by_dose <- numeric(ndose)
  y_by_dose <- numeric(ndose)
  nipde_by_dose <- numeric(ndose)
  unique_n_by_dose <- numeric(ndose)
  total_admin <- numeric(ntrial)
  total_unique <- numeric(ntrial)
  duration <- numeric(ntrial)
  final_mtd_by_trial <- rep(NA_integer_, ntrial)
  earlystop_by_trial <- integer(ntrial)

  raw_trials <- if (store_raw) vector("list", ntrial) else NULL

  for (itrial in seq_len(ntrial)) {
    fit <- simulate_AIDE_CRM_TITE_design(
      p_true = p.true,
      p_ipde = p.true_ipde,
      target = target,
      seed = seed + itrial - 1L,
      ...
    )

    if (store_raw) raw_trials[[itrial]] <- fit
    admin <- fit$admin

    if (nrow(admin) > 0L) {
      n_by_dose <- n_by_dose + tabulate(admin$dose, nbins = ndose)
      y_by_dose <- y_by_dose + tabulate(admin$dose[admin$y == 1L], nbins = ndose)
      nipde_by_dose <- nipde_by_dose +
        tabulate(admin$dose[admin$type == "retreat"], nbins = ndose)

      id_dose <- unique(admin[, c("id", "dose"), drop = FALSE])
      unique_n_by_dose <- unique_n_by_dose +
        tabulate(id_dose$dose, nbins = ndose)

      total_admin[itrial] <- nrow(admin)
      total_unique[itrial] <- length(unique(admin$id))
    }

    duration[itrial] <- fit$final$trial_time
    mtd <- fit$final$MTD
    final_mtd_by_trial[itrial] <- as.integer(mtd)
    earlystop_by_trial[itrial] <- as.integer(isTRUE(as.logical(fit$final$earlystop)))
    if (isTRUE(as.logical(fit$final$earlystop)) || identical(mtd, 99L)) {
      stop_count <- stop_count + 1L
    } else if (!is.na(mtd) && mtd >= 1L && mtd <= ndose) {
      sel_count[mtd] <- sel_count[mtd] + 1L
    } else {
      na_count <- na_count + 1L
    }
  }

  list(
    target = target,
    p.true = p.true,
    p.true_ipde = p.true_ipde,
    ntrial = ntrial,
    ndose = ndose,
    final_mtd_by_trial = final_mtd_by_trial,
    earlystop_by_trial = earlystop_by_trial,
    n_by_dose = n_by_dose / ntrial,
    unique_n_by_dose = unique_n_by_dose / ntrial,
    nipde_by_dose = nipde_by_dose / ntrial,
    total_admin_by_trial = total_admin,
    total_unique_by_trial = total_unique,
    duration_by_trial = duration,
    selpercent = 100 * sel_count / ntrial,
    npatients = n_by_dose / ntrial,
    ntox = y_by_dose / ntrial,
    nipde = nipde_by_dose / ntrial,
    nuniquepatients = unique_n_by_dose / ntrial,
    totaln = mean(total_admin),
    total_unique = mean(total_unique),
    percentstop = 100 * stop_count / ntrial,
    percentna = 100 * na_count / ntrial,
    duration = mean(duration, na.rm = TRUE),
    raw = raw_trials
  )
}
