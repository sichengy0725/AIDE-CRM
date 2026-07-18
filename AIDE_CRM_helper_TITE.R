## ============================================================
## AIDE-CRM TITE helper and simulation wrapper
##
## This file layers a TITE event-clock simulation over the updated
## likelihoods in AIDE_CRM_helper_modified.R.
##
## Cohort assignment rule (controlled by new_pat_first):
##   - new arrivals (1) or eligible IPDE/retreat patients (2) are prioritized
##   - ipde_design = 1 permits recycling from any lower prior dose; 2 permits
##     recycling only from the immediately lower prior dose
##   - recycled-patient toxicity uses the patient's actual prior dose
##   - remaining events wait until a cohort slot opens
## ============================================================

if (!exists("crm_fit", mode = "function") ||
    !exists("crm_move", mode = "function") ||
    !exists("select.mtd.crm", mode = "function")) {
  source("AIDE_CRM_helper_modified.R")
}

if (!exists("aide_ipde_toxicity_probability", mode = "function")) {
  aide_ipde_toxicity_probability <- function(p_true,
                                              previous_dose,
                                              current_dose,
                                              ipde_alpha = 0) {
    p_true <- as.numeric(p_true)
    ndose <- length(p_true)

    if (ndose < 1L || any(!is.finite(p_true)) || any(p_true < 0 | p_true > 1)) {
      stop("p_true must contain finite toxicity probabilities in [0, 1].")
    }
    if (length(previous_dose) != 1L || is.na(previous_dose) ||
        previous_dose != as.integer(previous_dose) ||
        previous_dose < 1L || previous_dose > ndose) {
      stop("previous_dose must be a valid dose index.")
    }
    if (length(current_dose) != 1L || is.na(current_dose) ||
        current_dose != as.integer(current_dose) ||
        current_dose < 1L || current_dose > ndose) {
      stop("current_dose must be a valid dose index.")
    }
    if (length(ipde_alpha) != 1L || !is.finite(ipde_alpha) || ipde_alpha < 0) {
      stop("ipde_alpha must be a single finite non-negative value.")
    }

    min(1, p_true[as.integer(current_dose)] +
          ipde_alpha * p_true[as.integer(previous_dose)])
  }
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
    # Retained for backward-compatible reporting; recycled-patient toxicity is
    # calculated from p_true, the patient's previous dose, and ipde_alpha.
    p_ipde = p_true,
    ipde_alpha = 0,
    target = 0.30,
    N_pat = 30L,
    Nmax_eff = 30L,
    cohortsize = 3L,
    T_assess = 28,
    cycle_max = 2L,
    arrival_rate = 0.2,
    t0 = 0,
    # 1 = prioritize new arrivals; 2 = prioritize eligible IPDE patients
    new_pat_first = 2L,
    # 1 = recycle from any lower dose; 2 = recycle from the adjacent lower dose
    ipde_design = 2L,
    startdose = 1L,
    cutoff = 0.95,
    dose_cap = 3L,
    # Escalation requires at least this many evaluated patients at the current dose.
    n_eval_escalate = 3L,
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
  if (length(ipde_alpha) != 1L || !is.finite(ipde_alpha) || ipde_alpha < 0) {
    stop("ipde_alpha must be a single finite non-negative value.")
  }
  ipde_alpha <- as.numeric(ipde_alpha)
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
  if (length(new_pat_first) != 1L || is.na(new_pat_first) ||
      !new_pat_first %in% c(1, 2)) {
    stop("new_pat_first must be 1 (new patients first) or 2 (IPDE patients first).")
  }
  new_pat_first <- as.integer(new_pat_first)
  if (length(ipde_design) != 1L || is.na(ipde_design) ||
      !ipde_design %in% c(1, 2)) {
    stop("ipde_design must be 1 (any lower dose) or 2 (adjacent lower dose).")
  }
  ipde_design <- as.integer(ipde_design)
  if (length(n_eval_escalate) != 1L || !is.finite(n_eval_escalate) ||
      n_eval_escalate < 0 || n_eval_escalate != as.integer(n_eval_escalate)) {
    stop("n_eval_escalate must be a single non-negative integer.")
  }
  n_eval_escalate <- as.integer(n_eval_escalate)

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
    assignment_suspended = logical(0),
    stringsAsFactors = FALSE
  )

  elimi <- rep(0L, ndose)
  earlystop <- 0L
  trial_stop <- FALSE
  trial_suspended <- FALSE
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

  has_scheduled_decision <- function(time) {
    any(
      future_events$type == "decision" &
        abs(future_events$time - time) <= 1e-10
    )
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

      ## Decision events only wake a suspended trial; they are not patients.
      if (ev$type == "decision") next

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

    eligible_prior_dose <- if (ipde_design == 1L) {
      last$dose < dose
    } else {
      last$dose == dose - 1L
    }

    isTRUE(
      last$y == 0L &&
        last$t_eval <= ev$time + 1e-10 &&
        last$ncycle < cycle_max &&
        eligible_prior_dose
    )
  }

  next_waiting_index <- function() {
    if (nrow(waiting) == 0L) return(integer(0))

    new_idx <- which(waiting$type == "new")
    retreat_idx <- which(waiting$type == "retreat")
    eligible_retreat_idx <- retreat_idx[vapply(
      retreat_idx,
      function(i) recycle_event_is_eligible(waiting[i, , drop = FALSE]),
      logical(1)
    )]

    if (new_pat_first == 1L) {
      if (length(new_idx) > 0L) return(new_idx[1L])
      if (length(eligible_retreat_idx) > 0L) return(eligible_retreat_idx[1L])
    } else {
      if (length(eligible_retreat_idx) > 0L) return(eligible_retreat_idx[1L])
      if (length(new_idx) > 0L) return(new_idx[1L])
    }

    integer(0)
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
      pi_val <- aide_ipde_toxicity_probability(
        p_true = p_true,
        previous_dose = last$dose,
        current_dose = current_dose,
        ipde_alpha = ipde_alpha
      )
      row_type <- "retreat"
    }

    dlt <- aide_crm_tite_draw_dlt(
      pi_val = pi_val,
      # A queued patient begins treatment only when a cohort slot is assigned.
      t_start = time,
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
        t_start = time,
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

    ## Reassess as soon as an outcome becomes evaluable, even if no new
    ## patient is waiting to open the following cohort.
    if (!has_scheduled_decision(dlt$t_eval)) {
      schedule_event(dlt$t_eval, "decision", NA_integer_)
    }

    if (dlt$y == 0L && row_cycle < cycle_max) {
      schedule_event(dlt$t_eval, "retreat", ev$id)
    }

    cohort_enrolled <<- cohort_enrolled + 1L
    if (cohort_enrolled >= cohortsize) {
      cohort_open <<- FALSE
    }

    if (verbose) {
      message(
        "t=", round(time, 2),
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

  open_next_cohort <- function(time,
                               advance_cohort = TRUE,
                               apply_move = TRUE) {
    if (nrow(admin) == 0L) {
      cohort_open <<- TRUE
      cohort_enrolled <<- 0L
      return(invisible(NULL))
    }

    dat_dec <- aide_crm_tite_make_decision_dat(admin, time, T_assess)
    n_trt_curr <- sum(admin$dose == current_dose & admin$t_start <= time)
    n_eval_curr <- sum(admin$dose == current_dose & admin$t_eval <= time)

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

    insufficient_evaluated <- n_eval_curr < n_eval_escalate
    suspend_for_evaluation <- isTRUE(apply_move) &&
      isTRUE(move$next_dose > current_dose) &&
      insufficient_evaluated

    if (suspend_for_evaluation) {
      ## Do not assign another cohort at the current dose. New arrivals remain
      ## in waiting while the trial waits for more information.
      move$action <- "suspend_insufficient_evaluated"
    } else if (!apply_move &&
               !isTRUE(move$stop_trial) &&
               !isTRUE(as.logical(move$earlystop))) {
      ## A partially enrolled cohort continues at its assigned dose, but the
      ## evaluation gate is recorded at every arrival and evaluation.
      move$action <- if (insufficient_evaluated) {
        "suspend_insufficient_evaluated"
      } else {
        "cohort_incomplete"
      }
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
        assignment_suspended = isTRUE(suspend_for_evaluation),
        stringsAsFactors = FALSE
      )
    )

    if (isTRUE(move$stop_trial) || isTRUE(as.logical(move$earlystop))) {
      earlystop <<- 1L
      trial_stop <<- TRUE
      return(invisible(NULL))
    }

    if (suspend_for_evaluation) {
      ## Any newly evaluable patient can update the CRM decision, whereas the
      ## gate itself remains based on evaluated patients at current_dose.
      next_eval <- admin$t_eval[admin$t_eval > time + 1e-10]
      next_eval <- next_eval[is.finite(next_eval)]
      if (length(next_eval) == 0L) {
        stop("Cannot suspend: no pending patient evaluations.")
      }

      next_eval <- min(next_eval)
      if (!has_scheduled_decision(next_eval)) {
        schedule_event(next_eval, "decision", NA_integer_)
      }
      trial_suspended <<- TRUE
      cohort_open <<- FALSE
      return(invisible(NULL))
    }

    if (!apply_move) {
      return(invisible(NULL))
    }

    trial_suspended <<- FALSE
    current_dose <<- max(1L, min(ndose, as.integer(move$next_dose)))
    if (advance_cohort) {
      cohort_id <<- cohort_id + 1L
      cohort_enrolled <<- 0L
    }
    cohort_open <<- TRUE
    invisible(NULL)
  }

  while (!trial_stop &&
         nrow(admin) < Nmax_eff &&
         (nrow(future_events) > 0L || nrow(waiting) > 0L)) {
    if (trial_suspended) {
      ## Continue accruing arrivals into waiting, but assign no dose until a
      ## pending evaluation or new arrival makes a new decision possible.
      if (nrow(future_events) == 0L) break
      future_events <- future_events[order(future_events$time, future_events$seq), , drop = FALSE]
      t_now <- future_events$time[1L]
      due_events <- pop_due_events(t_now)
      decision_due <- any(due_events$type == "decision")
      arrival_due <- any(due_events$type == "new")
      add_to_waiting(due_events)

      if (decision_due || arrival_due) {
        open_next_cohort(t_now)
        if (!trial_stop && !trial_suspended && cohort_open) {
          fill_open_cohort(t_now)
        }
      }
      next
    }

    decision_due <- FALSE
    arrival_due <- FALSE
    if (!has_assignable_waiting()) {
      if (nrow(future_events) == 0L) break
      future_events <- future_events[order(future_events$time, future_events$seq), , drop = FALSE]
      t_now <- future_events$time[1L]
      due_events <- pop_due_events(t_now)
      decision_due <- any(due_events$type == "decision")
      arrival_due <- any(due_events$type == "new")
      add_to_waiting(due_events)
    }

    if (!trial_suspended && cohort_open) {
      fill_open_cohort(t_now)
    }

    ## Fit the model after every new arrival and every newly evaluable outcome.
    ## A partially enrolled cohort retains its assigned dose until completion.
    ## The first arrival has no prior outcome information to update.
    check_due <- decision_due || (arrival_due && nrow(admin) > 1L)
    if (!trial_stop && !trial_suspended && check_due) {
      if (!cohort_open) {
        open_next_cohort(t_now)
      } else {
        open_next_cohort(
          t_now,
          advance_cohort = FALSE,
          apply_move = cohort_enrolled == 0L
        )
      }
    }

    repeat {
      if (trial_stop || trial_suspended || nrow(admin) >= Nmax_eff) break
      if (!cohort_open) {
        open_next_cohort(t_now)
        if (trial_stop || trial_suspended) break
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
      new_pat_first = new_pat_first,
      ipde_design = ipde_design,
      n_eval_escalate = n_eval_escalate,
      ipde_alpha = ipde_alpha,
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
    browser()
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
    new_pat_first = new_pat_first,
    ipde_design = ipde_design,
    n_eval_escalate = n_eval_escalate,
    ipde_alpha = ipde_alpha,
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
    ipde_alpha = 0,
    ntrial = 1000,
    seed = 1,
    new_pat_first = 2L,
    ipde_design = 2L,
    n_eval_escalate = 3L,
    store_raw = FALSE,
    ...
) {
  if (length(new_pat_first) != 1L || is.na(new_pat_first) ||
      !new_pat_first %in% c(1, 2)) {
    stop("new_pat_first must be 1 (new patients first) or 2 (IPDE patients first).")
  }
  new_pat_first <- as.integer(new_pat_first)
  if (length(ipde_design) != 1L || is.na(ipde_design) ||
      !ipde_design %in% c(1, 2)) {
    stop("ipde_design must be 1 (any lower dose) or 2 (adjacent lower dose).")
  }
  ipde_design <- as.integer(ipde_design)
  if (length(n_eval_escalate) != 1L || !is.finite(n_eval_escalate) ||
      n_eval_escalate < 0 || n_eval_escalate != as.integer(n_eval_escalate)) {
    stop("n_eval_escalate must be a single non-negative integer.")
  }
  n_eval_escalate <- as.integer(n_eval_escalate)
  if (length(ipde_alpha) != 1L || !is.finite(ipde_alpha) || ipde_alpha < 0) {
    stop("ipde_alpha must be a single finite non-negative value.")
  }
  ipde_alpha <- as.numeric(ipde_alpha)

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
      ipde_alpha = ipde_alpha,
      target = target,
      seed = seed + itrial - 1L,
      new_pat_first = new_pat_first,
      ipde_design = ipde_design,
      n_eval_escalate = n_eval_escalate,
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
    ipde_alpha = ipde_alpha,
    ntrial = ntrial,
    ndose = ndose,
    new_pat_first = new_pat_first,
    ipde_design = ipde_design,
    n_eval_escalate = n_eval_escalate,
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
