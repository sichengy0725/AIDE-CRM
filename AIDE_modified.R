## ============================================================
## AIDE simulation wrapper with continuous new-patient accrual,
## CFO, and five CRM backends:
##   fixed/r_fixed, random, level, alpha_crm, cumu_crm
##
## Source this after BOIN helpers and AIDE_CRM_helper_final.R.
## ============================================================

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

simulate_AIDE_design <- function(
    # --- trial size / timing ---
  N_pat       = 30L,
  Nmax_eff    = 30L,
  C           = 3L,
  T_assess    = 28,
  cycle_max   = 2L,

  # --- arrival process ---
  arrival_rate = 0.2,
  t0 = 0,
  continuous_enrollment = TRUE,
  # 1 = prioritize new arrivals; 2 = prioritize eligible IPDE patients
  new_pat_first = 2L,

  # --- DLT time generation ---
  dlt_dist  = 2,
  dlt_alpha = 0.5,

  # --- truth ---
  p_true,
  # Retained for backward-compatible reporting; recycled-patient toxicity is
  # calculated from p_true, the patient's previous dose, and ipde_alpha.
  p_ipde = p_true,
  ipde_alpha = 0,
  seed = NULL,
  verbose = FALSE,

  TARGET = 0.3,
  cutoff = 0.95,

  model = c("BOIN", "CRM", "CFO"),
  ipde_design = 2,
  d.cap = 100,
  day_obs = 0,

  # --- escalation gate ---
  dose_cap = 3L,
  # Escalation requires at least this many evaluated patients at the current dose.
  n_eval_escalate = 3L,

  # --- BOIN decision / final selection method ---
  decision_method = c("boin", "approx1", "approx2"),
  mtd_method = NULL,
  restrict_to_tried = TRUE,
  restrict_to_target = FALSE,
  r_carry = 0.1,
  r_estimator = c("r_fixed", "r_mle", "r_adaptive"),
  r_adaptive_prior = c(1, 9),
  p_adaptive_prior = c(1, 1),
  r_adaptive_max = 0.6,
  r_adaptive_plug_in = c("mean", "map"),
  r_adaptive_rel_tol = 1e-6,

  # --- CRM model options ---
  crm_r_model = c("fixed", "r_fixed", "random", "level", "alpha_crm", "cumu_crm"),
  crm_skeleton = NULL,
  crm_model_file = NULL,
  crm_fixed_model_file = "fix_CRM.bug",
  crm_random_model_file = "random_CRM.bug",
  crm_level_model_file = "random_CRM_level.bug",
  crm_alpha_sd = 2,
  crm_a_r = 1,
  crm_b_r = 9,

  # --- alpha-CRM options ---
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

  # --- cumulative CRM options ---
  crm_dose_scores = NULL,
  crm_cumu_model_file = NULL,
  crm_cumu_beta0_mean = stats::qlogis(TARGET),
  crm_cumu_beta0_prec = 4,
  crm_cumu_beta0_df = 1,
  crm_cumu_beta1_shape = 5.83,
  crm_cumu_beta1_rate = 1.21,
  crm_cumu_beta2_rate = 1,
  crm_cumu_include_current = FALSE,

  # --- CFO / PRIDE options ---
  cfo_method = c("empirical", "pride"),
  cfo_skeleton = NULL,
  cfo_mu = NULL,
  cfo_model_file = "PRIDE.bug",
  cfo_sigma2_beta = 10,
  cfo_eta = 1,
  cfo_pk_method = c("approx", "mc"),
  cfo_n_mc_w = 200,
  cfo_m_use = 1000,
  cfo_use_monotone_pair = FALSE,
  cfo_restrict_to_tried = TRUE,

  # --- MCMC options ---
  crm_n_chains = 2,
  crm_n_adapt = 500,
  crm_n_burnin = 500,
  crm_n_iter = 2000,
  crm_thin = 1,
  cfo_n_chains = 3,
  cfo_n_adapt = 1000,
  cfo_n_burnin = 2000,
  cfo_n_iter = 5000,
  cfo_thin = 2
) {
  if (!is.null(seed)) set.seed(seed)

  model <- match.arg(model)
  restrict_to_tried <- isTRUE(restrict_to_tried)
  restrict_to_target <- isTRUE(restrict_to_target)
  if (length(new_pat_first) != 1L || is.na(new_pat_first) ||
      !new_pat_first %in% c(1, 2)) {
    stop("new_pat_first must be 1 (new patients first) or 2 (IPDE patients first).")
  }
  new_pat_first <- as.integer(new_pat_first)
  if (length(n_eval_escalate) != 1L || !is.finite(n_eval_escalate) ||
      n_eval_escalate < 0 || n_eval_escalate != as.integer(n_eval_escalate)) {
    stop("n_eval_escalate must be a single non-negative integer.")
  }
  n_eval_escalate <- as.integer(n_eval_escalate)

  if (model == "BOIN") {
    decision_method <- match.arg(decision_method)
    r_estimator <- match.arg(r_estimator)
    r_adaptive_plug_in <- match.arg(r_adaptive_plug_in)

    if (is.null(mtd_method)) {
      mtd_method <- decision_method
    }

    mtd_method <- match.arg(
      mtd_method,
      choices = c("boin", "approx1", "approx2")
    )
  }

  if (model == "CRM") {
    crm_r_model <- match.arg(
      crm_r_model,
      choices = c("fixed", "r_fixed", "random", "level", "alpha_crm", "cumu_crm")
    )
    if (crm_r_model == "r_fixed") crm_r_model <- "fixed"

    if (is.null(crm_skeleton)) {
      stop("For model = 'CRM', provide crm_skeleton.")
    }

    if (length(crm_skeleton) != length(p_true)) {
      stop("crm_skeleton must have the same length as p_true.")
    }

    if (any(!is.finite(crm_skeleton)) || any(crm_skeleton <= 0 | crm_skeleton >= 1)) {
      stop("crm_skeleton must contain probabilities in (0, 1).")
    }

    if (!is.finite(crm_alpha_sd) || crm_alpha_sd <= 0) {
      stop("crm_alpha_sd must be positive.")
    }

    if (crm_r_model == "level" && length(p_true) != 5L) {
      stop("crm_r_model = 'level' requires exactly 5 dose levels.")
    }

    if (crm_r_model == "alpha_crm" &&
        !is.null(crm_dose_values) && length(crm_dose_values) != length(p_true)) {
      stop("crm_dose_values must have the same length as p_true.")
    }

    if (crm_r_model == "cumu_crm" &&
        !is.null(crm_dose_scores) && length(crm_dose_scores) != length(p_true)) {
      stop("crm_dose_scores must have the same length as p_true.")
    }

    decision_method <- paste0("crm_", crm_r_model)
    mtd_method <- decision_method
  }

  if (model == "CFO") {
    cfo_method <- match.arg(cfo_method)
    cfo_pk_method <- match.arg(cfo_pk_method)

    if (!file.exists("CFO_tox_utils.R")) {
      stop("Cannot find CFO utility file: CFO_tox_utils.R")
    }
    source("CFO_tox_utils.R")

    if (cfo_method == "pride") {
      if (is.null(cfo_mu)) {
        if (is.null(cfo_skeleton)) {
          stop("For PRIDE-backed CFO, provide cfo_mu or cfo_skeleton.")
        }
        cfo_mu <- stats::qlogis(cfo_skeleton)
      }
      if (length(cfo_mu) != length(p_true) || any(!is.finite(cfo_mu))) {
        stop("cfo_mu must be finite and have the same length as p_true.")
      }
      if (!file.exists("PRIDE.R")) {
        stop("Cannot find PRIDE helper file: PRIDE.R")
      }
      aide_select_mtd <- select.mtd
      source("PRIDE.R", local = TRUE)
      select.mtd <- aide_select_mtd
    } else if (!is.null(cfo_skeleton) && length(cfo_skeleton) != length(p_true)) {
      stop("cfo_skeleton must have the same length as p_true.")
    }

    decision_method <- paste0("cfo_", cfo_method)
    mtd_method <- decision_method
  }

  if (length(p_ipde) != length(p_true)) {
    stop("p_ipde must have the same length as p_true.")
  }
  if (length(ipde_alpha) != 1L || !is.finite(ipde_alpha) || ipde_alpha < 0) {
    stop("ipde_alpha must be a single finite non-negative value.")
  }
  ipde_alpha <- as.numeric(ipde_alpha)

  if (N_pat < C) {
    stop("Need N_pat >= C.")
  }

  if (length(arrival_rate) != 1L || is.na(arrival_rate) || arrival_rate <= 0) {
    stop("arrival_rate must be a single positive value.")
  }

  K <- length(p_true)
  ndose <- K
  window <- T_assess

  ## Continuous accrual: pre-generate all new-patient arrival times.
  ## This allows backlog to accumulate while active patients are being evaluated.
  arrival_times <- as.numeric(t0 + cumsum(stats::rexp(N_pat, rate = arrival_rate)))
  next_new_idx <- 1L

  draw_dlt_time <- function(pi_val, t_start) {
    ## Binary DLT is generated immediately, but evaluability is at t_start + window.
    ## dlt_dist/dlt_alpha are retained for API compatibility.
    y <- stats::rbinom(n = 1L, size = 1L, prob = pi_val)

    list(
      y      = as.integer(y),
      t_tox  = 0,
      t_eval = as.numeric(t_start + window)
    )
  }

  get_patient_state <- function(admin) {
    if (nrow(admin) == 0) return(admin)

    a <- admin[order(admin$id, admin$t_eval, admin$row_id), , drop = FALSE]
    a[!duplicated(a$id, fromLast = TRUE), , drop = FALSE]
  }

  eligible_ipde_ids <- function(admin, next_dose, cycle_max, design, t_available) {
    st <- get_patient_state(admin)
    if (nrow(st) == 0) return(integer(0))

    ## only patients whose most recent dose has completed evaluation can be recycled
    st <- st[st$t_eval <= t_available, , drop = FALSE]
    if (nrow(st) == 0) return(integer(0))

    if (design == 1) {
      ok <- (st$y == 0L) &
        (st$ncycle < cycle_max) &
        (st$dose < next_dose)
    } else if (design == 2) {
      ok <- (st$y == 0L) &
        (st$ncycle < cycle_max) &
        (st$dose == next_dose - 1L)
    } else {
      return(integer(0))
    }

    st <- st[ok, , drop = FALSE]
    st <- st[order(st$t_eval, st$id), , drop = FALSE]
    st$id
  }

  add_arrivals_to_waiting <- function(cut_time) {
    if (!continuous_enrollment) return(invisible(NULL))

    while (next_new_idx <= N_pat && arrival_times[next_new_idx] <= cut_time) {
      waiting <<- rbind(waiting, data.frame(
        id = next_new_idx,
        t_arrival = arrival_times[next_new_idx],
        stringsAsFactors = FALSE
      ))
      next_new_idx <<- next_new_idx + 1L
    }
    invisible(NULL)
  }

  take_waiting <- function(n_needed) {
    if (n_needed <= 0L || nrow(waiting) == 0L) {
      return(data.frame(id = integer(0), t_arrival = numeric(0), stringsAsFactors = FALSE))
    }

    n_take <- min(n_needed, nrow(waiting))
    out <- waiting[seq_len(n_take), , drop = FALSE]
    waiting_new <- waiting[-seq_len(n_take), , drop = FALSE]
    rownames(waiting_new) <- NULL
    waiting <<- waiting_new
    rownames(out) <- NULL
    out
  }

  wait_for_new_patients_from_queue <- function(n_needed, earliest_time) {
    if (n_needed <= 0L) {
      return(list(
        patients = data.frame(id = integer(0), t_arrival = numeric(0), stringsAsFactors = FALSE),
        t_start = as.numeric(earliest_time)
      ))
    }

    ## In continuous-enrollment runs, future arrivals enter the same waiting
    ## queue used for cohort assignment.  This avoids bypassing eligible IPDE
    ## patients merely to pre-fill a new-patient cohort.
    if (continuous_enrollment) {
      add_arrivals_to_waiting(earliest_time)
      patients <- take_waiting(n_needed)

      while (nrow(patients) < n_needed && next_new_idx <= N_pat) {
        next_arrival <- arrival_times[next_new_idx]
        add_arrivals_to_waiting(next_arrival)
        patients <- rbind(
          patients,
          take_waiting(n_needed - nrow(patients))
        )
      }

      rownames(patients) <- NULL
      return(list(
        patients = patients,
        t_start = if (nrow(patients) > 0L) {
          max(as.numeric(earliest_time), max(patients$t_arrival))
        } else {
          as.numeric(earliest_time)
        }
      ))
    }

    ## Without continuous enrollment, no waiting queue exists; retain the
    ## legacy just-in-time recruitment behavior.
    recruit_new_patients(n_needed = n_needed, earliest_time = earliest_time)
  }

  recruit_new_patients <- function(n_needed, earliest_time) {
    if (n_needed <= 0L) {
      return(list(
        patients = data.frame(id = integer(0), t_arrival = numeric(0), stringsAsFactors = FALSE),
        t_start = as.numeric(earliest_time)
      ))
    }

    if (continuous_enrollment) {
      add_arrivals_to_waiting(earliest_time)

      new_patients <- take_waiting(n_needed)
      need_more <- n_needed - nrow(new_patients)
      t_start <- as.numeric(earliest_time)

      if (need_more > 0L) {
        n_left <- N_pat - next_new_idx + 1L
        if (n_left <= 0L) {
          return(list(patients = new_patients, t_start = t_start))
        }

        n_take <- min(as.integer(need_more), as.integer(n_left))
        ids <- seq.int(next_new_idx, length.out = n_take)
        future <- data.frame(
          id = ids,
          t_arrival = arrival_times[ids],
          stringsAsFactors = FALSE
        )
        next_new_idx <<- next_new_idx + n_take
        new_patients <- rbind(new_patients, future)
        t_start <- max(t_start, max(future$t_arrival))
      }

      rownames(new_patients) <- NULL
      return(list(patients = new_patients, t_start = as.numeric(t_start)))
    }

    ## Backward-compatible just-in-time accrual: generate arrivals only when needed.
    n_left <- N_pat - next_new_idx + 1L
    if (n_left <= 0L) {
      return(list(
        patients = data.frame(id = integer(0), t_arrival = numeric(0), stringsAsFactors = FALSE),
        t_start = as.numeric(earliest_time)
      ))
    }

    n_take <- min(as.integer(n_needed), as.integer(n_left))
    ids <- seq.int(next_new_idx, length.out = n_take)
    inter <- stats::rexp(n_take, rate = arrival_rate)
    arr <- as.numeric(earliest_time + cumsum(inter))
    next_new_idx <<- next_new_idx + n_take

    list(
      patients = data.frame(id = ids, t_arrival = arr, stringsAsFactors = FALSE),
      t_start = max(as.numeric(earliest_time), max(arr))
    )
  }

  make_decision_dat <- function(admin, t_decision) {
    dat_dec <- admin[
      admin$t_eval <= t_decision,
      c("id", "dose", "y", "type", "t_arrival", "t_start", "t_tox", "t_eval", "ncycle", "cohort"),
      drop = FALSE
    ]
    if (nrow(dat_dec) > 0L) dat_dec$cycle <- dat_dec$ncycle
    dat_dec
  }

  ## ---------------- BOIN boundaries ----------------
  ncohort_boundary <- ceiling(Nmax_eff / C)

  temp <- get.boundary(
    target = TARGET,
    ncohort = ncohort_boundary,
    cohortsize = C,
    design = 3,
    cutoff.eli = cutoff
  )

  b.e <- temp[4, ]
  b.d <- temp[3, ]

  bd_lambda <- boin_boundary(TARGET)
  lambda_e <- bd_lambda$lambda_e
  lambda_d <- bd_lambda$lambda_d

  ## ---------------- storage ----------------
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

  waiting <- data.frame(
    id = integer(0),
    t_arrival = numeric(0),
    stringsAsFactors = FALSE
  )

  elimi <- rep(0L, K)
  earlystop <- 0L

  ## ---------------- 1) first cohort ----------------
  curr_dose <- 1L
  cohort_id <- 1L

  first_info <- recruit_new_patients(C, earliest_time = t0)
  first_new <- first_info$patients

  if (nrow(first_new) < C) {
    stop("Not enough unique patients to enroll cohort 1.")
  }

  first_ids <- first_new$id
  t_start1 <- first_info$t_start

  for (pid in first_ids) {
    dlt <- draw_dlt_time(
      pi_val = p_true[curr_dose],
      t_start = t_start1
    )

    admin <- rbind(admin, data.frame(
      row_id = nrow(admin) + 1L,
      id = pid,
      t_arrival = first_new$t_arrival[match(pid, first_new$id)],
      t_start = t_start1,
      t_tox = dlt$t_tox,
      t_eval = dlt$t_eval,
      dose = curr_dose,
      y = dlt$y,
      ncycle = 1L,
      cohort = cohort_id,
      type = "new",
      stringsAsFactors = FALSE
    ))
  }

  t_decision <- max(admin$t_eval[admin$cohort == cohort_id]) + day_obs
  stop_trial <- FALSE

  if (verbose) {
    message(
      "Cohort 1: dose = ", curr_dose,
      ", start = ", round(t_start1, 2),
      ", decision_time = ", round(t_decision, 2)
    )
  }

  ## ---------------- main loop ----------------
  repeat {
    if (stop_trial) break
    if (nrow(admin) >= Nmax_eff) break
    add_arrivals_to_waiting(t_decision)

    dat_dec <- make_decision_dat(admin, t_decision)

    n <- tabulate(dat_dec$dose, nbins = K)
    y <- tabulate(dat_dec$dose[dat_dec$y == 1L], nbins = K)

    n_curr <- n[curr_dose]
    y_curr <- y[curr_dose]

    ## BOIN empirical overdose elimination only.
    ## CRM dose-1 overdose stopping is handled inside crm_move().
    if (model == "BOIN") {
      post_overtox <- 1 - stats::pbeta(
        TARGET,
        1 + y_curr,
        1 + n_curr - y_curr
      )

      if (n_curr >= 3 && post_overtox > cutoff) {
        elimi[curr_dose:K] <- 1L

        if (curr_dose == 1L) {
          earlystop <- 1L
          break
        }
      }
    }

    n_trt_curr <- sum(
      admin$dose == curr_dose &
        admin$t_start <= t_decision
    )

    dat_curr <- dat_dec[dat_dec$dose == curr_dose, , drop = FALSE]

    NR_curr <- sum(dat_curr$type == "new")
    YR_curr <- sum(dat_curr$y[dat_curr$type == "new"] == 1L)

    NI_curr <- sum(dat_curr$type == "retreat")
    YI_curr <- sum(dat_curr$y[dat_curr$type == "retreat"] == 1L)

    n_regular_all_dec <- tabulate(
      dat_dec$dose[dat_dec$type == "new"],
      nbins = K
    )

    y_regular_all_dec <- tabulate(
      dat_dec$dose[dat_dec$type == "new" & dat_dec$y == 1L],
      nbins = K
    )

    post_dec_cfo <- NULL

    if (model == "CFO") {
      n_cfo_dec <- tabulate(dat_dec$dose, nbins = K)
      y_cfo_dec <- tabulate(dat_dec$dose[dat_dec$y == 1L], nbins = K)

      if (cfo_method == "pride") {
        dat_pride_dec <- data.frame(
          id = dat_dec$id,
          dose = dat_dec$dose,
          y = dat_dec$y,
          stringsAsFactors = FALSE
        )

        post_dec_cfo <- get_pride_posterior(
          tmp = dat_pride_dec,
          K = K,
          mu = cfo_mu,
          TARGET = TARGET,
          model_file = cfo_model_file,
          sigma2_beta = cfo_sigma2_beta,
          eta = cfo_eta,
          m_use = cfo_m_use,
          pk_method = cfo_pk_method,
          n_mc_w = cfo_n_mc_w,
          seed_subsample = seed,
          n.chains = cfo_n_chains,
          n.adapt = cfo_n_adapt,
          n.burn = cfo_n_burnin,
          n.iter = cfo_n_iter,
          thin = cfo_thin
        )

        post_overtox_curr <- post_dec_cfo$prob_overtox[curr_dose]
      } else {
        post_overtox_curr <- post.prob.fn(
          phi = TARGET,
          y = y_cfo_dec[curr_dose],
          n = n_cfo_dec[curr_dose],
          alp.prior = TARGET,
          bet.prior = 1 - TARGET
        )
      }

      if (n_cfo_dec[curr_dose] >= 3L && post_overtox_curr > cutoff) {
        elimi[curr_dose:K] <- 1L

        if (curr_dose == 1L) {
          earlystop <- 1L
          break
        }
      }
    }

    if (elimi[curr_dose] == 1L && curr_dose > 1L) {
      move <- list(
        next_dose = curr_dose - 1L,
        action = "de-escalate_elim",
        mu_hat = NA_real_,
        n_eff = n_curr,
        stop_trial = FALSE,
        earlystop = 0L,
        eliminated = elimi
      )
    } else {
      if (model == "BOIN") {
        move <- boin_move(
          current_dose = curr_dose,
          ndose = K,
          method = decision_method,

          y_curr = y_curr,
          n_curr = n_curr,
          b.e = b.e,
          b.d = b.d,
          C = C,

          YR = YR_curr,
          YI = YI_curr,
          NR_star = NR_curr,
          NI_star = NI_curr,
          lambda_e = lambda_e,
          lambda_d = lambda_d,
          phi = TARGET,
          r_carry = r_carry,
          r_estimator = r_estimator,
          r_adaptive_prior = r_adaptive_prior,
          p_adaptive_prior = p_adaptive_prior,
          r_adaptive_max = r_adaptive_max,
          r_adaptive_plug_in = r_adaptive_plug_in,
          r_adaptive_rel_tol = r_adaptive_rel_tol,

          y_regular_all = y_regular_all_dec,
          n_regular_all = n_regular_all_dec,
          y_recycle_all = tabulate(
            dat_dec$dose[dat_dec$type == "retreat" & dat_dec$y == 1L],
            nbins = K
          ),
          n_recycle_all = tabulate(
            dat_dec$dose[dat_dec$type == "retreat"],
            nbins = K
          ),

          elimi = elimi,
          n_trt_curr = n_trt_curr,
          dose_cap = dose_cap
        )

        if (is.null(move$stop_trial)) move$stop_trial <- FALSE
        if (is.null(move$earlystop)) move$earlystop <- 0L
        if (is.null(move$eliminated)) move$eliminated <- elimi

      } else if (model == "CRM") {
        move <- crm_move(
          current_dose = curr_dose,
          ndose = K,
          dat = dat_dec,

          target = TARGET,
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

          elimi = elimi,
          n_trt_curr = n_trt_curr,
          dose_cap = dose_cap,
          no_skip = TRUE,

          dose_values = crm_dose_values,
          dose_scores = crm_dose_scores,
          time_col = crm_time_col,
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
      } else if (model == "CFO") {
        if (cfo_method == "empirical") {
          if (K == 1L) {
            next_cfo <- 1L
          } else if (curr_dose == 1L) {
            cys <- c(NA, y_cfo_dec[1L], y_cfo_dec[2L])
            cns <- c(NA, n_cfo_dec[1L], n_cfo_dec[2L])
            cover.doses <- c(NA, elimi[1L], elimi[2L])
            next_cfo <- make.decision.CFO.fn(
              phi = TARGET,
              cys = cys,
              cns = cns,
              alp.prior = TARGET,
              bet.prior = 1 - TARGET,
              cover.doses = cover.doses
            ) - 2L + curr_dose
          } else {
            idx <- (curr_dose - 1L):(curr_dose + 1L)
            cys <- y_cfo_dec[idx]
            cns <- n_cfo_dec[idx]
            cover.doses <- elimi[idx]
            next_cfo <- make.decision.CFO.fn(
              phi = TARGET,
              cys = cys,
              cns = cns,
              alp.prior = TARGET,
              bet.prior = 1 - TARGET,
              cover.doses = cover.doses
            ) - 2L + curr_dose
          }

          next_cfo <- max(1L, min(K, as.integer(next_cfo)))
          action <- if (next_cfo > curr_dose) {
            "escalate"
          } else if (next_cfo < curr_dose) {
            "de-escalate"
          } else {
            "stay"
          }

          if (next_cfo > curr_dose &&
              (!is.finite(n_trt_curr) || n_trt_curr < dose_cap ||
               elimi[next_cfo] == 1L)) {
            next_cfo <- curr_dose
            action <- "stay_escalation_blocked"
          }

          move <- list(
            next_dose = next_cfo,
            action = action,
            mu_hat = if (n_cfo_dec[curr_dose] > 0) {
              y_cfo_dec[curr_dose] / n_cfo_dec[curr_dose]
            } else {
              NA_real_
            },
            n_eff = sum(n_cfo_dec),
            method = "cfo_empirical",
            stop_trial = FALSE,
            earlystop = 0L,
            eliminated = elimi
          )

          move$prob_overtox <- post.prob.fn(
            phi = TARGET,
            y = y_cfo_dec[curr_dose],
            n = n_cfo_dec[curr_dose],
            alp.prior = TARGET,
            bet.prior = 1 - TARGET
          )
          move$prob_p1_over_target <- post.prob.fn(
            phi = TARGET,
            y = y_cfo_dec[1L],
            n = n_cfo_dec[1L],
            alp.prior = TARGET,
            bet.prior = 1 - TARGET
          )
          move$model_file <- NA_character_
        } else {
          alp.prior <- TARGET
          bet.prior <- 1 - TARGET

          if (K == 1L) {
            gammaL <- Inf
            gammaR <- Inf
          } else if (curr_dose == 1L) {
            gammaL <- Inf
            gammaR <- optim.gamma.fn(
              n1 = n_cfo_dec[curr_dose],
              n2 = n_cfo_dec[curr_dose + 1L],
              phi = TARGET,
              type = "R",
              alp.prior = alp.prior,
              bet.prior = bet.prior
            )$gamma
          } else if (curr_dose == K) {
            gammaL <- optim.gamma.fn(
              n1 = n_cfo_dec[curr_dose - 1L],
              n2 = n_cfo_dec[curr_dose],
              phi = TARGET,
              type = "L",
              alp.prior = alp.prior,
              bet.prior = bet.prior
            )$gamma
            gammaR <- Inf
          } else {
            gammaL <- optim.gamma.fn(
              n1 = n_cfo_dec[curr_dose - 1L],
              n2 = n_cfo_dec[curr_dose],
              phi = TARGET,
              type = "L",
              alp.prior = alp.prior,
              bet.prior = bet.prior
            )$gamma
            gammaR <- optim.gamma.fn(
              n1 = n_cfo_dec[curr_dose],
              n2 = n_cfo_dec[curr_dose + 1L],
              phi = TARGET,
              type = "R",
              alp.prior = alp.prior,
              bet.prior = bet.prior
            )$gamma
          }

          next_cfo <- cfo_move_pride(
            cur_dose = curr_dose,
            pk_draws = post_dec_cfo$pk_draws,
            TARGET = TARGET,
            gammaL = gammaL,
            gammaR = gammaR,
            elim = elimi,
            use_monotone_pair = cfo_use_monotone_pair
          )

          if (next_cfo == 0L) {
            move <- list(
              next_dose = 99L,
              action = "earlystop_all_lower_eliminated",
              mu_hat = post_dec_cfo$posttox[curr_dose],
              n_eff = nrow(dat_dec),
              method = "cfo_pride",
              stop_trial = TRUE,
              earlystop = 1L,
              eliminated = elimi
            )
          } else {
            next_cfo <- max(1L, min(K, as.integer(next_cfo)))
            action <- if (next_cfo > curr_dose) {
              "escalate"
            } else if (next_cfo < curr_dose) {
              "de-escalate"
            } else {
              "stay"
            }

            if (next_cfo > curr_dose &&
                (!is.finite(n_trt_curr) || n_trt_curr < dose_cap ||
                 elimi[next_cfo] == 1L)) {
              next_cfo <- curr_dose
              action <- "stay_escalation_blocked"
            }

            move <- list(
              next_dose = next_cfo,
              action = action,
              mu_hat = post_dec_cfo$posttox[curr_dose],
              n_eff = nrow(dat_dec),
              method = "cfo_pride",
              stop_trial = FALSE,
              earlystop = 0L,
              eliminated = elimi
            )
          }

          move$p_hat <- post_dec_cfo$posttox
          move$prob_overtox <- post_dec_cfo$prob_overtox[curr_dose]
          move$prob_p1_over_target <- post_dec_cfo$prob_overtox[1L]
          move$model_file <- cfo_model_file
        }
      }
    }

    if (isTRUE(move$stop_trial) || isTRUE(as.logical(move$earlystop))) {
      earlystop <- as.integer(isTRUE(as.logical(move$earlystop)))
      if (!is.null(move$eliminated)) elimi <- move$eliminated
      break
    }

    if (!is.null(move$eliminated)) {
      elimi <- move$eliminated
    }

    next_dose <- move$next_dose
    if (next_dose > curr_dose && n_curr < n_eval_escalate) {
      next_dose <- curr_dose
      move$next_dose <- curr_dose
      move$action <- "stay_insufficient_evaluated"
    }

    if (n[curr_dose] >= d.cap && next_dose == curr_dose) {
      break
    }

    ## ---------------- form next cohort ----------------
    cohort_id <- cohort_id + 1L
    t_available <- t_decision
    t_start <- t_available
    new_ids <- integer(0)
    ret_ids <- integer(0)
    new_patients <- data.frame(
      id = integer(0),
      t_arrival = numeric(0),
      stringsAsFactors = FALSE
    )

    cand <- eligible_ipde_ids(
      admin = admin,
      next_dose = next_dose,
      cycle_max = cycle_max,
      design = ipde_design,
      t_available = t_available
    )

    if (new_pat_first == 1L) {
      ## Mirror the TITE priority rule: new patients who have already arrived
      ## and are in the waiting queue take precedence, but do not pull future
      ## arrivals ahead of currently eligible IPDE patients.
      new_patients <- take_waiting(C)
      new_ids <- new_patients$id

      need_ret <- C - length(new_ids)
      if (need_ret > 0L && length(cand) > 0L) {
        ret_ids <- head(cand, min(need_ret, length(cand)))
      }

      ## Only if the waiting queue and eligible IPDE patients cannot fill the
      ## cohort do we wait for additional new arrivals.  This preserves the
      ## non-TITE all-at-once cohort start while applying TITE's queue-first
      ## allocation priority.
      need_new <- C - length(new_ids) - length(ret_ids)
      if (need_new > 0L) {
        new_info <- wait_for_new_patients_from_queue(
          n_needed = need_new,
          earliest_time = t_available
        )
        if (nrow(new_info$patients) < need_new) {
          if (verbose) {
            message(
              "Could not fill cohort after using waiting new patients and eligible IPDE patients: need ", C,
              ", available ", length(new_ids) + length(ret_ids) + nrow(new_info$patients), "."
            )
          }
          break
        }
        new_patients <- rbind(new_patients, new_info$patients)
        rownames(new_patients) <- NULL
        new_ids <- new_patients$id
        t_start <- max(t_start, new_info$t_start)
      }
    } else {
      if (length(cand) > 0L) {
        ret_ids <- head(cand, min(C, length(cand)))
      }

      need_new <- C - length(ret_ids)

      if (need_new > 0L) {
        new_info <- recruit_new_patients(
          n_needed = need_new,
          earliest_time = t_available
        )
        new_patients <- new_info$patients

        if (nrow(new_patients) < need_new) {
          if (verbose) {
            message(
              "Ran out of new patients: need ", need_new,
              ", available ", nrow(new_patients), "."
            )
          }
          break
        }

        new_ids <- new_patients$id
        t_start <- max(t_start, new_info$t_start)
      }
    }

    ## append new-patient outcomes
    for (pid in new_ids) {
      if (nrow(admin) >= Nmax_eff) break

      dlt <- draw_dlt_time(
        pi_val = p_true[next_dose],
        t_start = t_start
      )

      admin <- rbind(admin, data.frame(
        row_id = nrow(admin) + 1L,
        id = pid,
        t_arrival = new_patients$t_arrival[match(pid, new_patients$id)],
        t_start = t_start,
        t_tox = dlt$t_tox,
        t_eval = dlt$t_eval,
        dose = next_dose,
        y = dlt$y,
        ncycle = 1L,
        cohort = cohort_id,
        type = "new",
        stringsAsFactors = FALSE
      ))
    }

    ## append IPDE / retreat outcomes
    if (length(ret_ids) > 0L) {
      st <- get_patient_state(admin)

      for (pid in ret_ids) {
        if (nrow(admin) >= Nmax_eff) break

        last <- st[st$id == pid, , drop = FALSE]
        if (nrow(last) != 1L) next
        
        dlt <- draw_dlt_time(
          pi_val = aide_ipde_toxicity_probability(
            p_true = p_true,
            previous_dose = last$dose,
            current_dose = next_dose,
            ipde_alpha = ipde_alpha
          ),
          t_start = t_start
        )

        admin <- rbind(admin, data.frame(
          row_id = nrow(admin) + 1L,
          id = pid,
          t_arrival = last$t_arrival,
          t_start = t_start,
          t_tox = dlt$t_tox,
          t_eval = dlt$t_eval,
          dose = next_dose,
          y = dlt$y,
          ncycle = as.integer(last$ncycle + 1L),
          cohort = cohort_id,
          type = "retreat",
          stringsAsFactors = FALSE
        ))
      }
    }

    t_decision <- max(admin$t_eval[admin$cohort == cohort_id]) + day_obs
    curr_dose <- next_dose

    if (verbose) {
      message(
        "Cohort ", cohort_id,
        ": dose = ", next_dose,
        ", action = ", move$action,
        ", mu_hat = ", round(move$mu_hat, 4),
        ", start = ", round(t_start, 2),
        ", decision_time = ", round(t_decision, 2),
        ", n_eff = ", nrow(admin),
        ", unique_n = ", length(unique(admin$id)),
        ", waiting = ", nrow(waiting),
        ", new_patients_used = ", length(new_ids),
        ", ipde_used = ", length(ret_ids)
      )
    }
  }

  ## ---------------- final MTD selection ----------------
  if (nrow(admin) == 0L) {
    return(list(
      admin = admin,
      waiting = waiting,
      arrival_times = arrival_times,
      new_pat_first = new_pat_first,
      n_eval_escalate = n_eval_escalate,
      ipde_alpha = ipde_alpha,
      final = list(
        t_end = NA_real_,
        MTD = 99L,
        trial_time = NA_real_,
        eliminated = elimi,
        earlystop = 1L,
        phat = rep(NA_real_, K),
        pj_iso = rep(NA_real_, K),
        decision_method = decision_method,
        mtd_method = mtd_method,
        r_carry = r_carry,
        r_estimator = if (model == "BOIN") r_estimator else NA_character_,
        r_model = if (model == "CRM") crm_r_model else NA_character_,
        cfo_method = if (model == "CFO") cfo_method else NA_character_,
        restrict_to_tried = restrict_to_tried,
        restrict_to_target = restrict_to_target,
        r_hat = rep(NA_real_, K),
        r_use = rep(NA_real_, K),
        r_cap = rep(NA_real_, K),
        model_file = if (model == "CRM") {
          crm_model_file
        } else if (model == "CFO") {
          if (cfo_method == "pride") cfo_model_file else NA_character_
        } else {
          NA_character_
        },
        model = model
      )
    ))
  }

  t_end <- max(admin$t_eval, na.rm = TRUE)
  add_arrivals_to_waiting(t_end)

  dat_final <- make_decision_dat(admin, t_end)

  n_final <- tabulate(dat_final$dose, nbins = K)
  y_final <- tabulate(dat_final$dose[dat_final$y == 1L], nbins = K)

  n_new_final <- tabulate(
    dat_final$dose[dat_final$type == "new"],
    nbins = K
  )

  y_new_final <- tabulate(
    dat_final$dose[dat_final$type == "new" & dat_final$y == 1L],
    nbins = K
  )

  n_recycle_final <- tabulate(
    dat_final$dose[dat_final$type == "retreat"],
    nbins = K
  )

  y_recycle_final <- tabulate(
    dat_final$dose[dat_final$type == "retreat" & dat_final$y == 1L],
    nbins = K
  )
  # browser()
  trial_time <- max(admin$t_eval, na.rm = TRUE) -
    min(admin$t_arrival, na.rm = TRUE)
  phat <- rep(NA_real_, K)
  r_hat <- if (model == "CRM" && crm_r_model == "fixed") r_carry else rep(NA_real_, K)
  r_use <- rep(NA_real_, K)
  r_cap <- rep(NA_real_, K)

  if (earlystop == 1L) {
    final_dose <- 99L
    final_fit <- list(
      phat = rep(NA_real_, K),
      pj_iso = rep(NA_real_, K),
      eliminated = elimi,
      r_hat = r_hat,
      r_use = r_use,
      r_cap = r_cap,
      p_hat = rep(NA_real_, K),
      earlystop = 1L,
      stop = 1L
    )
  } else {
    if (model == "BOIN") {
      final_fit <- select.mtd(
        target = TARGET,
        y = y_final,
        n = n_final,
        cutoff.eli = cutoff,
        approx = mtd_method,
        r_carry = r_carry,
        r_estimator = r_estimator,
        phi = TARGET,
        y_new = y_new_final,
        n_new = n_new_final,
        y_recycle = y_recycle_final,
        n_recycle = n_recycle_final,
        r_adaptive_prior = r_adaptive_prior,
        p_adaptive_prior = p_adaptive_prior,
        r_adaptive_max = r_adaptive_max,
        r_adaptive_plug_in = r_adaptive_plug_in,
        r_adaptive_rel_tol = r_adaptive_rel_tol,
        restrict_to_tried = restrict_to_tried,
        restrict_to_target = restrict_to_target
      )

    } else if (model == "CRM") {
      final_fit <- select.mtd.crm(
        target = TARGET,
        dat = dat_final,
        ndose = K,
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
        restrict_to_tried = restrict_to_tried,
        restrict_to_target = restrict_to_target,

        dose_values = crm_dose_values,
        dose_scores = crm_dose_scores,
        time_col = crm_time_col,
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
    } else if (model == "CFO") {
      final_fit <- select.mtd(
        target = TARGET,
        y = y_final,
        n = n_final,
        cutoff.eli = cutoff,
        approx = "boin",
        restrict_to_tried = restrict_to_tried
      )
      final_fit$approx <- paste0("cfo_", cfo_method)
      final_fit$model_file <- if (cfo_method == "pride") cfo_model_file else NA_character_
    }
    
    final_dose <- final_fit$MTD
    phat[seq_along(final_fit$phat)] <- final_fit$phat
    if (!is.null(final_fit$r_hat)) r_hat <- final_fit$r_hat
    if (!is.null(final_fit$r_use)) r_use <- final_fit$r_use
    if (!is.null(final_fit$r_cap)) r_cap <- final_fit$r_cap
  }

  final_elimi <- if (!is.null(final_fit$eliminated)) final_fit$eliminated else elimi
  list(
    admin = admin,
    waiting = waiting,
    arrival_times = arrival_times,
    new_pat_first = new_pat_first,
    n_eval_escalate = n_eval_escalate,
    ipde_alpha = ipde_alpha,
    final = list(
      t_end = t_end,
      MTD = final_dose,
      trial_time = trial_time,
      eliminated = final_elimi,
      earlystop = if (!is.null(final_fit$earlystop)) final_fit$earlystop else earlystop,
      phat = phat,
      pj_iso = final_fit$pj_iso,
      decision_method = decision_method,
      mtd_method = mtd_method,
      r_carry = r_carry,
      r_estimator = if (model == "BOIN") r_estimator else NA_character_,
      r_model = if (model == "CRM") crm_r_model else NA_character_,
      cfo_method = if (model == "CFO") cfo_method else NA_character_,
      restrict_to_tried = restrict_to_tried,
      restrict_to_target = restrict_to_target,
      r_hat = if (!is.null(final_fit$r_hat)) final_fit$r_hat else r_hat,
      r_use = if (model == "BOIN") r_use else rep(NA_real_, K),
      r_cap = if (model == "BOIN") r_cap else rep(NA_real_, K),
      prob_p1_over_target = if (!is.null(final_fit$prob_p1_over_target)) final_fit$prob_p1_over_target else NA_real_,
      model_file = if (model == "CRM") {
        if (!is.null(final_fit$model_file)) final_fit$model_file else crm_model_file
      } else if (model == "CFO") {
        if (cfo_method == "pride") {
          if (!is.null(final_fit$model_file)) final_fit$model_file else cfo_model_file
        } else {
          NA_character_
        }
      } else {
        NA_character_
      },
      model = model
    )
  )
  
}

## ------------------------------------------------------------
## Operating-characteristic simulation wrapper for AIDE-BOIN / AIDE-CRM / AIDE-CFO
## ------------------------------------------------------------
get_oc_sim_AIDE <- function(
    target,
    p.true,
    p.true_ipde = p.true,
    ipde_alpha = 0,
    ntrial = 1000,
    seed = 1,

    model = c("BOIN", "CRM", "CFO"),
    ipde_design = 2,

    N_pat = 30L,
    Nmax_eff = 30L,
    C = 3L,
    T_assess = 28,
    cycle_max = 2L,

    arrival_rate = 0.2,
    t0 = 0,
    continuous_enrollment = TRUE,
    new_pat_first = 2L,

    cutoff = 0.95,
    d.cap = 100,
    dose_cap = 3L,
    n_eval_escalate = 3L,
    day_obs = 0,

    dlt_dist = 2,
    dlt_alpha = 0.5,

    ## BOIN decision / final MTD method
    decision_method = c("boin", "approx1", "approx2"),
    mtd_method = NULL,
    restrict_to_tried = TRUE,
    restrict_to_target = FALSE,
    r_carry = 0.1,
    r_estimator = c("r_fixed", "r_mle", "r_adaptive"),
    r_adaptive_prior = c(1, 9),
    p_adaptive_prior = c(1, 1),
    r_adaptive_max = 0.6,
    r_adaptive_plug_in = c("mean", "map"),
    r_adaptive_rel_tol = 1e-6,

    ## CRM options
    crm_r_model = c("fixed", "r_fixed", "random", "level", "alpha_crm", "cumu_crm"),
    crm_skeleton = NULL,
    crm_model_file = NULL,
    crm_fixed_model_file = "fix_CRM.bug",
    crm_random_model_file = "random_CRM.bug",
    crm_level_model_file = "random_CRM_level.bug",
    crm_alpha_sd = 2,
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
    crm_cumu_model_file = NULL,
    crm_cumu_beta0_mean = stats::qlogis(target),
    crm_cumu_beta0_prec = 4,
    crm_cumu_beta0_df = 1,
    crm_cumu_beta1_shape = 5.83,
    crm_cumu_beta1_rate = 1.21,
    crm_cumu_beta2_rate = 1,
    crm_cumu_include_current = FALSE,

    ## CFO / PRIDE options
    cfo_method = c("empirical", "pride"),
    cfo_skeleton = NULL,
    cfo_mu = NULL,
    cfo_model_file = "PRIDE.bug",
    cfo_sigma2_beta = 10,
    cfo_eta = 1,
    cfo_pk_method = c("approx", "mc"),
    cfo_n_mc_w = 200,
    cfo_m_use = 1000,
    cfo_use_monotone_pair = FALSE,
    cfo_restrict_to_tried = TRUE,

    crm_n_chains = 2,
    crm_n_adapt = 500,
    crm_n_burnin = 500,
    crm_n_iter = 2000,
    crm_thin = 1,
    cfo_n_chains = 3,
    cfo_n_adapt = 1000,
    cfo_n_burnin = 2000,
    cfo_n_iter = 5000,
    cfo_thin = 2,

    store_raw = FALSE,
    verbose = FALSE
) {
  model <- match.arg(model)
  restrict_to_tried <- isTRUE(restrict_to_tried)
  restrict_to_target <- isTRUE(restrict_to_target)
  if (length(new_pat_first) != 1L || is.na(new_pat_first) ||
      !new_pat_first %in% c(1, 2)) {
    stop("new_pat_first must be 1 (new patients first) or 2 (IPDE patients first).")
  }
  new_pat_first <- as.integer(new_pat_first)
  if (length(n_eval_escalate) != 1L || !is.finite(n_eval_escalate) ||
      n_eval_escalate < 0 || n_eval_escalate != as.integer(n_eval_escalate)) {
    stop("n_eval_escalate must be a single non-negative integer.")
  }
  n_eval_escalate <- as.integer(n_eval_escalate)

  if (model == "BOIN") {
    decision_method <- match.arg(decision_method)
    r_estimator <- match.arg(r_estimator)
    r_adaptive_plug_in <- match.arg(r_adaptive_plug_in)

    if (is.null(mtd_method)) {
      mtd_method <- decision_method
    }

    mtd_method <- match.arg(
      mtd_method,
      choices = c("boin", "approx1", "approx2")
    )
  }

  if (model == "CRM") {
    crm_r_model <- match.arg(
      crm_r_model,
      choices = c("fixed", "r_fixed", "random", "level", "alpha_crm", "cumu_crm")
    )
    if (crm_r_model == "r_fixed") crm_r_model <- "fixed"

    if (is.null(crm_skeleton)) {
      stop("For model = 'CRM', provide crm_skeleton.")
    }

    if (length(crm_skeleton) != length(p.true)) {
      stop("crm_skeleton must have the same length as p.true.")
    }

    if (any(!is.finite(crm_skeleton)) || any(crm_skeleton <= 0 | crm_skeleton >= 1)) {
      stop("crm_skeleton must contain probabilities in (0, 1).")
    }

    if (!is.finite(crm_alpha_sd) || crm_alpha_sd <= 0) {
      stop("crm_alpha_sd must be positive.")
    }

    if (crm_r_model == "level" && length(p.true) != 5L) {
      stop("crm_r_model = 'level' requires exactly 5 dose levels.")
    }

    decision_method <- paste0("crm_", crm_r_model)
    mtd_method <- decision_method
  }

  if (model == "CFO") {
    cfo_method <- match.arg(cfo_method)
    cfo_pk_method <- match.arg(cfo_pk_method)

    if (cfo_method == "pride") {
      if (is.null(cfo_mu)) {
        if (is.null(cfo_skeleton)) {
          stop("For PRIDE-backed CFO, provide cfo_mu or cfo_skeleton.")
        }
        cfo_mu <- stats::qlogis(cfo_skeleton)
      }
      if (length(cfo_mu) != length(p.true) || any(!is.finite(cfo_mu))) {
        stop("cfo_mu must be finite and have the same length as p.true.")
      }
    } else if (!is.null(cfo_skeleton) && length(cfo_skeleton) != length(p.true)) {
      stop("cfo_skeleton must have the same length as p.true.")
    }

    decision_method <- paste0("cfo_", cfo_method)
    mtd_method <- decision_method
  }

  p.true <- as.numeric(p.true)
  p.true_ipde <- as.numeric(p.true_ipde)
  ndose <- length(p.true)

  if (length(p.true_ipde) != ndose) {
    stop("p.true_ipde must have the same length as p.true.")
  }
  if (length(ipde_alpha) != 1L || !is.finite(ipde_alpha) || ipde_alpha < 0) {
    stop("ipde_alpha must be a single finite non-negative value.")
  }
  ipde_alpha <- as.numeric(ipde_alpha)

  sel_count <- integer(ndose)
  stop_count <- 0L
  na_count <- 0L
  earlystop_vec <- integer(ntrial)

  n_by_dose <- numeric(ndose)
  unique_n_by_dose <- numeric(ndose)
  nipde_by_dose <- numeric(ndose)

  total_admin <- numeric(ntrial)
  total_unique <- numeric(ntrial)
  duration <- numeric(ntrial)

  final_MTD <- rep(NA_integer_, ntrial)
  pj_iso_mat <- matrix(NA_real_, nrow = ntrial, ncol = ndose)

  r_hat_mat <- matrix(NA_real_, nrow = ntrial, ncol = ndose)
  r_use_mat <- matrix(NA_real_, nrow = ntrial, ncol = ndose)
  r_cap_mat <- matrix(NA_real_, nrow = ntrial, ncol = ndose)

  raw_trials <- if (store_raw) vector("list", ntrial) else NULL

  for (itrial in seq_len(ntrial)) {
    trial_seed <- seed + itrial - 1L

    fit <- simulate_AIDE_design(
      N_pat = N_pat,
      Nmax_eff = Nmax_eff,
      C = C,
      T_assess = T_assess,
      cycle_max = cycle_max,
      arrival_rate = arrival_rate,
      t0 = t0,
      continuous_enrollment = continuous_enrollment,
      new_pat_first = new_pat_first,
      dlt_dist = dlt_dist,
      dlt_alpha = dlt_alpha,
      p_true = p.true,
      p_ipde = p.true_ipde,
      ipde_alpha = ipde_alpha,
      seed = trial_seed,
      verbose = verbose,
      TARGET = target,
      cutoff = cutoff,
      model = model,
      ipde_design = ipde_design,
      d.cap = d.cap,
      dose_cap = dose_cap,
      n_eval_escalate = n_eval_escalate,
      day_obs = day_obs,
      decision_method = decision_method,
      mtd_method = mtd_method,
      restrict_to_tried = restrict_to_tried,
      restrict_to_target = restrict_to_target,
      r_carry = r_carry,
      r_estimator = r_estimator,
      r_adaptive_prior = r_adaptive_prior,
      p_adaptive_prior = p_adaptive_prior,
      r_adaptive_max = r_adaptive_max,
      r_adaptive_plug_in = r_adaptive_plug_in,
      r_adaptive_rel_tol = r_adaptive_rel_tol,
      crm_r_model = crm_r_model,
      crm_skeleton = crm_skeleton,
      crm_model_file = crm_model_file,
      crm_fixed_model_file = crm_fixed_model_file,
      crm_random_model_file = crm_random_model_file,
      crm_level_model_file = crm_level_model_file,
      crm_alpha_sd = crm_alpha_sd,
      crm_a_r = crm_a_r,
      crm_b_r = crm_b_r,
      crm_dose_values = crm_dose_values,
      crm_time_col = crm_time_col,
      crm_alpha_grid = crm_alpha_grid,
      crm_alpha_T = crm_alpha_T,
      crm_theta_prior_mean = crm_theta_prior_mean,
      crm_theta_prior_sd = crm_theta_prior_sd,
      crm_alpha_L = crm_alpha_L,
      crm_alpha_rel_tol = crm_alpha_rel_tol,
      crm_alpha_eps = crm_alpha_eps,
      crm_alpha_n_draw_prior = crm_alpha_n_draw_prior,
      crm_dose_scores = crm_dose_scores,
      crm_cumu_model_file = crm_cumu_model_file,
      crm_cumu_beta0_mean = crm_cumu_beta0_mean,
      crm_cumu_beta0_prec = crm_cumu_beta0_prec,
      crm_cumu_beta0_df = crm_cumu_beta0_df,
      crm_cumu_beta1_shape = crm_cumu_beta1_shape,
      crm_cumu_beta1_rate = crm_cumu_beta1_rate,
      crm_cumu_beta2_rate = crm_cumu_beta2_rate,
      crm_cumu_include_current = crm_cumu_include_current,
      crm_n_chains = crm_n_chains,
      crm_n_adapt = crm_n_adapt,
      crm_n_burnin = crm_n_burnin,
      crm_n_iter = crm_n_iter,
      crm_thin = crm_thin,
      cfo_method = cfo_method,
      cfo_skeleton = cfo_skeleton,
      cfo_mu = cfo_mu,
      cfo_model_file = cfo_model_file,
      cfo_sigma2_beta = cfo_sigma2_beta,
      cfo_eta = cfo_eta,
      cfo_pk_method = cfo_pk_method,
      cfo_n_mc_w = cfo_n_mc_w,
      cfo_m_use = cfo_m_use,
      cfo_use_monotone_pair = cfo_use_monotone_pair,
      cfo_restrict_to_tried = cfo_restrict_to_tried,
      cfo_n_chains = cfo_n_chains,
      cfo_n_adapt = cfo_n_adapt,
      cfo_n_burnin = cfo_n_burnin,
      cfo_n_iter = cfo_n_iter,
      cfo_thin = cfo_thin
    )

    admin <- fit$admin
    mtd <- fit$final$MTD
    final_MTD[itrial] <- mtd

    earlystop_vec[itrial] <- if (!is.null(fit$final$earlystop)) {
      as.integer(isTRUE(as.logical(fit$final$earlystop)))
    } else {
      as.integer(!is.na(mtd) && mtd == 99L)
    }

    pj_i <- fit$final$pj_iso
    if (is.null(pj_i)) pj_i <- fit$final$phat

    if (!is.null(pj_i) && length(pj_i) > 0L) {
      pj_tmp <- rep(NA_real_, ndose)
      L <- min(length(pj_i), ndose)
      pj_tmp[seq_len(L)] <- as.numeric(pj_i[seq_len(L)])
      pj_iso_mat[itrial, ] <- pj_tmp
    }

    r_hat_i <- fit$final$r_hat
    if (!is.null(r_hat_i) && length(r_hat_i) > 0L) {
      r_tmp <- rep(NA_real_, ndose)
      L <- min(length(r_hat_i), ndose)
      r_tmp[seq_len(L)] <- as.numeric(r_hat_i[seq_len(L)])
      r_hat_mat[itrial, ] <- r_tmp
    }

    if (model == "BOIN") {
      r_use_i <- fit$final$r_use
      if (!is.null(r_use_i) && length(r_use_i) > 0L) {
        r_tmp <- rep(NA_real_, ndose)
        L <- min(length(r_use_i), ndose)
        r_tmp[seq_len(L)] <- as.numeric(r_use_i[seq_len(L)])
        r_use_mat[itrial, ] <- r_tmp
      }

      r_cap_i <- fit$final$r_cap
      if (!is.null(r_cap_i) && length(r_cap_i) > 0L) {
        r_tmp <- rep(NA_real_, ndose)
        L <- min(length(r_cap_i), ndose)
        r_tmp[seq_len(L)] <- as.numeric(r_cap_i[seq_len(L)])
        r_cap_mat[itrial, ] <- r_tmp
      }
    }

    if (length(mtd) == 0L || is.na(mtd)) {
      na_count <- na_count + 1L
    } else if (mtd == 99L) {
      stop_count <- stop_count + 1L
    } else if (mtd >= 1L && mtd <= ndose) {
      sel_count[mtd] <- sel_count[mtd] + 1L
    } else {
      na_count <- na_count + 1L
    }

    if (nrow(admin) > 0L) {
      n_by_dose <- n_by_dose + tabulate(admin$dose, nbins = ndose)

      nipde_by_dose <- nipde_by_dose +
        tabulate(admin$dose[admin$type == "retreat"], nbins = ndose)

      unique_by_dose_i <- sapply(seq_len(ndose), function(d) {
        length(unique(admin$id[admin$dose == d]))
      })

      unique_n_by_dose <- unique_n_by_dose + unique_by_dose_i

      total_admin[itrial] <- nrow(admin)
      total_unique[itrial] <- length(unique(admin$id))
      duration[itrial] <- fit$final$trial_time
    } else {
      total_admin[itrial] <- 0
      total_unique[itrial] <- 0
      duration[itrial] <- NA_real_
    }

    if (store_raw) {
      raw_trials[[itrial]] <- fit
    }
  }

  pj_iso_mean <- apply(pj_iso_mat, 2, function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  })

  r_hat_mean <- apply(r_hat_mat, 2, function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  })

  r_use_mean <- apply(r_use_mat, 2, function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  })

  r_cap_mean <- apply(r_cap_mat, 2, function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  })

  earlystop_count <- sum(earlystop_vec, na.rm = TRUE)

  out <- list(
    target = target,
    p.true = p.true,
    p.true_ipde = p.true_ipde,
    ipde_alpha = ipde_alpha,
    ntrial = ntrial,
    ndose = ndose,
    model = model,
    ipde_design = ipde_design,
    continuous_enrollment = continuous_enrollment,
    new_pat_first = new_pat_first,
    d.cap = d.cap,
    dose_cap = dose_cap,
    n_eval_escalate = n_eval_escalate,
    decision_method = decision_method,
    mtd_method = mtd_method,
    r_carry = r_carry,
    r_estimator = if (model == "BOIN") r_estimator else NA_character_,
    r_adaptive_prior = if (model == "BOIN" && r_estimator == "r_adaptive") r_adaptive_prior else c(NA_real_, NA_real_),
    p_adaptive_prior = if (model == "BOIN" && r_estimator == "r_adaptive") p_adaptive_prior else c(NA_real_, NA_real_),
    r_adaptive_max = if (model == "BOIN" && r_estimator == "r_adaptive") r_adaptive_max else NA_real_,
    r_adaptive_plug_in = if (model == "BOIN" && r_estimator == "r_adaptive") r_adaptive_plug_in else NA_character_,
    r_adaptive_rel_tol = if (model == "BOIN" && r_estimator == "r_adaptive") r_adaptive_rel_tol else NA_real_,
    crm_r_model = if (model == "CRM") crm_r_model else NA_character_,
    crm_skeleton = if (model == "CRM") crm_skeleton else NULL,
    crm_alpha_sd = if (model == "CRM") crm_alpha_sd else NA_real_,
    crm_a_r = if (model == "CRM") crm_a_r else NA_real_,
    crm_b_r = if (model == "CRM") crm_b_r else NA_real_,
    crm_model_file = if (model == "CRM") crm_model_file else NA_character_,
    crm_dose_values = if (model == "CRM") crm_dose_values else NULL,
    crm_dose_scores = if (model == "CRM") crm_dose_scores else NULL,
    cfo_method = if (model == "CFO") cfo_method else NA_character_,
    cfo_skeleton = if (model == "CFO") cfo_skeleton else NULL,
    cfo_mu = if (model == "CFO") cfo_mu else NULL,
    cfo_model_file = if (model == "CFO" && cfo_method == "pride") cfo_model_file else NA_character_,
    cfo_sigma2_beta = if (model == "CFO") cfo_sigma2_beta else NA_real_,
    cfo_eta = if (model == "CFO") cfo_eta else NA_real_,
    cfo_pk_method = if (model == "CFO") cfo_pk_method else NA_character_,
    restrict_to_tried = restrict_to_tried,
    restrict_to_target = restrict_to_target,
    sel_count = sel_count,
    stop_count = stop_count,
    na_count = na_count,
    selection_pct = 100 * sel_count / ntrial,
    early_stop_pct = 100 * earlystop_count / ntrial,
    na_pct = 100 * na_count / ntrial,
    earlystop = earlystop_vec,
    earlystop_count = earlystop_count,
    final_MTD = final_MTD,
    pj_iso_by_trial = pj_iso_mat,
    pj_iso_mean = pj_iso_mean,
    r_hat_by_trial = r_hat_mat,
    r_hat_mean = r_hat_mean,
    r_use_by_trial = if (model == "BOIN") r_use_mat else NULL,
    r_use_mean = if (model == "BOIN") r_use_mean else NULL,
    r_cap_by_trial = if (model == "BOIN") r_cap_mat else NULL,
    r_cap_mean = if (model == "BOIN") r_cap_mean else NULL,
    n_by_dose = n_by_dose / ntrial,
    unique_n_by_dose = unique_n_by_dose / ntrial,
    nipde_by_dose = nipde_by_dose / ntrial,
    total_admin_mean = mean(total_admin, na.rm = TRUE),
    total_unique_mean = mean(total_unique, na.rm = TRUE),
    duration_mean = mean(duration, na.rm = TRUE),
    total_admin = total_admin,
    total_unique = total_unique,
    duration = duration,
    raw_trials = raw_trials
  )

  class(out) <- "oc_AIDE"
  out
}
