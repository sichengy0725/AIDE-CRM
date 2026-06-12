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
  
  # --- DLT time generation ---
  dlt_dist  = 2,
  dlt_alpha = 0.5,
  
  # --- truth ---
  p_true,
  p_ipde = p_true,
  seed = NULL,
  verbose = FALSE,
  
  TARGET = 0.3,
  cutoff = 0.95,
  
  model = c("BOIN", "CRM"),
  ipde_design = 2,
  d.cap = 100,
  day_obs = 0,
  
  # --- escalation gate ---
  dose_cap = 3L,
  
  # --- BOIN decision / final selection method ---
  decision_method = c("boin", "approx1", "approx2"),
  mtd_method = NULL,
  r_carry = 0.1,
  
  # --- CRM discount-model options ---
  crm_r_model = c("fixed", "random"),
  crm_skeleton = NULL,
  crm_alpha_sd = 2,
  crm_a_r = 1,
  crm_b_r = 9,
  crm_fixed_model_file = "fix_CRM.bug",
  crm_random_model_file = "random_CRM.bug",
  crm_n_chains = 2,
  crm_n_adapt = 500,
  crm_n_burnin = 500,
  crm_n_iter = 2000,
  crm_thin = 1
) {
  if (!is.null(seed)) set.seed(seed)
  
  model <- match.arg(model)
  
  if (model == "BOIN") {
    decision_method <- match.arg(decision_method)
    
    if (is.null(mtd_method)) {
      mtd_method <- decision_method
    }
    
    mtd_method <- match.arg(
      mtd_method,
      choices = c("boin", "approx1", "approx2")
    )
  }
  
  if (model == "CRM") {
    crm_r_model <- match.arg(crm_r_model)
    
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
    
    decision_method <- paste0("crm_", crm_r_model)
    mtd_method <- decision_method
  }
  
  if (length(p_ipde) != length(p_true)) {
    stop("p_ipde must have the same length as p_true.")
  }
  
  if (N_pat < C) {
    stop("Need N_pat >= C.")
  }
  
  if (length(arrival_rate) != 1L || is.na(arrival_rate) || arrival_rate <= 0) {
    stop("arrival_rate must be a single positive value.")
  }
  
  K <- length(p_true)
  ndose <- K
  
  ## ---------------- helpers ----------------
  
  window <- T_assess
  
  draw_dlt_time <- function(pi_val, t_start) {
    ## Instant outcome setting: no time-to-DLT generation.
    ## Binary DLT is generated immediately, but evaluability is at t_start + window.
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
  
  eligible_ipde_ids <- function(admin, next_dose, cycle_max, design) {
    st <- get_patient_state(admin)
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
  
  recruit_new_patients <- function(n_needed, earliest_time) {
    if (n_needed <= 0L) {
      return(data.frame(
        id = integer(0),
        t_arrival = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    n_left <- N_pat - next_new_idx + 1L
    
    if (n_left <= 0L) {
      return(data.frame(
        id = integer(0),
        t_arrival = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    n_take <- min(as.integer(n_needed), as.integer(n_left))
    ids <- seq.int(next_new_idx, length.out = n_take)
    
    inter <- rexp(n_take, rate = arrival_rate)
    arr <- as.numeric(earliest_time + cumsum(inter))
    
    next_new_idx <<- next_new_idx + n_take
    
    data.frame(
      id = ids,
      t_arrival = arr,
      stringsAsFactors = FALSE
    )
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
  
  ## ---------------- recruitment state ----------------
  
  next_new_idx <- 1L
  
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
  
  first_new <- recruit_new_patients(C, earliest_time = t0)
  
  if (nrow(first_new) < C) {
    stop("Not enough unique patients to enroll cohort 1.")
  }
  
  first_ids <- first_new$id
  t_start1 <- max(t0, max(first_new$t_arrival))
  
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
      ", decision_time = ", round(t_decision, 2)
    )
  }
  
  ## ---------------- main loop ----------------
  
  repeat {
    if (stop_trial) break
    if (nrow(admin) >= Nmax_eff) break
    
    dat_dec <- admin[
      admin$t_eval <= t_decision,
      c("id", "dose", "y", "type"),
      drop = FALSE
    ]
    
    n <- tabulate(dat_dec$dose, nbins = K)
    y <- tabulate(dat_dec$dose[dat_dec$y == 1L], nbins = K)
    
    n_curr <- n[curr_dose]
    y_curr <- y[curr_dose]
    
    ## ---------------- BOIN elimination only ----------------
    ## CRM early stopping is NOT checked here.
    ## For CRM, dose-1 overdose stopping is handled inside crm_move().
    
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
          r_carry = r_carry,
          
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
          cutoff.eli = cutoff,
          skeleton = crm_skeleton,
          
          r_model = crm_r_model,
          r_carry = r_carry,
          a_r = crm_a_r,
          b_r = crm_b_r,
          
          alpha_sd = crm_alpha_sd,
          
          fixed_model_file = crm_fixed_model_file,
          random_model_file = crm_random_model_file,
          
          n_chains = crm_n_chains,
          n_adapt = crm_n_adapt,
          n_burnin = crm_n_burnin,
          n_iter = crm_n_iter,
          thin = crm_thin,
          seed = seed,
          
          elimi = elimi,
          n_trt_curr = n_trt_curr,
          dose_cap = dose_cap,
          no_skip = TRUE
        )
      }
    }
    
    ## CRM early stopping is consumed here only from crm_move().
    ## No separate posterior overdose calculation is performed in this trial code.
    if (isTRUE(move$stop_trial)) {
      earlystop <- as.integer(isTRUE(as.logical(move$earlystop)))
      if (!is.null(move$eliminated)) elimi <- move$eliminated
      break
    }
    
    if (!is.null(move$eliminated)) {
      elimi <- move$eliminated
    }
    
    next_dose <- move$next_dose
    
    if (n[curr_dose] >= d.cap && next_dose == curr_dose) {
      break
    }
    
    ## ---------------- form next cohort ----------------
    
    cohort_id <- cohort_id + 1L
    t_start <- t_decision
    
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
      design = ipde_design
    )
    
    if (length(cand) > 0L) {
      ret_ids <- head(cand, min(C, length(cand)))
    }
    
    need_new <- C - length(ret_ids)
    
    if (need_new > 0L) {
      new_patients <- recruit_new_patients(
        n_needed = need_new,
        earliest_time = t_decision
      )
      
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
      t_start <- max(t_start, max(new_patients$t_arrival))
    }
    
    ## ---------------- append new-patient outcomes ----------------
    
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
    
    ## ---------------- append IPDE / retreat outcomes ----------------
    
    if (length(ret_ids) > 0) {
      st <- get_patient_state(admin)
      
      for (pid in ret_ids) {
        if (nrow(admin) >= Nmax_eff) break
        
        last <- st[st$id == pid, , drop = FALSE]
        if (nrow(last) != 1L) next
        
        dlt <- draw_dlt_time(
          pi_val = p_ipde[next_dose],
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
        model = model
      )
    ))
  }
  
  t_end <- max(admin$t_eval, na.rm = TRUE)
  
  dat_final <- admin[
    admin$t_eval <= t_end,
    c("id", "dose", "y", "type"),
    drop = FALSE
  ]
  
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
  
  trial_time <- max(admin$t_eval, na.rm = TRUE) -
    min(admin$t_arrival, na.rm = TRUE)
  
  phat <- rep(NA_real_, K)
  r_hat <- if (model == "CRM" && crm_r_model == "fixed") r_carry else NA_real_
  
  if (earlystop == 1L) {
    final_dose <- 99L
    final_fit <- list(
      phat = rep(NA_real_, K),
      pj_iso = rep(NA_real_, K),
      eliminated = elimi,
      r_hat = r_hat
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
        y_new = y_new_final,
        n_new = n_new_final,
        y_recycle = y_recycle_final,
        n_recycle = n_recycle_final
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
        
        fixed_model_file = crm_fixed_model_file,
        random_model_file = crm_random_model_file,
        
        n_chains = crm_n_chains,
        n_adapt = crm_n_adapt,
        n_burnin = crm_n_burnin,
        n_iter = max(5000, crm_n_iter),
        thin = crm_thin,
        seed = seed,
        
        restrict_to_tried = TRUE
      )
    }
    
    final_dose <- final_fit$MTD
    phat[seq_along(final_fit$phat)] <- final_fit$phat
    if (!is.null(final_fit$r_hat)) r_hat <- final_fit$r_hat
  }
  
  final_elimi <- if (!is.null(final_fit$eliminated)) final_fit$eliminated else elimi
  
  list(
    admin = admin,
    waiting = waiting,
    final = list(
      t_end = t_end,
      MTD = final_dose,
      trial_time = trial_time,
      eliminated = final_elimi,
      earlystop = earlystop,
      phat = phat,
      pj_iso = final_fit$pj_iso,
      decision_method = decision_method,
      mtd_method = mtd_method,
      r_carry = r_carry,
      r_model = if (model == "CRM") crm_r_model else NA_character_,
      r_hat = r_hat,
      model = model
    )
  )
}