## ============================================================
## AIDE Phase I/II designs (non-TITE)
##
## Two allocation options are available:
##   1. "two_stage": toxicity-only AIDE until the cumulative Stage I sample
##      reaches N_s1 and the next toxicity recommendation is stay.  Stage II
##      then applies the same one-stage allocation rule after every cohort.
##   2. "one_stage": the response floor and the toxicity-model MTD define
##      the utility candidate interval after each cohort.
##
## Efficacy supports four one-cycle carryover models plus a shared arbitrary-
## cycle additive model. Toxicity uses the paired shared-r, shared-alpha, or
## arbitrary-cycle additive CRM model. Both binary toxicity and
## binary efficacy outcomes are generated with gen.tite;
## the IPDE true-outcome probabilities retain their selected carryover model
## for their respective endpoints. The non-conjugate efficacy posterior is
## fitted with JAGS MCMC draws.
## Patients accrue continuously, wait when a cohort is under evaluation, and
## are allocated only after the preceding cohort has completed assessment.
## This file intentionally does not use TITE weights or partial follow-up
## decisions.
## ============================================================

if (!exists("boin_move", mode = "function") ||
    !exists("select.mtd", mode = "function") ||
    !exists("gen.tite", mode = "function")) {
  source("AIDE_BOIN_helper.R")
}

if (!exists("crm_move", mode = "function") ||
    !exists("select.mtd.crm", mode = "function")) {
  source("AIDE_CRM_helper_modified.R")
}

aide_phase12_validate_beta_prior <- function(prior, name) {
  if (length(prior) != 2L || any(!is.finite(prior)) || any(prior <= 0)) {
    stop(name, " must contain two positive finite beta-prior parameters.")
  }
  as.numeric(prior)
}

aide_phase12_validate_logical_flag <- function(value, name) {
  if (length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.")
  }
  as.logical(value)
}

## The high-level option keeps toxicity and efficacy carryover assumptions
## aligned. Toxicity has shared carryover models; efficacy adds dose-specific
## one-cycle variants for the corresponding mechanisms.
aide_phase12_carryover_model_map <- function(carryover_model) {
  carryover_model <- match.arg(
    carryover_model,
    c(
      "discount_shared", "discount_dose_specific",
       "additive_shared", "additive_dose_specific", "multicycle_additive"
    )
  )
  switch(
    carryover_model,
    discount_shared = list(
      toxicity_model = "discount_shared",
      crm_r_model = "random",
      efficacy_model = "shared_carryover"
    ),
    discount_dose_specific = list(
      toxicity_model = "discount_shared",
      crm_r_model = "random",
      efficacy_model = "dose_specific_carryover"
    ),
    additive_shared = list(
      toxicity_model = "additive_shared",
      crm_r_model = "previous_dose",
      efficacy_model = "previous_dose_additive"
    ),
    additive_dose_specific = list(
      toxicity_model = "additive_shared",
      crm_r_model = "previous_dose",
      efficacy_model = "dose_specific_previous_dose_additive"
    ),
    multicycle_additive = list(
      toxicity_model = "multicycle_additive",
      crm_r_model = "multicycle_additive",
      efficacy_model = "multicycle_additive"
    )
  )
}

aide_phase12_is_additive_efficacy_model <- function(efficacy_model) {
  efficacy_model %in% c(
    "previous_dose_additive", "dose_specific_previous_dose_additive",
    "multicycle_additive"
  )
}

aide_phase12_is_multicycle_efficacy_model <- function(efficacy_model) {
  identical(efficacy_model, "multicycle_additive")
}

## The arbitrary-cycle state is deliberately uncapped. Only the probability
## used to generate or score the current administration is capped at one.
aide_phase12_multicycle_next_state <- function(regular_probability,
                                                alpha,
                                                previous_state = 0) {
  if (length(regular_probability) != 1L || !is.finite(regular_probability) ||
      regular_probability < 0 || regular_probability > 1 ||
      length(alpha) != 1L || !is.finite(alpha) || alpha < 0 ||
      length(previous_state) != 1L || !is.finite(previous_state) ||
      previous_state < 0) {
    stop("multicycle state inputs must be finite non-negative values, with regular_probability in [0, 1].")
  }
  state <- as.numeric(regular_probability) + as.numeric(alpha) * as.numeric(previous_state)
  list(state = state, probability = min(1, state))
}

## Build the JAGS-oriented full-history representation. The returned rows are
## ordered within patient so the previous state is always already defined.
## previous_state_index = 1 denotes the explicit zero-state sentinel; the
## state generated for row r resides at index r + 1.
aide_phase12_multicycle_history_data <- function(admin, ndose) {
  if (is.null(admin)) admin <- data.frame(dose = integer(0), eff = integer(0))
  if (!all(c("dose", "eff") %in% names(admin))) {
    stop("admin must contain dose and eff for the multicycle efficacy model.")
  }
  if (nrow(admin) == 0L) {
    return(data.frame(
      admin_id = integer(0), patient_id = integer(0), cycle = integer(0),
      dose = integer(0), efficacy = integer(0), ipde = integer(0),
      previous_admin_id = integer(0), previous_state_index = integer(0)
    ))
  }
  if (any(!is.finite(admin$dose)) || any(admin$dose < 1L | admin$dose > ndose) ||
      any(!admin$eff %in% c(0L, 1L))) {
    stop("admin contains invalid dose or efficacy values.")
  }
  patient_id <- if ("id" %in% names(admin)) admin$id else seq_len(nrow(admin))
  if (any(is.na(patient_id))) stop("admin$id cannot be missing for multicycle carryover.")
  cycle <- if ("ncycle" %in% names(admin)) admin$ncycle else if ("cycle" %in% names(admin)) admin$cycle else ave(seq_len(nrow(admin)), patient_id, FUN = seq_along)
  if (any(!is.finite(cycle)) || any(cycle < 1L)) stop("admin cycle values must be positive and finite.")
  time <- if ("t_start" %in% names(admin)) admin$t_start else if ("t_arrival" %in% names(admin)) admin$t_arrival else seq_len(nrow(admin))
  if (any(!is.finite(time))) stop("admin administration times must be finite.")
  admin_id <- if ("row_id" %in% names(admin)) admin$row_id else seq_len(nrow(admin))
  if (any(is.na(admin_id)) || anyDuplicated(admin_id)) stop("admin row identifiers must be unique and non-missing.")
  order_index <- order(patient_id, cycle, time, admin_id)
  patient_id <- patient_id[order_index]
  admin_id <- as.integer(admin_id[order_index])
  cycle <- as.integer(cycle[order_index])
  dose <- as.integer(admin$dose[order_index])
  efficacy <- as.integer(admin$eff[order_index])
  ipde <- as.integer(if ("type" %in% names(admin)) admin$type[order_index] == "retreat" else FALSE)
  previous_admin_id <- rep.int(0L, length(admin_id))
  previous_state_index <- rep.int(1L, length(admin_id))
  for (rows in split(seq_along(admin_id), as.character(patient_id))) {
    if (length(rows) > 1L) {
      previous_admin_id[rows[-1L]] <- admin_id[rows[-length(rows)]]
      previous_state_index[rows[-1L]] <- rows[-length(rows)] + 1L
    }
  }
  data.frame(
    admin_id = admin_id, patient_id = as.integer(patient_id), cycle = cycle,
    dose = dose, efficacy = efficacy, ipde = ipde,
    previous_admin_id = previous_admin_id,
    previous_state_index = as.integer(previous_state_index)
  )
}

## The multicycle likelihood uses only the terminal administration for each
## patient. The history is ordered by patient and cycle before this map is
## constructed, so the selected row is the current/latest cycle.
aide_phase12_multicycle_latest_admin_index <- function(history_data) {
  if (is.null(history_data) || nrow(history_data) == 0L) return(integer(0))
  if (!"patient_id" %in% names(history_data)) {
    stop("history_data must contain patient_id for multicycle likelihood mapping.")
  }
  as.integer(which(!duplicated(history_data$patient_id, fromLast = TRUE)))
}

aide_phase12_dose_threshold_reached <- function(n_by_dose, threshold) {
  n_by_dose <- as.numeric(n_by_dose)
  if (length(threshold) != 1L || !is.finite(threshold) || threshold < 1L) {
    stop("threshold must be a positive finite value.")
  }
  any(is.finite(n_by_dose) & n_by_dose >= threshold)
}

## Determine whether a patient most recently treated at previous_dose can be
## recycled at next_dose.  A same-dose recycle is permitted for both IPDE
## designs.  Without flexible_ipde, movement remains upward-only; with it,
## downward recycling follows the selected any-dose/adjacent-dose rule.
aide_phase12_recycle_dose_eligible <- function(previous_dose,
                                                next_dose,
                                                ipde_design,
                                                flexible_ipde) {
  previous_dose <- as.numeric(previous_dose)
  if (length(next_dose) != 1L || !is.finite(next_dose) ||
      next_dose < 1L || next_dose != as.integer(next_dose)) {
    stop("next_dose must be a positive integer.")
  }
  if (length(ipde_design) != 1L || !ipde_design %in% c(1L, 2L)) {
    stop("ipde_design must be 1 (any dose) or 2 (adjacent dose).")
  }
  flexible_ipde <- aide_phase12_validate_logical_flag(
    flexible_ipde, "flexible_ipde"
  )
  if (any(!is.finite(previous_dose)) || any(previous_dose < 1L) ||
      any(previous_dose != as.integer(previous_dose))) {
    stop("previous_dose must contain positive integer dose indices.")
  }
  previous_dose <- as.integer(previous_dose)
  next_dose <- as.integer(next_dose)
  if (flexible_ipde && ipde_design == 1L) {
    rep(TRUE, length(previous_dose))
  } else if (flexible_ipde && ipde_design == 2L) {
    abs(previous_dose - next_dose) <= 1L
  } else if (ipde_design == 1L) {
    previous_dose <= next_dose
  } else {
    previous_dose >= next_dose - 1L & previous_dose <= next_dose
  }
}

## True recycled-patient efficacy. At a proposed dose d2 for a patient most
## recently treated at d1, this is min(1, p2 + alpha * p1), where p2 is the
## recycled-efficacy base at d2 and p1 is regular efficacy at d1.
aide_phase12_ipde_efficacy_probability <- function(e_regular,
                                                    e_ipde_base,
                                                    previous_dose,
                                                    current_dose,
                                                    alpha = 0) {
  e_regular <- as.numeric(e_regular)
  e_ipde_base <- as.numeric(e_ipde_base)
  ndose <- length(e_regular)
  if (ndose < 1L || length(e_ipde_base) != ndose ||
      any(!is.finite(c(e_regular, e_ipde_base))) ||
      any(c(e_regular, e_ipde_base) < 0 | c(e_regular, e_ipde_base) > 1)) {
    stop("e_regular and e_ipde_base must be equal-length probability vectors in [0, 1].")
  }
  if (length(previous_dose) != 1L || !is.finite(previous_dose) ||
      previous_dose < 1L || previous_dose > ndose ||
      previous_dose != as.integer(previous_dose) ||
      length(current_dose) != 1L || !is.finite(current_dose) ||
      current_dose < 1L || current_dose > ndose ||
      current_dose != as.integer(current_dose)) {
    stop("previous_dose and current_dose must be valid dose indices.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
    stop("alpha must be a single finite non-negative value.")
  }
  min(
    1,
    e_ipde_base[as.integer(current_dose)] +
      as.numeric(alpha) * e_regular[as.integer(previous_dose)]
  )
}

## True recycled-patient toxicity. At a proposed dose d2 for a patient most
## recently treated at d1, this is min(1, p2 + alpha * p1), where p2 is the
## recycled-toxicity base at d2 and p1 is regular toxicity at d1.
aide_phase12_ipde_toxicity_probability <- function(p_regular,
                                                    p_ipde_base,
                                                    previous_dose,
                                                    current_dose,
                                                    alpha = 0) {
  p_regular <- as.numeric(p_regular)
  p_ipde_base <- as.numeric(p_ipde_base)
  ndose <- length(p_regular)
  if (ndose < 1L || length(p_ipde_base) != ndose ||
      any(!is.finite(c(p_regular, p_ipde_base))) ||
      any(c(p_regular, p_ipde_base) < 0 | c(p_regular, p_ipde_base) > 1)) {
    stop("p_regular and p_ipde_base must be equal-length probability vectors in [0, 1].")
  }
  if (length(previous_dose) != 1L || !is.finite(previous_dose) ||
      previous_dose < 1L || previous_dose > ndose ||
      previous_dose != as.integer(previous_dose) ||
      length(current_dose) != 1L || !is.finite(current_dose) ||
      current_dose < 1L || current_dose > ndose ||
      current_dose != as.integer(current_dose)) {
    stop("previous_dose and current_dose must be valid dose indices.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
    stop("alpha must be a single finite non-negative value.")
  }
  min(
    1,
    p_ipde_base[as.integer(current_dose)] +
      as.numeric(alpha) * p_regular[as.integer(previous_dose)]
  )
}

## Posterior safety gate for a recycled assignment under the random-r CRM.
## The CRM JAGS model supplies p_j and r, so this is evaluated draw by draw
## as Pr{r + (1-r) p_j > phi}.
aide_phase12_random_crm_recycle_toxicity_gate <- function(post,
                                                           next_dose,
                                                           phi,
                                                           cutoff) {
  post <- as.matrix(post)
  if (nrow(post) < 1L || !"r" %in% colnames(post)) {
    stop("Random-CRM posterior samples must contain at least one draw of r.")
  }
  if (length(next_dose) != 1L || !is.finite(next_dose) ||
      next_dose < 1L || next_dose != as.integer(next_dose)) {
    stop("next_dose must be a positive integer.")
  }
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0 || phi >= 1) {
    stop("phi must be a scalar in (0, 1).")
  }
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff >= 1) {
    stop("cutoff must be a scalar in (0, 1).")
  }

  p_column <- paste0("p[", as.integer(next_dose), "]")
  if (!p_column %in% colnames(post)) {
    stop("Random-CRM posterior samples do not contain ", p_column, ".")
  }
  r_draw <- as.numeric(post[, "r"])
  p_draw <- as.numeric(post[, p_column])
  if (any(!is.finite(r_draw) | r_draw < 0 | r_draw > 1) ||
      any(!is.finite(p_draw) | p_draw < 0 | p_draw > 1)) {
    stop("Random-CRM posterior samples contain invalid probability draws.")
  }
  theta_draw <- r_draw + (1 - r_draw) * p_draw
  probability_over_phi <- mean(theta_draw > phi)
  list(
    allowed = probability_over_phi < cutoff,
    probability_over_phi = probability_over_phi,
    theta_posterior_mean = mean(theta_draw),
    phi = phi,
    cutoff = cutoff,
    next_dose = as.integer(next_dose)
  )
}

## Expand a shared Beta(a, b) prior or a dose-specific ndose-by-2 prior
## matrix. The latter implements the prespecified independent priors in the
## dose-specific IPDE efficacy model.
aide_phase12_expand_beta_prior <- function(prior, ndose, name) {
  if (length(ndose) != 1L || !is.finite(ndose) || ndose < 1L ||
      ndose != as.integer(ndose)) {
    stop("ndose must be a positive integer.")
  }
  if (is.matrix(prior) || is.data.frame(prior)) {
    out <- as.matrix(prior)
    if (!identical(dim(out), c(as.integer(ndose), 2L))) {
      stop(name, " must be length 2 or an ndose-by-2 matrix.")
    }
  } else {
    prior <- as.numeric(prior)
    if (length(prior) == 2L) {
      out <- matrix(rep(prior, each = as.integer(ndose)), ncol = 2L)
    } else if (length(prior) == 2L * as.integer(ndose)) {
      out <- matrix(prior, ncol = 2L, byrow = TRUE)
    } else {
      stop(name, " must be length 2 or contain two values per dose.")
    }
  }
  storage.mode(out) <- "double"
  if (any(!is.finite(out)) || any(out <= 0)) {
    stop(name, " must contain positive finite Beta-prior parameters.")
  }
  colnames(out) <- c("a", "b")
  out
}

## Load the common JAGS implementation lazily so sourcing this trial-design
## file does not require rjags until an efficacy model is actually fitted.
aide_phase12_load_efficacy_jags <- function(require_standard_fitter = TRUE) {
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop(
      "The phase I/II efficacy model requires rjags and a working JAGS installation. ",
      "Install them before running simulate_AIDE_phase_I_II()."
    )
  }
  if (isTRUE(require_standard_fitter) &&
      !exists("fit_beta_binomial_efficacy", mode = "function")) {
    source("generate_dose_specific_beta_binomial_efficacy_comparison.R")
  }
  if (isTRUE(require_standard_fitter) &&
      !exists("fit_beta_binomial_efficacy", mode = "function")) {
    stop("Could not load fit_beta_binomial_efficacy().")
  }
  invisible(NULL)
}

## Produce the individual-level data needed when IPDE efficacy depends on a
## patient's immediately preceding dose. The input order is used as the final
## within-patient tie breaker if no cycle/time information is supplied.
aide_phase12_previous_dose_efficacy_data <- function(admin, ndose) {
  if (is.null(admin)) admin <- data.frame(dose = integer(0), eff = integer(0))
  if (!all(c("dose", "eff") %in% names(admin))) {
    stop("admin must contain dose and eff for the efficacy model.")
  }
  if (nrow(admin) == 0L) {
    return(data.frame(
      dose = integer(0), efficacy = integer(0), ipde = integer(0),
      previous_dose = integer(0)
    ))
  }
  if (any(!is.finite(admin$dose)) || any(admin$dose < 1L | admin$dose > ndose) ||
      any(!admin$eff %in% c(0L, 1L))) {
    stop("admin contains invalid dose or efficacy values.")
  }

  ipde <- if ("type" %in% names(admin)) admin$type == "retreat" else rep(FALSE, nrow(admin))
  if (any(is.na(ipde))) stop("admin$type cannot be missing.")
  if (any(ipde) && !"id" %in% names(admin)) {
    stop("The previous_dose_additive efficacy model requires admin$id for recycled/IPDE administrations.")
  }
  id <- if ("id" %in% names(admin)) admin$id else seq_len(nrow(admin))
  if (any(is.na(id))) stop("admin$id cannot be missing.")
  cycle <- if ("ncycle" %in% names(admin)) {
    admin$ncycle
  } else if ("cycle" %in% names(admin)) {
    admin$cycle
  } else {
    ave(seq_len(nrow(admin)), id, FUN = seq_along)
  }
  if (any(!is.finite(cycle))) stop("admin cycle values must be finite.")
  time <- if ("t_start" %in% names(admin)) {
    admin$t_start
  } else if ("t_arrival" %in% names(admin)) {
    admin$t_arrival
  } else {
    seq_len(nrow(admin))
  }
  if (any(!is.finite(time))) stop("admin time values must be finite.")

  order_index <- order(id, cycle, time, seq_len(nrow(admin)))
  dose <- as.integer(admin$dose[order_index])
  efficacy <- as.integer(admin$eff[order_index])
  ipde <- as.integer(ipde[order_index])
  id <- id[order_index]
  previous_dose <- rep.int(1L, length(dose))
  for (rows in split(seq_along(dose), as.character(id))) {
    if (ipde[rows[1L]] == 1L) {
      stop("Each recycled/IPDE efficacy administration must follow an earlier administration for the same patient.")
    }
    if (length(rows) > 1L) {
      for (k in 2:length(rows)) {
        if (ipde[rows[k]] == 1L) previous_dose[rows[k]] <- dose[rows[k - 1L]]
      }
    }
  }
  data.frame(
    dose = dose,
    efficacy = efficacy,
    ipde = ipde,
    previous_dose = as.integer(previous_dose)
  )
}

## Fit independent dose-specific beta-binomial efficacy probabilities. The
## default model uses a dose-specific IPDE carryover probability. The optional
## previous_dose_additive model uses p[current] + alpha * p[previous] for each
## recycled administration, requiring individual-level data.
aide_phase12_efficacy_posterior <- function(admin,
                                             ndose,
                                             efficacy_prior = c(1, 1),
                                             efficacy_carryover_prior = c(1, 9),
                                             efficacy_model = c(
                                                "shared_carryover",
                                                "dose_specific_carryover",
                                                "previous_dose_additive",
                                                "dose_specific_previous_dose_additive",
                                                "multicycle_additive"
                                             ),
                                             efficacy_additive_alpha_prior = c(1, 9),
                                             model_file = NULL,
                                             n_chains = 3L,
                                             n_adapt = 1000L,
                                             n_burnin = 1000L,
                                             n_iter = 4000L,
                                             thin = 2L,
                                             cache = NULL,
                                             threshold = NULL) {
  efficacy_model <- match.arg(efficacy_model)
  efficacy_prior <- aide_phase12_expand_beta_prior(
    efficacy_prior, ndose, "efficacy_prior"
  )
  efficacy_carryover_prior <- aide_phase12_expand_beta_prior(
    efficacy_carryover_prior, ndose, "efficacy_carryover_prior"
  )
  efficacy_additive_alpha_prior <- aide_phase12_expand_beta_prior(
    efficacy_additive_alpha_prior, ndose, "efficacy_additive_alpha_prior"
  )
  is_additive <- aide_phase12_is_additive_efficacy_model(efficacy_model)
  is_dose_specific_additive <- identical(
    efficacy_model, "dose_specific_previous_dose_additive"
  )
  is_multicycle_additive <- aide_phase12_is_multicycle_efficacy_model(efficacy_model)
  if (is_additive && !is_dose_specific_additive &&
      (any(efficacy_additive_alpha_prior[, "a"] !=
           efficacy_additive_alpha_prior[1L, "a"]) ||
       any(efficacy_additive_alpha_prior[, "b"] !=
           efficacy_additive_alpha_prior[1L, "b"]))) {
    stop(
      "For previous_dose_additive, efficacy_additive_alpha_prior must be ",
      "shared (a length-two vector or identical values at every dose)."
    )
  }
  if (!is.null(threshold) &&
      (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 || threshold >= 1)) {
    stop("threshold must be NULL or a scalar in (0, 1).")
  }
  if (!is.null(cache) && !is.environment(cache)) {
    stop("cache must be NULL or an environment.")
  }

  if (is.null(admin)) admin <- data.frame(dose = integer(0), eff = integer(0))
  if (nrow(admin) > 0L) {
    if (!all(c("dose", "eff") %in% names(admin))) {
      stop("admin must contain dose and eff for the efficacy model.")
    }
    if (any(!is.finite(admin$dose)) || any(admin$dose < 1L | admin$dose > ndose) ||
        any(!admin$eff %in% c(0L, 1L))) {
      stop("admin contains invalid dose or efficacy values.")
    }
    is_ipde <- if ("type" %in% names(admin)) {
      admin$type == "retreat"
    } else {
      rep(FALSE, nrow(admin))
    }
    if (any(is.na(is_ipde))) stop("admin$type cannot be missing.")
  } else {
    is_ipde <- logical(0)
  }

  dose_data <- data.frame(
    dose = as.integer(admin$dose),
    group = ifelse(is_ipde, "ipde", "regular"),
    efficacy = as.integer(admin$eff),
    stringsAsFactors = FALSE
  )
  previous_dose_data <- if (is_multicycle_additive) {
    aide_phase12_multicycle_history_data(admin, ndose)
  } else if (is_additive) {
    aide_phase12_previous_dose_efficacy_data(admin, ndose)
  } else {
    NULL
  }
  latest_admin_index <- if (is_multicycle_additive) {
    aide_phase12_multicycle_latest_admin_index(previous_dose_data)
  } else {
    NULL
  }
  model_key <- if (is.null(model_file)) "<default>" else as.character(model_file)
  cache_key <- paste(
    efficacy_model, ndose,
    paste(c(efficacy_prior), collapse = ","),
    paste(c(efficacy_carryover_prior), collapse = ","),
    paste(c(efficacy_additive_alpha_prior), collapse = ","),
    model_key, n_chains, n_adapt, n_burnin, n_iter, thin,
    if (is_additive) {
      if (is_multicycle_additive) {
        paste(previous_dose_data$admin_id, previous_dose_data$patient_id,
              previous_dose_data$cycle, previous_dose_data$dose,
              previous_dose_data$efficacy,
              previous_dose_data$previous_state_index,
              sep = ":", collapse = ";")
      } else {
        paste(previous_dose_data$dose, previous_dose_data$efficacy,
              previous_dose_data$ipde, previous_dose_data$previous_dose,
              sep = ":", collapse = ";")
      }
    } else {
      paste(dose_data$dose, dose_data$group, dose_data$efficacy,
            sep = ":", collapse = ";")
    },
    sep = "|"
  )
  cache_name <- paste0("efficacy_", cache_key)
  base <- if (!is.null(cache) && exists(cache_name, envir = cache, inherits = FALSE)) {
    get(cache_name, envir = cache, inherits = FALSE)
  } else {
    if (efficacy_model %in% c("dose_specific_carryover", "shared_carryover")) {
      aide_phase12_load_efficacy_jags(require_standard_fitter = TRUE)
      fit <- fit_beta_binomial_efficacy(
        dose_data = dose_data,
        a_r = efficacy_prior[, "a"],
        b_r = efficacy_prior[, "b"],
        a_carry = efficacy_carryover_prior[, "a"],
        b_carry = efficacy_carryover_prior[, "b"],
        prior_type = if (efficacy_model == "dose_specific_carryover") {
          "dose_specific"
        } else {
          "shared_carryover"
        },
        model_file = model_file,
        n_chains = n_chains,
        n_adapt = n_adapt,
        n_burnin = n_burnin,
        n_iter = n_iter,
        thin = thin,
        ndose = ndose
      )
      draws <- as.matrix(fit$samples)
      parameter_draws <- function(parameter) {
        columns <- paste0(parameter, "[", seq_len(ndose), "]")
        if (!all(columns %in% colnames(draws))) {
          stop("The JAGS efficacy fit did not return all ", parameter, " posterior draws.")
        }
        draws[, columns, drop = FALSE]
      }
      carryover_draws <- if (efficacy_model == "dose_specific_carryover") {
        parameter_draws("r_e")
      } else {
        if (!"r_e" %in% colnames(draws)) {
          stop("The shared-carryover efficacy fit did not return r_e posterior draws.")
        }
        matrix(draws[, "r_e"], nrow = nrow(draws), ncol = ndose)
      }
      result <- list(
        n_regular = fit$counts$n_regular,
        y_regular = fit$counts$y_regular,
        n_ipde = fit$counts$n_ipde,
        y_ipde = fit$counts$y_ipde,
        regular_draws = parameter_draws("p_regular"),
        carryover_draws = carryover_draws,
        ipde_draws = parameter_draws("p_ipde"),
        alpha_draws = NULL,
        history_data = NULL,
        model_file = model_file
      )
    } else {
      aide_phase12_load_efficacy_jags(require_standard_fitter = FALSE)
      model_path <- if (is.null(model_file)) {
        if (is_multicycle_additive) {
          "multicycle_additive_beta_binomial_efficacy.jags"
        } else if (is_dose_specific_additive) {
          "previous_dose_additive_dose_specific_beta_binomial_efficacy.jags"
        } else {
          "previous_dose_additive_beta_binomial_efficacy.jags"
        }
      } else {
        model_file
      }
      if (!file.exists(model_path)) {
        stop("Cannot find JAGS efficacy model file: ", model_path)
      }
      fit_data <- previous_dose_data
      ## JAGS cannot define an empty indexed loop; an unobserved sentinel row
      ## yields the beta priors unchanged when no administration is available.
      if (nrow(fit_data) == 0L) {
        fit_data <- if (is_multicycle_additive) {
          data.frame(
            admin_id = 1L, patient_id = 1L, cycle = 1L, dose = 1L,
            efficacy = NA_integer_, ipde = 0L, previous_admin_id = 0L,
            previous_state_index = 1L
          )
        } else {
          data.frame(
            dose = 1L, efficacy = NA_integer_, ipde = 0L, previous_dose = 1L
          )
        }
      }
      if (is_multicycle_additive) {
        fit_latest_admin_index <- if (nrow(previous_dose_data) == 0L) {
          1L
        } else {
          latest_admin_index
        }
        jags_data <- list(
          N_admin = nrow(fit_data),
          N_patient = length(fit_latest_admin_index),
          eff_patient = as.integer(fit_data$efficacy[fit_latest_admin_index]),
          dose = as.integer(fit_data$dose),
          previous_state_index = as.integer(fit_data$previous_state_index),
          latest_admin_index = as.integer(fit_latest_admin_index),
          ndose = ndose,
          a_regular = efficacy_prior[, "a"],
          b_regular = efficacy_prior[, "b"],
          a_alpha = as.numeric(efficacy_additive_alpha_prior[1L, "a"]),
          b_alpha = as.numeric(efficacy_additive_alpha_prior[1L, "b"])
        )
      } else {
        jags_data <- list(
          N = nrow(fit_data),
          ndose = ndose,
          eff = as.integer(fit_data$efficacy),
          dose = as.integer(fit_data$dose),
          ipde = as.integer(fit_data$ipde),
          a_regular = efficacy_prior[, "a"],
          b_regular = efficacy_prior[, "b"],
          a_alpha = if (is_dose_specific_additive) {
            efficacy_additive_alpha_prior[, "a"]
          } else {
            as.numeric(efficacy_additive_alpha_prior[1L, "a"])
          },
          b_alpha = if (is_dose_specific_additive) {
            efficacy_additive_alpha_prior[, "b"]
          } else {
            as.numeric(efficacy_additive_alpha_prior[1L, "b"])
          }
        )
        jags_data$previous_dose <- as.integer(fit_data$previous_dose)
      }
      jm <- rjags::jags.model(
        file = model_path, data = jags_data, n.chains = n_chains,
        n.adapt = n_adapt, quiet = TRUE
      )
      stats::update(jm, n.iter = n_burnin, progress.bar = "none")
      draws <- as.matrix(rjags::coda.samples(
        model = jm, variable.names = c("alpha", "p_regular"),
        n.iter = n_iter, thin = thin, progress.bar = "none"
      ))
      p_columns <- paste0("p_regular[", seq_len(ndose), "]")
      alpha_columns <- if (is_dose_specific_additive) {
        paste0("alpha[", seq_len(ndose), "]")
      } else {
        "alpha"
      }
      if (!all(c(alpha_columns, p_columns) %in% colnames(draws))) {
        stop("The additive previous-dose efficacy fit did not return alpha and all p_regular draws.")
      }
      regular_draws <- draws[, p_columns, drop = FALSE]
      alpha_draws <- if (is_dose_specific_additive) {
        draws[, alpha_columns, drop = FALSE]
      } else {
        as.numeric(draws[, "alpha"])
      }
      count_by_dose <- function(rows) {
        out <- numeric(ndose)
        if (length(rows) > 0L) {
          out <- tabulate(previous_dose_data$dose[rows], nbins = ndose)
        }
        out
      }
      sum_by_dose <- function(rows) {
        out <- numeric(ndose)
        if (length(rows) > 0L) {
          sums <- tapply(
            previous_dose_data$efficacy[rows], previous_dose_data$dose[rows], sum
          )
          out[as.integer(names(sums))] <- as.numeric(sums)
        }
        out
      }
      regular_rows <- which(previous_dose_data$ipde == 0L)
      ipde_rows <- which(previous_dose_data$ipde == 1L)
      n_regular <- count_by_dose(regular_rows)
      y_regular <- sum_by_dose(regular_rows)
      n_ipde <- count_by_dose(ipde_rows)
      y_ipde <- sum_by_dose(ipde_rows)
      ipde_draws <- if (is_dose_specific_additive) {
        regular_draws + alpha_draws * regular_draws
      } else {
        regular_draws + sweep(regular_draws, 1L, alpha_draws, "*")
      }
      ipde_draws[] <- pmin(1, ipde_draws)
      result <- list(
        n_regular = n_regular,
        y_regular = y_regular,
        n_ipde = n_ipde,
        y_ipde = y_ipde,
        regular_draws = regular_draws,
        carryover_draws = if (is_dose_specific_additive) {
          alpha_draws
        } else {
          matrix(alpha_draws, nrow = length(alpha_draws), ncol = ndose)
        },
        ## This is a same-dose reference only. Recycled posterior predictions
        ## use the actual current/previous dose pair in the recycle gate.
        ipde_draws = ipde_draws,
        alpha_draws = alpha_draws,
        history_data = if (is_multicycle_additive) previous_dose_data else NULL,
        latest_admin_index = if (is_multicycle_additive) latest_admin_index else NULL,
        model_file = model_path
      )
    }
    if (!is.null(cache)) assign(cache_name, result, envir = cache)
    result
  }

  prob_below_threshold <- if (is.null(threshold)) {
    rep(NA_real_, ndose)
  } else {
    colMeans(base$regular_draws < threshold)
  }

  list(
    n_regular = base$n_regular,
    y_regular = base$y_regular,
    n_ipde = base$n_ipde,
    y_ipde = base$y_ipde,
    n = base$n_regular + base$n_ipde,
    y = base$y_regular + base$y_ipde,
    regular_posterior_mean = colMeans(base$regular_draws),
    carryover_posterior_mean = colMeans(base$carryover_draws),
    ipde_posterior_mean = colMeans(base$ipde_draws),
    prob_regular_below_threshold = prob_below_threshold,
    efficacy_model = efficacy_model,
    efficacy_prior = efficacy_prior,
    efficacy_carryover_prior = efficacy_carryover_prior,
    efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
    regular_draws = base$regular_draws,
    alpha_draws = base$alpha_draws,
    history_data = base$history_data,
    latest_admin_index = base$latest_admin_index,
    mcmc = list(
      model_file = base$model_file,
      n_chains = as.integer(n_chains),
      n_adapt = as.integer(n_adapt),
      n_burnin = as.integer(n_burnin),
      n_iter = as.integer(n_iter),
      thin = as.integer(thin)
    )
  )
}

## Reconstruct posterior uncapped states from the shared-alpha multicycle
## draws. This is the patient-specific prediction interface for a proposed
## next cycle: q = min(1, p[next_dose] + alpha * current_state).
aide_phase12_multicycle_posterior_states <- function(posterior) {
  if (!identical(posterior$efficacy_model, "multicycle_additive")) {
    stop("posterior must come from the multicycle_additive efficacy model.")
  }
  history <- posterior$history_data
  if (is.null(history) || nrow(history) == 0L) {
    return(list(history = history, state_draws = matrix(0, nrow(posterior$regular_draws), 1L)))
  }
  alpha_draws <- as.numeric(posterior$alpha_draws)
  if (length(alpha_draws) != nrow(posterior$regular_draws)) {
    stop("Multicycle efficacy posterior alpha draws have an unexpected length.")
  }
  state_draws <- matrix(0, nrow = nrow(posterior$regular_draws),
                        ncol = nrow(history) + 1L)
  for (row in seq_len(nrow(history))) {
    state_draws[, row + 1L] <-
      posterior$regular_draws[, history$dose[row]] +
      alpha_draws * state_draws[, history$previous_state_index[row]]
  }
  list(history = history, state_draws = state_draws)
}

## Dose-specific IPDE efficacy gate using pi*_E,d2 - pi_E,d1 > delta.
aide_phase12_dose_specific_efficacy_recycle_gate <- function(admin,
                                                              current_dose,
                                                              next_dose,
                                                               ndose,
                                                               efficacy_prior = c(1, 1),
                                                               efficacy_carryover_prior = c(1, 9),
                                                               efficacy_model = c(
                                                                 "shared_carryover",
                                                                 "dose_specific_carryover",
                                                                 "previous_dose_additive",
                                                                 "dose_specific_previous_dose_additive",
                                                                 "multicycle_additive"
                                                               ),
                                                               efficacy_additive_alpha_prior = c(1, 9),
                                                               efficacy_model_file = NULL,
                                                              efficacy_n_chains = 3L,
                                                              efficacy_n_adapt = 1000L,
                                                              efficacy_n_burnin = 1000L,
                                                               efficacy_n_iter = 4000L,
                                                               efficacy_thin = 2L,
                                                               cache = NULL,
                                                               patient_id = NULL,
                                                               delta = 0.20) {
  efficacy_model <- match.arg(efficacy_model)
  if (length(current_dose) != 1L || !is.finite(current_dose) ||
      current_dose < 1L || current_dose > ndose ||
      current_dose != as.integer(current_dose)) {
    stop("current_dose must be an integer in 1:ndose.")
  }
  if (length(next_dose) != 1L || !is.finite(next_dose) ||
      next_dose < 1L || next_dose > ndose ||
      next_dose != as.integer(next_dose)) {
    stop("next_dose must be an integer in 1:ndose.")
  }
  if (length(delta) != 1L || !is.finite(delta) || delta < 0 || delta > 1) {
    stop("delta must be a scalar in [0, 1].")
  }
  posterior <- aide_phase12_efficacy_posterior(
    admin = admin,
    ndose = ndose,
    efficacy_prior = efficacy_prior,
    efficacy_carryover_prior = efficacy_carryover_prior,
    efficacy_model = efficacy_model,
    efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
    model_file = efficacy_model_file,
    n_chains = efficacy_n_chains,
    n_adapt = efficacy_n_adapt,
    n_burnin = efficacy_n_burnin,
    n_iter = efficacy_n_iter,
    thin = efficacy_thin,
    cache = cache
  )
  if (identical(efficacy_model, "multicycle_additive")) {
    if (length(patient_id) != 1L || is.na(patient_id)) {
      stop("patient_id is required for multicycle_additive next-cycle predictions.")
    }
    state <- aide_phase12_multicycle_posterior_states(posterior)
    history_rows <- which(state$history$patient_id == as.integer(patient_id))
    if (length(history_rows) == 0L) {
      stop("patient_id has no administration history in the multicycle posterior.")
    }
    current_row <- tail(history_rows, 1L)
    if (state$history$dose[current_row] != as.integer(current_dose)) {
      stop("current_dose does not match the patient's most recent administration.")
    }
    alpha_draws <- as.numeric(posterior$alpha_draws)
    theta_ipde_next_draws <- pmin(
      1,
      posterior$regular_draws[, as.integer(next_dose)] +
        alpha_draws * state$state_draws[, current_row + 1L]
    )
    theta_ipde_next <- mean(theta_ipde_next_draws)
    r_next <- mean(alpha_draws)
  } else if (aide_phase12_is_additive_efficacy_model(efficacy_model)) {
    alpha_draws <- if (is.matrix(posterior$alpha_draws)) {
      posterior$alpha_draws[, as.integer(next_dose)]
    } else {
      posterior$alpha_draws
    }
    theta_ipde_next_draws <- pmin(
      1,
      posterior$regular_draws[, as.integer(next_dose)] +
        alpha_draws * posterior$regular_draws[, as.integer(current_dose)]
    )
    theta_ipde_next <- mean(theta_ipde_next_draws)
    r_next <- mean(alpha_draws)
  } else {
    theta_ipde_next <- posterior$ipde_posterior_mean[as.integer(next_dose)]
    r_next <- posterior$carryover_posterior_mean[as.integer(next_dose)]
  }
  increment <- theta_ipde_next - posterior$regular_posterior_mean[as.integer(current_dose)]
  list(
    allowed = increment > delta,
    posterior_mean_increment = increment,
    p_regular_current = posterior$regular_posterior_mean[as.integer(current_dose)],
    r_next = r_next,
    theta_ipde_next = theta_ipde_next,
    delta = delta,
    current_dose = as.integer(current_dose),
    next_dose = as.integer(next_dose),
    n_regular = posterior$n_regular,
    y_regular = posterior$y_regular,
    n_ipde = posterior$n_ipde,
    y_ipde = posterior$y_ipde
  )
}

aide_phase12_beta_binomial_futility <- function(admin,
                                                  ndose,
                                                  efficacy_threshold = 0.20,
                                                  futility_cutoff = 0.95,
                                                  min_eff_n_for_futility = 0L) {
  if (length(ndose) != 1L || !is.finite(ndose) || ndose < 1L ||
      ndose != as.integer(ndose)) {
    stop("ndose must be a positive integer.")
  }
  ndose <- as.integer(ndose)
  if (length(efficacy_threshold) != 1L || !is.finite(efficacy_threshold) ||
      efficacy_threshold <= 0 || efficacy_threshold >= 1) {
    stop("efficacy_threshold must be a scalar in (0, 1).")
  }
  if (length(futility_cutoff) != 1L || !is.finite(futility_cutoff) ||
      futility_cutoff <= 0 || futility_cutoff >= 1) {
    stop("futility_cutoff must be a scalar in (0, 1).")
  }
  if (length(min_eff_n_for_futility) != 1L ||
      !is.finite(min_eff_n_for_futility) || min_eff_n_for_futility < 0 ||
      min_eff_n_for_futility != as.integer(min_eff_n_for_futility)) {
    stop("min_eff_n_for_futility must be a non-negative integer.")
  }
  if (is.null(admin)) admin <- data.frame(dose = integer(0), eff = integer(0))
  if (!all(c("dose", "eff") %in% names(admin))) {
    stop("admin must contain dose and eff for the efficacy-futility rule.")
  }
  if (nrow(admin) > 0L &&
      (any(!is.finite(admin$dose)) || any(admin$dose < 1L | admin$dose > ndose) ||
       any(!admin$eff %in% c(0L, 1L)))) {
    stop("admin contains invalid dose or efficacy values.")
  }

  n <- tabulate(as.integer(admin$dose), nbins = ndose)
  y <- tabulate(as.integer(admin$dose[admin$eff == 1L]), nbins = ndose)

  ## With the Beta(1, 1) prior in the supplied example,
  ## Pr(p_E,j < eta | y_j, n_j) = pbeta(eta, y_j + 1, n_j - y_j + 1).
  ## Equivalently, the posterior chance of efficacy exceeding eta is
  ## 1 - pbeta(...). A dose is futile when that latter chance is below
  ## 1 - futility_cutoff.
  prob_below_threshold <- stats::pbeta(
    efficacy_threshold,
    y + 1L,
    n - y + 1L
  )
  prob_above_threshold <- 1 - prob_below_threshold
  futility_eliminated <- as.integer(
    n >= as.integer(min_eff_n_for_futility) &
      prob_below_threshold > futility_cutoff
  )

  list(
    n = as.integer(n),
    y = as.integer(y),
    prob_below_threshold = prob_below_threshold,
    prob_above_threshold = prob_above_threshold,
    futility_eliminated = futility_eliminated,
    efficacy_threshold = efficacy_threshold,
    futility_cutoff = futility_cutoff,
    min_eff_n_for_futility = as.integer(min_eff_n_for_futility)
  )
}

aide_phase12_efficacy_summary <- function(admin,
                                           ndose,
                                           efficacy_prior = c(1, 1),
                                           efficacy_carryover_prior = c(1, 9),
                                           efficacy_model = c(
                                             "shared_carryover",
                                             "dose_specific_carryover",
                                             "previous_dose_additive",
                                             "dose_specific_previous_dose_additive",
                                             "multicycle_additive"
                                           ),
                                           efficacy_additive_alpha_prior = c(1, 9),
                                           efficacy_model_file = NULL,
                                           efficacy_n_chains = 3L,
                                           efficacy_n_adapt = 1000L,
                                           efficacy_n_burnin = 1000L,
                                           efficacy_n_iter = 4000L,
                                           efficacy_thin = 2L,
                                           cache = NULL,
                                           efficacy_threshold = 0.20,
                                           futility_cutoff = 0.95,
                                           min_eff_n_for_futility = 0L,
                                           futility_eliminated = NULL) {
  efficacy_model <- match.arg(efficacy_model)
  if (length(efficacy_threshold) != 1L || !is.finite(efficacy_threshold) ||
      efficacy_threshold <= 0 || efficacy_threshold >= 1) {
    stop("efficacy_threshold must be a scalar in (0, 1).")
  }
  if (length(futility_cutoff) != 1L || !is.finite(futility_cutoff) ||
      futility_cutoff <= 0 || futility_cutoff >= 1) {
    stop("futility_cutoff must be a scalar in (0, 1).")
  }
  if (length(min_eff_n_for_futility) != 1L ||
      !is.finite(min_eff_n_for_futility) || min_eff_n_for_futility < 0 ||
      min_eff_n_for_futility != as.integer(min_eff_n_for_futility)) {
    stop("min_eff_n_for_futility must be a non-negative integer.")
  }
  posterior <- aide_phase12_efficacy_posterior(
    admin = admin,
    ndose = ndose,
    efficacy_prior = efficacy_prior,
    efficacy_carryover_prior = efficacy_carryover_prior,
    efficacy_model = efficacy_model,
    efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
    model_file = efficacy_model_file,
    n_chains = efficacy_n_chains,
    n_adapt = efficacy_n_adapt,
    n_burnin = efficacy_n_burnin,
    n_iter = efficacy_n_iter,
    thin = efficacy_thin,
    cache = cache,
    threshold = efficacy_threshold
  )
  observed_futility <- if (is.null(futility_eliminated)) {
    aide_phase12_beta_binomial_futility(
      admin = admin,
      ndose = ndose,
      efficacy_threshold = efficacy_threshold,
      futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = min_eff_n_for_futility
    )
  } else {
    if (length(futility_eliminated) != ndose ||
        any(!futility_eliminated %in% c(0L, 1L))) {
      stop("futility_eliminated must contain one 0/1 value per dose.")
    }
    ## The simulator supplies its persistent state here.  This prevents an
    ## interim summary from reassessing unobserved or non-current doses.
    list(
      n = posterior$n,
      y = posterior$y,
      prob_below_threshold = rep(NA_real_, ndose),
      prob_above_threshold = rep(NA_real_, ndose),
      futility_eliminated = as.integer(futility_eliminated),
      efficacy_threshold = efficacy_threshold,
      futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = as.integer(min_eff_n_for_futility)
    )
  }
  c(
    posterior,
    list(
      posterior_mean = posterior$regular_posterior_mean,
      prob_below_threshold = observed_futility$prob_below_threshold,
      prob_above_threshold = observed_futility$prob_above_threshold,
      futility_eliminated = observed_futility$futility_eliminated,
      efficacy_futility_n = observed_futility$n,
      efficacy_futility_y = observed_futility$y,
      efficacy_threshold = observed_futility$efficacy_threshold,
      futility_cutoff = observed_futility$futility_cutoff,
      min_eff_n_for_futility = observed_futility$min_eff_n_for_futility
    )
  )
}

aide_phase12_compute_utility <- function(efficacy,
                                          toxicity,
                                          utility_type = 2L,
                                          lambda_T = 1,
                                          utility_scores = c(u00 = 0, u01 = 1,
                                                             u10 = 0, u11 = 0)) {
  efficacy <- as.numeric(efficacy)
  toxicity <- as.numeric(toxicity)
  if (length(efficacy) < 1L || length(toxicity) != length(efficacy) ||
      any(!is.na(efficacy) & (!is.finite(efficacy) | efficacy < 0 | efficacy > 1)) ||
      any(!is.na(toxicity) & (!is.finite(toxicity) | toxicity < 0 | toxicity > 1))) {
    stop("efficacy and toxicity must be equal-length probabilities in [0, 1] or NA.")
  }
  if (length(utility_type) != 1L || is.na(utility_type) ||
      !utility_type %in% 1:3) {
    stop("utility_type must be 1, 2, or 3.")
  }
  utility_type <- as.integer(utility_type)
  if (length(lambda_T) != 1L || !is.finite(lambda_T) || lambda_T < 0) {
    stop("lambda_T must be a single non-negative finite value.")
  }
  utility_scores <- as.numeric(utility_scores)
  if (length(utility_scores) != 4L || any(!is.finite(utility_scores))) {
    stop("utility_scores must contain four finite values: u00, u01, u10, u11.")
  }
  names(utility_scores) <- c("u00", "u01", "u10", "u11")
  switch(
    as.character(utility_type),
    "1" = efficacy,
    "2" = efficacy - lambda_T * toxicity,
    "3" = utility_scores[["u00"]] * (1 - toxicity) * (1 - efficacy) +
      utility_scores[["u01"]] * (1 - toxicity) * efficacy +
      utility_scores[["u10"]] * toxicity * (1 - efficacy) +
      utility_scores[["u11"]] * toxicity * efficacy
  )
}

aide_phase12_utility <- function(admin,
                                  ndose,
                                  efficacy_prior = c(1, 1),
                                  efficacy_carryover_prior = c(1, 9),
                                  efficacy_model = c(
                                    "shared_carryover",
                                    "dose_specific_carryover",
                                    "previous_dose_additive",
                                    "dose_specific_previous_dose_additive",
                                    "multicycle_additive"
                                  ),
                                  efficacy_additive_alpha_prior = c(1, 9),
                                  efficacy_model_file = NULL,
                                  efficacy_n_chains = 3L,
                                  efficacy_n_adapt = 1000L,
                                  efficacy_n_burnin = 1000L,
                                  efficacy_n_iter = 4000L,
                                  efficacy_thin = 2L,
                                  cache = NULL,
                                  toxicity_estimate,
                                  utility_type = 2L,
                                  lambda_T = 1,
                                  utility_scores = c(u00 = 0, u01 = 1,
                                                     u10 = 0, u11 = 0),
                                  utility_weight = NULL) {
  efficacy_model <- match.arg(efficacy_model)
  if (length(utility_type) != 1L || is.na(utility_type) ||
      !utility_type %in% 1:3) {
    stop("utility_type must be 1, 2, or 3.")
  }
  utility_type <- as.integer(utility_type)
  if (!is.null(utility_weight)) {
    if (length(utility_weight) != 1L || !is.finite(utility_weight) ||
        utility_weight < 0) {
      stop("utility_weight must be NULL or a single non-negative finite value.")
    }
    lambda_T <- as.numeric(utility_weight)
  }
  if (length(lambda_T) != 1L || !is.finite(lambda_T) || lambda_T < 0) {
    stop("lambda_T must be a single non-negative finite value.")
  }
  utility_scores <- as.numeric(utility_scores)
  if (length(utility_scores) != 4L || any(!is.finite(utility_scores))) {
    stop("utility_scores must contain four finite values: u00, u01, u10, u11.")
  }
  names(utility_scores) <- c("u00", "u01", "u10", "u11")
  posterior <- aide_phase12_efficacy_posterior(
    admin = admin,
    ndose = ndose,
    efficacy_prior = efficacy_prior,
    efficacy_carryover_prior = efficacy_carryover_prior,
    efficacy_model = efficacy_model,
    efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
    model_file = efficacy_model_file,
    n_chains = efficacy_n_chains,
    n_adapt = efficacy_n_adapt,
    n_burnin = efficacy_n_burnin,
    n_iter = efficacy_n_iter,
    thin = efficacy_thin,
    cache = cache
  )
  tox <- as.numeric(toxicity_estimate)
  if (length(tox) != ndose ||
      any(!is.na(tox) & (!is.finite(tox) | tox < 0 | tox > 1))) {
    stop("toxicity_estimate must contain probabilities in [0, 1] or NA for excluded doses.")
  }
  efficacy <- posterior$regular_posterior_mean
  utility <- aide_phase12_compute_utility(
    efficacy = efficacy,
    toxicity = tox,
    utility_type = utility_type,
    lambda_T = lambda_T,
    utility_scores = utility_scores
  )
  list(
    efficacy = efficacy,
    toxicity = tox,
    utility = utility,
    utility_type = utility_type,
    lambda_T = as.numeric(lambda_T),
    utility_scores = utility_scores,
    efficacy_posterior = posterior
  )
}

aide_phase12_select_lowest_tied <- function(score, candidates) {
  candidates <- as.integer(candidates)
  candidates <- candidates[is.finite(candidates)]
  if (length(candidates) == 0L) return(NA_integer_)
  candidates <- sort(unique(candidates))
  candidates[which.max(score[candidates])]
}

## Apply the common one-stage allocation rule used by the one-stage design
## and by Stage II of the two-stage design.  Safety and futility exclusions
## are supplied through admissible_doses; an observed efficacy response sets
## the lower end of the utility interval, while the toxicity-model MTD is the
## fallback target when no such response exists.
aide_phase12_one_stage_allocation <- function(current_dose,
                                               mtd,
                                               admissible_doses,
                                               utility,
                                               response_observed,
                                               tried_doses = integer(0)) {
  current_dose <- as.integer(current_dose)
  mtd <- as.integer(mtd)
  admissible_doses <- sort(unique(as.integer(admissible_doses)))
  tried_doses <- unique(as.integer(tried_doses))
  ndose <- length(utility)
  if (length(current_dose) != 1L || is.na(current_dose) ||
      current_dose < 1L || current_dose > ndose ||
      length(mtd) != 1L || is.na(mtd) || mtd < 1L || mtd > ndose ||
      length(response_observed) != ndose) {
    stop("Invalid one-stage allocation inputs.")
  }

  below_mtd <- admissible_doses[admissible_doses <= mtd]
  if (!length(below_mtd)) {
    return(list(dose = 99L, target = NA_integer_, candidates = integer(0),
                action = "no_admissible_one_stage_dose",
                response_floor = NA_integer_))
  }

  response_doses <- below_mtd[as.logical(response_observed[below_mtd])]
  if (length(response_doses)) {
    response_floor <- min(response_doses)
    candidates <- below_mtd[below_mtd >= response_floor]
    target <- aide_phase12_select_lowest_tied(utility, candidates)
    target_action <- "response_floor_utility_target"
  } else {
    response_floor <- NA_integer_
    candidates <- below_mtd
    ## The MTD is the specified no-response target.  If a standing
    ## elimination removes it, use the highest remaining admissible dose
    ## below that ceiling rather than assigning an eliminated dose.
    target <- if (mtd %in% candidates) mtd else max(candidates)
    target_action <- "mtd_no_response_fallback"
  }

  if (target <= current_dose) {
    return(list(dose = target, target = target, candidates = candidates,
                action = target_action, response_floor = response_floor))
  }

  ## Upward movement is governed by the path from the current dose, not by
  ## the highest dose ever tried.  This prevents a previously tried high dose
  ## from allowing an intervening untried dose to be skipped.
  path <- seq.int(current_dose + 1L, target)
  first_untried <- path[
    path %in% admissible_doses & !(path %in% tried_doses)
  ]
  assigned <- if (length(first_untried)) min(first_untried) else target
  list(
    dose = assigned,
    target = target,
    candidates = candidates,
    action = if (assigned == target) target_action else "no_skip_untried_dose",
    response_floor = response_floor
  )
}

simulate_AIDE_phase_I_II <- function(
    p_true,
    e_true,
    ## Recycled toxicity at d2 is p_ipde[d2] + toxicity_ipde_alpha * p_true[d1],
    ## capped at one, where d1 is the patient's immediately previous dose.
    p_ipde = p_true,
    toxicity_ipde_alpha = 0,
    ## Recycled efficacy at d2 is e_ipde[d2] + efficacy_ipde_alpha * e_true[d1],
    ## capped at one, where d1 is the patient's immediately previous dose.
    e_ipde = e_true,
    efficacy_ipde_alpha = 0,
    allocation = c("two_stage", "one_stage"),

    ## Stage I ends when the cumulative Stage I sample reaches N_s1 and the
    ## toxicity recommendation is stay. Stage II continues to Nmax using the
    ## same one-stage allocation rule.
    Nmax = 30L,
    N_s1 = 15L,
    ## Retained for existing calls; the revised design has no per-dose Stage
    ## II cap, so this value no longer affects allocation.
    N_s2 = Nmax,
    N_pat = Nmax,
    C = 3L,
    cycle_max = 1L,
    ## One-stage threshold for switching from {d-1,d,d+1} to {d-1,d}.
    m_U = 6L,
    ## ipde_design = 1 allows any eligible dose; 2 allows only an adjacent
    ## dose. flexible_ipde additionally permits recycling downward.
    ipde_design = 2L,
    flexible_ipde = FALSE,
    ## "continuous" keeps enrollment open and gives new patients priority.
    ## "ipde_first" closes recruitment while eligible IPDE patients fill the
    ## next cohort, enrolling only the remaining number of new patients.
    enrollment_scheme = c("continuous", "ipde_first"),
    ## Deprecated compatibility mapping: "new_first" -> continuous and
    ## "recycle_first" -> ipde_first.
    enrollment_priority = NULL,

    ## Non-TITE decisions use fully observed binary endpoints.  gen.tite is
    ## used to generate the event time and binary outcome for each endpoint.
    arrival_rate = 0.2,
    t0 = 0,
    T_assess = 28,
    dlt_dist = 2,
    dlt_alpha = 0.5,
    efficacy_dist = dlt_dist,
    efficacy_alpha = dlt_alpha,

    ## Toxicity model and AIDE controls.
    target = 0.30,
    cutoff = 0.95,
    model = c("BOIN", "CRM"),
    decision_method = c("boin", "approx1", "approx2"),
    mtd_method = NULL,
    r_carry = 0.10,
    dose_cap = 3L,

    ## CRM settings, used only when model = "CRM".
    crm_r_model = c(
      "fixed", "random", "level", "alpha_crm", "cumu_crm", "ipcrm",
      "previous_dose", "previous_dose_additive", "multicycle_additive"
    ),
    crm_skeleton = NULL,
    crm_alpha_sd = 2,
    crm_a_r = 1,
    crm_b_r = 9,
    crm_model_file = NULL,
    crm_fixed_model_file = "fix_CRM.bug",
    crm_random_model_file = "random_CRM.bug",
    crm_level_model_file = "random_CRM_level.bug",
    crm_previous_dose_model_file = "previous_dose_additive_CRM.bug",
    crm_n_chains = 2,
    crm_n_adapt = 500,
    crm_n_burnin = 500,
    crm_n_iter = 2000,
    crm_thin = 1,

    ## Efficacy model and allocation utility.
    ## This paired option controls both endpoint models.  It maps the four
    ## efficacy choices to the two shared toxicity choices.  Leave NULL to
    ## retain direct control through crm_r_model and efficacy_model.
    carryover_model = NULL,
    efficacy_prior = c(1, 1),
    efficacy_carryover_prior = c(1, 9),
    efficacy_model = c(
      "shared_carryover", "dose_specific_carryover",
      "previous_dose_additive", "dose_specific_previous_dose_additive",
      "multicycle_additive"
    ),
    efficacy_additive_alpha_prior = c(1, 9),
    efficacy_model_file = NULL,
    efficacy_n_chains = 3L,
    efficacy_n_adapt = 1000L,
    efficacy_n_burnin = 1000L,
    efficacy_n_iter = 4000L,
    efficacy_thin = 2L,
    efficacy_threshold = 0.20,
    futility_cutoff = 0.95,
    min_eff_n_for_futility = 0L,
    ## Utility 1 = efficacy; Utility 2 = efficacy - lambda_T * toxicity;
    ## Utility 3 = user-specified four-outcome joint score.
    utility_type = 2L,
    lambda_T = 1,
    utility_scores = c(u00 = 0, u01 = 1, u10 = 0, u11 = 0),
    ## Deprecated alias for lambda_T, retained for existing scripts.
    utility_weight = NULL,

    ## Recycling gates. They are used only for model = "CRM" with
    ## crm_r_model = "random"; the switches make each gate optional.
    apply_random_crm_recycle_toxicity_rule = TRUE,
    random_crm_recycle_toxicity_cutoff = cutoff,
    apply_random_crm_recycle_efficacy_rule = TRUE,
    random_crm_recycle_efficacy_delta = 0.05,
    seed = NULL,
    verbose = FALSE
) {
  if (!is.null(seed)) set.seed(seed)

  allocation <- match.arg(allocation)
  model <- match.arg(model)
  enrollment_scheme <- match.arg(enrollment_scheme)
  if (!is.null(enrollment_priority)) {
    enrollment_priority <- match.arg(
      enrollment_priority, c("recycle_first", "new_first")
    )
    enrollment_scheme <- if (enrollment_priority == "recycle_first") {
      "ipde_first"
    } else {
      "continuous"
    }
  }
  decision_method <- match.arg(decision_method)
  if (!is.null(carryover_model)) {
    if (model != "CRM") {
      stop("carryover_model is available only when model = 'CRM'.")
    }
    carryover_model <- match.arg(
      carryover_model,
      c(
        "discount_shared", "discount_dose_specific",
        "additive_shared", "additive_dose_specific", "multicycle_additive"
      )
    )
    carryover_settings <- aide_phase12_carryover_model_map(carryover_model)
    crm_r_model <- carryover_settings$crm_r_model
    efficacy_model <- carryover_settings$efficacy_model
  }
  crm_r_model <- crm_normalize_r_model(match.arg(crm_r_model))
  efficacy_model <- match.arg(efficacy_model)
  multicycle_toxicity_active <- model == "CRM" &&
    identical(crm_r_model, "multicycle_additive")
  multicycle_efficacy_active <- aide_phase12_is_multicycle_efficacy_model(
    efficacy_model
  )
  if (is.null(mtd_method)) mtd_method <- decision_method
  mtd_method <- match.arg(mtd_method, c("boin", "approx1", "approx2"))

  p_true <- as.numeric(p_true)
  e_true <- as.numeric(e_true)
  p_ipde <- as.numeric(p_ipde)
  e_ipde <- as.numeric(e_ipde)
  ndose <- length(p_true)
  if (ndose < 1L || length(e_true) != ndose || length(p_ipde) != ndose ||
      length(e_ipde) != ndose) {
    stop("p_true, e_true, p_ipde, and e_ipde must have the same positive length.")
  }
  if (any(!is.finite(c(p_true, e_true, p_ipde, e_ipde))) ||
      any(c(p_true, e_true, p_ipde, e_ipde) < 0 | c(p_true, e_true, p_ipde, e_ipde) > 1)) {
    stop("All true toxicity and efficacy probabilities must lie in [0, 1].")
  }
  if (length(toxicity_ipde_alpha) != 1L || !is.finite(toxicity_ipde_alpha) ||
      toxicity_ipde_alpha < 0) {
    stop("toxicity_ipde_alpha must be a single finite non-negative value.")
  }
  toxicity_ipde_alpha <- as.numeric(toxicity_ipde_alpha)
  if (length(efficacy_ipde_alpha) != 1L || !is.finite(efficacy_ipde_alpha) ||
      efficacy_ipde_alpha < 0) {
    stop("efficacy_ipde_alpha must be a single finite non-negative value.")
  }
  efficacy_ipde_alpha <- as.numeric(efficacy_ipde_alpha)
  if (length(N_s1) != 1L || N_s1 < 1L || N_s1 != as.integer(N_s1)) {
    stop("N_s1 must be a positive integer.")
  }
  if (length(N_s2) != 1L || N_s2 < 1L || N_s2 != as.integer(N_s2)) {
    stop("N_s2 must be a positive integer.")
  }
  if (length(C) != 1L || C < 1L || C != as.integer(C)) {
    stop("C must be a positive integer.")
  }
  if (length(Nmax) != 1L || Nmax < C || Nmax != as.integer(Nmax)) {
    stop("Nmax must be an integer at least as large as C.")
  }
  if (allocation == "two_stage" && N_s1 > Nmax) {
    stop("For allocation = 'two_stage', N_s1 must not exceed Nmax.")
  }
  if (length(N_pat) != 1L || N_pat < 1L || N_pat != as.integer(N_pat)) {
    stop("N_pat must be a positive integer.")
  }
  if (length(cycle_max) != 1L || cycle_max < 1L || cycle_max != as.integer(cycle_max)) {
    stop("cycle_max must be a positive integer.")
  }
  if (length(m_U) != 1L || !is.finite(m_U) || m_U < 1L ||
      m_U != as.integer(m_U)) {
    stop("m_U must be a positive integer.")
  }
  m_U <- as.integer(m_U)
  if (!ipde_design %in% c(1L, 2L)) {
    stop("ipde_design must be 1 (any dose) or 2 (adjacent dose).")
  }
  flexible_ipde <- aide_phase12_validate_logical_flag(
    flexible_ipde,
    "flexible_ipde"
  )
  if (length(arrival_rate) != 1L || !is.finite(arrival_rate) || arrival_rate <= 0) {
    stop("arrival_rate must be a positive finite value.")
  }
  if (length(T_assess) != 1L || !is.finite(T_assess) || T_assess < 0) {
    stop("T_assess must be a single non-negative finite value.")
  }
  validate_tite_controls <- function(dist, alpha, endpoint) {
    if (length(dist) != 1L || is.na(dist) || !dist %in% 1:4) {
      stop(endpoint, "_dist must be one of 1, 2, 3, or 4.")
    }
    if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
      stop(endpoint, "_alpha must be a single value in [0, 1].")
    }
    if (dist %in% c(2, 3) && (alpha <= 0 || alpha >= 1)) {
      stop(endpoint, "_alpha must lie strictly between 0 and 1 for ", endpoint, "_dist 2 or 3.")
    }
    list(dist = as.integer(dist), alpha = as.numeric(alpha))
  }
  dlt_tite <- validate_tite_controls(dlt_dist, dlt_alpha, "dlt")
  efficacy_tite <- validate_tite_controls(
    efficacy_dist, efficacy_alpha, "efficacy"
  )
  if (length(target) != 1L || !is.finite(target) || target <= 0 || target >= 1) {
    stop("target must be a scalar in (0, 1).")
  }
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff >= 1) {
    stop("cutoff must be a scalar in (0, 1).")
  }
  if (!is.null(utility_weight)) {
    if (length(utility_weight) != 1L || !is.finite(utility_weight) ||
        utility_weight < 0) {
      stop("utility_weight must be NULL or a single non-negative finite value.")
    }
    lambda_T <- as.numeric(utility_weight)
  }
  if (length(utility_type) != 1L || is.na(utility_type) ||
      !utility_type %in% 1:3) {
    stop("utility_type must be 1, 2, or 3.")
  }
  utility_type <- as.integer(utility_type)
  if (length(lambda_T) != 1L || !is.finite(lambda_T) || lambda_T < 0) {
    stop("lambda_T must be a single non-negative finite value.")
  }
  lambda_T <- as.numeric(lambda_T)
  utility_scores <- as.numeric(utility_scores)
  if (length(utility_scores) != 4L || any(!is.finite(utility_scores))) {
    stop("utility_scores must contain four finite values: u00, u01, u10, u11.")
  }
  names(utility_scores) <- c("u00", "u01", "u10", "u11")
  efficacy_prior <- aide_phase12_expand_beta_prior(
    efficacy_prior, ndose, "efficacy_prior"
  )
  efficacy_carryover_prior <- aide_phase12_expand_beta_prior(
    efficacy_carryover_prior, ndose, "efficacy_carryover_prior"
  )
  efficacy_additive_alpha_prior <- aide_phase12_expand_beta_prior(
    efficacy_additive_alpha_prior, ndose, "efficacy_additive_alpha_prior"
  )
  efficacy_mcmc <- as.integer(c(
    efficacy_n_chains, efficacy_n_adapt, efficacy_n_burnin,
    efficacy_n_iter, efficacy_thin
  ))
  if (any(is.na(efficacy_mcmc)) || efficacy_mcmc[1L] < 1L ||
      efficacy_mcmc[2L] < 0L || efficacy_mcmc[3L] < 0L ||
      efficacy_mcmc[4L] < 1L || efficacy_mcmc[5L] < 1L) {
    stop("The efficacy JAGS MCMC controls must be valid positive integers (or zero for adaptation/burn-in).")
  }
  efficacy_n_chains <- efficacy_mcmc[1L]
  efficacy_n_adapt <- efficacy_mcmc[2L]
  efficacy_n_burnin <- efficacy_mcmc[3L]
  efficacy_n_iter <- efficacy_mcmc[4L]
  efficacy_thin <- efficacy_mcmc[5L]
  aide_phase12_load_efficacy_jags(
    require_standard_fitter = efficacy_model %in% c(
      "dose_specific_carryover", "shared_carryover"
    )
  )
  apply_random_crm_recycle_toxicity_rule <- aide_phase12_validate_logical_flag(
    apply_random_crm_recycle_toxicity_rule,
    "apply_random_crm_recycle_toxicity_rule"
  )
  apply_random_crm_recycle_efficacy_rule <- aide_phase12_validate_logical_flag(
    apply_random_crm_recycle_efficacy_rule,
    "apply_random_crm_recycle_efficacy_rule"
  )
  if (length(random_crm_recycle_toxicity_cutoff) != 1L ||
      !is.finite(random_crm_recycle_toxicity_cutoff) ||
      random_crm_recycle_toxicity_cutoff <= 0 ||
      random_crm_recycle_toxicity_cutoff >= 1) {
    stop("random_crm_recycle_toxicity_cutoff must be a scalar in (0, 1).")
  }
  if (length(random_crm_recycle_efficacy_delta) != 1L ||
      !is.finite(random_crm_recycle_efficacy_delta) ||
      random_crm_recycle_efficacy_delta < 0 ||
      random_crm_recycle_efficacy_delta > 1) {
    stop("random_crm_recycle_efficacy_delta must be a scalar in [0, 1].")
  }
  if (model == "CRM" && is.null(crm_skeleton)) {
    stop("For model = 'CRM', provide crm_skeleton.")
  }
  if (model == "CRM") crm_validate_skeleton(crm_skeleton, ndose)

  ## BOIN boundaries are prepared once and used for every toxicity decision.
  boin_bounds <- get.boundary(
    target = target,
    ncohort = ceiling(Nmax / C),
    cohortsize = C,
    design = 3,
    cutoff.eli = cutoff
  )
  b_e <- boin_bounds[4L, ]
  b_d <- boin_bounds[3L, ]
  boin_lambdas <- boin_boundary(target)

  admin <- data.frame(
    row_id = integer(0),
    id = integer(0),
    t_arrival = numeric(0),
    t_start = numeric(0),
    t_tox = numeric(0),
    t_eff = numeric(0),
    t_eval = numeric(0),
    dose = integer(0),
    y = integer(0),
    eff = integer(0),
    true_toxicity_probability = numeric(0),
    true_efficacy_probability = numeric(0),
    ## These state fields are used only by the multicycle model. They retain
    ## the uncapped deterministic state needed for the next administration.
    previous_admin_id = integer(0),
    true_toxicity_state = numeric(0),
    true_efficacy_state = numeric(0),
    ncycle = integer(0),
    cohort = integer(0),
    stage = character(0),
    type = character(0),
    stringsAsFactors = FALSE
  )
  decision_log <- data.frame(
    cohort = integer(0),
    stage = character(0),
    current_dose = integer(0),
    toxicity_next_dose = integer(0),
    allocated_dose = integer(0),
    toxicity_action = character(0),
    n_current = integer(0),
    efficacy_futility_eliminated = integer(0),
    n_new_enrolled = integer(0),
    n_recycled_enrolled = integer(0),
    waiting_queue_size = integer(0),
    stringsAsFactors = FALSE
  )
  ## One row per patient considered for a recycled/IPDE administration.  This
  ## preserves the exact posterior quantities used by the two random-CRM
  ## recycling gates, including candidates that are rejected.
  recycling_decision_log <- data.frame(
    cohort = integer(0),
    stage = character(0),
    patient_id = integer(0),
    current_dose = integer(0),
    next_dose = integer(0),
    toxicity_rule_applied = logical(0),
    toxicity_probability_over_phi = numeric(0),
    toxicity_phi = numeric(0),
    toxicity_cutoff = numeric(0),
    toxicity_allowed = logical(0),
    efficacy_rule_applied = logical(0),
    efficacy_p_regular_current = numeric(0),
    efficacy_r_next = numeric(0),
    efficacy_theta_ipde_next = numeric(0),
    efficacy_posterior_mean_increment = numeric(0),
    efficacy_delta = numeric(0),
    efficacy_allowed = logical(0),
    eligible_for_recycling = logical(0),
    selected_for_recycling = logical(0),
    stringsAsFactors = FALSE
  )

  ## Multiple interim quantities use the same efficacy posterior. Cache the
  ## JAGS draws by the current administration data and MCMC settings so a
  ## cohort decision is based on one common posterior fit.
  efficacy_posterior_cache <- new.env(parent = emptyenv())
  efficacy_summary_current <- function() {
    summary <- aide_phase12_efficacy_summary(
      admin = admin,
      ndose = ndose,
      efficacy_prior = efficacy_prior,
      efficacy_carryover_prior = efficacy_carryover_prior,
      efficacy_model = efficacy_model,
      efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
      efficacy_model_file = efficacy_model_file,
      efficacy_n_chains = efficacy_n_chains,
      efficacy_n_adapt = efficacy_n_adapt,
      efficacy_n_burnin = efficacy_n_burnin,
      efficacy_n_iter = efficacy_n_iter,
      efficacy_thin = efficacy_thin,
      cache = efficacy_posterior_cache,
      efficacy_threshold = efficacy_threshold,
      futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = min_eff_n_for_futility,
      futility_eliminated = efficacy_futility_eliminated
    )
    summary
  }
  efficacy_utility_current <- function(toxicity_estimate) {
    aide_phase12_utility(
      admin = admin,
      ndose = ndose,
      efficacy_prior = efficacy_prior,
      efficacy_carryover_prior = efficacy_carryover_prior,
      efficacy_model = efficacy_model,
      efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
      efficacy_model_file = efficacy_model_file,
      efficacy_n_chains = efficacy_n_chains,
      efficacy_n_adapt = efficacy_n_adapt,
      efficacy_n_burnin = efficacy_n_burnin,
      efficacy_n_iter = efficacy_n_iter,
      efficacy_thin = efficacy_thin,
      cache = efficacy_posterior_cache,
      toxicity_estimate = toxicity_estimate,
      utility_type = utility_type,
      lambda_T = lambda_T,
      utility_scores = utility_scores
    )
  }
  efficacy_recycle_gate_current <- function(current_dose, next_dose) {
    aide_phase12_dose_specific_efficacy_recycle_gate(
      admin = admin,
      current_dose = current_dose,
      next_dose = next_dose,
      ndose = ndose,
      efficacy_prior = efficacy_prior,
      efficacy_carryover_prior = efficacy_carryover_prior,
      efficacy_model = efficacy_model,
      efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
      efficacy_model_file = efficacy_model_file,
      efficacy_n_chains = efficacy_n_chains,
      efficacy_n_adapt = efficacy_n_adapt,
      efficacy_n_burnin = efficacy_n_burnin,
      efficacy_n_iter = efficacy_n_iter,
      efficacy_thin = efficacy_thin,
      cache = efficacy_posterior_cache,
      delta = random_crm_recycle_efficacy_delta
    )
  }

  toxicity_eliminated <- rep(0L, ndose)
  efficacy_futility_eliminated <- rep(0L, ndose)
  futility_state <- function() {
    list(futility_eliminated = as.integer(efficacy_futility_eliminated))
  }
  evaluate_current_dose_futility <- function(dose) {
    if (length(dose) != 1L || !is.finite(dose) || dose < 1L || dose > ndose) {
      stop("dose must identify one valid dose level.")
    }
    dose <- as.integer(dose)
    observed <- aide_phase12_beta_binomial_futility(
      admin = admin,
      ndose = ndose,
      efficacy_threshold = efficacy_threshold,
      futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = min_eff_n_for_futility
    )
    ## Only the just-treated dose is assessed.  Once a dose is futile, it
    ## remains eliminated for every later allocation.
    efficacy_futility_eliminated[dose] <<- max(
      efficacy_futility_eliminated[dose],
      observed$futility_eliminated[dose]
    )
    observed$futility_eliminated <- as.integer(efficacy_futility_eliminated)
    observed
  }
  apply_early_stop_rule <- function(stage) {
    if (toxicity_eliminated[1L] == 1L) {
      earlystop <<- TRUE
      stop_reason <<- paste0("dose1_toxicity_eliminated_", stage)
      return(TRUE)
    }
    if (all(efficacy_futility_eliminated == 1L)) {
      earlystop <<- TRUE
      stop_reason <<- paste0("all_doses_efficacy_futile_", stage)
      return(TRUE)
    }
    FALSE
  }
  ## Under continuous enrollment, all potential new-patient arrivals are
  ## generated up front and accumulate in the waiting queue. Under IPDE-first
  ## recruitment, no new patient is recruited until an open cohort still has
  ## unfilled positions after eligible IPDE patients are assigned.
  arrival_schedule <- if (enrollment_scheme == "continuous") {
    data.frame(
      id = seq_len(N_pat),
      t_arrival = as.numeric(t0 + cumsum(stats::rexp(N_pat, rate = arrival_rate))),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(id = integer(0), t_arrival = numeric(0), stringsAsFactors = FALSE)
  }
  next_new_id <- 1L
  cohort_id <- 0L
  decision_time <- as.numeric(t0)
  current_dose <- 1L
  earlystop <- FALSE
  stop_reason <- NA_character_
  random_crm_recycle_rules_active <- model == "CRM" && crm_r_model == "random"
  random_crm_recycle_post_cache <- list(
    n_admin = NA_integer_,
    post = NULL
  )

  latest_patient_state <- function() {
    if (nrow(admin) == 0L) return(admin)
    z <- admin[order(admin$id, admin$t_eval, admin$row_id), , drop = FALSE]
    z[!duplicated(z$id, fromLast = TRUE), , drop = FALSE]
  }

  cache_random_crm_recycle_post <- function(model_fit) {
    if (!random_crm_recycle_rules_active || is.null(model_fit)) return(invisible(NULL))
    post <- if (!is.null(model_fit$fit$post)) model_fit$fit$post else model_fit$post
    if (is.null(post)) return(invisible(NULL))
    post <- as.matrix(post)
    p_columns <- paste0("p[", seq_len(ndose), "]")
    if (nrow(post) > 0L && "r" %in% colnames(post) &&
        all(p_columns %in% colnames(post))) {
      random_crm_recycle_post_cache <<- list(
        n_admin = nrow(admin),
        post = post
      )
    }
    invisible(NULL)
  }

  random_crm_recycle_post <- function() {
    if (!random_crm_recycle_rules_active) return(NULL)
    if (identical(random_crm_recycle_post_cache$n_admin, nrow(admin)) &&
        !is.null(random_crm_recycle_post_cache$post)) {
      return(random_crm_recycle_post_cache$post)
    }

    ## Reuse the same user-selected random CRM and MCMC controls as AIDE's
    ## toxicity decisions whenever a new posterior is required for recycling.
    fitted <- toxicity_utility_fit()
    cache_random_crm_recycle_post(fitted)
    if (is.null(random_crm_recycle_post_cache$post)) {
      stop("The random CRM fit did not return posterior samples needed for the recycling safety rule.")
    }
    random_crm_recycle_post_cache$post
  }

  enrolled_new_ids <- function() {
    if (nrow(admin) == 0L) return(integer(0))
    as.integer(unique(admin$id[admin$type == "new"]))
  }

  waiting_new_ids <- function(exclude = integer(0)) {
    enrolled <- enrolled_new_ids()
    available <- arrival_schedule$id[
      arrival_schedule$t_arrival <= decision_time &
        !arrival_schedule$id %in% enrolled &
        !arrival_schedule$id %in% as.integer(exclude)
    ]
    as.integer(available)
  }

  wait_for_new_ids <- function(n_needed, exclude = integer(0)) {
    if (n_needed <= 0L) return(integer(0))
    if (enrollment_scheme != "continuous") {
      stop("wait_for_new_ids() is only used for continuous enrollment.")
    }
    available <- waiting_new_ids(exclude)
    if (length(available) >= n_needed) return(head(available, n_needed))

    enrolled <- enrolled_new_ids()
    future <- arrival_schedule[
      !arrival_schedule$id %in% enrolled &
        !arrival_schedule$id %in% as.integer(exclude) &
        arrival_schedule$t_arrival > decision_time,
      , drop = FALSE
    ]
    if (length(available) + nrow(future) < n_needed) return(integer(0))
    n_from_future <- n_needed - length(available)
    decision_time <<- future$t_arrival[n_from_future]
    c(available, head(waiting_new_ids(c(exclude, available)), n_from_future))
  }

  recruit_ipde_first_new_ids <- function(n_needed) {
    if (n_needed <= 0L) return(integer(0))
    if (enrollment_scheme != "ipde_first") {
      stop("recruit_ipde_first_new_ids() is only used for IPDE-first recruitment.")
    }
    available <- as.integer(N_pat - next_new_id + 1L)
    if (available < n_needed) return(integer(0))

    ids <- seq.int(next_new_id, length.out = n_needed)
    arrival_times <- decision_time + cumsum(stats::rexp(n_needed, rate = arrival_rate))
    arrival_schedule <<- rbind(
      arrival_schedule,
      data.frame(
        id = as.integer(ids),
        t_arrival = as.numeric(arrival_times),
        stringsAsFactors = FALSE
      )
    )
    next_new_id <<- max(ids) + 1L
    ## The cohort waits only until the final required new patient arrives.
    decision_time <<- max(arrival_times)
    as.integer(ids)
  }

  eligible_ipde_ids <- function(next_dose, stage, futility = NULL) {
    if (cycle_max <= 1L) {
      return(integer(0))
    }
    if (is.null(futility)) futility <- futility_state()
    if (toxicity_eliminated[next_dose] == 1L ||
        futility$futility_eliminated[next_dose] == 1L) {
      ## A dose eliminated for toxicity or efficacy cannot be an IPDE
      ## destination, including during Stage I.
      return(integer(0))
    }
    state <- latest_patient_state()
    if (nrow(state) == 0L) return(integer(0))
    dose_ok <- aide_phase12_recycle_dose_eligible(
      previous_dose = state$dose,
      next_dose = next_dose,
      ipde_design = ipde_design,
      flexible_ipde = flexible_ipde
    )
    ## An IPDE administration requires completion of a cycle with neither a
    ## DLT nor an efficacy response, as specified for the efficacy-enabled
    ## design.  t_eval is set only after both endpoint times are observed.
    state <- state[
      state$y == 0L & state$eff == 0L & state$ncycle < cycle_max & dose_ok,
      , drop = FALSE
    ]
    if (nrow(state) == 0L) return(integer(0))

    ## These two constraints are deliberately limited to the random-r CRM.
    ## They protect the proposed recycled administration, not the regular
    ## allocation decision or the selected AIDE toxicity model itself.
    toxicity_rule_applied <- random_crm_recycle_rules_active &&
      apply_random_crm_recycle_toxicity_rule
    efficacy_rule_applied <- random_crm_recycle_rules_active &&
      apply_random_crm_recycle_efficacy_rule
    toxicity_gate <- list(
      allowed = TRUE,
      probability_over_phi = NA_real_,
      phi = target,
      cutoff = random_crm_recycle_toxicity_cutoff
    )
    if (random_crm_recycle_rules_active &&
        apply_random_crm_recycle_toxicity_rule) {
      toxicity_gate <- aide_phase12_random_crm_recycle_toxicity_gate(
        post = random_crm_recycle_post(),
        next_dose = next_dose,
        phi = target,
        cutoff = random_crm_recycle_toxicity_cutoff
      )
    }
    efficacy_gate_by_current_dose <- list()
    if (efficacy_rule_applied) {
      efficacy_gate_by_current_dose <- lapply(
        unique(as.integer(state$dose)),
        function(current_dose_for_patient) {
          efficacy_recycle_gate_current(current_dose_for_patient, next_dose)
        }
      )
      names(efficacy_gate_by_current_dose) <- as.character(unique(as.integer(state$dose)))
    }
    efficacy_allowed <- if (efficacy_rule_applied) {
      vapply(
        as.character(as.integer(state$dose)),
        function(current_dose_for_patient) {
          efficacy_gate_by_current_dose[[current_dose_for_patient]]$allowed
        },
        logical(1)
      )
    } else {
      rep(TRUE, nrow(state))
    }
    toxicity_allowed <- rep(isTRUE(toxicity_gate$allowed), nrow(state))
    eligible_for_recycling <- toxicity_allowed & efficacy_allowed

    efficacy_gate_for_patient <- if (efficacy_rule_applied) {
      lapply(as.character(as.integer(state$dose)), function(current_dose_for_patient) {
        efficacy_gate_by_current_dose[[current_dose_for_patient]]
      })
    } else {
      rep(list(NULL), nrow(state))
    }
    recycling_decision_log <<- rbind(
      recycling_decision_log,
      data.frame(
        cohort = rep.int(cohort_id + 1L, nrow(state)),
        stage = rep(stage, nrow(state)),
        patient_id = as.integer(state$id),
        current_dose = as.integer(state$dose),
        next_dose = rep.int(as.integer(next_dose), nrow(state)),
        toxicity_rule_applied = rep(toxicity_rule_applied, nrow(state)),
        toxicity_probability_over_phi = rep(toxicity_gate$probability_over_phi, nrow(state)),
        toxicity_phi = rep(toxicity_gate$phi, nrow(state)),
        toxicity_cutoff = rep(toxicity_gate$cutoff, nrow(state)),
        toxicity_allowed = toxicity_allowed,
        efficacy_rule_applied = rep(efficacy_rule_applied, nrow(state)),
        efficacy_p_regular_current = vapply(
          efficacy_gate_for_patient,
          function(x) if (is.null(x)) NA_real_ else x$p_regular_current,
          numeric(1)
        ),
        efficacy_r_next = vapply(
          efficacy_gate_for_patient,
          function(x) if (is.null(x)) NA_real_ else x$r_next,
          numeric(1)
        ),
        efficacy_theta_ipde_next = vapply(
          efficacy_gate_for_patient,
          function(x) if (is.null(x)) NA_real_ else x$theta_ipde_next,
          numeric(1)
        ),
        efficacy_posterior_mean_increment = vapply(
          efficacy_gate_for_patient,
          function(x) if (is.null(x)) NA_real_ else x$posterior_mean_increment,
          numeric(1)
        ),
        efficacy_delta = vapply(
          efficacy_gate_for_patient,
          function(x) if (is.null(x)) NA_real_ else x$delta,
          numeric(1)
        ),
        efficacy_allowed = efficacy_allowed,
        eligible_for_recycling = eligible_for_recycling,
        selected_for_recycling = rep(FALSE, nrow(state)),
        stringsAsFactors = FALSE
      )
    )
    state <- state[eligible_for_recycling, , drop = FALSE]
    state <- state[order(state$t_eval, state$id), , drop = FALSE]
    as.integer(state$id)
  }

  draw_tite_endpoint <- function(probability, start_time, dist, alpha) {
    z <- gen.tite(
      dist = dist,
      n = 1L,
      pi = probability,
      alpha = alpha,
      Tobs = T_assess
    )
    list(
      outcome = as.integer(z$tox[1L]),
      time = as.numeric(z$t.tox[1L]),
      t_observed = as.numeric(start_time + z$t.tox[1L])
    )
  }

  enroll_cohort <- function(dose, stage, n_to_enroll, futility = NULL) {
    if (!is.finite(dose) || dose < 1L || dose > ndose || n_to_enroll < 1L) {
      return(FALSE)
    }
    if (is.null(futility)) futility <- futility_state()
    if (toxicity_eliminated[dose] == 1L ||
        futility$futility_eliminated[dose] == 1L) {
      return(FALSE)
    }
    n_to_enroll <- as.integer(n_to_enroll)
    recycle_candidates <- eligible_ipde_ids(dose, stage, futility = futility)
    if (enrollment_scheme == "continuous") {
      ## Continuous enrollment gives already-waiting new patients priority;
      ## eligible IPDE patients fill only the remaining cohort positions.
      queued_new <- waiting_new_ids()
      new_ids <- head(queued_new, n_to_enroll)
      ret_ids <- head(recycle_candidates, n_to_enroll - length(new_ids))
      still_needed <- n_to_enroll - length(ret_ids) - length(new_ids)
      if (still_needed > 0L) {
        later_new_ids <- wait_for_new_ids(still_needed, exclude = new_ids)
        if (length(later_new_ids) != still_needed) return(FALSE)
        new_ids <- c(new_ids, later_new_ids)
      }
    } else {
      ## IPDE-first recruitment does not retain a continuous new-patient
      ## queue. Recruit exactly the remaining slots after recycling.
      ret_ids <- head(recycle_candidates, n_to_enroll)
      new_ids <- recruit_ipde_first_new_ids(n_to_enroll - length(ret_ids))
      if (length(new_ids) != n_to_enroll - length(ret_ids)) return(FALSE)
    }
    new_arrivals <- arrival_schedule$t_arrival[match(new_ids, arrival_schedule$id)]

    cohort_id <<- cohort_id + 1L
    if (length(ret_ids) > 0L) {
      recycle_log_rows <- recycling_decision_log$cohort == cohort_id &
        recycling_decision_log$stage == stage &
        recycling_decision_log$patient_id %in% ret_ids
      recycling_decision_log$selected_for_recycling[recycle_log_rows] <<- TRUE
    }
    start_time <- decision_time
    state <- latest_patient_state()

    append_row <- function(id, arrival, toxicity_probability, efficacy_probability,
                           toxicity_state, efficacy_state, previous_admin_id,
                           cycle, type) {
      tox <- draw_tite_endpoint(
        probability = toxicity_probability,
        start_time = start_time,
        dist = dlt_tite$dist,
        alpha = dlt_tite$alpha
      )
      eff <- draw_tite_endpoint(
        probability = efficacy_probability,
        start_time = start_time,
        dist = efficacy_tite$dist,
        alpha = efficacy_tite$alpha
      )
      admin <<- rbind(
        admin,
        data.frame(
          row_id = nrow(admin) + 1L,
          id = as.integer(id),
          t_arrival = as.numeric(arrival),
          t_start = start_time,
          t_tox = tox$time,
          t_eff = eff$time,
          ## This non-TITE implementation updates models only after both
          ## binary endpoint outcomes have become observable.
          t_eval = max(tox$t_observed, eff$t_observed),
          dose = as.integer(dose),
          y = tox$outcome,
          eff = eff$outcome,
          true_toxicity_probability = as.numeric(toxicity_probability),
          true_efficacy_probability = as.numeric(efficacy_probability),
          previous_admin_id = as.integer(previous_admin_id),
          true_toxicity_state = as.numeric(toxicity_state),
          true_efficacy_state = as.numeric(efficacy_state),
          ncycle = as.integer(cycle),
          cohort = cohort_id,
          stage = stage,
          type = type,
          stringsAsFactors = FALSE
        )
      )
    }

    for (i in seq_along(new_ids)) {
      tox_state <- aide_phase12_multicycle_next_state(
        p_true[dose], toxicity_ipde_alpha, previous_state = 0
      )
      eff_state <- aide_phase12_multicycle_next_state(
        e_true[dose], efficacy_ipde_alpha, previous_state = 0
      )
      append_row(
        id = new_ids[i],
        arrival = new_arrivals[i],
        toxicity_probability = p_true[dose],
        efficacy_probability = e_true[dose],
        toxicity_state = tox_state$state,
        efficacy_state = eff_state$state,
        previous_admin_id = 0L,
        cycle = 1L,
        type = "new"
      )
    }
    for (id in ret_ids) {
      last <- state[state$id == id, , drop = FALSE]
      tox_state <- aide_phase12_multicycle_next_state(
        p_true[dose], toxicity_ipde_alpha, last$true_toxicity_state
      )
      eff_state <- aide_phase12_multicycle_next_state(
        e_true[dose], efficacy_ipde_alpha, last$true_efficacy_state
      )
      append_row(
        id = id,
        arrival = last$t_arrival,
        toxicity_probability = if (multicycle_toxicity_active) {
          tox_state$probability
        } else {
          aide_phase12_ipde_toxicity_probability(
            p_regular = p_true, p_ipde_base = p_ipde,
            previous_dose = last$dose, current_dose = dose,
            alpha = toxicity_ipde_alpha
          )
        },
        efficacy_probability = if (multicycle_efficacy_active) {
          eff_state$probability
        } else {
          aide_phase12_ipde_efficacy_probability(
            e_regular = e_true, e_ipde_base = e_ipde,
            previous_dose = last$dose, current_dose = dose,
            alpha = efficacy_ipde_alpha
          )
        },
        toxicity_state = tox_state$state,
        efficacy_state = eff_state$state,
        previous_admin_id = last$row_id,
        cycle = last$ncycle + 1L,
        type = "retreat"
      )
    }

    decision_time <<- max(admin$t_eval[admin$cohort == cohort_id])
    if (verbose) {
      message(
        "Cohort ", cohort_id,
        " (", stage, "): dose ", dose,
        ", administrations = ", n_to_enroll,
        ", total = ", nrow(admin)
      )
    }
    TRUE
  }

  toxicity_counts <- function() {
    list(
      n = tabulate(as.integer(admin$dose), nbins = ndose),
      y = tabulate(as.integer(admin$dose[admin$y == 1L]), nbins = ndose),
      n_new = tabulate(as.integer(admin$dose[admin$type == "new"]), nbins = ndose),
      y_new = tabulate(
        as.integer(admin$dose[admin$type == "new" & admin$y == 1L]),
        nbins = ndose
      ),
      n_recycle = tabulate(as.integer(admin$dose[admin$type == "retreat"]), nbins = ndose),
      y_recycle = tabulate(
        as.integer(admin$dose[admin$type == "retreat" & admin$y == 1L]),
        nbins = ndose
      )
    )
  }

  dose_threshold_reached <- function(threshold) {
    aide_phase12_dose_threshold_reached(toxicity_counts()$n, threshold)
  }

  absorb_toxicity_model_elimination <- function(model_fit) {
    if (!is.null(model_fit$eliminated) &&
        length(model_fit$eliminated) == ndose) {
      toxicity_eliminated <<- pmax(
        toxicity_eliminated,
        as.integer(model_fit$eliminated)
      )
    }
    invisible(NULL)
  }

  ## Toxicity exclusions are obtained only from the selected BOIN/CRM model.
  ## In particular, do not add a separate beta-binomial toxicity rule here.
  update_toxicity_elimination <- function() {
    model_fit <- toxicity_utility_fit()
    absorb_toxicity_model_elimination(model_fit)
    cache_random_crm_recycle_post(model_fit)
    model_fit
  }

  toxicity_move <- function(dose) {
    counts <- toxicity_counts()
    if (toxicity_eliminated[dose] == 1L) {
      return(list(
        next_dose = if (dose == 1L) 99L else dose - 1L,
        action = if (dose == 1L) "earlystop_dose1_eliminated" else "de-escalate_elim",
        mu_hat = NA_real_
      ))
    }

    if (model == "BOIN") {
      current_rows <- admin$dose == dose
      return(boin_move(
        current_dose = dose,
        ndose = ndose,
        method = decision_method,
        y_curr = counts$y[dose],
        n_curr = counts$n[dose],
        b.e = b_e,
        b.d = b_d,
        C = C,
        YR = counts$y_new[dose],
        YI = counts$y_recycle[dose],
        NR_star = counts$n_new[dose],
        NI_star = counts$n_recycle[dose],
        lambda_e = boin_lambdas$lambda_e,
        lambda_d = boin_lambdas$lambda_d,
        r_carry = r_carry,
        elimi = toxicity_eliminated,
        n_trt_curr = sum(current_rows),
        dose_cap = dose_cap
      ))
    }

    crm_move(
      current_dose = dose,
      ndose = ndose,
      dat = toxicity_model_data(),
      target = target,
      cutoff = cutoff,
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
      previous_dose_model_file = crm_previous_dose_model_file,
      n_chains = crm_n_chains,
      n_adapt = crm_n_adapt,
      n_burnin = crm_n_burnin,
      n_iter = crm_n_iter,
      thin = crm_thin,
      elimi = toxicity_eliminated,
      n_trt_curr = counts$n[dose],
      dose_cap = dose_cap,
      no_skip = TRUE
    )
  }

  final_toxicity_fit <- function() {
    if (nrow(admin) == 0L) {
      return(list(
        MTD = 99L,
        phat = rep(NA_real_, ndose),
        pj_iso = rep(NA_real_, ndose),
        eliminated = toxicity_eliminated
      ))
    }
    counts <- toxicity_counts()
    if (model == "BOIN") {
      return(select.mtd(
        target = target,
        y = counts$y,
        n = counts$n,
        cutoff.eli = cutoff,
        approx = mtd_method,
        r_carry = r_carry,
        y_new = counts$y_new,
        n_new = counts$n_new,
        y_recycle = counts$y_recycle,
        n_recycle = counts$n_recycle,
        restrict_to_tried = TRUE
      ))
    }
    select.mtd.crm(
      target = target,
      dat = toxicity_model_data(),
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
      previous_dose_model_file = crm_previous_dose_model_file,
      n_chains = crm_n_chains,
      n_adapt = crm_n_adapt,
      n_burnin = crm_n_burnin,
      n_iter = max(5000L, crm_n_iter),
      thin = crm_thin,
      restrict_to_tried = TRUE
    )
  }

  ## Refit the user-selected toxicity model without restricting the estimate
  ## vector to previously treated doses. This supplies model-based toxicity
  ## probabilities for every dose that can enter the one-stage utility rule.
  toxicity_utility_fit <- function() {
    if (nrow(admin) == 0L) {
      return(list(
        phat = rep(NA_real_, ndose),
        eliminated = toxicity_eliminated
      ))
    }
    counts <- toxicity_counts()
    if (model == "BOIN") {
      return(select.mtd(
        target = target,
        y = counts$y,
        n = counts$n,
        cutoff.eli = cutoff,
        approx = mtd_method,
        r_carry = r_carry,
        y_new = counts$y_new,
        n_new = counts$n_new,
        y_recycle = counts$y_recycle,
        n_recycle = counts$n_recycle,
        restrict_to_tried = FALSE
      ))
    }
    select.mtd.crm(
      target = target,
      dat = toxicity_model_data(),
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
      previous_dose_model_file = crm_previous_dose_model_file,
      n_chains = crm_n_chains,
      n_adapt = crm_n_adapt,
      n_burnin = crm_n_burnin,
      n_iter = max(5000L, crm_n_iter),
      thin = crm_thin,
      restrict_to_tried = FALSE
    )
  }

  toxicity_estimate_from_fit <- function(model_fit) {
    estimate <- NULL
    for (name in c("p_hat", "phat", "pj_iso")) {
      if (!is.null(model_fit[[name]])) {
        estimate <- as.numeric(model_fit[[name]])
        break
      }
    }
    if (is.null(estimate) || length(estimate) != ndose ||
        any(!is.na(estimate) & (!is.finite(estimate) |
                                 estimate < 0 | estimate > 1))) {
      stop("The selected toxicity model did not return valid per-dose estimates for utility allocation.")
    }
    estimate
  }

  record_decision <- function(stage, dose, tox_move_out, allocated_dose, futility) {
    cohort_rows <- admin$cohort == cohort_id
    decision_log <<- rbind(
      decision_log,
      data.frame(
        cohort = cohort_id,
        stage = stage,
        current_dose = as.integer(dose),
        toxicity_next_dose = as.integer(tox_move_out$next_dose),
        allocated_dose = as.integer(allocated_dose),
        toxicity_action = as.character(tox_move_out$action),
        n_current = as.integer(toxicity_counts()$n[dose]),
        efficacy_futility_eliminated = sum(futility$futility_eliminated),
        n_new_enrolled = sum(cohort_rows & admin$type == "new"),
        n_recycled_enrolled = sum(cohort_rows & admin$type == "retreat"),
        waiting_queue_size = length(waiting_new_ids()),
        stringsAsFactors = FALSE
      )
    )
  }

  ## Existing CRM backends receive their original four-column view. The new
  ## arbitrary-cycle backend receives the administration history needed to
  ## reconstruct same-patient state pointers.
  toxicity_model_data <- function() {
    if (multicycle_toxicity_active) {
      admin
    } else {
      admin[, c("id", "dose", "y", "type"), drop = FALSE]
    }
  }

  observed_efficacy_response <- function() {
    tabulate(
      as.integer(admin$dose[admin$eff == 1L]),
      nbins = ndose
    ) > 0L
  }

  current_utility_allocation <- function(toxicity_fit, futility, current_dose) {
    MTD <- as.integer(toxicity_fit$MTD)
    admissible <- seq_len(ndose) <= MTD &
      toxicity_eliminated == 0L &
      futility$futility_eliminated == 0L
    utility <- efficacy_utility_current(
      toxicity_estimate_from_fit(toxicity_fit)
    )$utility
    one_stage <- aide_phase12_one_stage_allocation(
      current_dose = current_dose,
      mtd = MTD,
      admissible_doses = which(admissible),
      utility = utility,
      response_observed = observed_efficacy_response(),
      tried_doses = unique(as.integer(admin$dose))
    )
    list(
      MTD = MTD,
      admissible = which(admissible),
      utility = utility,
      provisional_OBD = one_stage$target,
      one_stage = one_stage
    )
  }

  choose_one_stage_dose <- function(allocation_state) {
    allocation_state$one_stage
  }

  ## Stage I remains toxicity-directed.  Futility cannot redirect an AIDE
  ## recommendation to a higher dose, but a futile recommended dose is
  ## removed before allocation and the closest lower admissible dose is used.
  choose_stage1_dose <- function(dose, tox_move_out, futility) {
    if (isTRUE(tox_move_out$stop_trial) || tox_move_out$next_dose == 99L) {
      return(99L)
    }
    admissible <- toxicity_eliminated == 0L &
      futility$futility_eliminated == 0L
    recommended <- as.integer(tox_move_out$next_dose)
    if (recommended >= 1L && recommended <= ndose && admissible[recommended]) {
      return(recommended)
    }
    if (recommended > dose) {
      candidates <- seq.int(dose, min(recommended, ndose))
    } else {
      candidates <- seq_len(max(1L, min(recommended, ndose)))
    }
    candidates <- candidates[admissible[candidates]]
    if (length(candidates) == 0L) return(99L)
    max(candidates)
  }

  stage1_fit <- NULL
  stage1_futility <- NULL
  stage1_admissible <- rep(FALSE, ndose)
  stage1_transition_dose <- NA_integer_

  if (allocation == "two_stage") {
    ## Stage I: standard AIDE toxicity allocation. Efficacy does not rank
    ## doses.  A dose is assessed for efficacy futility only after it has just
    ## been treated; once eliminated, it is excluded from later allocation and
    ## IPDE screening. Stage I ends only after the current dose has
    ## accumulated N_s1 administrations and the allocation decision remains
    ## at that dose; reaching N_s1 alone does not trigger the transition.
    while (nrow(admin) < Nmax && !earlystop) {
      futility_before <- futility_state()
      if (apply_early_stop_rule("stage1")) break
      if (toxicity_eliminated[current_dose] == 1L ||
          futility_before$futility_eliminated[current_dose] == 1L) {
        stop_reason <- "no_safe_nonfutile_stage1_current_dose"
        break
      }
      n_to_enroll <- min(
        C,
        Nmax - nrow(admin)
      )
      if (!enroll_cohort(
        current_dose, "stage1", n_to_enroll, futility = futility_before
      )) {
        stop_reason <- "ran_out_of_new_patients_stage1"
        break
      }
      update_toxicity_elimination()
      futility_after <- evaluate_current_dose_futility(current_dose)
      if (apply_early_stop_rule("stage1")) break

      tox_move_out <- toxicity_move(current_dose)
      absorb_toxicity_model_elimination(tox_move_out)
      if (isTRUE(tox_move_out$stop_trial) || tox_move_out$next_dose == 99L) {
        if (apply_early_stop_rule("stage1")) break
        stop_reason <- "no_safe_nonfutile_stage1_candidate"
        break
      }
      allocated_dose <- choose_stage1_dose(
        current_dose, tox_move_out, futility_after
      )
      record_decision(
        "stage1",
        current_dose,
        tox_move_out,
        allocated_dose,
        futility_after
      )
      if (allocated_dose == 99L) {
        if (apply_early_stop_rule("stage1")) break
        stop_reason <- "no_safe_nonfutile_stage1_candidate"
        break
      }
      if (sum(admin$stage == "stage1") >= N_s1 &&
          isTRUE(tox_move_out$action == "stay") &&
          tox_move_out$next_dose == current_dose &&
          allocated_dose == current_dose) {
        stage1_transition_dose <- as.integer(current_dose)
        break
      }
      current_dose <- allocated_dose
    }

    stage1_fit <- final_toxicity_fit()
    stage1_futility <- efficacy_summary_current()
    stage1_mtd <- stage1_fit$MTD
    stage1_tox_elim <- if (!is.null(stage1_fit$eliminated)) {
      as.integer(stage1_fit$eliminated)
    } else {
      toxicity_eliminated
    }
    if (!earlystop && !is.na(stage1_mtd) && stage1_mtd >= 1L && stage1_mtd <= ndose) {
      stage1_admissible <- seq_len(ndose) <= stage1_mtd &
        stage1_futility$futility_eliminated == 0L &
        stage1_tox_elim == 0L
    }

    ## Stage II reuses the one-stage allocation rule after every evaluated
    ## cohort. There is no separate Stage II allocation mode or per-dose cap.
    stage2_admissible <- stage1_admissible
    while (!earlystop && !is.na(stage1_transition_dose) && nrow(admin) < Nmax) {
      toxicity_fit_now <- update_toxicity_elimination()
      eff_now <- efficacy_summary_current()
      if (apply_early_stop_rule("stage2")) break
      allocation_state <- current_utility_allocation(
        toxicity_fit_now, eff_now, current_dose
      )
      stage2_admissible <- allocation_state$admissible
      if (!length(stage2_admissible)) {
        stop_reason <- "no_stage2_admissible_dose"
        break
      }
      stage2_choice <- choose_one_stage_dose(allocation_state)
      allocated_dose <- stage2_choice$dose
      if (allocated_dose == 99L) {
        stop_reason <- "no_safe_nonfutile_stage2_candidate"
        break
      }
      n_to_enroll <- min(C, Nmax - nrow(admin))
      if (!enroll_cohort(
        allocated_dose, "stage2", n_to_enroll, futility = eff_now
      )) {
        stop_reason <- "ran_out_of_new_patients_stage2"
        break
      }
      update_toxicity_elimination()
      futility_after <- evaluate_current_dose_futility(allocated_dose)
      record_decision(
        "stage2",
        current_dose,
        list(next_dose = allocated_dose, action = stage2_choice$action),
        allocated_dose,
        futility_after
      )
      if (apply_early_stop_rule("stage2")) break
      current_dose <- allocated_dose
    }
    if (!earlystop && nrow(admin) < Nmax && !any(stage2_admissible)) {
      stop_reason <- "no_stage2_admissible_dose"
    }
  } else {
    ## One-stage design: every completed cohort first obtains an AIDE toxicity
    ## decision, then utility selects from the prescribed local candidate set.
    ## Efficacy futility is assessed only at the dose just treated.
    while (nrow(admin) < Nmax && !earlystop) {
      futility_before <- futility_state()
      if (apply_early_stop_rule("one_stage")) break
      if (toxicity_eliminated[current_dose] == 1L ||
          futility_before$futility_eliminated[current_dose] == 1L) {
        stop_reason <- "no_safe_nonfutile_one_stage_current_dose"
        break
      }
      n_to_enroll <- min(C, Nmax - nrow(admin))
      if (!enroll_cohort(
        current_dose, "one_stage", n_to_enroll, futility = futility_before
      )) {
        stop_reason <- "ran_out_of_new_patients_one_stage"
        break
      }
      toxicity_fit_now <- update_toxicity_elimination()
      futility_now <- evaluate_current_dose_futility(current_dose)
      if (apply_early_stop_rule("one_stage")) break
      if (nrow(admin) >= Nmax) break

      tox_move_out <- toxicity_move(current_dose)
      absorb_toxicity_model_elimination(tox_move_out)
      if (isTRUE(tox_move_out$stop_trial) || tox_move_out$next_dose == 99L) {
        if (apply_early_stop_rule("one_stage")) break
        stop_reason <- "no_safe_nonfutile_one_stage_candidate"
        break
      }
      allocation_state <- current_utility_allocation(
        toxicity_fit_now, futility_now, current_dose
      )
      one_stage_choice <- choose_one_stage_dose(allocation_state)
      allocated_dose <- one_stage_choice$dose
      record_decision(
        "one_stage",
        current_dose,
        list(next_dose = allocated_dose, action = one_stage_choice$action),
        allocated_dose,
        futility_now
      )
      if (allocated_dose == 99L) {
        if (apply_early_stop_rule("one_stage")) break
        stop_reason <- "no_safe_nonfutile_one_stage_candidate"
        break
      }
      current_dose <- allocated_dose
    }
  }

  final_fit <- final_toxicity_fit()
  final_mtd <- final_fit$MTD
  final_futility <- efficacy_summary_current()
  ## These dose-wise beta-binomial probabilities are calculated only for the
  ## completed-trial report.  They do not alter the persistent interim
  ## elimination state used for allocation or early stopping.
  final_observed_futility <- aide_phase12_beta_binomial_futility(
    admin = admin,
    ndose = ndose,
    efficacy_threshold = efficacy_threshold,
    futility_cutoff = futility_cutoff,
    min_eff_n_for_futility = min_eff_n_for_futility
  )
  final_futility$prob_below_threshold <- final_observed_futility$prob_below_threshold
  final_futility$prob_above_threshold <- final_observed_futility$prob_above_threshold
  final_futility$efficacy_futility_n <- final_observed_futility$n
  final_futility$efficacy_futility_y <- final_observed_futility$y
  final_utility_fit <- toxicity_utility_fit()
  final_utility <- efficacy_utility_current(
    toxicity_estimate_from_fit(final_utility_fit)
  )
  final_tox_elim <- if (!is.null(final_fit$eliminated)) {
    as.integer(final_fit$eliminated)
  } else {
    toxicity_eliminated
  }

  final_admissible <- rep(FALSE, ndose)
  tried_dose <- tabulate(as.integer(admin$dose), nbins = ndose) > 0L
  response_observed_final <- observed_efficacy_response()
  if (!is.na(final_mtd) && final_mtd >= 1L && final_mtd <= ndose && !earlystop) {
    final_admissible <- tried_dose & seq_len(ndose) <= final_mtd &
      final_tox_elim == 0L &
      final_futility$futility_eliminated == 0L &
      response_observed_final
  }
  final_candidates <- which(final_admissible & is.finite(final_utility$utility))
  final_obd <- if (length(final_candidates) == 0L) {
    99L
  } else {
    aide_phase12_select_lowest_tied(
      final_utility$utility,
      final_candidates
    )
  }

  trial_time <- if (nrow(admin) == 0L) {
    NA_real_
  } else {
    max(admin$t_eval) - min(admin$t_arrival)
  }
  waiting_queue <- arrival_schedule[
    arrival_schedule$t_arrival <= decision_time &
      !arrival_schedule$id %in% enrolled_new_ids(),
    , drop = FALSE
  ]
  pending_arrivals <- arrival_schedule[
    arrival_schedule$t_arrival > decision_time,
    , drop = FALSE
  ]
  list(
    admin = admin,
    decision_log = decision_log,
    recycling_decision_log = recycling_decision_log,
    arrival_schedule = arrival_schedule,
    waiting_queue = waiting_queue,
    pending_arrivals = pending_arrivals,
    stage1 = list(
      N_s1 = as.integer(N_s1),
      n_by_dose = tabulate(
        as.integer(admin$dose[admin$stage == "stage1"]),
        nbins = ndose
      ),
      threshold_reached = !is.na(stage1_transition_dose),
      transition_dose = stage1_transition_dose,
      MTD = if (is.null(stage1_fit)) NA_integer_ else stage1_fit$MTD,
      toxicity_fit = stage1_fit,
      efficacy = stage1_futility,
      admissible = stage1_admissible
    ),
    stage2 = list(
      N_s2 = NA_integer_,
      n_by_dose = tabulate(
        as.integer(admin$dose[admin$stage == "stage2"]),
        nbins = ndose
      ),
      threshold_reached = nrow(admin) >= Nmax
    ),
    final = list(
      allocation = allocation,
      MTD = as.integer(final_mtd),
      OBD = as.integer(final_obd),
      admissible = final_admissible,
      tried_dose = tried_dose,
      response_observed = response_observed_final,
      toxicity_fit = final_fit,
      toxicity_utility_fit = final_utility_fit,
      efficacy = final_futility,
      utility = final_utility,
      toxicity_eliminated = final_tox_elim,
      efficacy_futility_eliminated = final_futility$futility_eliminated,
      earlystop = as.integer(earlystop),
      stop_reason = stop_reason,
      n_admin = nrow(admin),
      n_unique_patients = length(unique(admin$id)),
      n_waiting_new_patients = nrow(waiting_queue),
      n_pending_new_arrivals = nrow(pending_arrivals),
      trial_time = trial_time,
      model = model,
      carryover_model = carryover_model,
      toxicity_carryover_model = if (crm_r_model == "random") {
        "discount_shared"
      } else if (crm_r_model == "previous_dose") {
        "additive_shared"
      } else if (crm_r_model == "multicycle_additive") {
        "multicycle_additive"
      } else {
        NA_character_
      },
      efficacy_carryover_model = efficacy_model,
      decision_method = if (model == "BOIN") decision_method else paste0("crm_", crm_r_model),
      N_s1 = as.integer(N_s1),
      N_s2 = as.integer(N_s2),
      m_U = m_U,
      utility_type = utility_type,
      lambda_T = lambda_T,
      utility_scores = utility_scores,
      recycling_rules = list(
        ipde_design = as.integer(ipde_design),
        flexible_ipde = flexible_ipde,
        enrollment_scheme = enrollment_scheme,
        enrollment_priority = enrollment_priority,
        active = random_crm_recycle_rules_active,
        toxicity = list(
          enabled = apply_random_crm_recycle_toxicity_rule,
          phi = target,
          cutoff = random_crm_recycle_toxicity_cutoff,
          true_generation = if (multicycle_toxicity_active) {
            "min(1, p_true[d_current] + toxicity_ipde_alpha * prior_uncapped_toxicity_state)"
          } else {
            "min(1, p_ipde[d2] + toxicity_ipde_alpha * p_true[d1])"
          },
          ipde_alpha = toxicity_ipde_alpha,
          ipde_base = p_ipde,
          event_generator = "gen.tite",
          dist = dlt_tite$dist,
          alpha = dlt_tite$alpha,
          assessment_window = T_assess
        ),
        efficacy = list(
          enabled = apply_random_crm_recycle_efficacy_rule,
          delta = random_crm_recycle_efficacy_delta,
          model = efficacy_model,
          regular_prior = efficacy_prior,
          carryover_prior = efficacy_carryover_prior,
          additive_alpha_prior = efficacy_additive_alpha_prior,
          jags_model_file = efficacy_model_file,
          jags_n_chains = efficacy_n_chains,
          jags_n_adapt = efficacy_n_adapt,
          jags_n_burnin = efficacy_n_burnin,
          jags_n_iter = efficacy_n_iter,
          jags_thin = efficacy_thin,
          true_generation = if (multicycle_efficacy_active) {
            "min(1, e_true[d_current] + efficacy_ipde_alpha * prior_uncapped_efficacy_state)"
          } else {
            "min(1, e_ipde[d2] + efficacy_ipde_alpha * e_true[d1])"
          },
          ipde_alpha = efficacy_ipde_alpha,
          ipde_base = e_ipde,
          event_generator = "gen.tite",
          dist = efficacy_tite$dist,
          alpha = efficacy_tite$alpha,
          assessment_window = T_assess
        )
      )
    )
  )
}

get_oc_sim_AIDE_phase_I_II <- function(p_true,
                                        e_true,
                                        ntrial = 1000L,
                                        seed = 1L,
                                        store_raw = FALSE,
                                        ...) {
  p_true <- as.numeric(p_true)
  e_true <- as.numeric(e_true)
  if (length(p_true) != length(e_true) || length(p_true) == 0L) {
    stop("p_true and e_true must have the same positive length.")
  }
  if (length(ntrial) != 1L || ntrial < 1L || ntrial != as.integer(ntrial)) {
    stop("ntrial must be a positive integer.")
  }
  ndose <- length(p_true)
  mtd <- rep(NA_integer_, ntrial)
  obd <- rep(NA_integer_, ntrial)
  stopped <- integer(ntrial)
  n_admin <- numeric(ntrial)
  n_unique <- numeric(ntrial)
  duration <- rep(NA_real_, ntrial)
  n_by_dose_by_trial <- matrix(0, nrow = ntrial, ncol = ndose)
  unique_n_by_dose_by_trial <- matrix(0, nrow = ntrial, ncol = ndose)
  nipde_by_dose_by_trial <- matrix(0, nrow = ntrial, ncol = ndose)
  crm_pj_by_trial <- matrix(NA_real_, nrow = ntrial, ncol = ndose)
  efficacy_pj_by_trial <- matrix(NA_real_, nrow = ntrial, ncol = ndose)
  utility_by_trial <- matrix(NA_real_, nrow = ntrial, ncol = ndose)
  r_hat_by_trial <- matrix(NA_real_, nrow = ntrial, ncol = ndose)
  obd_selection <- integer(ndose)
  raw <- if (isTRUE(store_raw)) vector("list", ntrial) else NULL

  for (i in seq_len(ntrial)) {
    fit <- simulate_AIDE_phase_I_II(
      p_true = p_true,
      e_true = e_true,
      seed = seed + i - 1L,
      ...
    )
    mtd[i] <- fit$final$MTD
    obd[i] <- fit$final$OBD
    stopped[i] <- fit$final$earlystop
    n_admin[i] <- fit$final$n_admin
    n_unique[i] <- fit$final$n_unique_patients
    duration[i] <- fit$final$trial_time
    if (nrow(fit$admin) > 0L) {
      n_by_dose_by_trial[i, ] <- tabulate(fit$admin$dose, nbins = ndose)
      unique_n_by_dose_by_trial[i, ] <- tabulate(
        unique(fit$admin[, c("id", "dose"), drop = FALSE])$dose,
        nbins = ndose
      )
      ipde_rows <- fit$admin$type == "retreat"
      if (any(ipde_rows)) {
        nipde_by_dose_by_trial[i, ] <- tabulate(
          fit$admin$dose[ipde_rows], nbins = ndose
        )
      }
    }
    crm_pj_by_trial[i, ] <- as.numeric(fit$final$utility$toxicity)
    efficacy_pj_by_trial[i, ] <- as.numeric(fit$final$utility$efficacy)
    utility_by_trial[i, ] <- as.numeric(fit$final$utility$utility)
    r_hat <- as.numeric(fit$final$toxicity_utility_fit$r_hat)
    if (length(r_hat) == 1L && is.finite(r_hat)) {
      r_hat_by_trial[i, ] <- rep(r_hat, ndose)
    } else if (length(r_hat) == ndose) {
      r_hat_by_trial[i, ] <- r_hat
    }
    if (!is.na(obd[i]) && obd[i] >= 1L && obd[i] <= ndose) {
      obd_selection[obd[i]] <- obd_selection[obd[i]] + 1L
    }
    if (isTRUE(store_raw)) raw[[i]] <- fit
  }

  ## Keep the OC vector consistent with the trial-level no-OBD sentinel.
  obd[is.na(obd)] <- 99L

  list(
    p_true = p_true,
    e_true = e_true,
    ntrial = as.integer(ntrial),
    MTD_by_trial = mtd,
    OBD_by_trial = obd,
    early_stop_by_trial = stopped,
    n_administrations_by_trial = n_admin,
    n_unique_patients_by_trial = n_unique,
    duration_by_trial = duration,
    n_by_dose_by_trial = n_by_dose_by_trial,
    unique_n_by_dose_by_trial = unique_n_by_dose_by_trial,
    nipde_by_dose_by_trial = nipde_by_dose_by_trial,
    crm_pj_by_trial = crm_pj_by_trial,
    efficacy_pj_by_trial = efficacy_pj_by_trial,
    utility_by_trial = utility_by_trial,
    r_hat_by_trial = r_hat_by_trial,
    MTD_selection_percent = 100 * tabulate(
      mtd[!is.na(mtd) & mtd >= 1L & mtd <= ndose], nbins = ndose
    ) / ntrial,
    OBD_selection_percent = 100 * obd_selection / ntrial,
    No_OBD_selection_percent = 100 * mean(obd == 99L),
    early_stop_percent = 100 * mean(stopped),
    mean_administrations = mean(n_admin),
    mean_administrations_by_dose = colMeans(n_by_dose_by_trial),
    mean_unique_patients = mean(n_unique),
    mean_unique_patients_by_dose = colMeans(unique_n_by_dose_by_trial),
    mean_ipde_doses_by_dose = colMeans(nipde_by_dose_by_trial),
    mean_duration = mean(duration, na.rm = TRUE),
    raw = raw
  )
}
