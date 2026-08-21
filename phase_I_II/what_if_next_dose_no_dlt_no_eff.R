## ============================================================
## One-step non-TITE AIDE next-dose calculation
##
## The fixed history has three completed new patients at every dose level,
## all with no DLT (y = 0) and no efficacy response (eff = 0). The script
## refits the same JAGS CRM and efficacy models used in AIDE_phase_I_II.R
## 100 times and demonstrates the common no-response rule: both the one-stage
## design and Stage II of the two-stage design target the current MTD, subject
## to the universal no-skipping constraint.
##
## Run from the project root:
##   Rscript phase_I_II/what_if_next_dose_no_dlt_no_eff.R
## ============================================================

source("phase_I_II/AIDE_phase_I_II.R")
`%||%` <- function(x, y) if (is.null(x)) y else x

if (!requireNamespace("rjags", quietly = TRUE) ||
    !requireNamespace("coda", quietly = TRUE)) {
  stop("This script requires rjags, coda, and a working JAGS installation.")
}

## ------------------------------------------------------------
## Settings: match run_oc_AIDE_phase_I_II.R where applicable.
## ------------------------------------------------------------

ndose <- 5L
patients_per_dose <- 3L
nrep <- 100L
seed_base <- 9100L
target <- .30
cutoff <- .95
crm_r_model <- "previous_dose"
crm_skeleton <- c(.15, .20, .30, .35, .45)
crm_alpha_sd <- sqrt(2)
crm_a_r <- .15
crm_b_r <- .85
crm_n_chains <- 2L
crm_n_adapt <- 500L
crm_n_burnin <- 500L
crm_n_iter <- 2000L
crm_thin <- 1L
efficacy_prior <- c(.5, .5)
efficacy_carryover_prior <- c(.15, .85)
efficacy_model <- "previous_dose_additive"
efficacy_additive_alpha_prior <- c(.15, .85)
efficacy_n_chains <- 3L
efficacy_n_adapt <- 1000L
efficacy_n_burnin <- 1000L
efficacy_n_iter <- 4000L
efficacy_thin <- 2L
efficacy_threshold <- .20
futility_cutoff <- .85
min_eff_n_for_futility <- 0L
utility_type <- 3L
lambda_T <- 1
utility_scores <- c(u00 = 40, u01 = 100, u10 = 0, u11 = 60)

## Three no-event patients at every dose: D1, D1, D1, D2, D2, D2, ..., D5.
admin <- data.frame(
  id = seq_len(ndose * patients_per_dose),
  dose = rep(seq_len(ndose), each = patients_per_dose),
  y = 0L,
  eff = 0L,
  type = "new",
  stringsAsFactors = FALSE
)

extract_toxicity_estimate <- function(fit) {
  for (name in c("p_hat", "phat", "pj_iso")) {
    if (!is.null(fit[[name]])) return(as.numeric(fit[[name]]))
  }
  stop("The CRM fit did not return dose-specific toxicity estimates.")
}

one_step_models <- function(seed) {
  set.seed(seed)
  toxicity_fit <- select.mtd.crm(
    target = target,
    dat = admin[, c("id", "dose", "y", "type"), drop = FALSE],
    ndose = ndose,
    skeleton = crm_skeleton,
    cutoff.eli = cutoff,
    r_model = crm_r_model,
    a_r = crm_a_r,
    b_r = crm_b_r,
    alpha_sd = crm_alpha_sd,
    previous_dose_model_file = "previous_dose_additive_CRM.bug",
    n_chains = crm_n_chains,
    n_adapt = crm_n_adapt,
    n_burnin = crm_n_burnin,
    n_iter = crm_n_iter,
    thin = crm_thin,
    restrict_to_tried = FALSE
  )
  efficacy_fit <- aide_phase12_efficacy_posterior(
    admin = admin,
    ndose = ndose,
    efficacy_prior = efficacy_prior,
    efficacy_carryover_prior = efficacy_carryover_prior,
    efficacy_model = efficacy_model,
    efficacy_additive_alpha_prior = efficacy_additive_alpha_prior,
    n_chains = efficacy_n_chains,
    n_adapt = efficacy_n_adapt,
    n_burnin = efficacy_n_burnin,
    n_iter = efficacy_n_iter,
    thin = efficacy_thin
  )
  futility <- aide_phase12_beta_binomial_futility(
    admin = admin,
    ndose = ndose,
    efficacy_threshold = efficacy_threshold,
    futility_cutoff = futility_cutoff,
    min_eff_n_for_futility = min_eff_n_for_futility
  )
  toxicity <- extract_toxicity_estimate(toxicity_fit)
  utility <- aide_phase12_compute_utility(
    efficacy = efficacy_fit$regular_posterior_mean,
    toxicity = toxicity,
    utility_type = utility_type,
    lambda_T = lambda_T,
    utility_scores = utility_scores
  )
  eliminated <- toxicity_fit$eliminated %||% rep.int(0L, ndose)
  MTD <- as.integer(toxicity_fit$MTD)
  admissible <- seq_len(ndose) <= MTD &
    as.integer(eliminated) == 0L &
    futility$futility_eliminated == 0L
  list(MTD = MTD, utility = utility, admissible = which(admissible),
       efficacy_futility = futility$futility_eliminated)
}

next_one_stage <- function(model_state) {
  response_observed <- tabulate(
    as.integer(admin$dose[admin$eff == 1L]), nbins = ndose
  ) > 0L
  aide_phase12_one_stage_allocation(
    current_dose = max(admin$dose),
    mtd = model_state$MTD,
    admissible_doses = model_state$admissible,
    utility = model_state$utility,
    response_observed = response_observed,
    tried_doses = unique(admin$dose)
  )$dose
}

next_two_stage <- function(model_state) {
  ## After Stage I, the two-stage design invokes exactly the same rule.
  next_one_stage(model_state)
}

results <- lapply(seq_len(nrep), function(i) {
  state <- one_step_models(seed_base + i)
  data.frame(
    replicate = i,
    MTD = state$MTD,
    one_stage_next_dose = next_one_stage(state),
    two_stage_next_dose = next_two_stage(state),
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, results)

summary <- data.frame(
  design = c("one_stage", "two_stage"),
  average_next_dose = c(
    mean(results$one_stage_next_dose),
    mean(results$two_stage_next_dose)
  ),
  minimum_next_dose = c(
    min(results$one_stage_next_dose),
    min(results$two_stage_next_dose)
  ),
  maximum_next_dose = c(
    max(results$one_stage_next_dose),
    max(results$two_stage_next_dose)
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Non-TITE AIDE one-step next-dose check\n")
cat("Fixed data: 3 patients at every dose, all with no DLT and no efficacy\n")
cat("Replicates:", nrep, "\n")
cat("============================================================\n\n")
print(summary, row.names = FALSE)

cat("\nReplicate-level next moves\n")
print(results, row.names = FALSE)
