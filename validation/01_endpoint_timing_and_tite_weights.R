args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
validation_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg[1L], winslash = "/"))
} else {
  source_files <- vapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value)) "" else as.character(value)[1L]
  }, character(1))
  source_file <- tail(source_files[nzchar(source_files)], 1L)
  if (length(source_file)) dirname(normalizePath(source_file, winslash = "/")) else normalizePath(getwd(), winslash = "/")
}

source(file.path(dirname(validation_dir), "TITE-AIDE-Rebuild", "TITE-AIDE.R"))
source(file.path(validation_dir, "validation_models.R"))

T_assess <- 28
admin <- data.frame(
  admin_id = 1:3,
  patient_id = 1:3,
  cycle = 1L,
  t_start = c(0, 0, 0),
  dose = c(1L, 1L, 1L),
  previous_dose = c(1L, 1L, 1L),
  assignment_type = c("new", "new", "new"),
  t_dlt = c(Inf, 10, Inf),
  t_response = c(Inf, Inf, 12),
  assessment_end = c(T_assess, T_assess, T_assess),
  dlt_final = c(0L, 1L, 0L),
  eff_final = c(0L, 0L, 1L),
  stringsAsFactors = FALSE
)

status_times <- c(7, 10, 12, 14, 28)
endpoint_statuses <- do.call(rbind, lapply(status_times, function(t_now) {
  toxicity <- aide_phase12_tite_toxicity_data(admin, t_now, T_assess)
  efficacy <- aide_phase12_tite_efficacy_data(admin, t_now, T_assess)
  rbind(
    data.frame(t_now = t_now, endpoint = "toxicity", toxicity[, c("admin_id", "y", "weight", "ascertained")]),
    data.frame(t_now = t_now, endpoint = "efficacy", efficacy[, c("admin_id", "y", "weight", "ascertained")])
  )
}))

weight_schedule <- endpoint_statuses[
  endpoint_statuses$admin_id == 1L & endpoint_statuses$t_now %in% c(7, 14, 28),
  c("t_now", "endpoint", "weight")
]

theta <- 0.30
likelihood_contributions <- data.frame(
  endpoint = c("toxicity", "efficacy"),
  outcome = c("pending non-DLT", "pending nonresponse"),
  t_now = c(14, 14),
  weight = c(0.50, 0.50),
  theta = c(theta, theta),
  jags_likelihood = c(1 - 0.50 * theta, (1 - theta)^0.50)
)

timing_config <- aide_phase12_config(
  Nmax = 3L,
  cohort_size = 3L,
  n_eval = 1L,
  time = list(arrival_rate = 0.2, t0 = 0, T_assess = T_assess),
  toxicity = list(
    model = "discount_r",
    target = 0.30,
    cutoff = 0.95,
    skeleton = 0.30,
    beta_prior_mean = 0,
    beta_prior_sd = sqrt(2),
    carryover_prior = c(1, 9),
    n_chains = 2L,
    n_adapt = 500L,
    n_burnin = 500L,
    n_iter = 4000L,
    thin = 1L
  ),
  efficacy = list(
    model = "shared_carryover",
    prior = c(1, 1),
    carryover_prior = c(1, 9),
    threshold = 0.20,
    futility_cutoff = 0.90,
    min_eff_n_for_futility = 3L,
    n_chains = 2L,
    n_adapt = 500L,
    n_burnin = 500L,
    n_iter = 4000L,
    thin = 1L
  )
)

fit_times <- c(7, 14, 28)
timing_fits <- lapply(fit_times, function(t_now) {
  toxicity_data <- aide_phase12_tite_toxicity_data(admin, t_now, T_assess)
  efficacy_data <- aide_phase12_tite_efficacy_data(admin, t_now, T_assess)
  set.seed(20260808L + t_now)
  toxicity_fit <- aide_fit_toxicity(toxicity_data, timing_config, ndose = 1L)
  set.seed(20260908L + t_now)
  efficacy_fit <- aide_fit_efficacy(efficacy_data, timing_config, ndose = 1L)
  list(
    t_now = t_now,
    toxicity_data = toxicity_data,
    efficacy_data = efficacy_data,
    toxicity_fit = toxicity_fit,
    efficacy_fit = efficacy_fit
  )
})

timing_posterior_means <- do.call(rbind, lapply(timing_fits, function(x) {
  rbind(
    data.frame(t_now = x$t_now, endpoint = "toxicity", p_regular_mean = x$toxicity_fit$p_regular_mean,
               p_ipde_mean = x$toxicity_fit$p_ipde_mean),
    data.frame(t_now = x$t_now, endpoint = "efficacy", p_regular_mean = x$efficacy_fit$p_regular_mean,
               p_ipde_mean = x$efficacy_fit$p_ipde_mean)
  )
}))

results_dir <- file.path(validation_dir, "results")
dir.create(results_dir, showWarnings = FALSE)
write.csv(endpoint_statuses, file.path(results_dir, "endpoint_statuses.csv"), row.names = FALSE)
write.csv(weight_schedule, file.path(results_dir, "weight_schedule.csv"), row.names = FALSE)
write.csv(likelihood_contributions, file.path(results_dir, "pending_outcome_likelihoods.csv"), row.names = FALSE)
write.csv(timing_posterior_means, file.path(results_dir, "endpoint_timing_posterior_means.csv"), row.names = FALSE)
saveRDS(timing_fits, file.path(results_dir, "endpoint_timing_fits.rds"))

print(weight_schedule)
print(likelihood_contributions)
print(timing_posterior_means)
