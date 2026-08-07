aide_tite_fmt <- function(x, digits = 4L) {
  out <- format(round(as.numeric(x), digits), scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out); out <- sub("\\.$", "", out)
  gsub("\\.", "p", gsub("-", "m", out))
}

aide_tite_cross_join <- function(x, y) merge(x, y, by = NULL, sort = FALSE)

aide_tite_oc_settings <- function() list(
  scenario_file = file.path("..", "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"),
  scenario_ids = 1:37,
  seed_base = 1L, ntrial = 1L, cohort_size = 3L, cycle_max = 2L,
  T_assess = 28, n_eval = 3L, m_U = 6L,
  skeleton = c(.15, .20, .30, .35, .45),
  utility_scores = c(u00 = 40, u01 = 100, u10 = 0, u11 = 60),
  design_grid = rbind(
    data.frame(allocation = "two_stage", Nmax = 30L, s1_Max = 6L, N_s2 = 30L),
    data.frame(allocation = "one_stage", Nmax = 30L, s1_Max = 30L, N_s2 = 30L)
  ),
  model_grid = data.frame(
    model_id = c("discount_shared", "discount_dose_specific", "additive_shared", "additive_dose_specific"),
    toxicity_model = c("discount_r", "discount_r", "additive_alpha", "additive_alpha"),
    efficacy_model = c("shared_carryover", "dose_specific_carryover", "previous_dose_additive", "dose_specific_previous_dose_additive")
  ),
  toxicity_prior_grid = data.frame(toxicity_a = .15, toxicity_b = .85),
  efficacy_prior_grid = data.frame(efficacy_a = .5, efficacy_b = .5, carry_a = .15, carry_b = .85,
    additive_a = c(.15, .3, .5, 1), additive_b = c(.85, .7, .5, 1)),
  utility_grid = rbind(data.frame(utility_type = 2L, lambda_T = .3), data.frame(utility_type = 3L, lambda_T = 1)),
  arrival_grid = data.frame(arrival_rate = 1 / 56),
  alpha_grid = data.frame(toxicity_ipde_alpha = c(0, .3, .6, .9), efficacy_ipde_alpha = c(0, .3, .6, .9))
)

aide_tite_read_truth <- function(file, scenario_ids) {
  if (!file.exists(file)) stop("Scenario summary does not exist: ", file)
  x <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("Scenario", "Source_Scenario", "Scenario_Group", "Attempt", "True_MTD_Level", "Target_Toxicity",
            paste0("Tox_Dose", 1:5), paste0("Eff_Dose", 1:5))
  if (length(setdiff(need, names(x)))) stop("Scenario summary is missing: ", paste(setdiff(need, names(x)), collapse = ", "))
  x <- x[x$Scenario %in% scenario_ids, , drop = FALSE]
  if (nrow(x) != length(scenario_ids)) stop("Requested scenario IDs are missing from the scenario file.")
  x
}

aide_tite_setting_grid <- function(settings) {
  out <- Reduce(aide_tite_cross_join, list(settings$design_grid, settings$model_grid, settings$toxicity_prior_grid,
    settings$efficacy_prior_grid, settings$utility_grid, settings$arrival_grid, settings$alpha_grid))
  out$setting_id <- seq_len(nrow(out)); out
}

aide_tite_tag <- function(task) paste0(
  "P12TITE-a", if (task$allocation == "two_stage") "2s" else "1s",
  "-N", task$Nmax, "-s1", task$s1_Max, "-s2", task$N_s2,
  "-u", task$utility_type, "-l", aide_tite_fmt(task$lambda_T),
  "-tm", task$toxicity_model, "-em", task$efficacy_model,
  "-rp", aide_tite_fmt(task$toxicity_a), "x", aide_tite_fmt(task$toxicity_b),
  "-ep", aide_tite_fmt(task$efficacy_a), "x", aide_tite_fmt(task$efficacy_b),
  "-cp", aide_tite_fmt(task$carry_a), "x", aide_tite_fmt(task$carry_b),
  "-ap", aide_tite_fmt(task$additive_a), "x", aide_tite_fmt(task$additive_b),
  "-w", aide_tite_fmt(task$T_assess), "-c", task$cohort_size, "-cyc", task$cycle_max,
  "-rate", aide_tite_fmt(task$arrival_rate), "-neval", task$n_eval,
  "-ta", aide_tite_fmt(task$toxicity_ipde_alpha), "-ea", aide_tite_fmt(task$efficacy_ipde_alpha)
)

aide_tite_make_tasks <- function(settings, truth, job_i, ntrial, root) {
  grid <- aide_tite_setting_grid(settings); tasks <- vector("list", nrow(grid) * nrow(truth)); k <- 0L
  for (i in seq_len(nrow(grid))) for (j in seq_len(nrow(truth))) {
    k <- k + 1L; g <- as.list(grid[i, , drop = FALSE]); tr <- truth[j, , drop = FALSE]
    tasks[[k]] <- c(g, list(task_id = k, job_i = job_i, ntrial = ntrial, seed = settings$seed_base + (job_i - 1L) * ntrial,
      Scenario = as.integer(tr$Scenario), Source_Scenario = as.integer(tr$Source_Scenario), Scenario_Group = as.character(tr$Scenario_Group),
      Attempt = as.integer(tr$Attempt), true_mtd = as.integer(tr$True_MTD_Level), target = as.numeric(tr$Target_Toxicity),
      p_true = as.numeric(unlist(tr[paste0("Tox_Dose", 1:5)], use.names = FALSE)),
      e_true = as.numeric(unlist(tr[paste0("Eff_Dose", 1:5)], use.names = FALSE)),
      cohort_size = settings$cohort_size, cycle_max = settings$cycle_max, T_assess = settings$T_assess,
      n_eval = settings$n_eval, m_U = settings$m_U, skeleton = settings$skeleton, utility_scores = settings$utility_scores, root = root))
  }
  tasks
}

aide_tite_task_config <- function(task) aide_phase12_config(
  allocation = task$allocation, cohort_size = task$cohort_size, Nmax = task$Nmax, s1_Max = task$s1_Max,
  N_s2 = task$N_s2, n_eval = task$n_eval, m_U = task$m_U, cycle_max = task$cycle_max,
  time = list(arrival_rate = task$arrival_rate, t0 = 0, T_assess = task$T_assess),
  toxicity = list(model = task$toxicity_model, target = task$target, cutoff = .95, skeleton = task$skeleton,
    beta_prior_mean = 0, beta_prior_sd = sqrt(2), carryover_prior = c(task$toxicity_a, task$toxicity_b)),
  efficacy = list(model = task$efficacy_model, prior = c(task$efficacy_a, task$efficacy_b),
    carryover_prior = c(task$carry_a, task$carry_b), threshold = .20, futility_cutoff = .95, min_eff_n_for_futility = 0L),
  utility = list(type = task$utility_type, lambda_T = task$lambda_T, scores = task$utility_scores)
)
