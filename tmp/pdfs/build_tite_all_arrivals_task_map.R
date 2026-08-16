setwd("TITE-AIDE-Rebuild")
source("TITE-AIDE.R")

# Recreate the task ordering in the full three-accrual-rate TITE result file.
settings <- aide_tite_oc_settings()
settings$arrival_grid <- data.frame(arrival_rate = c(1 / 56, 1 / 28, 1 / 14))
truth <- aide_tite_read_truth(settings$scenario_file, settings$scenario_ids)
tasks <- aide_tite_make_tasks(settings, truth, job_i = 1L, ntrial = 1L, root = getwd())

task_map <- do.call(rbind, lapply(tasks, function(task) {
  data.frame(
    Task_ID = task$task_id,
    Scenario = task$Scenario,
    Allocation = task$allocation,
    Model_ID = task$model_id,
    Utility_Type = task$utility_type,
    Lambda_T = task$lambda_T,
    Toxicity_IPDE_Alpha = task$toxicity_ipde_alpha,
    Efficacy_IPDE_Alpha = task$efficacy_ipde_alpha,
    Arrival_Rate = task$arrival_rate
  )
}))

out <- "../tmp/pdfs/tite_all_arrivals_task_map.csv"
utils::write.csv(task_map, out, row.names = FALSE)
cat("Wrote", nrow(task_map), "task-map rows to", out, "\n")
