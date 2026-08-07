## ============================================================
## run_oc_TITE_AIDE_phase_I_II.R
## JAGS TITE-AIDE Phase I/II cluster-style OC runner.
## Mirrors run_oc_AIDE_phase_I_II.R: one IDX, one deterministic seed block,
## all 37 scenarios and the same design/model/prior/utility/DGM grids.
## Usage: Rscript run_oc_TITE_AIDE_phase_I_II.R [workers] [ntrial] [IDX]
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
if (length(script_file) == 1L && file.exists(script_file)) setwd(dirname(normalizePath(script_file)))
source("TITE-AIDE.R")
if (!requireNamespace("rjags", quietly = TRUE) || !requireNamespace("coda", quietly = TRUE))
  stop("This runner requires rjags and coda in the active R library.")

positive_integer <- function(x, name) {
  x <- suppressWarnings(as.integer(x)); if (length(x) != 1L || is.na(x) || x < 1L) stop(name, " must be a positive integer.")
  x
}
job_index <- function(args) {
  if (length(args) >= 3L && !is.na(suppressWarnings(as.integer(args[3L])))) return(positive_integer(args[3L], "IDX"))
  for (name in c("JOB_I", "LSB_JOBINDEX")) if (nzchar(Sys.getenv(name, ""))) return(positive_integer(Sys.getenv(name), name))
  1L
}

settings <- aide_tite_oc_settings()
workers <- if (length(args) >= 1L) positive_integer(args[1L], "workers") else 1L
ntrial <- if (length(args) >= 2L) positive_integer(args[2L], "ntrial") else settings$ntrial
job_i <- job_index(args)
if (workers != 1L) warning("workers is retained for run_oc_AIDE_phase_I_II compatibility; execution is sequential.")
truth <- aide_tite_read_truth(settings$scenario_file, settings$scenario_ids)
scenario_set <- tools::file_path_sans_ext(basename(settings$scenario_file))
results_root <- paste0("oc_results_cluster_TITE_AIDE_phase_I_II_", scenario_set)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
tasks <- aide_tite_make_tasks(settings, truth, job_i, ntrial, getwd())

cat("LSF job/block index:", job_i, "\nTrials per setting:", ntrial, "\nNumber of tasks:", length(tasks), "\n")

run_task <- function(task) {
  tag <- aide_tite_tag(task)
  folder <- file.path(results_root, tag); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0("P12TITE-SC", task$Scenario, "-", task$model_id, "-j", task$job_i,
                     "-b1-s", task$seed, "-n", task$ntrial, ".rds")
  outfile <- file.path(folder, filename)
  cat("Task", task$task_id, "/", length(tasks), "scenario", task$Scenario, "\n")
  config <- aide_tite_task_config(task)
  scenario <- aide_phase12_scenario(task$p_true, task$e_true,
    toxicity_ipde_dgm = list(alpha_true = task$toxicity_ipde_alpha),
    efficacy_ipde_dgm = list(alpha_true = task$efficacy_ipde_alpha))
  result <- get_oc_sim_AIDE_phase_I_II(scenario, config, ntrial = task$ntrial, seed = task$seed, store_raw_tables = FALSE)
  result$task <- task
  result$truth <- list(p_true = task$p_true, e_true = task$e_true, target = task$target, true_mtd = task$true_mtd)
  result$runner_settings <- list(n_eval_rule = "stay_and_escalate_blocked", new_first = TRUE,
    seed = task$seed, ntrial = task$ntrial, T_assess = task$T_assess, cycle_max = task$cycle_max)
  saveRDS(result, outfile)
  combined <- file.path(folder, paste0("P12TITE-SC", task$Scenario, "_", tag, "-job-", task$job_i, "-combined.rds"))
  saveRDS(result, combined)
  cat("Saved:", outfile, "\n")
  outfile
}

files <- vapply(tasks, run_task, character(1))
cat("All TITE-AIDE Phase I/II OC tasks completed.", "\n")
