## Run one non-TITE AIDE trial and one TITE AIDE-CRM trial for debugging.
## From this folder:
##   Rscript debug_one_trial_AIDE_TITE.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(file_arg[1L]))
} else {
  getwd()
}
setwd(script_dir)

source("AIDE_BOIN_helper.R")
source("AIDE_modified.R")
source("AIDE_CRM_helper_TITE.R")

## ---- Debug settings: edit these as needed ----
p_true <- c(0.10,0.20,0.30,0.40,0.50)
crm_skeleton <- c(0.05, 0.10, 0.20, 0.35, 0.50)
target <- 0.30
ipde_alpha <- 0.50
seed <- 20260716L

## More patients/cycles make recycled patients easier to observe.
n_patients <- 30L
cohort_size <- 3L
cycle_max <- 2L

show_recycled_patients <- function(admin, p_true, ipde_alpha) {
  by_id <- split(admin, admin$id)
  out <- lapply(by_id, function(patient_history) {
    patient_history <- patient_history[
      order(patient_history$ncycle, patient_history$row_id),
      , drop = FALSE
    ]
    retreat_rows <- which(patient_history$type == "retreat")
    if (length(retreat_rows) == 0L) return(NULL)

    do.call(rbind, lapply(retreat_rows, function(i) {
      previous <- patient_history[i - 1L, , drop = FALSE]
      current <- patient_history[i, , drop = FALSE]
      data.frame(
        id = current$id,
        cycle = current$ncycle,
        previous_dose = previous$dose,
        current_dose = current$dose,
        ipde_probability = aide_ipde_toxicity_probability(
          p_true = p_true,
          previous_dose = previous$dose,
          current_dose = current$dose,
          ipde_alpha = ipde_alpha
        ),
        observed_dlt = current$y,
        stringsAsFactors = FALSE
      )
    }))
  })

  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) {
    message("No recycled patients occurred in this run. Change seed or increase n_patients.")
    return(invisible(NULL))
  }
  do.call(rbind, out)
}

cat("\n========== One non-TITE AIDE trial ==========\n")
aide_fit <- simulate_AIDE_design(
  p_true = p_true,
  ipde_alpha = ipde_alpha,
  seed = seed,
  verbose = TRUE,
  model = "CRM",
  ## 1 permits recycling from any lower dose; use 2 for adjacent-dose only.
  ipde_design = 1L,
  arrival_rate = 1/56,
  N_pat = n_patients,
  Nmax_eff = n_patients,
  crm_skeleton = c(0.1,0.2,0.3,0.4,0.5),
  crm_r_model = 'r_fixed',
  C = cohort_size,
  cycle_max = cycle_max,
  T_assess = 28,
  new_pat_first = 1L,
  n_eval_escalate = 3L,
  dose_cap = 3,
  TARGET = target
)

print(aide_fit$admin)
print(show_recycled_patients(aide_fit$admin, p_true, ipde_alpha))
print(aide_fit$final)

# cat("\n========== One TITE AIDE-CRM trial ==========\n")
# tite_fit <- simulate_AIDE_CRM_TITE_design(
#   p_true = p_true,
#   ipde_alpha = ipde_alpha,
#   target = target,
#   seed = seed,
#   verbose = TRUE,
#   N_pat = n_patients,
#   Nmax_eff = n_patients,
#   cohortsize = cohort_size,
#   cycle_max = cycle_max,
#   T_assess = 28,
#   arrival_rate = 1/56,
#   new_pat_first = 2L,
#   n_eval_escalate = 3L,
#   dose_cap = 3,
#   ## alpha_crm is deterministic and does not require JAGS/rjags.
#   ## Change to "fixed" or "random" when debugging a JAGS CRM backend.
#   crm_r_model = "r_fixed",
#   r_carry = 0,
#   crm_skeleton = crm_skeleton,
#   ## A smaller alpha grid keeps this debugging run quick.
#   crm_alpha_grid = seq(0.05, 0.95, length.out = 11L)
# )
# 
# print(tite_fit$admin)
# print(tite_fit$decision_log)
# print(show_recycled_patients(tite_fit$admin, p_true, ipde_alpha))
# print(tite_fit$final)
