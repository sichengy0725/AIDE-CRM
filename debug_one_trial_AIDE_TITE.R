## Run one or more non-TITE AIDE-CRM and TITE AIDE-CRM trials for debugging.
## From this folder:
##   Rscript debug_one_trial_AIDE_TITE.R
##   Rscript debug_one_trial_AIDE_TITE.R 20
##   Rscript debug_one_trial_AIDE_TITE.R --ntrial=20 --seed=20260716
##
## Useful flags:
##   --print-trials=TRUE    print admin/decision/final output for every trial
##   --verbose=TRUE         pass verbose=TRUE into the simulator
##   --keep-browser=TRUE    keep browser() calls in sourced helper code

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

cli_args <- commandArgs(trailingOnly = TRUE)
pos_args <- cli_args[!grepl("^--", cli_args)]

get_named_arg <- function(name, default = NULL) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, cli_args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(pat, "", hit[1L])
}

as_bool <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

## ---- Debug settings: edit these as needed ----
p_true <- c(0.3, 0.35, 0.40, 0.45, 0.50)
crm_skeleton <- c(0.05, 0.10, 0.20, 0.35, 0.50)
target <- 0.30
ipde_alpha <- 0.50
seed_default <- 20260716L

## More patients/cycles make recycled patients easier to observe.
n_patients <- 30L
cohort_size <- 3L
cycle_max <- 1L
t_assess <- 28
arrival_rate <- 1 / 56
new_pat_first <- 1L
ipde_design <- 1L
n_eval_escalate <- 3L
dose_cap <- 3L
crm_r_model <- "r_fixed"
r_carry <- 0
restrict_to_tried <- TRUE
restrict_to_target <- FALSE
crm_alpha_grid <- seq(0.05, 0.95, length.out = 11L)

n_trials <- as.integer(get_named_arg(
  "ntrial",
  if (length(pos_args) >= 1L) pos_args[1L] else 10L
))
seed <- as.integer(get_named_arg(
  "seed",
  if (length(pos_args) >= 2L) pos_args[2L] else seed_default
))

if (is.na(n_trials) || n_trials < 1L) {
  stop("ntrial must be a positive integer.")
}
if (is.na(seed)) {
  stop("seed must be an integer.")
}

print_trials <- as_bool(get_named_arg("print-trials", n_trials <= 5L))
verbose_trials <- as_bool(get_named_arg("verbose", n_trials == 1L))
keep_browser <- as_bool(get_named_arg("keep-browser", n_trials == 1L))

if (!keep_browser) {
  ## The helper currently contains a browser() before final TITE CRM selection.
  ## In multi-trial command-line runs that would stop the batch, so mask it here.
  browser <- function(...) invisible(NULL)
}

show_recycled_patients <- function(admin, p_true, ipde_alpha, quiet = FALSE) {
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
    if (!quiet) {
      message("No recycled patients occurred in this run. Change seed or increase n_patients.")
    }
    return(invisible(NULL))
  }
  do.call(rbind, out)
}

summarize_debug_trial <- function(fit, itrial, trial_seed, design) {
  if (!is.null(fit$error)) {
    return(data.frame(
      design = design,
      trial = itrial,
      seed = trial_seed,
      error = fit$error,
      MTD = NA_integer_,
      earlystop = NA_integer_,
      total_admin = NA_integer_,
      total_unique = NA_integer_,
      ipde_admin = NA_integer_,
      dlt = NA_integer_,
      max_dose = NA_integer_,
      trial_time = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  admin <- fit$admin
  data.frame(
    design = design,
    trial = itrial,
    seed = trial_seed,
    error = NA_character_,
    MTD = as.integer(fit$final$MTD),
    earlystop = as.integer(isTRUE(as.logical(fit$final$earlystop))),
    total_admin = nrow(admin),
    total_unique = if (nrow(admin) > 0L) length(unique(admin$id)) else 0L,
    ipde_admin = if (nrow(admin) > 0L) sum(admin$type == "retreat") else 0L,
    dlt = if (nrow(admin) > 0L) sum(admin$y == 1L) else 0L,
    max_dose = if (nrow(admin) > 0L) max(admin$dose, na.rm = TRUE) else NA_integer_,
    trial_time = fit$final$trial_time,
    stringsAsFactors = FALSE
  )
}

patients_by_dose <- function(fit, ndose) {
  if (!is.null(fit$error) || nrow(fit$admin) == 0L) {
    return(rep(NA_integer_, ndose))
  }
  tabulate(fit$admin$dose, nbins = ndose)
}

print_batch_summary <- function(label, fits, summary_df, p_true) {
  cat("\n==========", label, "trial summary ==========\n")
  print(summary_df)

  ok <- is.na(summary_df$error)
  if (any(!ok)) {
    cat("\n==========", label, "failed trials ==========\n")
    print(summary_df[!ok, c("design", "trial", "seed", "error"), drop = FALSE])
  }

  if (!any(ok)) {
    return(invisible(NULL))
  }

  ndose <- length(p_true)
  dose_mat <- do.call(rbind, lapply(fits[ok], patients_by_dose, ndose = ndose))
  colnames(dose_mat) <- paste0("D", seq_len(ndose))

  cat("\n==========", label, "MTD selection among completed trials ==========\n")
  mtd_levels <- c(as.character(seq_len(ndose)), "99")
  mtd_tab <- table(factor(as.character(summary_df$MTD[ok]), levels = mtd_levels))
  print(data.frame(
    MTD = names(mtd_tab),
    count = as.integer(mtd_tab),
    percent = 100 * as.integer(mtd_tab) / sum(ok),
    row.names = NULL
  ))

  cat("\n==========", label, "average administered patients by dose ==========\n")
  print(colMeans(dose_mat, na.rm = TRUE))

  cat("\n==========", label, "overall averages ==========\n")
  print(data.frame(
    total_admin = mean(summary_df$total_admin[ok], na.rm = TRUE),
    total_unique = mean(summary_df$total_unique[ok], na.rm = TRUE),
    ipde_admin = mean(summary_df$ipde_admin[ok], na.rm = TRUE),
    dlt = mean(summary_df$dlt[ok], na.rm = TRUE),
    max_dose = mean(summary_df$max_dose[ok], na.rm = TRUE),
    trial_time = mean(summary_df$trial_time[ok], na.rm = TRUE)
  ))

  invisible(NULL)
}

cat("\n========== Non-TITE AIDE-CRM debug run ==========\n")
cat("n_trials:", n_trials, "\n")
cat("base seed:", seed, "\n")
cat("crm_r_model:", crm_r_model, "\n")
cat("ipde_design:", ipde_design, "\n")
cat("new_pat_first:", new_pat_first, "\n")
cat("restrict_to_tried:", restrict_to_tried, "\n")
cat("restrict_to_target:", restrict_to_target, "\n")

aide_fits <- vector("list", n_trials)

for (itrial in seq_len(n_trials)) {
  trial_seed <- seed + itrial - 1L
  cat("\n----- AIDE trial", itrial, "seed", trial_seed, "-----\n")

  aide_fits[[itrial]] <- tryCatch(
    simulate_AIDE_design(
      p_true = p_true,
      ipde_alpha = ipde_alpha,
      seed = trial_seed,
      verbose = verbose_trials,
      model = "CRM",
      ipde_design = ipde_design,
      arrival_rate = arrival_rate,
      N_pat = n_patients,
      Nmax_eff = n_patients,
      crm_skeleton = crm_skeleton,
      crm_r_model = crm_r_model,
      r_carry = r_carry,
      crm_alpha_grid = crm_alpha_grid,
      C = cohort_size,
      cycle_max = cycle_max,
      T_assess = t_assess,
      new_pat_first = new_pat_first,
      n_eval_escalate = n_eval_escalate,
      dose_cap = dose_cap,
      restrict_to_tried = restrict_to_tried,
      restrict_to_target = restrict_to_target,
      TARGET = target
    ),
    error = function(e) {
      list(error = conditionMessage(e), seed = trial_seed)
    }
  )

  if (print_trials) {
    if (!is.null(aide_fits[[itrial]]$error)) {
      cat("ERROR:", aide_fits[[itrial]]$error, "\n")
    } else {
      print(aide_fits[[itrial]]$admin)
      print(show_recycled_patients(
        aide_fits[[itrial]]$admin,
        p_true,
        ipde_alpha
      ))
      print(aide_fits[[itrial]]$final)
    }
  }
}

cat("\n========== TITE AIDE-CRM debug run ==========\n")
cat("n_trials:", n_trials, "\n")
cat("base seed:", seed, "\n")
cat("crm_r_model:", crm_r_model, "\n")
cat("ipde_design:", ipde_design, "\n")
cat("new_pat_first:", new_pat_first, "\n")
cat("restrict_to_tried:", restrict_to_tried, "\n")
cat("restrict_to_target:", restrict_to_target, "\n")
cat("keep_browser:", keep_browser, "\n")

tite_fits <- vector("list", n_trials)

for (itrial in seq_len(n_trials)) {
  trial_seed <- seed + itrial - 1L
  cat("\n----- TITE trial", itrial, "seed", trial_seed, "-----\n")

  tite_fits[[itrial]] <- tryCatch(
    simulate_AIDE_CRM_TITE_design(
      p_true = p_true,
      ipde_alpha = ipde_alpha,
      target = target,
      seed = trial_seed,
      verbose = verbose_trials,
      N_pat = n_patients,
      Nmax_eff = n_patients,
      cohortsize = cohort_size,
      cycle_max = cycle_max,
      T_assess = t_assess,
      arrival_rate = arrival_rate,
      new_pat_first = new_pat_first,
      ipde_design = ipde_design,
      n_eval_escalate = n_eval_escalate,
      dose_cap = dose_cap,
      restrict_to_tried = restrict_to_tried,
      restrict_to_target = restrict_to_target,
      ## alpha_crm is deterministic and does not require JAGS/rjags.
      ## Change to "fixed" or "random" when debugging a JAGS CRM backend.
      crm_r_model = crm_r_model,
      r_carry = r_carry,
      crm_skeleton = crm_skeleton,
      ## A smaller alpha grid keeps alpha-CRM debugging runs quick.
      crm_alpha_grid = crm_alpha_grid
    ),
    error = function(e) {
      list(error = conditionMessage(e), seed = trial_seed)
    }
  )

  if (print_trials) {
    if (!is.null(tite_fits[[itrial]]$error)) {
      cat("ERROR:", tite_fits[[itrial]]$error, "\n")
    } else {
      print(tite_fits[[itrial]]$admin)
      print(tite_fits[[itrial]]$decision_log)
      print(show_recycled_patients(
        tite_fits[[itrial]]$admin,
        p_true,
        ipde_alpha
      ))
      print(tite_fits[[itrial]]$final)
    }
  }
}

aide_summary <- do.call(
  rbind,
  lapply(seq_len(n_trials), function(i) {
    summarize_debug_trial(aide_fits[[i]], i, seed + i - 1L, "AIDE")
  })
)

tite_summary <- do.call(
  rbind,
  lapply(seq_len(n_trials), function(i) {
    summarize_debug_trial(tite_fits[[i]], i, seed + i - 1L, "TITE")
  })
)

print_batch_summary("AIDE", aide_fits, aide_summary, p_true)
print_batch_summary("TITE", tite_fits, tite_summary, p_true)

cat("\n========== Combined AIDE vs TITE trial summary ==========\n")
print(rbind(aide_summary, tite_summary))

## Convenience alias for old one-trial debugging snippets.
aide_fit <- aide_fits[[1L]]
tite_fit <- tite_fits[[1L]]
