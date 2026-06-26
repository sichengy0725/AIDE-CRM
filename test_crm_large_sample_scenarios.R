## ============================================================
## Large-sample estimation check for fixed-r AIDE-CRM
## Data are generated dose-by-dose from the user-specified true
## scenario probabilities.
##
## Required files in the same folder:
##   - AIDE_CRM_helper.R
##   - fix_CRM.bug
## ============================================================

rm(list = ls())

## ---- paths --------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else {
  getwd()
}

helper_file <- file.path(script_dir, "AIDE_CRM_helper.R")
fixed_model_file <- file.path(script_dir, "fix_CRM.bug")

if (!file.exists(helper_file)) stop("Cannot find: ", helper_file)
if (!file.exists(fixed_model_file)) stop("Cannot find: ", fixed_model_file)

source(helper_file)

## ---- user scenarios -----------------------------------------
scenarios <- rbind(
  c(0.07, 0.12, 0.17, 0.22, 0.30),
  c(0.05, 0.10, 0.18, 0.30, 0.40),
  c(0.15, 0.20, 0.30, 0.35, 0.45),
  c(0.15, 0.30, 0.38, 0.45, 0.55),
  c(0.30, 0.35, 0.40, 0.45, 0.50),
  c(0.50, 0.55, 0.60, 0.65, 0.70)
)

TARGET <- 0.30
K <- ncol(scenarios)
scenario_id <- seq_len(nrow(scenarios))

## Fixed carryover/discount parameter used in fix_CRM.bug.
r_true <- 0

## Large sample per dose.
## Increase these for a stricter asymptotic check.
n_new_per_dose <- 100L
n_ipde_per_dose <- 100L

## alpha_sd is a standard deviation, not variance.
## AIDE_CRM_helper.R internally sends tau_alpha = 1 / alpha_sd^2.
alpha_sd <- 2

## JAGS/MCMC settings. Increase n_iter for final checks.
n_chains <- 3
n_adapt <- 1000
n_burnin <- 1000
n_iter <- 5000
thin <- 2

## ------------------------------------------------------------
## Skeleton choice
## ------------------------------------------------------------
## Use skeleton_mode = "scenario" if the goal is to verify that the
## fixed CRM model can recover the true probabilities in large samples.
## In this mode, each scenario uses skeleton = p_true, so the true alpha is 0.
##
## Use skeleton_mode = "fixed" if the goal is to evaluate your actual CRM
## with one prespecified skeleton under the six true scenarios. In that case,
## because the CRM has only one alpha parameter, p_hat generally converges to
## the best-fitting CRM curve, not exactly to every p_true vector.

skeleton_mode <- "scenario"   # change to "fixed" if desired
fixed_skeleton <- scenarios[1, ]

## ---- helper functions ---------------------------------------

sim_crm_scenario_data <- function(p_true,
                                  r_true,
                                  n_new_per_dose,
                                  n_ipde_per_dose = 0L,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  K <- length(p_true)
  theta_ipde_true <- r_true + (1 - r_true) * p_true

  out <- list()

  if (n_new_per_dose > 0L) {
    dose_new <- rep(seq_len(K), each = n_new_per_dose)
    y_new <- stats::rbinom(
      n = length(dose_new),
      size = 1L,
      prob = p_true[dose_new]
    )
    out[[length(out) + 1L]] <- data.frame(
      dose = as.integer(dose_new),
      y = as.integer(y_new),
      type = "new"
    )
  }

  if (n_ipde_per_dose > 0L) {
    dose_ipde <- rep(seq_len(K), each = n_ipde_per_dose)
    y_ipde <- stats::rbinom(
      n = length(dose_ipde),
      size = 1L,
      prob = theta_ipde_true[dose_ipde]
    )
    out[[length(out) + 1L]] <- data.frame(
      dose = as.integer(dose_ipde),
      y = as.integer(y_ipde),
      type = "retreat"
    )
  }

  dat <- do.call(rbind, out)
  rownames(dat) <- NULL
  dat
}

safe_prob <- function(x, eps = 1e-12) {
  pmin(pmax(x, eps), 1 - eps)
}

crm_curve <- function(alpha, skeleton) {
  safe_prob(skeleton ^ exp(alpha))
}

## Large-sample pseudo-true CRM curve for a fixed skeleton.
## This is useful when skeleton != p_true. It computes the alpha that maximizes
## the expected Bernoulli log-likelihood under the true scenario.
crm_pseudotruth <- function(p_true,
                            skeleton,
                            r_true,
                            n_new_per_dose,
                            n_ipde_per_dose) {
  theta_ipde_true <- r_true + (1 - r_true) * p_true

  neg_expected_loglik <- function(alpha) {
    p <- crm_curve(alpha, skeleton)
    theta_ipde <- r_true + (1 - r_true) * p

    ll_new <- 0
    if (n_new_per_dose > 0L) {
      ll_new <- sum(
        n_new_per_dose *
          (p_true * log(safe_prob(p)) +
             (1 - p_true) * log(safe_prob(1 - p)))
      )
    }

    ll_ipde <- 0
    if (n_ipde_per_dose > 0L) {
      ll_ipde <- sum(
        n_ipde_per_dose *
          (theta_ipde_true * log(safe_prob(theta_ipde)) +
             (1 - theta_ipde_true) * log(safe_prob(1 - theta_ipde)))
      )
    }

    -(ll_new + ll_ipde)
  }

  opt <- stats::optimize(
    f = neg_expected_loglik,
    interval = c(-10, 10),
    tol = 1e-10
  )

  alpha_star <- opt$minimum
  p_star <- crm_curve(alpha_star, skeleton)
  theta_star <- r_true + (1 - r_true) * p_star

  list(
    alpha_star = alpha_star,
    p_star = p_star,
    theta_ipde_star = theta_star
  )
}

fit_one_scenario <- function(s, seed_base = 20260611L) {
  p_true <- as.numeric(scenarios[s, ])
  theta_ipde_true <- r_true + (1 - r_true) * p_true
  true_mtd <- which.min(abs(p_true - TARGET))

  skeleton <- switch(
    skeleton_mode,
    scenario = p_true,
    fixed = fixed_skeleton,
    stop("Unknown skeleton_mode: ", skeleton_mode)
  )

  dat <- sim_crm_scenario_data(
    p_true = p_true,
    r_true = r_true,
    n_new_per_dose = n_new_per_dose,
    n_ipde_per_dose = n_ipde_per_dose,
    seed = seed_base + s
  )

  obs_new <- with(subset(dat, type == "new"), tapply(y, dose, mean))
  obs_ipde <- if (n_ipde_per_dose > 0L) {
    with(subset(dat, type == "retreat"), tapply(y, dose, mean))
  } else {
    rep(NA_real_, K)
  }

  pseudo <- crm_pseudotruth(
    p_true = p_true,
    skeleton = skeleton,
    r_true = r_true,
    n_new_per_dose = n_new_per_dose,
    n_ipde_per_dose = n_ipde_per_dose
  )

  fit <- crm_fit_discount(
    dat = dat,
    ndose = K,
    skeleton = skeleton,
    target = TARGET,
    r_model = "fixed",
    r_carry = r_true,
    alpha_sd = alpha_sd,
    fixed_model_file = fixed_model_file,
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_burnin = n_burnin,
    n_iter = n_iter,
    thin = thin,
    seed = seed_base + 1000L + s
  )

  est_mtd <- which.min(abs(fit$p_hat - TARGET))
  pseudo_mtd <- which.min(abs(pseudo$p_star - TARGET))

  dose_tab <- data.frame(
    scenario = s,
    dose = seq_len(K),
    skeleton = skeleton,
    p_true = p_true,
    obs_new_rate = as.numeric(obs_new),
    p_pseudotrue = pseudo$p_star,
    p_hat = fit$p_hat,
    p_hat_minus_true = fit$p_hat - p_true,
    p_hat_minus_pseudotrue = fit$p_hat - pseudo$p_star,
    theta_ipde_true = theta_ipde_true,
    obs_ipde_rate = as.numeric(obs_ipde),
    theta_ipde_pseudotrue = pseudo$theta_ipde_star,
    theta_ipde_hat = fit$theta_ipde_hat,
    theta_ipde_hat_minus_true = fit$theta_ipde_hat - theta_ipde_true,
    theta_ipde_hat_minus_pseudotrue = fit$theta_ipde_hat - pseudo$theta_ipde_star
  )

  scenario_tab <- data.frame(
    scenario = s,
    skeleton_mode = skeleton_mode,
    total_n = nrow(dat),
    n_new_per_dose = n_new_per_dose,
    n_ipde_per_dose = n_ipde_per_dose,
    r_true = r_true,
    r_hat = fit$r_hat,
    true_mtd = true_mtd,
    pseudo_mtd = pseudo_mtd,
    est_mtd = est_mtd,
    alpha_pseudotrue = pseudo$alpha_star,
    max_abs_p_hat_minus_true = max(abs(fit$p_hat - p_true)),
    max_abs_p_hat_minus_pseudotrue = max(abs(fit$p_hat - pseudo$p_star)),
    max_abs_theta_hat_minus_true = max(abs(fit$theta_ipde_hat - theta_ipde_true)),
    max_abs_theta_hat_minus_pseudotrue = max(abs(fit$theta_ipde_hat - pseudo$theta_ipde_star))
  )

  list(
    scenario_summary = scenario_tab,
    dose_summary = dose_tab,
    fit = fit
  )
}

## ---- run all scenarios --------------------------------------
results <- vector("list", length(scenario_id))

for (s in scenario_id) {
  cat("\n============================================================\n")
  cat("Fitting scenario", s, "\n")
  cat("True probabilities:", paste(scenarios[s, ], collapse = ", "), "\n")
  cat("Skeleton mode:", skeleton_mode, "\n")
  cat("============================================================\n")

  results[[s]] <- fit_one_scenario(s)

  print(round(results[[s]]$dose_summary[, c(
    "dose",
    "skeleton",
    "p_true",
    "obs_new_rate",
    "p_pseudotrue",
    "p_hat",
    "p_hat_minus_true",
    "p_hat_minus_pseudotrue",
    "theta_ipde_true",
    "obs_ipde_rate",
    "theta_ipde_hat"
  )], 4))

  print(results[[s]]$scenario_summary)
}

scenario_summary <- do.call(rbind, lapply(results, `[[`, "scenario_summary"))
dose_summary <- do.call(rbind, lapply(results, `[[`, "dose_summary"))

cat("\n===== Overall scenario summary =====\n")
print(scenario_summary)

## ---- write output -------------------------------------------
out_prefix <- paste0(
  "crm_large_sample_",
  skeleton_mode,
  "_r", gsub("\\.", "p", as.character(r_true)),
  "_new", n_new_per_dose,
  "_ipde", n_ipde_per_dose
)

write.csv(
  scenario_summary,
  file = file.path(script_dir, paste0(out_prefix, "_scenario_summary.csv")),
  row.names = FALSE
)

write.csv(
  dose_summary,
  file = file.path(script_dir, paste0(out_prefix, "_dose_summary.csv")),
  row.names = FALSE
)

saveRDS(
  list(
    settings = list(
      target = TARGET,
      r_true = r_true,
      alpha_sd = alpha_sd,
      skeleton_mode = skeleton_mode,
      fixed_skeleton = fixed_skeleton,
      n_new_per_dose = n_new_per_dose,
      n_ipde_per_dose = n_ipde_per_dose,
      n_chains = n_chains,
      n_adapt = n_adapt,
      n_burnin = n_burnin,
      n_iter = n_iter,
      thin = thin
    ),
    scenarios = scenarios,
    scenario_summary = scenario_summary,
    dose_summary = dose_summary,
    results = results
  ),
  file = file.path(script_dir, paste0(out_prefix, ".rds"))
)

cat("\nSaved files with prefix:", out_prefix, "\n")

## ---- interpretation note ------------------------------------
cat("\nInterpretation note:\n")
cat("- If skeleton_mode = 'scenario', the CRM working model is correctly specified,\n")
cat("  because skeleton = p_true and alpha_true = 0. In large samples, p_hat\n")
cat("  should be close to p_true.\n")
cat("- If skeleton_mode = 'fixed', the model is generally misspecified for these\n")
cat("  arbitrary scenario vectors. Then p_hat should be compared primarily with\n")
cat("  p_pseudotrue, the best-fitting CRM curve, not necessarily p_true.\n")
