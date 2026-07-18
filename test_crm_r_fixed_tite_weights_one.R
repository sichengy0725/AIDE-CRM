## Regression check: non-TITE CRM versus TITE CRM with complete follow-up.
##
## Both fits use fixed-r CRM.  Because every TITE weight is one, they have
## exactly the same likelihood and should yield the same fitted dose-toxicity
## curve (apart from negligible MCMC variation).
##
## Run from this directory with:
##   & 'C:/Program Files/R/R-4.3.2/bin/Rscript.exe' test_crm_r_fixed_tite_weights_one.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1L])) else getwd()
setwd(script_dir)

source("AIDE_CRM_helper_modified.R")

if (!requireNamespace("rjags", quietly = TRUE) ||
    !requireNamespace("coda", quietly = TRUE)) {
  stop("This test requires the rjags and coda packages (and a working JAGS installation).")
}

## Patient-level data.  Edit this block to use another completed dataset.
## type = "retreat" identifies a recycled/IPDE administration; all other
## rows are treated as ordinary/new administrations.
patient_data <- data.frame(
  id    = 1:15,
  dose  = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5),
  y     = c(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1),
  type  = "new",
  cycle = 1L,
  stringsAsFactors = FALSE
)

## Design/model inputs.
ndose <- 5L
skeleton <- c(0.05, 0.10, 0.20, 0.35, 0.50)
target <- 0.30
r_carry <- 0.10
mcmc_seed <- 20260718L

## Standard (non-TITE) CRM: no TITE-weight column is supplied.
non_tite_dat <- patient_data

## TITE CRM: every patient has complete follow-up, so each weight is one.
tite_dat <- transform(patient_data, tite_weight = rep(1, nrow(patient_data)))

## First confirm that both data sets reduce to the identical likelihood input.
prepared_non_tite <- crm_prepare_dat(non_tite_dat, ndose = ndose)
prepared_tite <- crm_prepare_dat(
  tite_dat,
  ndose = ndose,
  weight_col = "tite_weight"
)
stopifnot(
  all(prepared_non_tite$tite_weight == 1),
  all(prepared_tite$tite_weight == 1),
  identical(prepared_non_tite, prepared_tite)
)

fit_args <- list(
  ndose = ndose,
  skeleton = skeleton,
  target = target,
  r_model = "r_fixed",
  r_carry = r_carry,
  n_chains = 2L,
  n_adapt = 500L,
  n_burnin = 1000L,
  n_iter = 10000L,
  thin = 1L,
  seed = mcmc_seed
)

non_tite_fit <- do.call(crm_fit, c(list(dat = non_tite_dat), fit_args))
tite_fit <- do.call(
  crm_fit,
  c(list(dat = tite_dat, weight_col = "tite_weight"), fit_args)
)

comparison <- data.frame(
  dose = seq_len(ndose),
  non_tite_p_hat = non_tite_fit$p_hat,
  tite_p_hat = tite_fit$p_hat,
  absolute_difference = abs(non_tite_fit$p_hat - tite_fit$p_hat)
)

print(patient_data)
cat("\nFitted toxicity probabilities\n")
print(round(comparison, 6), row.names = FALSE)
cat("\nFixed r:", non_tite_fit$r_hat, "\n")
cat("Effective sample sizes (non-TITE, TITE):",
    non_tite_fit$n_eff, tite_fit$n_eff, "\n")

## Repeating the fit with the same seed makes the two posterior draws
## deterministic in the usual rjags setup.  The tolerance also permits tiny
## platform-specific MCMC differences while still catching a weighting bug.
tolerance <- 0.01
max_difference <- max(comparison$absolute_difference)
if (max_difference > tolerance) {
  stop(
    sprintf(
      "Complete-follow-up TITE CRM disagrees with non-TITE CRM: max |difference| = %.6f (tolerance %.3f).",
      max_difference, tolerance
    )
  )
}

cat(sprintf("\nPASS: max fitted-probability difference = %.6f (tolerance %.3f).\n",
            max_difference, tolerance))
