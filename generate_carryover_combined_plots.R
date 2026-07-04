txt <- readLines("carryover_pk_r_estimators_demo.R")
end_idx <- grep("^## Example test run", txt)[1] - 1L
eval(parse(text = txt[seq_len(end_idx)]))

source("AIDE_CRM_helper_final.R")

## User settings
NR <- 9L
NI <- 9L
N_TRIALS <- 100L
SIM_SEED <- 20260704L
alpha_values <- c(0, 0.3, 0.6, 0.9)
include_crm_random <- TRUE

crm_random_settings <- list(
  target = 0.30,
  skeleton = NULL,
  a_r = 0.15,
  b_r = 0.85,
  alpha_sd = sqrt(2),
  model_file = "random_CRM.bug",
  n_chains = 2L,
  n_adapt = 500L,
  n_burnin = 500L,
  n_iter = 2000L,
  thin = 1L,
  seed = 20260704L
)

if (isTRUE(include_crm_random) &&
    (!requireNamespace("rjags", quietly = TRUE) ||
     !requireNamespace("coda", quietly = TRUE))) {
  warning(
    "Random CRM comparator requires rjags and coda; continuing without it.",
    call. = FALSE
  )
  include_crm_random <- FALSE
}

if (length(N_TRIALS) != 1L || !is.finite(N_TRIALS) || N_TRIALS < 1L) {
  stop("N_TRIALS must be a positive integer.")
}

N_TRIALS <- as.integer(N_TRIALS)

summary_source <- read.csv("carryover_pk_r_summary_results.csv")
summary_source <- summary_source[summary_source$scenario %in% paste0("sce", 1:5), ]

if (is.null(alpha_values)) {
  alpha_values <- sort(unique(summary_source$alpha_true))
}

make_aide_alpha_q <- function(p, alpha) {
  q <- p
  if (length(p) >= 2L) {
    q[-1] <- pmin(1, p[-1] + alpha * p[-length(p)])
  }
  q
}

fit_binomial_trial_summary <- function(scenario_name,
                                       p,
                                       alpha,
                                       n_regular,
                                       n_ipde,
                                       n_trials,
                                       seed) {
  all_trials <- vector("list", n_trials)
  
  for (trial_id in seq_len(n_trials)) {
    trial_seed <- seed + trial_id
    
    dat <- simulate_regular_ipde_counts(
      p_true = p,
      n_regular = n_regular,
      n_ipde = n_ipde,
      alpha_true = alpha,
      seed = trial_seed,
      ipde_model = "aide_alpha"
    )
    
    fit <- estimate_pk_with_global_r(
      dat = dat,
      r_prior = c(0.15, 0.85),
      p_prior = c(0.15, 0.85),
      r_max = 0.99,
      plug_in = "mean",
      debug_on_error = FALSE
    )
    
    est <- fit$estimates
    crm_p_hat <- rep(NA_real_, length(p))
    crm_r_hat <- NA_real_
    
    if (isTRUE(include_crm_random)) {
      crm_args <- crm_random_settings
      crm_args$skeleton <- if (is.null(crm_args$skeleton)) {
        p
      } else {
        crm_args$skeleton
      }
      if (!is.null(crm_args$seed)) {
        crm_args$seed <- crm_args$seed + trial_seed
      }
      
      crm_fit_random <- tryCatch(
        {
          do.call(
            fit_random_crm_comparator,
            c(list(dat = dat), crm_args)
          )
        },
        error = function(e) {
          warning(
            "Random CRM comparator failed for ",
            scenario_name,
            ", alpha=",
            alpha,
            ", trial=",
            trial_id,
            ": ",
            conditionMessage(e),
            call. = FALSE
          )
          NULL
        }
      )
      
      if (!is.null(crm_fit_random)) {
        crm_p_hat <- crm_fit_random$p_hat
        crm_r_hat <- if (length(crm_fit_random$r_hat) == 1L) {
          crm_fit_random$r_hat
        } else {
          mean(crm_fit_random$r_hat, na.rm = TRUE)
        }
      }
    }
    
    all_trials[[trial_id]] <- data.frame(
      scenario = scenario_name,
      alpha_true = alpha,
      replicate = trial_id,
      r_hat_used = fit$r_hat,
      r_hat_mean = fit$r$r_mean,
      r_hat_map = fit$r$r_map,
      r_hat_sd = fit$r$r_sd,
      dose = est$dose,
      nR = est$nR,
      yR = est$yR,
      nI = est$nI,
      yI = est$yI,
      p_true = est$p_true,
      q_true = est$q_true,
      obs_pooled = est$obs_pooled,
      p_pool_plugin = est$p_pool_plugin,
      p_mle_plugin = est$p_mle_plugin,
      p_crm_random = crm_p_hat,
      r_hat_crm_random = crm_r_hat,
      stringsAsFactors = FALSE
    )
  }
  
  trial_results <- do.call(rbind, all_trials)
  rownames(trial_results) <- NULL
  summarize_replicate_estimates(trial_results)
}

make_base_scenarios <- function(summary_rows) {
  scenario_names <- paste0("sce", 1:5)
  out <- setNames(vector("list", length(scenario_names)), scenario_names)
  
  for (scenario_name in scenario_names) {
    one <- summary_rows[summary_rows$scenario == scenario_name, , drop = FALSE]
    one <- one[order(one$dose), ]
    one <- one[!duplicated(one$dose), ]
    out[[scenario_name]] <- one$p_true
  }
  
  out
}

scenarios <- make_base_scenarios(summary_source)
scenarios$sce6 <- c(0.50, 0.55, 0.60, 0.65, 0.70)
scenario_names <- names(scenarios)

all_rows <- list()
idx <- 1L

for (scenario_name in scenario_names) {
  for (alpha in alpha_values) {
    cat(
      "Simulating ",
      N_TRIALS,
      " trials for ",
      scenario_name,
      ", alpha=",
      alpha,
      "\n",
      sep = ""
    )
    
    all_rows[[idx]] <- fit_binomial_trial_summary(
      scenario_name = scenario_name,
      p = scenarios[[scenario_name]],
      alpha = alpha,
      n_regular = NR,
      n_ipde = NI,
      n_trials = N_TRIALS,
      seed = SIM_SEED + idx * 100000L
    )
    idx <- idx + 1L
  }
}

combined_rows <- do.call(rbind, all_rows)
rownames(combined_rows) <- NULL

write.csv(
  combined_rows,
  "carryover_pk_r_summary_results_sce1_to_sce6_for_plots.csv",
  row.names = FALSE
)

files <- plot_carryover_average_alpha_files(
  summary_rows = combined_rows,
  alpha_values = alpha_values,
  scenario_names = scenario_names,
  plot_dir = "carryover_pk_r_demo_plots",
  file_suffix = paste0("_NR", NR, "_NI", NI)
)

cat(paste(files, collapse = "\n"), "\n")
