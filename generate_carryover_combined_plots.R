txt <- readLines("carryover_pk_r_estimators_demo.R")
end_idx <- grep("^## Example test run", txt)[1] - 1L
eval(parse(text = txt[seq_len(end_idx)]))

summary_rows <- read.csv("carryover_pk_r_summary_results.csv")
summary_rows <- summary_rows[summary_rows$scenario %in% paste0("sce", 1:5), ]

make_expected_sce6 <- function(alpha_values) {
  out <- list()
  p <- c(0.50, 0.55, 0.60, 0.65, 0.70)
  nR <- rep(100L, 5L)
  nI <- rep(100L, 5L)
  
  for (idx in seq_along(alpha_values)) {
    a <- alpha_values[idx]
    q <- pmin(1, a + (1 - a) * p)
    
    dat <- data.frame(
      dose = seq_along(p),
      p_true = p,
      q_true = q,
      nR = nR,
      yR = round(nR * p),
      nI = nI,
      yI = round(nI * q)
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
    
    out[[idx]] <- data.frame(
      scenario = "sce6",
      alpha_true = a,
      dose = est$dose,
      n_reps = NA_integer_,
      p_true = est$p_true,
      q_true = est$q_true,
      mean_yR = est$yR,
      mean_yI = est$yI,
      mean_obs_pooled = est$obs_pooled,
      sd_obs_pooled = NA_real_,
      mean_p_pool_plugin = est$p_pool_plugin,
      sd_p_pool_plugin = NA_real_,
      bias_p_pool_plugin = est$p_pool_plugin - est$p_true,
      mean_p_mle_plugin = est$p_mle_plugin,
      sd_p_mle_plugin = NA_real_,
      bias_p_mle_plugin = est$p_mle_plugin - est$p_true,
      mean_r_hat = fit$r_hat,
      sd_r_hat = NA_real_,
      mean_r_hat_posterior_mean = fit$r$r_mean,
      mean_r_hat_map = fit$r$r_map,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

alpha_values <- c(0, 0.3, 0.6, 0.9)
scenario_names <- paste0("sce", 1:6)

sce6_rows <- make_expected_sce6(alpha_values)
combined_rows <- rbind(summary_rows, sce6_rows)

write.csv(
  combined_rows,
  "carryover_pk_r_summary_results_sce1_to_sce6_for_plots.csv",
  row.names = FALSE
)

files <- plot_carryover_average_alpha_files(
  summary_rows = combined_rows,
  alpha_values = alpha_values,
  scenario_names = scenario_names,
  plot_dir = "carryover_pk_r_demo_plots"
)

cat(paste(files, collapse = "\n"), "\n")
