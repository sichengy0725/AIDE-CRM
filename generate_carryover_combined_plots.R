txt <- readLines("carryover_pk_r_estimators_demo.R")
end_idx <- grep("^## Example test run", txt)[1] - 1L
eval(parse(text = txt[seq_len(end_idx)]))

source("AIDE_CRM_helper_final.R")
source("AIDE_BOIN_helper.R")

parse_integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  out <- suppressWarnings(as.integer(value))
  if (length(out) != 1L || is.na(out)) {
    stop("Environment variable ", name, " must be a single integer.")
  }

  out
}

parse_numeric_vector_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  out <- suppressWarnings(as.numeric(parts))
  if (length(out) == 0L || any(!is.finite(out))) {
    stop("Environment variable ", name, " must be a comma-separated numeric vector.")
  }

  out
}

parse_logical_env <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(isTRUE(default))
  }

  if (value %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }

  stop("Environment variable ", name, " must be true/false or 1/0.")
}

## User settings.  Environment overrides are useful for quick validation runs.
DEFAULT_NR <- parse_integer_env("CARRYOVER_NR", 9L)
DEFAULT_NI <- parse_integer_env("CARRYOVER_NI", 9L)
DEFAULT_N_TRIALS <- parse_integer_env("CARRYOVER_N_TRIALS", 100L)
DEFAULT_SIM_SEED <- parse_integer_env("CARRYOVER_SIM_SEED", 20260704L)
DEFAULT_ALPHA_VALUES <- parse_numeric_vector_env(
  "CARRYOVER_ALPHA_VALUES",
  c(0, 0.3, 0.6, 0.9)
)
DEFAULT_INCLUDE_CRM_RANDOM <- parse_logical_env("CARRYOVER_INCLUDE_CRM_RANDOM", TRUE)
DEFAULT_PLOT_DIR <- Sys.getenv("CARRYOVER_PLOT_DIR", unset = "carryover_pk_r_demo_plots")

TOX_TARGET <- 0.30
BOIN12_UTILITY <- c(
  no_tox_eff = 100,
  no_tox_no_eff = 40,
  tox_eff = 60,
  tox_no_eff = 0
)
BOIN12_UTILITY_PRIOR <- c(1, 1)

BOIN12_SUPPLEMENT_SCENARIOS <- list(
  sce1 = list(
    p_tox = c(0.05, 0.12, 0.27, 0.35, 0.50),
    p_eff = c(0.20, 0.35, 0.36, 0.37, 0.38),
    utility = c(50, 56.2, 50.8, 48.2, 42.8) / 100
  ),
  sce2 = list(
    p_tox = c(0.05, 0.12, 0.27, 0.35, 0.50),
    p_eff = c(0.01, 0.05, 0.30, 0.60, 0.60),
    utility = c(38.6, 38.2, 47.2, 62, 56) / 100
  ),
  sce3 = list(
    p_tox = c(0.03, 0.05, 0.20, 0.22, 0.45),
    p_eff = c(0.05, 0.10, 0.50, 0.68, 0.70),
    utility = c(41.8, 44, 62, 72, 64) / 100
  ),
  sce4 = list(
    p_tox = c(0.05, 0.10, 0.15, 0.20, 0.27),
    p_eff = c(0.01, 0.05, 0.10, 0.20, 0.40),
    utility = c(38.6, 39, 40, 44, 53.2) / 100
  ),
  sce5 = list(
    p_tox = c(0.03, 0.06, 0.10, 0.30, 0.45),
    p_eff = c(0.10, 0.20, 0.40, 0.45, 0.50),
    utility = c(44.8, 49.6, 60, 55, 52) / 100
  )
)

## Adaptive/global r estimator settings.  These match the OC run defaults.
R_ADAPTIVE_PRIOR <- c(TOX_TARGET / 2, 1 - TOX_TARGET / 2)
P_ADAPTIVE_PRIOR <- c(TOX_TARGET / 2, 1 - TOX_TARGET / 2)
R_ADAPTIVE_MAX <- 0.99
R_ADAPTIVE_PLUG_IN <- "mean"
R_ADAPTIVE_REL_TOL <- 1e-6

default_crm_random_settings <- list(
  target = TOX_TARGET,
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

make_aide_alpha_q <- function(p, alpha) {
  q <- p
  if (length(p) >= 2L) {
    q[-1] <- pmin(1, p[-1] + alpha * p[-length(p)])
  }
  q
}

boin12_expected_utility <- function(p_tox,
                                    p_eff,
                                    utility = BOIN12_UTILITY) {
  if (length(p_tox) != length(p_eff)) {
    stop("p_tox and p_eff must have the same length.")
  }
  if (length(utility) != 4L) {
    stop("utility must contain four scores.")
  }

  u1 <- utility[["no_tox_eff"]]
  u2 <- utility[["no_tox_no_eff"]]
  u3 <- utility[["tox_eff"]]
  u4 <- utility[["tox_no_eff"]]

  (
    u1 * (1 - p_tox) * p_eff +
      u2 * (1 - p_tox) * (1 - p_eff) +
      u3 * p_tox * p_eff +
      u4 * p_tox * (1 - p_eff)
  ) / 100
}

simulate_regular_ipde_boin12_counts <- function(p_tox_true,
                                                p_eff_true,
                                                true_utility = NULL,
                                                n_regular = 20,
                                                n_ipde = 20,
                                                alpha_true = 0.1,
                                                seed = NULL,
                                                ipde_model = c("aide_alpha", "document")) {
  ipde_model <- match.arg(ipde_model)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  p_tox_true <- as.numeric(p_tox_true)
  p_eff_true <- as.numeric(p_eff_true)
  K <- length(p_tox_true)

  if (length(p_eff_true) != K) {
    stop("p_eff_true must have the same length as p_tox_true.")
  }
  if (any(!is.finite(p_tox_true)) || any(p_tox_true < 0 | p_tox_true > 1)) {
    stop("p_tox_true must contain probabilities in [0, 1].")
  }
  if (any(!is.finite(p_eff_true)) || any(p_eff_true < 0 | p_eff_true > 1)) {
    stop("p_eff_true must contain probabilities in [0, 1].")
  }
  if (!is.null(true_utility)) {
    true_utility <- as.numeric(true_utility)
    if (length(true_utility) != K ||
        any(!is.finite(true_utility)) ||
        any(true_utility < 0 | true_utility > 1)) {
      stop("true_utility must be NULL or a length-K vector in [0, 1].")
    }
  }
  if (length(alpha_true) != 1L || !is.finite(alpha_true) || alpha_true < 0) {
    stop("alpha_true must be a finite nonnegative scalar.")
  }

  nR <- expand_n_by_dose(n_regular, K, "n_regular")
  nI <- expand_n_by_dose(n_ipde, K, "n_ipde")

  if (ipde_model == "aide_alpha") {
    q_tox_true <- p_tox_true
    if (K >= 2L) {
      q_tox_true[-1] <- pmin(1, p_tox_true[-1] + alpha_true * p_tox_true[-K])
    }
  } else {
    q_tox_true <- pmin(1, alpha_true + (1 - alpha_true) * p_tox_true)
  }

  yR <- stats::rbinom(K, nR, p_tox_true)
  yI <- stats::rbinom(K, nI, q_tox_true)

  r_t_eff <- stats::rbinom(K, yR, p_eff_true)
  r_t_noeff <- yR - r_t_eff
  r_nt_eff <- stats::rbinom(K, nR - yR, p_eff_true)
  r_nt_noeff <- nR - yR - r_nt_eff

  i_t_eff <- stats::rbinom(K, yI, p_eff_true)
  i_t_noeff <- yI - i_t_eff
  i_nt_eff <- stats::rbinom(K, nI - yI, p_eff_true)
  i_nt_noeff <- nI - yI - i_nt_eff

  data.frame(
    dose = seq_len(K),
    p_true = p_tox_true,
    q_true = q_tox_true,
    e_true = p_eff_true,
    true_boin12_utility = if (is.null(true_utility)) {
      boin12_expected_utility(p_tox_true, p_eff_true)
    } else {
      true_utility
    },
    nR = nR,
    yR = yR,
    eR = r_nt_eff + r_t_eff,
    nI = nI,
    yI = yI,
    eI = i_nt_eff + i_t_eff,
    r_nt_eff = r_nt_eff,
    r_nt_noeff = r_nt_noeff,
    r_t_eff = r_t_eff,
    r_t_noeff = r_t_noeff,
    i_nt_eff = i_nt_eff,
    i_nt_noeff = i_nt_noeff,
    i_t_eff = i_t_eff,
    i_t_noeff = i_t_noeff
  )
}

compute_modified_boin12_estimates <- function(dat,
                                              utility = BOIN12_UTILITY,
                                              utility_prior = BOIN12_UTILITY_PRIOR,
                                              r_prior = R_ADAPTIVE_PRIOR,
                                              p_prior = P_ADAPTIVE_PRIOR,
                                              r_max = R_ADAPTIVE_MAX,
                                              r_plug_in = R_ADAPTIVE_PLUG_IN,
                                              r_rel_tol = R_ADAPTIVE_REL_TOL) {
  required <- c(
    "nR", "yR", "eR", "nI", "yI", "eI",
    "r_nt_eff", "r_nt_noeff", "r_t_eff", "r_t_noeff",
    "i_nt_eff", "i_nt_noeff", "i_t_eff", "i_t_noeff"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("dat is missing columns: ", paste(missing, collapse = ", "))
  }

  if (length(utility) != 4L || length(utility_prior) != 2L) {
    stop("utility must have length 4 and utility_prior must have length 2.")
  }

  r_fit <- estimate_global_r_adaptive(
    yR = dat$yR,
    nR = dat$nR,
    yI = dat$yI,
    nI = dat$nI,
    r_prior = r_prior,
    p_prior = p_prior,
    r_max = r_max,
    plug_in = r_plug_in,
    rel.tol = r_rel_tol
  )

  r_use <- clamp(r_fit$r_hat, 0, min(r_max, 1 - 1e-8))

  p2_hat <- vapply(
    seq_len(nrow(dat)),
    function(i) {
      exact_mixed_mle_model2(
        yR = dat$yR[i],
        yI = dat$yI[i],
        N = dat$nR[i] + dat$nI[i],
        NR = dat$nR[i],
        NI = dat$nI[i],
        r_carry = r_use
      )
    },
    numeric(1)
  )

  denom <- r_use + (1 - r_use) * p2_hat
  attribution_weight <- ifelse(
    is.finite(denom) & denom > 0,
    (1 - r_use) * p2_hat / denom,
    0
  )
  attribution_weight <- clamp(attribution_weight, 0, 1)

  u1 <- utility[["no_tox_eff"]]
  u2 <- utility[["no_tox_no_eff"]]
  u3 <- utility[["tox_eff"]]
  u4 <- utility[["tox_no_eff"]]

  regular_utility <- u1 * dat$r_nt_eff +
    u2 * dat$r_nt_noeff +
    u3 * dat$r_t_eff +
    u4 * dat$r_t_noeff

  ipde_raw_utility <- u1 * dat$i_nt_eff +
    u2 * dat$i_nt_noeff +
    u3 * dat$i_t_eff +
    u4 * dat$i_t_noeff

  ipde_modified_utility <- u1 * dat$i_nt_eff +
    u2 * dat$i_nt_noeff +
    (attribution_weight * u3 + (1 - attribution_weight) * u1) * dat$i_t_eff +
    (attribution_weight * u4 + (1 - attribution_weight) * u2) * dat$i_t_noeff

  n_total <- dat$nR + dat$nI
  x_raw <- (regular_utility + ipde_raw_utility) / 100
  x_modified <- (regular_utility + ipde_modified_utility) / 100

  a_u <- utility_prior[1]
  b_u <- utility_prior[2]

  raw_utility_mean <- (x_raw + a_u) / (n_total + a_u + b_u)
  modified_utility_mean <- (x_modified + a_u) / (n_total + a_u + b_u)

  data.frame(
    p2_hat_approx1 = p2_hat,
    r_hat_boin12_adaptive = r_use,
    r_hat_boin12_mean = r_fit$r_mean,
    r_hat_boin12_map = r_fit$r_map,
    r_hat_boin12_sd = r_fit$r_sd,
    attribution_weight = attribution_weight,
    eff_obs = ifelse(n_total > 0, (dat$eR + dat$eI) / n_total, NA_real_),
    boin12_raw_utility = raw_utility_mean,
    boin12_modified_utility = modified_utility_mean,
    boin12_raw_x = x_raw,
    boin12_modified_x = x_modified,
    stringsAsFactors = FALSE
  )
}

summarize_boin12_replicate_estimates <- function(replicate_results) {
  if (is.null(replicate_results) || nrow(replicate_results) == 0L) {
    stop("replicate_results must contain at least one row.")
  }

  keys <- interaction(
    replicate_results$scenario,
    replicate_results$alpha_true,
    replicate_results$dose,
    drop = TRUE
  )

  finite_mean <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0L) {
      return(NA_real_)
    }
    mean(z)
  }

  finite_sd <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) <= 1L) {
      return(NA_real_)
    }
    stats::sd(z)
  }

  out <- do.call(
    rbind,
    lapply(
      split(replicate_results, keys),
      function(x) {
        unique_r_by_rep <- unique(x[, c("replicate", "r_hat_boin12_adaptive")])

        data.frame(
          scenario = x$scenario[1],
          alpha_true = x$alpha_true[1],
          dose = x$dose[1],
          e_true = x$e_true[1],
          true_boin12_utility = x$true_boin12_utility[1],
          mean_eR = mean(x$eR),
          mean_eI = mean(x$eI),
          mean_eff_obs = finite_mean(x$eff_obs),
          sd_eff_obs = finite_sd(x$eff_obs),
          mean_p2_hat_approx1 = finite_mean(x$p2_hat_approx1),
          sd_p2_hat_approx1 = finite_sd(x$p2_hat_approx1),
          mean_attribution_weight = finite_mean(x$attribution_weight),
          sd_attribution_weight = finite_sd(x$attribution_weight),
          mean_boin12_raw_utility = finite_mean(x$boin12_raw_utility),
          sd_boin12_raw_utility = finite_sd(x$boin12_raw_utility),
          bias_boin12_raw_utility = finite_mean(
            x$boin12_raw_utility - x$true_boin12_utility
          ),
          mean_boin12_modified_utility = finite_mean(x$boin12_modified_utility),
          sd_boin12_modified_utility = finite_sd(x$boin12_modified_utility),
          bias_boin12_modified_utility = finite_mean(
            x$boin12_modified_utility - x$true_boin12_utility
          ),
          mean_r_hat_boin12_adaptive = finite_mean(unique_r_by_rep$r_hat_boin12_adaptive),
          sd_r_hat_boin12_adaptive = finite_sd(unique_r_by_rep$r_hat_boin12_adaptive),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(out) <- NULL
  out[order(out$scenario, out$alpha_true, out$dose), ]
}

plot_modified_boin12_average_by_alpha <- function(summary_rows,
                                                  alpha_true,
                                                  scenario_names = NULL,
                                                  file = NULL,
                                                  ylim = c(0, 1),
                                                  panel_ncol = 3L) {
  dat <- summary_rows[summary_rows$alpha_true == alpha_true, , drop = FALSE]

  if (!is.null(scenario_names)) {
    dat <- dat[dat$scenario %in% scenario_names, , drop = FALSE]
  } else {
    scenario_names <- unique(dat$scenario)
  }

  if (nrow(dat) == 0L) {
    stop("No rows to plot for alpha_true = ", alpha_true, ".")
  }

  required <- c(
    "true_boin12_utility",
    "mean_boin12_raw_utility",
    "mean_boin12_modified_utility",
    "e_true",
    "mean_eff_obs"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("summary_rows is missing columns: ", paste(missing, collapse = ", "))
  }

  scenario_names <- scenario_names[scenario_names %in% unique(dat$scenario)]
  panel_ncol <- as.integer(panel_ncol)
  panel_nrow <- ceiling(length(scenario_names) / panel_ncol)

  if (!is.null(file)) {
    grDevices::png(file, width = 1500, height = 900, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(
    mfrow = c(panel_nrow, panel_ncol),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 3, 0)
  )

  for (i in seq_along(scenario_names)) {
    scenario_name <- scenario_names[i]
    one <- dat[dat$scenario == scenario_name, , drop = FALSE]
    one <- one[order(one$dose), ]

    main <- paste0(
      scenario_name,
      ", mean r_hat=",
      sprintf("%.3f", unique(one$mean_r_hat_boin12_adaptive)[1])
    )

    graphics::plot(
      one$dose,
      one$true_boin12_utility,
      type = "b",
      pch = 16,
      lwd = 2,
      ylim = ylim,
      xlab = "Dose level",
      ylab = "Standardized utility / efficacy probability",
      main = main
    )

    graphics::lines(
      one$dose,
      one$mean_boin12_raw_utility,
      type = "b",
      pch = 1,
      lwd = 2,
      lty = 3,
      col = "gray35"
    )

    graphics::lines(
      one$dose,
      one$mean_boin12_modified_utility,
      type = "b",
      pch = 15,
      lwd = 2,
      col = "#CC79A7"
    )

    graphics::lines(
      one$dose,
      one$e_true,
      type = "b",
      pch = 17,
      lwd = 2,
      lty = 2,
      col = "#009E73"
    )

    graphics::lines(
      one$dose,
      one$mean_eff_obs,
      type = "b",
      pch = 2,
      lwd = 2,
      lty = 3,
      col = "#0072B2"
    )

    if (i == 1L) {
      graphics::legend(
        "topleft",
        bty = "n",
        cex = 0.72,
        lwd = 2,
        pch = c(16, 1, 15, 17, 2),
        lty = c(1, 3, 1, 2, 3),
        col = c("black", "gray35", "#CC79A7", "#009E73", "#0072B2"),
        legend = c(
          "True BOIN12 utility",
          "Raw BOIN12 utility",
          "Modified BOIN12 utility",
          "True efficacy p_E",
          "Observed efficacy p_E"
        )
      )
    }
  }

  graphics::mtext(
    paste0("Modified BOIN12 utility average curves, alpha = ", alpha_true),
    outer = TRUE,
    cex = 1.2,
    font = 2
  )

  invisible(file)
}

plot_modified_boin12_average_alpha_files <- function(summary_rows,
                                                     alpha_values,
                                                     scenario_names,
                                                     plot_dir = "carryover_pk_r_demo_plots",
                                                     file_suffix = "") {
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  if (length(file_suffix) != 1L || is.na(file_suffix)) {
    stop("file_suffix must be a single non-NA string.")
  }

  file_suffix <- as.character(file_suffix)

  files <- vapply(
    alpha_values,
    function(a) {
      plot_file <- file.path(
        plot_dir,
        paste0(
          "modified_boin12_supplement_s4_sce1_to_sce5_alpha_",
          gsub("\\.", "p", a),
          file_suffix,
          "_avg.png"
        )
      )

      plot_modified_boin12_average_by_alpha(
        summary_rows = summary_rows,
        alpha_true = a,
        scenario_names = scenario_names,
        file = plot_file
      )

      plot_file
    },
    character(1)
  )

  names(files) <- paste0("alpha_", alpha_values)
  files
}

plot_toxicity_supplement_average_alpha_files <- function(summary_rows,
                                                         alpha_values,
                                                         scenario_names,
                                                         plot_dir = "carryover_pk_r_demo_plots",
                                                         file_suffix = "") {
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  if (length(file_suffix) != 1L || is.na(file_suffix)) {
    stop("file_suffix must be a single non-NA string.")
  }

  file_suffix <- as.character(file_suffix)

  files <- vapply(
    alpha_values,
    function(a) {
      plot_file <- file.path(
        plot_dir,
        paste0(
          "toxicity_supplement_s4_sce1_to_sce5_alpha_",
          gsub("\\.", "p", a),
          file_suffix,
          "_avg.png"
        )
      )

      plot_carryover_average_by_alpha(
        summary_rows = summary_rows,
        alpha_true = a,
        scenario_names = scenario_names,
        file = plot_file
      )

      plot_file
    },
    character(1)
  )

  names(files) <- paste0("alpha_", alpha_values)
  files
}

fit_binomial_trial_summary <- function(scenario_name,
                                       p,
                                       p_eff,
                                       true_utility,
                                       alpha,
                                       n_regular,
                                       n_ipde,
                                       n_trials,
                                       seed,
                                       include_crm_random = FALSE,
                                       crm_random_settings = default_crm_random_settings,
                                       r_adaptive_prior = R_ADAPTIVE_PRIOR,
                                       p_adaptive_prior = P_ADAPTIVE_PRIOR,
                                       r_adaptive_max = R_ADAPTIVE_MAX,
                                       r_adaptive_plug_in = R_ADAPTIVE_PLUG_IN,
                                       r_adaptive_rel_tol = R_ADAPTIVE_REL_TOL,
                                       utility = BOIN12_UTILITY,
                                       utility_prior = BOIN12_UTILITY_PRIOR) {
  all_trials <- vector("list", n_trials)

  for (trial_id in seq_len(n_trials)) {
    trial_seed <- seed + trial_id

    dat <- simulate_regular_ipde_boin12_counts(
      p_tox_true = p,
      p_eff_true = p_eff,
      true_utility = true_utility,
      n_regular = n_regular,
      n_ipde = n_ipde,
      alpha_true = alpha,
      seed = trial_seed,
      ipde_model = "aide_alpha"
    )

    fit <- estimate_pk_with_global_r(
      dat = dat,
      r_prior = r_adaptive_prior,
      p_prior = p_adaptive_prior,
      r_max = r_adaptive_max,
      plug_in = r_adaptive_plug_in,
      debug_on_error = FALSE
    )

    boin12_fit <- compute_modified_boin12_estimates(
      dat = dat,
      utility = utility,
      utility_prior = utility_prior,
      r_prior = r_adaptive_prior,
      p_prior = p_adaptive_prior,
      r_max = r_adaptive_max,
      r_plug_in = r_adaptive_plug_in,
      r_rel_tol = r_adaptive_rel_tol
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
      eR = dat$eR,
      nI = est$nI,
      yI = est$yI,
      eI = dat$eI,
      p_true = est$p_true,
      q_true = est$q_true,
      e_true = dat$e_true,
      true_boin12_utility = dat$true_boin12_utility,
      obs_pooled = est$obs_pooled,
      p_pool_plugin = est$p_pool_plugin,
      p_mle_plugin = est$p_mle_plugin,
      p_crm_random = crm_p_hat,
      r_hat_crm_random = crm_r_hat,
      r_nt_eff = dat$r_nt_eff,
      r_nt_noeff = dat$r_nt_noeff,
      r_t_eff = dat$r_t_eff,
      r_t_noeff = dat$r_t_noeff,
      i_nt_eff = dat$i_nt_eff,
      i_nt_noeff = dat$i_nt_noeff,
      i_t_eff = dat$i_t_eff,
      i_t_noeff = dat$i_t_noeff,
      p2_hat_approx1 = boin12_fit$p2_hat_approx1,
      r_hat_boin12_adaptive = boin12_fit$r_hat_boin12_adaptive,
      r_hat_boin12_mean = boin12_fit$r_hat_boin12_mean,
      r_hat_boin12_map = boin12_fit$r_hat_boin12_map,
      r_hat_boin12_sd = boin12_fit$r_hat_boin12_sd,
      attribution_weight = boin12_fit$attribution_weight,
      eff_obs = boin12_fit$eff_obs,
      boin12_raw_utility = boin12_fit$boin12_raw_utility,
      boin12_modified_utility = boin12_fit$boin12_modified_utility,
      boin12_raw_x = boin12_fit$boin12_raw_x,
      boin12_modified_x = boin12_fit$boin12_modified_x,
      stringsAsFactors = FALSE
    )
  }

  trial_results <- do.call(rbind, all_trials)
  rownames(trial_results) <- NULL

  tox_summary <- summarize_replicate_estimates(trial_results)
  boin12_summary <- summarize_boin12_replicate_estimates(trial_results)

  out <- merge(
    tox_summary,
    boin12_summary,
    by = c("scenario", "alpha_true", "dose"),
    all.x = TRUE,
    sort = FALSE
  )

  out[order(out$scenario, out$alpha_true, out$dose), ]
}

make_supplement_scenario_vectors <- function(supplement_scenarios = BOIN12_SUPPLEMENT_SCENARIOS) {
  if (length(supplement_scenarios) == 0L) {
    stop("supplement_scenarios must contain at least one scenario.")
  }

  scenario_names <- names(supplement_scenarios)
  if (is.null(scenario_names) || any(!nzchar(scenario_names))) {
    scenario_names <- paste0("sce", seq_along(supplement_scenarios))
  }

  p_tox <- setNames(vector("list", length(supplement_scenarios)), scenario_names)
  p_eff <- setNames(vector("list", length(supplement_scenarios)), scenario_names)
  utility <- setNames(vector("list", length(supplement_scenarios)), scenario_names)

  for (scenario_name in scenario_names) {
    one <- supplement_scenarios[[scenario_name]]
    required <- c("p_tox", "p_eff", "utility")
    missing <- setdiff(required, names(one))
    if (length(missing) > 0L) {
      stop(
        "Supplement scenario ",
        scenario_name,
        " is missing: ",
        paste(missing, collapse = ", ")
      )
    }

    p_tox[[scenario_name]] <- one$p_tox
    p_eff[[scenario_name]] <- one$p_eff
    utility[[scenario_name]] <- one$utility
  }

  list(
    p_tox = p_tox,
    p_eff = p_eff,
    utility = utility,
    scenario_names = scenario_names
  )
}

generate_carryover_combined_plot <- function(NR = DEFAULT_NR,
                                             NI = DEFAULT_NI,
                                             N_TRIALS = DEFAULT_N_TRIALS,
                                             SIM_SEED = DEFAULT_SIM_SEED,
                                             alpha_values = DEFAULT_ALPHA_VALUES,
                                             include_crm_random = DEFAULT_INCLUDE_CRM_RANDOM,
                                             crm_random_settings = default_crm_random_settings,
                                             plot_dir = DEFAULT_PLOT_DIR,
                                             supplement_scenarios = BOIN12_SUPPLEMENT_SCENARIOS,
                                             summary_output = "carryover_pk_r_summary_results_supplement_s4_sce1_to_sce5_for_plots.csv",
                                             r_adaptive_prior = R_ADAPTIVE_PRIOR,
                                             p_adaptive_prior = P_ADAPTIVE_PRIOR,
                                             r_adaptive_max = R_ADAPTIVE_MAX,
                                             r_adaptive_plug_in = R_ADAPTIVE_PLUG_IN,
                                             r_adaptive_rel_tol = R_ADAPTIVE_REL_TOL,
                                             utility = BOIN12_UTILITY,
                                             utility_prior = BOIN12_UTILITY_PRIOR) {
  if (length(N_TRIALS) != 1L || !is.finite(N_TRIALS) || N_TRIALS < 1L) {
    stop("N_TRIALS must be a positive integer.")
  }

  N_TRIALS <- as.integer(N_TRIALS)

  if (isTRUE(include_crm_random) &&
      (!requireNamespace("rjags", quietly = TRUE) ||
       !requireNamespace("coda", quietly = TRUE))) {
    warning(
      "Random CRM comparator requires rjags and coda; continuing without it.",
      call. = FALSE
    )
    include_crm_random <- FALSE
  }

  if (is.null(alpha_values)) {
    alpha_values <- DEFAULT_ALPHA_VALUES
  }

  scenario_inputs <- make_supplement_scenario_vectors(supplement_scenarios)
  scenarios <- scenario_inputs$p_tox
  efficacy_scenarios <- scenario_inputs$p_eff
  utility_scenarios <- scenario_inputs$utility
  scenario_names <- scenario_inputs$scenario_names

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
        p_eff = efficacy_scenarios[[scenario_name]],
        true_utility = utility_scenarios[[scenario_name]],
        alpha = alpha,
        n_regular = NR,
        n_ipde = NI,
        n_trials = N_TRIALS,
        seed = SIM_SEED + idx * 100000L,
        include_crm_random = include_crm_random,
        crm_random_settings = crm_random_settings,
        r_adaptive_prior = r_adaptive_prior,
        p_adaptive_prior = p_adaptive_prior,
        r_adaptive_max = r_adaptive_max,
        r_adaptive_plug_in = r_adaptive_plug_in,
        r_adaptive_rel_tol = r_adaptive_rel_tol,
        utility = utility,
        utility_prior = utility_prior
      )
      idx <- idx + 1L
    }
  }

  combined_rows <- do.call(rbind, all_rows)
  rownames(combined_rows) <- NULL

  write.csv(
    combined_rows,
    summary_output,
    row.names = FALSE
  )

  file_suffix <- paste0("_NR", NR, "_NI", NI)

  tox_files <- plot_toxicity_supplement_average_alpha_files(
    summary_rows = combined_rows,
    alpha_values = alpha_values,
    scenario_names = scenario_names,
    plot_dir = plot_dir,
    file_suffix = file_suffix
  )

  boin12_files <- plot_modified_boin12_average_alpha_files(
    summary_rows = combined_rows,
    alpha_values = alpha_values,
    scenario_names = scenario_names,
    plot_dir = plot_dir,
    file_suffix = file_suffix
  )

  files <- c(toxicity = tox_files, modified_boin12 = boin12_files)
  cat(paste(unname(files), collapse = "\n"), "\n")

  invisible(list(
    summary_rows = combined_rows,
    toxicity_files = tox_files,
    modified_boin12_files = boin12_files
  ))
}

if (parse_logical_env("CARRYOVER_COMBINED_AUTORUN", TRUE)) {
  generate_carryover_combined_plot()
}
