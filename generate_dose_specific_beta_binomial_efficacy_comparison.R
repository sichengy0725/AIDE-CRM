## Beta-binomial efficacy simulation and estimation
##
## This file defines three functions only:
##   1. generate_efficacy_data()
##   2. fit_beta_binomial_efficacy()
##   3. run_efficacy_trial()
##
## The two JAGS models are stored beside this file:
##   beta_binomial_dose_specific.jags
##   beta_binomial_hierarchical.jags

#' Generate patient-level regular and IPDE efficacy data.
#'
#' @param p_regular True regular-patient efficacy probabilities by dose.
#' @param p_ipde True IPDE efficacy probabilities by dose.
#' @param ndose Number of dose levels.
#' @param n_ipde Number of IPDE patients per dose (a scalar or length-ndose vector).
#' @param n_regular Number of regular patients per dose (a scalar or length-ndose vector).
#' @param seed Optional random-number seed.
#' @return A patient-level data frame with dose, group, efficacy, and true
#'   probabilities.
generate_efficacy_data <- function(p_regular,
                                   p_ipde,
                                   ndose,
                                   n_ipde,
                                   n_regular,
                                   seed = NULL) {
  ndose <- as.integer(ndose)
  if (length(ndose) != 1L || is.na(ndose) || ndose < 2L) {
    stop("ndose must be one integer of at least 2.")
  }

  p_regular <- as.numeric(p_regular)
  p_ipde <- as.numeric(p_ipde)
  if (length(p_regular) != ndose || length(p_ipde) != ndose ||
      any(!is.finite(p_regular)) || any(!is.finite(p_ipde)) ||
      any(p_regular <= 0 | p_regular >= 1) ||
      any(p_ipde <= 0 | p_ipde >= 1)) {
    stop("p_regular and p_ipde must be length-ndose probability vectors in (0, 1).")
  }

  expand_count <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, ndose)
    }
    if (length(x) != ndose || any(!is.finite(x)) || any(x < 0) ||
        any(abs(x - round(x)) > 1e-8)) {
      stop(name, " must be a nonnegative integer scalar or length-ndose vector.")
    }
    as.integer(round(x))
  }

  n_ipde <- expand_count(n_ipde, "n_ipde")
  n_regular <- expand_count(n_regular, "n_regular")
  if (sum(n_ipde) + sum(n_regular) == 0L) {
    stop("At least one patient must be generated.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  regular_rows <- lapply(seq_len(ndose), function(j) {
    data.frame(
      dose = rep.int(j, n_regular[j]),
      group = rep("regular", n_regular[j]),
      efficacy = stats::rbinom(n_regular[j], size = 1L, prob = p_regular[j]),
      p_regular_true = rep(p_regular[j], n_regular[j]),
      p_ipde_true = rep(p_ipde[j], n_regular[j]),
      stringsAsFactors = FALSE
    )
  })
  ipde_rows <- lapply(seq_len(ndose), function(j) {
    data.frame(
      dose = rep.int(j, n_ipde[j]),
      group = rep("ipde", n_ipde[j]),
      efficacy = stats::rbinom(n_ipde[j], size = 1L, prob = p_ipde[j]),
      p_regular_true = rep(p_regular[j], n_ipde[j]),
      p_ipde_true = rep(p_ipde[j], n_ipde[j]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, c(regular_rows, ipde_rows))
  out <- out[order(out$dose, out$group), , drop = FALSE]
  rownames(out) <- NULL
  out$patient_id <- seq_len(nrow(out))
  out[, c("patient_id", "dose", "group", "efficacy", "p_regular_true", "p_ipde_true")]
}

#' Fit a shared-carryover, dose-specific-carryover, or hierarchical
#' beta-binomial efficacy model in JAGS.
#'
#' For the dose-specific model, a_r and b_r define independent priors for
#' regular efficacy pi_Ej, while a_carry and b_carry define independent priors
#' for the IPDE carryover. The dose-specific model uses r_Ej at dose j, while
#' the shared_carryover model uses one r_E across all dose levels. Both use
#' pi*_Ej = r_E + (1-r_E)pi_Ej. In the hierarchical model, a_carry/b_carry
#' are unused.
#'
#' @param dose_data Output from generate_efficacy_data().
#' @param a_r,b_r Positive Beta shapes for regular efficacy.
#' @param a_carry,b_carry Positive Beta shapes for IPDE carryover. For
#'   shared_carryover, these must be scalars or constant across dose levels.
#' @param prior_type "dose_specific", "shared_carryover", or "hierarchical".
#' @param model_file Optional path to the corresponding JAGS model file.
#' @param n_chains,n_adapt,n_burnin,n_iter,thin JAGS MCMC settings.
#' @param seed Optional JAGS random-number seed.
#' @param ndose Optional total number of doses. This permits unobserved dose
#'   levels to remain in the JAGS model with their prior only.
#' @return A list containing posterior summaries by dose, posterior samples,
#'   dose-level counts, and the selected model type.
fit_beta_binomial_efficacy <- function(dose_data,
                                       a_r,
                                       b_r,
                                        prior_type = c("dose_specific", "shared_carryover", "hierarchical"),
                                       a_carry = 1,
                                       b_carry = 9,
                                       model_file = NULL,
                                       n_chains = 3L,
                                       n_adapt = 1000L,
                                       n_burnin = 1000L,
                                       n_iter = 4000L,
                                       thin = 2L,
                                       seed = NULL,
                                       ndose = NULL) {
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("This function requires the rjags package and a working JAGS installation.")
  }
  prior_type <- match.arg(prior_type)
  required <- c("dose", "group", "efficacy")
  missing <- setdiff(required, names(dose_data))
  if (length(missing) > 0L) {
    stop("dose_data is missing: ", paste(missing, collapse = ", "))
  }

  dose_data <- dose_data[, required, drop = FALSE]
  dose_data$dose <- as.integer(dose_data$dose)
  dose_data$group <- tolower(as.character(dose_data$group))
  dose_data$efficacy <- as.numeric(dose_data$efficacy)
  if (any(is.na(dose_data$dose)) || any(dose_data$dose < 1L) ||
      any(!dose_data$group %in% c("regular", "ipde")) ||
      any(!dose_data$efficacy %in% c(0, 1))) {
    stop("dose_data must contain positive doses, regular/ipde groups, and binary efficacy outcomes.")
  }

  observed_ndose <- if (nrow(dose_data) == 0L) 0L else max(dose_data$dose)
  if (is.null(ndose)) {
    ndose <- observed_ndose
  }
  ndose <- as.integer(ndose)
  if (length(ndose) != 1L || is.na(ndose) || ndose < 2L || observed_ndose > ndose) {
    stop("ndose must be an integer of at least 2 and no smaller than the largest observed dose.")
  }

  expand_shape <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, ndose)
    }
    if (length(x) != ndose || any(!is.finite(x)) || any(x <= 0)) {
      stop(name, " must be a positive scalar or length-ndose vector.")
    }
    x
  }
  a_r <- expand_shape(a_r, "a_r")
  b_r <- expand_shape(b_r, "b_r")
  a_carry <- expand_shape(a_carry, "a_carry")
  b_carry <- expand_shape(b_carry, "b_carry")
  if (prior_type == "shared_carryover" &&
      (any(a_carry != a_carry[1L]) || any(b_carry != b_carry[1L]))) {
    stop("For shared_carryover, a_carry and b_carry must be scalar or constant across dose levels.")
  }

  count_by_dose <- function(group_name) {
    vapply(seq_len(ndose), function(j) {
      sum(dose_data$group == group_name & dose_data$dose == j)
    }, numeric(1))
  }
  response_by_dose <- function(group_name) {
    vapply(seq_len(ndose), function(j) {
      sum(dose_data$efficacy[dose_data$group == group_name & dose_data$dose == j])
    }, numeric(1))
  }
  n_regular <- count_by_dose("regular")
  y_regular <- response_by_dose("regular")
  n_ipde <- count_by_dose("ipde")
  y_ipde <- response_by_dose("ipde")

  if (is.null(model_file)) {
    model_file <- switch(
      prior_type,
      dose_specific = "beta_binomial_dose_specific.jags",
      shared_carryover = "beta_binomial_shared_carryover.jags",
      hierarchical = "beta_binomial_hierarchical.jags"
    )
  }
  if (!file.exists(model_file)) {
    stop("JAGS model file not found: ", model_file)
  }

  n_chains <- as.integer(n_chains)
  n_adapt <- as.integer(n_adapt)
  n_burnin <- as.integer(n_burnin)
  n_iter <- as.integer(n_iter)
  thin <- as.integer(thin)
  if (any(is.na(c(n_chains, n_adapt, n_burnin, n_iter, thin))) ||
      n_chains < 1L || n_adapt < 0L || n_burnin < 0L || n_iter < 1L || thin < 1L) {
    stop("JAGS MCMC settings must be valid positive integers (or zero for adaptation/burn-in).")
  }

  jags_data <- list(
    ndose = ndose,
    n_regular = n_regular,
    y_regular = y_regular,
    n_ipde = n_ipde,
    y_ipde = y_ipde
  )
  if (prior_type %in% c("dose_specific", "shared_carryover")) {
    jags_data$a_r <- a_r
    jags_data$b_r <- b_r
    if (prior_type == "dose_specific") {
      jags_data$a_carry <- a_carry
      jags_data$b_carry <- b_carry
    } else {
      jags_data$a_carry <- a_carry[1L]
      jags_data$b_carry <- b_carry[1L]
    }
  } else {
    jags_data$a0 <- mean(a_r)
    jags_data$b0 <- mean(b_r)
    jags_data$kappa0 <- mean(a_r + b_r)
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (is.na(seed) || seed < 1L) {
      stop("seed must be a positive integer.")
    }
    inits <- lapply(seq_len(n_chains), function(chain) {
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + chain)
    })
  } else {
    inits <- NULL
  }

  jags_model <- rjags::jags.model(
    file = model_file,
    data = jags_data,
    inits = inits,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  if (n_burnin > 0L) {
    stats::update(jags_model, n.iter = n_burnin, progress.bar = "none")
  }
  monitored_parameters <- if (prior_type == "dose_specific") {
    c("p_regular", "p_ipde", "r_e")
  } else if (prior_type == "shared_carryover") {
    c("p_regular", "p_ipde", "r_e")
  } else {
    c("p_regular", "p_ipde")
  }
  samples <- rjags::coda.samples(
    model = jags_model,
    variable.names = monitored_parameters,
    n.iter = n_iter,
    thin = thin,
    progress.bar = "none"
  )
  draws <- as.matrix(samples)
  summarize_probability <- function(parameter_name, dose) {
    draw_column <- draws[, paste0(parameter_name, "[", dose, "]")]
    c(
      mean = mean(draw_column),
      lcl = stats::quantile(draw_column, 0.025, names = FALSE),
      ucl = stats::quantile(draw_column, 0.975, names = FALSE)
    )
  }
  regular_summary <- t(vapply(
    seq_len(ndose),
    function(j) summarize_probability("p_regular", j),
    numeric(3)
  ))
  ipde_summary <- t(vapply(
    seq_len(ndose),
    function(j) summarize_probability("p_ipde", j),
    numeric(3)
  ))
  carryover_summary <- if (prior_type == "dose_specific") {
    t(vapply(
      seq_len(ndose),
      function(j) summarize_probability("r_e", j),
      numeric(3)
    ))
  } else if (prior_type == "shared_carryover") {
    shared_summary <- c(
      mean = mean(draws[, "r_e"]),
      lcl = stats::quantile(draws[, "r_e"], 0.025, names = FALSE),
      ucl = stats::quantile(draws[, "r_e"], 0.975, names = FALSE)
    )
    matrix(
      rep(shared_summary, each = ndose), nrow = ndose,
      dimnames = list(NULL, c("mean", "lcl", "ucl"))
    )
  } else {
    matrix(NA_real_, nrow = ndose, ncol = 3L,
           dimnames = list(NULL, c("mean", "lcl", "ucl")))
  }

  list(
    posterior = data.frame(
      dose = seq_len(ndose),
      p_regular_hat = regular_summary[, "mean"],
      p_regular_lcl = regular_summary[, "lcl"],
      p_regular_ucl = regular_summary[, "ucl"],
      p_ipde_hat = ipde_summary[, "mean"],
      p_ipde_lcl = ipde_summary[, "lcl"],
      p_ipde_ucl = ipde_summary[, "ucl"],
      r_e_hat = carryover_summary[, "mean"],
      r_e_lcl = carryover_summary[, "lcl"],
      r_e_ucl = carryover_summary[, "ucl"],
      stringsAsFactors = FALSE
    ),
    samples = samples,
    counts = data.frame(
      dose = seq_len(ndose),
      n_regular = n_regular,
      y_regular = y_regular,
      n_ipde = n_ipde,
      y_ipde = y_ipde
    ),
    prior_type = prior_type
  )
}

#' Run repeated efficacy trials and plot average regular-patient estimates.
#'
#' @param p_regular,p_ipde,ndose,n_ipde,n_regular Data-generating inputs passed
#'   to generate_efficacy_data().
#' @param ntrial Number of simulated trials.
#' @param a_r,b_r,a_carry,b_carry Beta prior shapes passed to fit_beta_binomial_efficacy().
#' @param prior_type JAGS model to fit: "dose_specific", "shared_carryover",
#'   or "hierarchical".
#' @param plot_file Optional PNG output path.  Use NULL to display the plot on
#'   the active graphics device.
#' @return A list containing the average posterior estimates and all trial-level
#'   posterior estimates.
run_efficacy_trial <- function(p_regular,
                               p_ipde,
                               ndose,
                               n_ipde,
                               n_regular,
                               ntrial,
                               a_r,
                               b_r,
                                prior_type = c("dose_specific", "shared_carryover", "hierarchical"),
                               a_carry = 1,
                               b_carry = 9,
                               model_file = NULL,
                               n_chains = 3L,
                               n_adapt = 1000L,
                               n_burnin = 1000L,
                               n_iter = 4000L,
                               thin = 2L,
                               seed = NULL,
                               plot_file = NULL,
                               make_plot = TRUE) {
  prior_type <- match.arg(prior_type)
  ntrial <- as.integer(ntrial)
  if (length(ntrial) != 1L || is.na(ntrial) || ntrial < 1L) {
    stop("ntrial must be a positive integer.")
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (is.na(seed) || seed < 1L) {
      stop("seed must be a positive integer.")
    }
  }

  trial_results <- vector("list", ntrial)
  for (trial in seq_len(ntrial)) {
    trial_seed <- if (is.null(seed)) NULL else seed + trial * 10L
    dose_data <- generate_efficacy_data(
      p_regular = p_regular,
      p_ipde = p_ipde,
      ndose = ndose,
      n_ipde = n_ipde,
      n_regular = n_regular,
      seed = trial_seed
    )
    fit <- fit_beta_binomial_efficacy(
      dose_data = dose_data,
      a_r = a_r,
      b_r = b_r,
      a_carry = a_carry,
      b_carry = b_carry,
      prior_type = prior_type,
      model_file = model_file,
      n_chains = n_chains,
      n_adapt = n_adapt,
      n_burnin = n_burnin,
      n_iter = n_iter,
      thin = thin,
      seed = if (is.null(trial_seed)) NULL else trial_seed + 1L,
      ndose = ndose
    )
    one <- fit$posterior
    one$trial <- trial
    trial_results[[trial]] <- one
  }

  trial_estimates <- do.call(rbind, trial_results)
  p_regular <- as.numeric(p_regular)
  if (length(p_regular) != as.integer(ndose)) {
    stop("p_regular must have length ndose.")
  }
  average_estimate <- do.call(rbind, lapply(split(trial_estimates, trial_estimates$dose), function(one) {
    data.frame(
      dose = one$dose[1L],
      true_regular_probability = p_regular[one$dose[1L]],
      average_regular_estimate = mean(one$p_regular_hat),
      average_ipde_estimate = mean(one$p_ipde_hat),
      stringsAsFactors = FALSE
    )
  }))
  average_estimate <- average_estimate[order(average_estimate$dose), ]
  rownames(average_estimate) <- NULL

  if (isTRUE(make_plot)) {
    if (!is.null(plot_file)) {
      grDevices::png(plot_file, width = 1000, height = 700, res = 120)
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit({
        graphics::par(old_par)
        grDevices::dev.off()
      }, add = TRUE)
    }
    graphics::plot(
      average_estimate$dose,
      average_estimate$true_regular_probability,
      type = "b",
      pch = 16,
      lwd = 2,
      ylim = c(0, 1),
      xlab = "Dose level",
      ylab = "Regular-patient efficacy probability",
      main = paste("Average beta-binomial estimate:", prior_type, "prior")
    )
    graphics::lines(
      average_estimate$dose,
      average_estimate$average_regular_estimate,
      type = "b",
      pch = 1,
      lwd = 2,
      lty = 2,
      col = "#0072B2"
    )
    graphics::legend(
      "topleft",
      bty = "n",
      legend = c("True regular probability", "Average posterior estimate"),
      col = c("black", "#0072B2"),
      lty = c(1, 2),
      lwd = 2,
      pch = c(16, 1)
    )
  }

  invisible(list(
    average_estimate = average_estimate,
    trial_estimates = trial_estimates,
    prior_type = prior_type,
    plot_file = plot_file
  ))
}
## Interactive prior-sensitivity experiment.  It uses only the efficacy
## scenarios and creates two graph sets: n_regular = n_ipde = 3 and 6.
## Set this to TRUE only after rjags is installed and you are ready to run the
## 16,000 JAGS fits implied by 5 scenarios x 4 alphas x 4 priors x 2 sample
## sizes x 100 trials.
RUN_BOIN12_PRIOR_SENSITIVITY_EXPERIMENT <- FALSE

if (isTRUE(RUN_BOIN12_PRIOR_SENSITIVITY_EXPERIMENT)) {
alpha_values <- c(0, 0.3, 0.6, 0.9)
ntrial <- 100L
prior_type <- "dose_specific"  # Set to "hierarchical" to fit that JAGS model.

BOIN12_EFFICACY_SCENARIOS <- list(
  sce1 = c(0.20, 0.35, 0.36, 0.37, 0.38),
  sce2 = c(0.01, 0.05, 0.30, 0.60, 0.60),
  sce3 = c(0.05, 0.10, 0.50, 0.68, 0.70),
  sce4 = c(0.01, 0.05, 0.10, 0.20, 0.40),
  sce5 = c(0.10, 0.20, 0.40, 0.45, 0.50)
)

prior_settings <- data.frame(
  prior_label = c("Beta(0.15, 0.85)", "Beta(0.05, 0.95)", "Beta(0.30, 0.70)", "Beta(0.50, 0.50)"),
  a_r = c(0.15, 0.05, 0.30, 0.50),
  b_r = c(0.85, 0.95, 0.70, 0.50),
  stringsAsFactors = FALSE
)

patient_settings <- data.frame(
  patient_label = c("n_regular_3_n_ipde_3", "n_regular_6_n_ipde_6"),
  n_regular = c(3L, 6L),
  n_ipde = c(3L, 6L),
  stringsAsFactors = FALSE
)

output_dir <- "boin12_efficacy_prior_sensitivity"
dir.create(output_dir, showWarnings = FALSE)

all_results <- list()
plot_files <- character(0)
result_id <- 1L
prior_colours <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")

for (patient_index in seq_len(nrow(patient_settings))) {
  patient_setting <- patient_settings[patient_index, ]
  patient_dir <- file.path(output_dir, patient_setting$patient_label)
  dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)

  for (alpha_eff in alpha_values) {
    plot_rows <- list()
    plot_row_id <- 1L

    for (scenario_name in names(BOIN12_EFFICACY_SCENARIOS)) {
      p_regular <- BOIN12_EFFICACY_SCENARIOS[[scenario_name]]
      ndose <- length(p_regular)
      p_ipde <- pmin(1 - 1e-8, p_regular + alpha_eff * c(0, p_regular[-ndose]))
      p_ipde[1L] <- p_regular[1L]

      for (prior_index in seq_len(nrow(prior_settings))) {
        prior_setting <- prior_settings[prior_index, ]
        result <- run_efficacy_trial(
          p_regular = p_regular,
          p_ipde = p_ipde,
          ndose = ndose,
          n_regular = patient_setting$n_regular,
          n_ipde = patient_setting$n_ipde,
          ntrial = ntrial,
          a_r = prior_setting$a_r,
          b_r = prior_setting$b_r,
          prior_type = prior_type,
          seed = 20260718L + result_id * 1000L,
          make_plot = FALSE
        )

        one <- result$average_estimate
        one$scenario <- scenario_name
        one$alpha_eff <- alpha_eff
        one$prior_type <- prior_type
        one$prior_label <- prior_setting$prior_label
        one$a_r <- prior_setting$a_r
        one$b_r <- prior_setting$b_r
        one$n_regular <- patient_setting$n_regular
        one$n_ipde <- patient_setting$n_ipde
        one$true_ipde_probability <- p_ipde

        all_results[[result_id]] <- one
        plot_rows[[plot_row_id]] <- one
        result_id <- result_id + 1L
        plot_row_id <- plot_row_id + 1L
      }
    }

    one_plot_data <- do.call(rbind, plot_rows)
    alpha_label <- gsub("\\.", "p", alpha_eff)
    plot_file <- file.path(patient_dir, paste0("average_regular_efficacy_alpha_", alpha_label, ".png"))
    grDevices::png(plot_file, width = 1600, height = 950, res = 130)
    old_par <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

    for (scenario_name in names(BOIN12_EFFICACY_SCENARIOS)) {
      scenario_data <- one_plot_data[one_plot_data$scenario == scenario_name, ]
      truth <- scenario_data[scenario_data$prior_label == prior_settings$prior_label[1L], ]
      truth <- truth[order(truth$dose), ]
      graphics::plot(
        truth$dose,
        truth$true_regular_probability,
        type = "b",
        pch = 16,
        lwd = 2,
        ylim = c(0, 1),
        xlab = "Dose level",
        ylab = "Regular efficacy probability",
        main = scenario_name
      )

      for (prior_index in seq_len(nrow(prior_settings))) {
        method_data <- scenario_data[
          scenario_data$prior_label == prior_settings$prior_label[prior_index],
        ]
        method_data <- method_data[order(method_data$dose), ]
        graphics::lines(
          method_data$dose,
          method_data$average_regular_estimate,
          type = "b",
          pch = 1,
          lwd = 2,
          lty = 2,
          col = prior_colours[prior_index]
        )
      }

      if (scenario_name == names(BOIN12_EFFICACY_SCENARIOS)[1L]) {
        graphics::legend(
          "topleft",
          bty = "n",
          cex = 0.68,
          legend = c("True regular efficacy", prior_settings$prior_label),
          col = c("black", prior_colours),
          lty = c(1, rep(2, nrow(prior_settings))),
          lwd = 2,
          pch = c(16, rep(1, nrow(prior_settings)))
        )
      }
    }

    graphics::mtext(
      paste0(
        "Average regular-efficacy estimates, alpha = ", alpha_eff,
        ", ", patient_setting$patient_label
      ),
      outer = TRUE,
      cex = 1.15,
      font = 2
    )
    graphics::par(old_par)
    grDevices::dev.off()
    plot_files <- c(plot_files, plot_file)
  }
}

average_results <- do.call(rbind, all_results)
average_results <- average_results[
  order(
    average_results$n_regular,
    average_results$alpha_eff,
    average_results$scenario,
    average_results$a_r,
    average_results$dose
  ),
]
rownames(average_results) <- NULL

utils::write.csv(
  average_results,
  file.path(output_dir, "boin12_efficacy_prior_sensitivity_average_estimates.csv"),
  row.names = FALSE
)

plot_files
}
