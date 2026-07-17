## ============================================================
## TITE operating-characteristic plots, matched to the
## Presentation 7-6-2026 OC Ncycle2 metric definitions.
##
## Creates:
##   1. TITE (56-day arrivals) versus non-TITE CRM comparisons
##      for baseline scenarios 1-6.
##   2. TITE-only (14-day arrivals) plots for baseline, larger
##      gaps below MTD (26-31), and larger gaps above MTD (32-37).
##
## Trial duration is intentionally omitted from every figure.
## ============================================================

target <- 0.30
ncycle_use <- 2L
nmax_use <- 30L
restrict_to_target_use <- 0L
alpha_values <- c(0, 0.3, 0.6, 0.9)
comparison_tite_arrival_rates <- c(0.02, 0.04)

tite_file <- file.path(
  "Results", "TITE",
  "All_AIDE_TITE_OC_Set_5dose_adaptive_r_37_modelsCRM_crmfixed_random_alpha_crm_ipcrm_target0p3_w28_c3_cyc2_rate0p02_0p04_0p07_Nmax30_dose_summary.csv"
)
non_tite_file <- file.path(
  "Presentation 7-5",
  "All_AIDE_OC_Set_5dose_adaptive_r_37_modelsCRM_BOIN_CFO_boinapprox1_approx2_restr_fixed_crmr_fixed_random_alpha_crm_cumu_crm_cfoempirical_pride_dose_summary.csv"
)

out_dir <- file.path("Results", "TITE", "Plots", "OC Ncycle2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tite <- read.csv(tite_file, stringsAsFactors = FALSE, check.names = FALSE)
non_tite <- read.csv(non_tite_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c(
  "Scenario", "Scenario_Group", "True_MTD", "Model", "CRM_r_model",
  "Alpha_true", "r_carry", "Accrual", "Nmax_eff", "Cycle_Max",
  "Restrict_To_Target", "Dose", "True_DLT_rate", "MTD_Selection_pct",
  "Pts_Treated", "IPDE_Doses", "Total_Unique_Patients",
  "Early_Stopping_pct"
)
for (source_name in c("tite", "non_tite")) {
  source_data <- get(source_name)
  missing_cols <- setdiff(required_cols, names(source_data))
  if (length(missing_cols) > 0L) {
    stop(
      source_name, " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

scenario_sets <- list(
  baseline = list(
    label = "Baseline Scenarios",
    scenarios = 1:6
  ),
  larger_gaps_below_MTD_s26_31 = list(
    label = "Larger Gaps Below MTD",
    scenarios = 26:31
  ),
  larger_gaps_above_MTD_s32_37 = list(
    label = "Larger Gaps Above MTD",
    scenarios = 32:37
  )
)

comparison_specs <- data.frame(
  key = c(
    "r_fixed_non_tite", "r_fixed_tite",
    "ipcrm_non_tite", "ipcrm_tite",
    "alpha_crm_non_tite", "alpha_crm_tite",
    "random_non_tite", "random_tite"
  ),
  label = c(
    "r_fixed CRM", "r_fixed CRM TITE",
    "IP-CRM (cumu_crm)", "IP-CRM TITE",
    "alpha-CRM", "alpha-CRM TITE",
    "random CRM", "random CRM TITE"
  ),
  color = c(
    "#1f77b4", "#1f77b4",
    "#ff7f0e", "#ff7f0e",
    "#2ca02c", "#2ca02c",
    "#e377c2", "#e377c2"
  ),
  pch = c(16, 1, 15, 0, 17, 2, 18, 5),
  lty = c(1, 2, 1, 2, 1, 2, 1, 2),
  lwd = rep(2, 8),
  stringsAsFactors = FALSE
)

tite_specs <- data.frame(
  key = c("r_fixed", "ipcrm", "alpha_crm", "random", "oracle"),
  label = c(
    "r_fixed CRM TITE", "IP-CRM TITE",
    "alpha-CRM TITE", "random CRM TITE",
    "r_fixed CRM oracle (alpha = 0)"
  ),
  color = c("#1f77b4", "#ff7f0e", "#2ca02c", "#e377c2", "#17becf"),
  pch = c(16, 15, 17, 18, 4),
  lty = c(1, 1, 1, 1, 3),
  lwd = rep(2, 5),
  stringsAsFactors = FALSE
)

alpha_tag <- function(alpha_value) {
  value <- format(alpha_value, trim = TRUE, scientific = FALSE)
  paste0("alpha", gsub("\\.", "p", value))
}

arrival_tag <- function(arrival_rate) {
  value <- format(arrival_rate, trim = TRUE, scientific = FALSE)
  paste0("rate", gsub("\\.", "p", value))
}

arrival_setting_label <- function(arrival_rate) {
  if (abs(arrival_rate - 0.02) < 1e-10) {
    return("56-day arrival setting")
  }
  if (abs(arrival_rate - 0.04) < 1e-10) {
    return("28-day arrival setting")
  }
  if (abs(arrival_rate - 0.0714) < 1e-10) {
    return("14-day arrival setting")
  }
  paste0("arrival-rate setting ", format(arrival_rate, trim = TRUE, scientific = FALSE))
}

first_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  x[[1]]
}

filter_base <- function(data, scenario_ids, alpha_value, arrival_rate) {
  data[
    data$Model == "CRM" &
      data$Scenario %in% scenario_ids &
      data$Nmax_eff == nmax_use &
      data$Cycle_Max == ncycle_use &
      data$Restrict_To_Target == restrict_to_target_use &
      abs(data$Alpha_true - alpha_value) < 1e-12 &
      abs(data$Accrual - arrival_rate) < 1e-10,
    ,
    drop = FALSE
  ]
}

add_design <- function(data, key) {
  data$Design <- key
  data
}

select_comparison_rows <- function(scenario_ids, alpha_value, tite_arrival_rate) {
  tite_base <- filter_base(tite, scenario_ids, alpha_value, arrival_rate = tite_arrival_rate)
  ## The non-TITE source only has the 14-day arrival setting (0.0714).
  non_tite_base <- filter_base(non_tite, scenario_ids, alpha_value, arrival_rate = 0.0714)

  pieces <- list(
    add_design(
      non_tite_base[
        non_tite_base$CRM_r_model == "r_fixed" &
          abs(non_tite_base$r_carry) < 1e-12,
        , drop = FALSE
      ],
      "r_fixed_non_tite"
    ),
    add_design(
      tite_base[tite_base$CRM_r_model == "fixed", , drop = FALSE],
      "r_fixed_tite"
    ),
    add_design(
      non_tite_base[non_tite_base$CRM_r_model == "cumu_crm", , drop = FALSE],
      "ipcrm_non_tite"
    ),
    add_design(
      tite_base[tite_base$CRM_r_model == "ipcrm", , drop = FALSE],
      "ipcrm_tite"
    ),
    add_design(
      non_tite_base[non_tite_base$CRM_r_model == "alpha_crm", , drop = FALSE],
      "alpha_crm_non_tite"
    ),
    add_design(
      tite_base[tite_base$CRM_r_model == "alpha_crm", , drop = FALSE],
      "alpha_crm_tite"
    ),
    add_design(
      non_tite_base[non_tite_base$CRM_r_model == "random", , drop = FALSE],
      "random_non_tite"
    ),
    add_design(
      tite_base[tite_base$CRM_r_model == "random", , drop = FALSE],
      "random_tite"
    )
  )
  do.call(rbind, pieces)
}

select_tite_rows <- function(scenario_ids, alpha_value) {
  ## 0.0714 corresponds to a 14-day inter-arrival time for TITE simulations.
  tite_base <- filter_base(tite, scenario_ids, alpha_value, arrival_rate = 0.0714)
  pieces <- list(
    add_design(tite_base[tite_base$CRM_r_model == "fixed", , drop = FALSE], "r_fixed"),
    add_design(tite_base[tite_base$CRM_r_model == "ipcrm", , drop = FALSE], "ipcrm"),
    add_design(tite_base[tite_base$CRM_r_model == "alpha_crm", , drop = FALSE], "alpha_crm"),
    add_design(tite_base[tite_base$CRM_r_model == "random", , drop = FALSE], "random")
  )

  if (abs(alpha_value) > 1e-12) {
    oracle_base <- filter_base(tite, scenario_ids, alpha_value = 0, arrival_rate = 0.0714)
    pieces <- c(
      pieces,
      list(add_design(oracle_base[oracle_base$CRM_r_model == "fixed", , drop = FALSE], "oracle"))
    )
  }
  do.call(rbind, pieces)
}

validate_plot_rows <- function(rows, design_specs, scenario_ids) {
  expected_designs <- design_specs$key
  found_designs <- unique(rows$Design)
  missing_designs <- setdiff(expected_designs, found_designs)
  if (length(missing_designs) > 0L) {
    stop("Missing requested designs: ", paste(missing_designs, collapse = ", "))
  }

  observed <- table(
    factor(rows$Design, levels = expected_designs),
    factor(rows$Scenario, levels = scenario_ids)
  )
  if (any(observed != 5L)) {
    bad <- which(observed != 5L, arr.ind = TRUE)
    bad_labels <- apply(
      bad,
      1,
      function(index) {
        paste0(
          rownames(observed)[index[[1]]], "/scenario ",
          colnames(observed)[index[[2]]], " (", observed[index[[1]], index[[2]]], " rows)"
        )
      }
    )
    stop(
      "Each design/scenario combination must have five dose rows. Problems: ",
      paste(bad_labels, collapse = "; ")
    )
  }
  invisible(rows)
}

summarize_group <- function(d) {
  d <- d[order(d$Dose), , drop = FALSE]

  true_rates <- d$True_DLT_rate
  total_treated <- sum(d$Pts_Treated, na.rm = TRUE)
  if (!is.finite(total_treated) || total_treated <= 0) {
    total_treated <- NA_real_
  }

  true_mtd <- first_value(d$True_MTD)
  selected_toxic <- sum(d$MTD_Selection_pct[true_rates > target + 1e-12], na.rm = TRUE)
  treated_toxic <- sum(d$Pts_Treated[true_rates > target + 1e-12], na.rm = TRUE)
  treated_subtherapeutic <- sum(d$Pts_Treated[true_rates < target - 1e-12], na.rm = TRUE)

  if (is.na(true_mtd) || true_mtd < 1L || true_mtd > nrow(d)) {
    correct_selection <- first_value(d$Early_Stopping_pct)
    treated_mtd <- 0
  } else {
    correct_selection <- d$MTD_Selection_pct[d$Dose == true_mtd][1]
    treated_mtd <- d$Pts_Treated[d$Dose == true_mtd][1]
  }

  data.frame(
    Scenario = first_value(d$Scenario),
    Scenario_Group = first_value(d$Scenario_Group),
    True_MTD = true_mtd,
    Design = first_value(d$Design),
    Correct_MTD_Selection = correct_selection,
    Overdose_Selection = selected_toxic,
    Patients_Allocated_MTD = 100 * treated_mtd / total_treated,
    Patients_Overdosed = 100 * treated_toxic / total_treated,
    Patients_Underdosed = 100 * treated_subtherapeutic / total_treated,
    Total_IPDE = sum(d$IPDE_Doses, na.rm = TRUE),
    Sample_Size = first_value(d$Total_Unique_Patients),
    stringsAsFactors = FALSE
  )
}

summarize_plot_data <- function(rows, design_specs, scenario_ids) {
  validate_plot_rows(rows, design_specs, scenario_ids)
  split_rows <- split(rows, paste(rows$Design, rows$Scenario, sep = "__"))
  out <- do.call(rbind, lapply(split_rows, summarize_group))
  out$Scenario <- factor(out$Scenario, levels = scenario_ids)
  out$Scenario_Index <- as.integer(out$Scenario)
  out$Design <- factor(out$Design, levels = design_specs$key)
  out <- out[order(out$Design, out$Scenario_Index), , drop = FALSE]
  row.names(out) <- NULL
  out
}

metric_limits <- function(plot_data, metric) {
  values <- plot_data[[metric]]
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(0, 1))
  }
  if (metric %in% c("Correct_MTD_Selection", "Patients_Overdosed", "Patients_Underdosed")) {
    return(c(0, 100))
  }
  if (metric == "Overdose_Selection") {
    return(c(0, max(70, ceiling(max(values) / 10) * 10)))
  }
  if (metric == "Patients_Allocated_MTD") {
    return(c(0, max(80, ceiling(max(values) / 10) * 10)))
  }
  if (metric == "Total_IPDE") {
    return(c(0, max(10, ceiling(max(values) + 1))))
  }
  padding <- max(1, diff(range(values)) * 0.08)
  c(max(0, floor(min(values) - padding)), ceiling(max(values) + padding))
}

plot_metric <- function(plot_data, design_specs, metric, panel_title, ylab, ylim, scenario_ids) {
  xvals <- seq_along(scenario_ids)
  plot(
    range(xvals), ylim,
    type = "n", axes = FALSE,
    xlab = "Scenario", ylab = ylab, main = panel_title,
    font.main = 2, cex.main = 1.12, cex.lab = 1.02
  )
  axis(1, at = xvals, labels = scenario_ids)
  axis(2)
  box()
  abline(v = xvals, col = "#d9d9d9", lwd = 0.8)
  abline(h = pretty(ylim), col = "#e6e6e6", lwd = 0.8)

  for (ii in seq_len(nrow(design_specs))) {
    spec <- design_specs[ii, , drop = FALSE]
    d <- plot_data[plot_data$Design == spec$key, , drop = FALSE]
    d <- d[order(d$Scenario_Index), , drop = FALSE]
    keep <- is.finite(d[[metric]])
    if (sum(keep) == 0L) {
      next
    }
    lines(
      d$Scenario_Index[keep], d[[metric]][keep],
      type = "b", col = spec$color, pch = spec$pch,
      lty = spec$lty, lwd = spec$lwd, cex = 1
    )
  }
}

draw_plot <- function(plot_data, design_specs, set_info, alpha_value, title, source_note) {
  metrics <- list(
    list("Correct_MTD_Selection", "A. Percentage of correct MTD selection", "Percentage"),
    list("Overdose_Selection", "B. Percentage of overdosing selection", "Percentage"),
    list("Patients_Allocated_MTD", "C. Percentage of patients allocated to the MTD", "Percentage"),
    list("Patients_Overdosed", "D. Percentage of patients overdosed", "Percentage"),
    list("Patients_Underdosed", "E. Percentage of patients underdosed", "Percentage"),
    list("Total_IPDE", "F. Average total number of IPDE", "Average"),
    list("Sample_Size", "G. Average sample size", "Average")
  )
  scenario_ids <- set_info$scenarios

  ## Metric G spans both columns; H (trial duration) is intentionally omitted.
  layout(
    rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 7), c(8, 8), c(9, 9)),
    heights = c(1, 1, 1, 1, 0.34, 0.24)
  )
  par(oma = c(0, 0, 3.0, 0), mar = c(3.4, 4.2, 2.1, 1.4), mgp = c(2.2, 0.7, 0))
  for (metric in metrics) {
    metric_name <- metric[[1]]
    plot_metric(
      plot_data, design_specs, metric_name, metric[[2]], metric[[3]],
      metric_limits(plot_data, metric_name), scenario_ids
    )
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center", legend = design_specs$label,
    col = design_specs$color, pch = design_specs$pch,
    lty = design_specs$lty, lwd = design_specs$lwd,
    ncol = if (nrow(design_specs) > 5L) 4L else nrow(design_specs),
    bty = "n", cex = 0.96, x.intersp = 0.80, y.intersp = 1.15
  )

  plot.new()
  text(
    0.5, 0.72,
    paste0("Scenarios shown: ", paste(scenario_ids, collapse = ", "), ". True MTD labels come from the dose-summary CSV."),
    cex = 0.82
  )
  text(0.5, 0.42, source_note, cex = 0.82)
  text(
    0.5, 0.14,
    "Patient-overdose percentages use treated counts above target divided by total treated counts. Trial duration is omitted.",
    cex = 0.82
  )

  heading <- paste0(
    title, " (alpha = ", format(alpha_value, trim = TRUE, scientific = FALSE),
    ", ncycle = ", ncycle_use, ", Nmax = ", nmax_use,
    ", restrict_to_target = ", restrict_to_target_use, ")"
  )
  mtext(
    heading,
    outer = TRUE, side = 3, line = 0.8,
    cex = if (nchar(heading) > 100L) 1.12 else 1.40,
    font = 2
  )
}

write_plot <- function(plot_data, design_specs, set_info, alpha_value, filename, title, source_note) {
  output_file <- file.path(out_dir, filename)
  pdf(output_file, width = 13.5, height = 10.5, onefile = FALSE, useDingbats = FALSE)
  draw_plot(plot_data, design_specs, set_info, alpha_value, title, source_note)
  dev.off()
  output_file
}

out_files <- character(0)

## Part 1: baseline TITE versus non-TITE comparisons.
for (tite_arrival_rate in comparison_tite_arrival_rates) {
  for (alpha_value in alpha_values) {
    set_info <- scenario_sets$baseline
    rows <- select_comparison_rows(set_info$scenarios, alpha_value, tite_arrival_rate)
    plot_data <- summarize_plot_data(rows, comparison_specs, set_info$scenarios)
    out_files <- c(
      out_files,
      write_plot(
        plot_data, comparison_specs, set_info, alpha_value,
        paste0(
          "Set_5dose_baseline_TITE_vs_nonTITE_OC_", alpha_tag(alpha_value),
          "_ncycle2_Nmax30_TITE", arrival_tag(tite_arrival_rate), ".pdf"
        ),
        "Baseline Scenarios: TITE vs non-TITE CRM",
        paste0(
          "TITE uses arrival rate ", format(tite_arrival_rate, trim = TRUE, scientific = FALSE),
          " (", arrival_setting_label(tite_arrival_rate),
          "); non-TITE uses the Presentation 7-5 Nmax = 30 summary (arrival rate 0.0714)."
        )
      )
    )
  }
}

## Part 2: TITE-only 14-day arrival plots, including the alpha = 0 r_fixed oracle for alpha > 0.
for (set_key in names(scenario_sets)) {
  set_info <- scenario_sets[[set_key]]
  for (alpha_value in alpha_values) {
    use_specs <- if (abs(alpha_value) < 1e-12) {
      tite_specs[tite_specs$key != "oracle", , drop = FALSE]
    } else {
      tite_specs
    }
    rows <- select_tite_rows(set_info$scenarios, alpha_value)
    plot_data <- summarize_plot_data(rows, use_specs, set_info$scenarios)
    out_files <- c(
      out_files,
      write_plot(
        plot_data, use_specs, set_info, alpha_value,
        paste0(
          "Set_5dose_", set_key, "_TITE_OC_", alpha_tag(alpha_value),
          "_ncycle2_Nmax30_rate0p0714.pdf"
        ),
        paste0(set_info$label, ": TITE CRM (14-day inter-arrival)"),
        if (abs(alpha_value) < 1e-12) {
          "All four TITE CRM models use arrival rate 0.0714 (14-day inter-arrival)."
        } else {
          "All TITE models use arrival rate 0.0714 (14-day inter-arrival); the oracle is r_fixed CRM with alpha = 0."
        }
      )
    )
  }
}

cat("Wrote ", length(out_files), " PDF files:\n", sep = "")
cat(paste0("  ", out_files), sep = "\n")
cat("\n")
