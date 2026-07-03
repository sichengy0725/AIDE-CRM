## ============================================================
## Selected S20_8d scenarios: toxicity curves and OC plots
## ============================================================

target <- 0.30

input_file <- file.path(
  "Results",
  "All_AIDE_OC_S20_8d_modelsCRM_CFO_crmfixed_random_alpha_crm_cumu_crm_cfoempirical_pride_target0p3_w28_c3_cyc2_3_rate0p07_Nmax30_45_60_dosecap3_cont1_tried1_jobs1to2000_dose_summary.csv"
)

out_dir <- file.path("Results", "S20_8d_selected_mtd_scenarios")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

selected_input <- data.frame(
  List_Order = seq_len(7L),
  Original_Scenario = c(12L, 3L, 7L, 9L, 1L, 14L, 16L),
  Given_MTD = c(3L, 4L, 4L, 4L, 5L, 5L, 6L)
)

selected_map <- selected_input[order(selected_input$Given_MTD, selected_input$List_Order), ]
selected_map$New_Scenario <- seq_len(nrow(selected_map))
selected_map <- selected_map[
  ,
  c("New_Scenario", "Original_Scenario", "Given_MTD", "List_Order")
]
rownames(selected_map) <- NULL

## Confirm the requested MTD labels match the data file.
scenario_mtd <- unique(dat[, c("Scenario", "True_MTD")])
selected_map$File_True_MTD <- scenario_mtd$True_MTD[
  match(selected_map$Original_Scenario, scenario_mtd$Scenario)
]
if (any(selected_map$Given_MTD != selected_map$File_True_MTD)) {
  stop("At least one requested MTD does not match True_MTD in the summary file.")
}

write.csv(
  selected_map,
  file.path(out_dir, "selected_scenario_relabeling.csv"),
  row.names = FALSE
)

alpha_values <- sort(unique(dat$Alpha_true))
cycle_values <- c(2L, 3L)
nmax_values <- sort(unique(dat$Nmax_eff))

design_specs <- data.frame(
  key = c("crm_fixed", "crm_random", "crm_alpha", "crm_ip", "cfo_pride", "oracle"),
  label = c(
    "AIDE-CRM fixed, r = 0",
    "AIDE-CRM random",
    "AIDE-alpha-CRM",
    "AIDE-IP-CRM",
    "CFO PRIDE",
    "Oracle CRM"
  ),
  color = c("#1f77b4", "#e377c2", "#2ca02c", "#ff7f0e", "#9467bd", "#17becf"),
  pch = c(16, 18, 17, 15, 8, 4),
  lty = c(1, 1, 1, 1, 4, 3),
  lwd = c(2, 2, 2, 2, 2, 2),
  stringsAsFactors = FALSE
)

alpha_tag <- function(alpha_value) {
  tag <- format(alpha_value, trim = TRUE, scientific = FALSE)
  paste0("alpha", gsub("\\.", "p", tag))
}

first_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  x[[1]]
}

make_selected_rows <- function(rows) {
  out <- merge(
    selected_map,
    rows,
    by.x = "Original_Scenario",
    by.y = "Scenario",
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  out[order(out$New_Scenario, out$Dose), ]
}

label_crm_design <- function(r_model) {
  out <- rep(NA_character_, length(r_model))
  out[r_model == "fixed"] <- "crm_fixed"
  out[r_model == "random"] <- "crm_random"
  out[r_model == "alpha_crm"] <- "crm_alpha"
  out[r_model == "cumu_crm"] <- "crm_ip"
  out
}

## -------------------------------
## Toxicity curve plot
## -------------------------------

tox_source <- unique(dat[, c("Scenario", "Dose", "True_DLT_rate")])
tox_rows <- make_selected_rows(tox_source)

tox_pdf <- file.path(out_dir, "S20_8d_selected_scenarios_toxicity_curves.pdf")
pdf(tox_pdf, width = 9.5, height = 6.5, onefile = FALSE, useDingbats = FALSE)
par(mar = c(4.3, 4.5, 3.0, 1.0), xpd = NA)
plot(
  c(1, 8), c(0, max(1, max(tox_rows$True_DLT_rate, na.rm = TRUE))),
  type = "n",
  xlab = "Dose Level",
  ylab = "True DLT probability",
  main = "Selected S20_8d Toxicity Curves",
  axes = FALSE
)
axis(1, at = 1:8)
axis(2)
box()
abline(h = target, col = "#8c8c8c", lty = 2, lwd = 1.5)
grid(nx = NA, ny = NULL, col = "#e6e6e6")

curve_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")
curve_lty <- c(1, 1, 1, 1, 1, 3, 1)
for (ii in seq_len(nrow(selected_map))) {
  d <- tox_rows[tox_rows$New_Scenario == ii, ]
  d <- d[order(d$Dose), ]
  lines(d$Dose, d$True_DLT_rate, type = "b", col = curve_cols[ii], pch = ii, lty = curve_lty[ii], lwd = 2)
}
legend(
  "topleft",
  legend = paste0(
    "Scenario ", selected_map$New_Scenario,
    " (orig ", selected_map$Original_Scenario,
    ", MTD ", selected_map$Given_MTD, ")"
  ),
  col = curve_cols,
  pch = seq_len(nrow(selected_map)),
  lty = curve_lty,
  lwd = 2,
  bty = "n",
  cex = 0.82
)
legend("bottomright", legend = "Target = 0.30", col = "#8c8c8c", lty = 2, lwd = 1.5, bty = "n")
dev.off()

## -------------------------------
## OC plot helpers
## -------------------------------

make_plot_rows <- function(alpha_value, ncycle, nmax_eff) {
  crm_use <- subset(
    dat,
    Model == "CRM" &
      abs(Alpha_true - alpha_value) < 1e-12 &
      Cycle_Max == ncycle &
      Nmax_eff == nmax_eff &
      CRM_r_model %in% c("fixed", "random", "alpha_crm", "cumu_crm")
  )
  crm_use$Design <- label_crm_design(crm_use$CRM_r_model)
  crm_use$Is_Oracle <- FALSE

  cfo_pride <- subset(
    dat,
    Model == "CFO" &
      Method == "cfo_pride" &
      CFO_Method == "pride" &
      abs(Alpha_true - alpha_value) < 1e-12 &
      Cycle_Max == ncycle &
      Nmax_eff == nmax_eff
  )
  cfo_pride$Design <- "cfo_pride"
  cfo_pride$Is_Oracle <- FALSE

  oracle <- subset(
    dat,
    Model == "CRM" &
      CRM_r_model == "fixed" &
      abs(Alpha_true - 0) < 1e-12 &
      Cycle_Max == ncycle &
      Nmax_eff == nmax_eff
  )
  oracle$Design <- "oracle"
  oracle$Is_Oracle <- TRUE

  common_cols <- Reduce(intersect, list(names(crm_use), names(cfo_pride), names(oracle)))
  combined <- rbind(crm_use[common_cols], cfo_pride[common_cols], oracle[common_cols])
  combined <- combined[!is.na(combined$Design), ]
  make_selected_rows(combined)
}

summarize_group <- function(d) {
  d <- d[order(d$Dose), ]

  true_rates <- d$True_DLT_rate
  true_mtd <- first_value(d$Given_MTD)
  design <- first_value(d$Design)
  is_oracle <- isTRUE(first_value(d$Is_Oracle))

  total_treated <- sum(d$Pts_Treated, na.rm = TRUE)
  if (!is.finite(total_treated) || total_treated <= 0) {
    total_treated <- NA_real_
  }

  selected_toxic <- sum(d$MTD_Selection_pct[true_rates > target + 1e-12], na.rm = TRUE)
  treated_toxic <- sum(d$Pts_Treated[true_rates > target + 1e-12], na.rm = TRUE)
  treated_subtherapeutic <- sum(d$Pts_Treated[true_rates < target - 1e-12], na.rm = TRUE)
  correct_selection <- d$MTD_Selection_pct[d$Dose == true_mtd][1]
  treated_mtd <- d$Pts_Treated[d$Dose == true_mtd][1]

  sample_size <- if (is_oracle) total_treated else first_value(d$Total_Unique_Patients)
  duration <- if (is_oracle) NA_real_ else first_value(d$Duration)

  data.frame(
    Scenario = first_value(d$New_Scenario),
    Original_Scenario = first_value(d$Original_Scenario),
    True_MTD = true_mtd,
    Design = design,
    Correct_MTD_Selection = correct_selection,
    Overdose_Selection = selected_toxic,
    Patients_Allocated_MTD = 100 * treated_mtd / total_treated,
    Patients_Overdosed = 100 * treated_toxic / total_treated,
    Patients_Underdosed = 100 * treated_subtherapeutic / total_treated,
    Total_IPDE = sum(d$IPDE_Doses, na.rm = TRUE),
    Sample_Size = sample_size,
    Duration = duration,
    stringsAsFactors = FALSE
  )
}

summarize_plot_data <- function(rows) {
  split_rows <- split(rows, paste(rows$Design, rows$New_Scenario, sep = "__"))
  out <- do.call(rbind, lapply(split_rows, summarize_group))
  out$Design <- factor(out$Design, levels = design_specs$key)
  out <- out[order(out$Design, out$Scenario), ]
  row.names(out) <- NULL
  out
}

metric_limits <- function(plot_data, metric) {
  if (metric %in% c("Correct_MTD_Selection", "Patients_Overdosed", "Patients_Underdosed")) {
    return(c(0, 100))
  }
  if (metric == "Overdose_Selection") {
    return(c(0, max(70, ceiling(max(plot_data[[metric]], na.rm = TRUE) / 10) * 10)))
  }
  if (metric == "Patients_Allocated_MTD") {
    return(c(0, max(80, ceiling(max(plot_data[[metric]], na.rm = TRUE) / 10) * 10)))
  }
  if (metric == "Total_IPDE") {
    max_val <- max(plot_data[[metric]], na.rm = TRUE)
    return(c(0, max(10, ceiling(max_val + 1))))
  }

  vals <- plot_data[[metric]]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(c(0, 1))
  pad <- max(1, diff(range(vals)) * 0.08)
  c(max(0, floor(min(vals) - pad)), ceiling(max(vals) + pad))
}

plot_metric <- function(plot_data, metric, panel_title, ylab, ylim) {
  scenarios <- sort(unique(plot_data$Scenario))
  plot(
    range(scenarios), ylim,
    type = "n",
    axes = FALSE,
    xlab = "Relabeled Scenario",
    ylab = ylab,
    main = panel_title,
    font.main = 2,
    cex.main = 1.12,
    cex.lab = 1.02
  )
  axis(1, at = scenarios, labels = scenarios)
  axis(2)
  box()
  abline(v = scenarios, col = "#d9d9d9", lwd = 0.8)
  abline(h = pretty(ylim), col = "#e6e6e6", lwd = 0.8)

  for (ii in seq_len(nrow(design_specs))) {
    spec <- design_specs[ii, ]
    d <- plot_data[plot_data$Design == spec$key, ]
    d <- d[order(d$Scenario), ]
    y <- d[[metric]]
    x <- d$Scenario
    keep <- is.finite(y)
    if (sum(keep) == 0L) next
    lines(
      x[keep], y[keep],
      type = "b",
      col = spec$color,
      pch = spec$pch,
      lty = spec$lty,
      lwd = spec$lwd,
      cex = 1
    )
  }
}

plot_one_setting <- function(alpha_value, ncycle, nmax_eff) {
  plot_data <- summarize_plot_data(make_plot_rows(alpha_value, ncycle, nmax_eff))

  metrics <- list(
    list("Correct_MTD_Selection", "A. Percentage of correct MTD selection", "Percentage"),
    list("Overdose_Selection", "B. Percentage of overdosing selection", "Percentage"),
    list("Patients_Allocated_MTD", "C. Percentage of patients allocated to the MTD", "Percentage"),
    list("Patients_Overdosed", "D. Percentage of patients overdosed", "Percentage"),
    list("Patients_Underdosed", "E. Percentage of patients underdosed", "Percentage"),
    list("Total_IPDE", "F. Average total number of IPDE", "Average"),
    list("Sample_Size", "G. Average sample size", "Average"),
    list("Duration", "H. Average trial duration", "Days")
  )

  out_file <- file.path(
    out_dir,
    paste0(
      "S20_8d_selected_CRM_OC_alpha_ipcrm_random_cfo_pride_",
      alpha_tag(alpha_value),
      "_ncycle", ncycle,
      "_Nmax", nmax_eff,
      ".pdf"
    )
  )

  pdf(out_file, width = 13.5, height = 10.5, onefile = FALSE, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  layout(
    rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 9), c(10, 10)),
    heights = c(1, 1, 1, 1, 0.34, 0.2)
  )
  par(oma = c(0, 0, 2.8, 0), mar = c(3.4, 4.2, 2.1, 1.4), mgp = c(2.2, 0.7, 0))

  for (metric in metrics) {
    plot_metric(
      plot_data,
      metric[[1]],
      metric[[2]],
      metric[[3]],
      metric_limits(plot_data, metric[[1]])
    )
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center",
    legend = design_specs$label,
    col = design_specs$color,
    pch = design_specs$pch,
    lty = design_specs$lty,
    lwd = design_specs$lwd,
    ncol = 3,
    bty = "n",
    cex = 1.02,
    x.intersp = 0.85,
    y.intersp = 1.15
  )

  plot.new()
  text(
    0.5, 0.78,
    "Relabeled scenarios are ordered by requested true MTD from low to high; ties preserve input order.",
    cex = 0.82
  )
  text(
    0.5, 0.48,
    "Oracle CRM uses fixed CRM with r = 0, carryover alpha = 0, the same ncycle, and the same Nmax.",
    cex = 0.82
  )
  text(
    0.5, 0.18,
    "Patient-overdose percentages use dose-level treated counts above the target divided by total treated counts. Oracle duration is not plotted.",
    cex = 0.82
  )

  mtext(
    paste0(
      "Selected S20_8d Operating Characteristics (alpha = ",
      format(alpha_value, trim = TRUE, scientific = FALSE),
      ", ncycle = ", ncycle,
      ", Nmax = ", nmax_eff,
      ")"
    ),
    outer = TRUE,
    side = 3,
    line = 0.8,
    cex = 1.5,
    font = 2
  )

  invisible(out_file)
}

out_files <- c(tox_pdf)
for (nmax_eff in nmax_values) {
  for (alpha_value in alpha_values) {
    for (ncycle in cycle_values) {
      out_files <- c(out_files, plot_one_setting(alpha_value, ncycle, nmax_eff))
    }
  }
}

cat("Scenario relabeling:\n")
print(selected_map)
cat("\nWrote:\n")
cat(paste0("  ", out_files, collapse = "\n"), "\n")
