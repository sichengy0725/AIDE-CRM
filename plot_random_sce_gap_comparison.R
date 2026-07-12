## ============================================================
## Random-scenario gap comparison plots
##
## Uses the three table summaries for target gaps 0.05, 0.10,
## and 0.15.  It produces one six-panel figure per alpha value.
## Allocation metrics are expressed as a percentage of Nmax, while
## sample size is the recorded average number of unique trial patients.
## ============================================================

input_dir <- file.path("Results", "Random Sce", "Random Sce")
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  "0.05" = "randomsce_targetgap0p05_table_summary.csv",
  "0.10" = "randomsce_targetgap0p10_table_summary.csv",
  "0.15" = "randomsce_targetgap0p15_table_summary.csv"
)

gap_levels <- names(input_files)
alpha_values <- c(0, 0.3, 0.6, 0.9)

method_labels <- c(
  r_fixed = "CRM",
  alpha_crm = "Alpha-CRM",
  cumu_crm = "IPCRM",
  random = "Discount CRM"
)

method_colours <- c(
  "CRM" = "#F8766D",
  "Alpha-CRM" = "#A3A500",
  "IPCRM" = "#00BA38",
  "Discount CRM" = "#00B0F6",
  "Oracle (CRM, alpha = 0)" = "#E76BF3"
)

metric_specs <- data.frame(
  Metric = c(
    "Average MTD selection %",
    "Average MTD allocation",
    "Average overdose selection %",
    "Average overdose allocation",
    "Average total unique patients",
    "Average trial duration"
  ),
  Panel = LETTERS[1:6],
  Title = c(
    "MTD selection",
    "MTD allocation",
    "Overdose selection",
    "Overdose allocation",
    "Sample size",
    "Average trial duration"
  ),
  Y_Label = c(
    "Percentage (%)",
    "Percentage (%)",
    "Percentage (%)",
    "Percentage (%)",
    "Patients",
    "Time (in weeks)"
  ),
  stringsAsFactors = FALSE
)

read_summary <- function(gap, filename) {
  path <- file.path(input_dir, filename)
  if (!file.exists(path)) {
    stop("Required summary file is missing: ", path)
  }

  out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  out$Gap <- gap
  out
}

summary_data <- do.call(
  rbind,
  Map(read_summary, names(input_files), unname(input_files))
)

required_columns <- c(
  "Model", "CRM_r_model", "Alpha_true", "Metric", "Total", "Duration",
  "Nmax_eff", "Gap"
)
missing_columns <- setdiff(required_columns, names(summary_data))
if (length(missing_columns) > 0L) {
  stop("The summary data are missing columns: ", paste(missing_columns, collapse = ", "))
}

summary_data$Alpha_true <- as.numeric(summary_data$Alpha_true)
summary_data$Nmax_eff <- as.numeric(summary_data$Nmax_eff)
summary_data$Total <- as.numeric(summary_data$Total)
summary_data$Duration <- as.numeric(summary_data$Duration)

make_plot_value <- function(data) {
  data$Raw_Value <- ifelse(
    data$Metric == "Average trial duration",
    data$Duration,
    data$Total
  )
  data$Plot_Value <- data$Raw_Value

  allocation_metrics <- c("Average MTD allocation", "Average overdose allocation")
  allocation_rows <- data$Metric %in% allocation_metrics
  data$Plot_Value[allocation_rows] <-
    100 * data$Raw_Value[allocation_rows] / data$Nmax_eff[allocation_rows]

  duration_rows <- data$Metric == "Average trial duration"
  data$Plot_Value[duration_rows] <- data$Raw_Value[duration_rows] / 7
  data
}

validate_plot_data <- function(data, alpha_value, display_methods) {
  expected <- expand.grid(
    Gap = gap_levels,
    Display_Method = display_methods,
    Metric = metric_specs$Metric,
    stringsAsFactors = FALSE
  )

  observed_key <- paste(data$Gap, data$Display_Method, data$Metric, sep = "__")
  expected_key <- paste(expected$Gap, expected$Display_Method, expected$Metric, sep = "__")
  key_count <- table(factor(observed_key, levels = expected_key))

  if (any(key_count != 1L)) {
    bad_keys <- names(key_count)[key_count != 1L]
    stop(
      "Expected exactly one value for every alpha/gap/method/metric combination at alpha = ",
      alpha_value,
      ". Problem keys: ", paste(bad_keys, collapse = "; ")
    )
  }
  if (any(!is.finite(data$Raw_Value)) || any(!is.finite(data$Plot_Value))) {
    stop("Non-finite value found in the extracted plot data for alpha = ", alpha_value)
  }
}

extract_alpha_data <- function(alpha_value) {
  core <- summary_data[
    summary_data$Model == "CRM" &
      summary_data$CRM_r_model %in% names(method_labels) &
      abs(summary_data$Alpha_true - alpha_value) < 1e-12 &
      summary_data$Metric %in% metric_specs$Metric,
    , drop = FALSE
  ]
  core$Display_Method <- unname(method_labels[core$CRM_r_model])
  core$Source_Alpha <- core$Alpha_true

  if (alpha_value > 0) {
    oracle <- summary_data[
      summary_data$Model == "CRM" &
        summary_data$CRM_r_model == "r_fixed" &
        abs(summary_data$Alpha_true) < 1e-12 &
        summary_data$Metric %in% metric_specs$Metric,
      , drop = FALSE
    ]
    oracle$Display_Method <- "Oracle (CRM, alpha = 0)"
    oracle$Source_Alpha <- oracle$Alpha_true
    core <- rbind(core, oracle)
  }

  core$Plot_Alpha <- alpha_value
  core$Gap <- factor(core$Gap, levels = gap_levels)
  core <- make_plot_value(core)

  display_methods <- if (alpha_value == 0) {
    unname(method_labels)
  } else {
    c(unname(method_labels), "Oracle (CRM, alpha = 0)")
  }
  core$Display_Method <- factor(core$Display_Method, levels = display_methods)
  core$Metric <- factor(core$Metric, levels = metric_specs$Metric)
  core <- core[order(core$Metric, core$Gap, core$Display_Method), , drop = FALSE]

  validate_plot_data(core, alpha_value, display_methods)
  core
}

pretty_ymax <- function(metric, values) {
  max_value <- max(values, na.rm = TRUE)

  if (metric == "Average MTD selection %") return(max(70, ceiling(max_value / 10) * 10))
  if (metric == "Average MTD allocation") return(max(50, ceiling(max_value / 10) * 10))
  if (metric == "Average overdose selection %") return(max(35, ceiling(max_value / 5) * 5))
  if (metric == "Average overdose allocation") return(max(30, ceiling(max_value / 5) * 5))
  if (metric == "Average total unique patients") return(max(30, ceiling(max_value / 5) * 5))
  if (metric == "Average trial duration") return(max(60, ceiling(max_value / 5) * 5))

  ceiling(max_value)
}

plot_one_panel <- function(data, metric_spec, display_methods) {
  panel_data <- data[data$Metric == metric_spec$Metric, , drop = FALSE]
  values <- sapply(gap_levels, function(gap) {
    gap_data <- panel_data[panel_data$Gap == gap, , drop = FALSE]
    gap_data$Plot_Value[match(display_methods, as.character(gap_data$Display_Method))]
  })
  rownames(values) <- display_methods

  ymax <- pretty_ymax(metric_spec$Metric, values)
  bar_midpoints <- barplot(values, beside = TRUE, plot = FALSE)
  xlim <- range(bar_midpoints) + c(-0.8, 0.8)
  y_ticks <- pretty(c(0, ymax), n = 5)
  y_ticks <- y_ticks[y_ticks >= 0 & y_ticks <= ymax]

  plot(
    NA,
    xlim = xlim,
    ylim = c(0, ymax),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  abline(h = y_ticks, col = "white", lwd = 1.1)
  barplot(
    values,
    beside = TRUE,
    add = TRUE,
    axes = FALSE,
    axisnames = FALSE,
    col = unname(method_colours[display_methods]),
    border = NA
  )
  axis(1, at = colMeans(bar_midpoints), labels = gap_levels, tck = -0.02)
  axis(2, at = y_ticks, labels = y_ticks, las = 1, tck = -0.02)
  box(col = "white")
  mtext(metric_spec$Y_Label, side = 2, line = 2.8, cex = 0.95)
  title(
    main = paste0("(", metric_spec$Panel, ") ", metric_spec$Title),
    cex.main = 1.05,
    font.main = 1
  )
}

draw_figure <- function(alpha_value, plot_data) {
  display_methods <- levels(plot_data$Display_Method)
  layout(
    matrix(1:12, nrow = 3, byrow = TRUE),
    widths = c(3.6, 1.55, 3.6, 1.55),
    heights = c(1, 1, 1)
  )
  par(oma = c(0, 0, 2.1, 0), bg = "#EBEBEB", fg = "#4D4D4D")

  for (ii in seq_len(nrow(metric_specs))) {
    metric_spec <- metric_specs[ii, ]
    par(
      mar = c(3.6, 4.0, 2.5, 0.5),
      mgp = c(2.2, 0.55, 0),
      tcl = -0.25,
      col.axis = "#4D4D4D",
      col.lab = "#333333",
      col.main = "#111111"
    )
    plot_one_panel(plot_data, metric_spec, display_methods)

    par(mar = c(0, 0, 0, 0), xpd = NA)
    plot.new()
    legend(
      "center",
      legend = display_methods,
      fill = unname(method_colours[display_methods]),
      border = NA,
      bty = "n",
      cex = if (alpha_value == 0) 0.82 else 0.70,
      y.intersp = 0.9
    )
  }

  mtext(
    paste0("Random scenarios (alpha = ", format(alpha_value, trim = TRUE, scientific = FALSE), ")"),
    outer = TRUE,
    side = 3,
    line = 0.35,
    cex = 1.28,
    font = 1
  )
}

alpha_tag <- function(alpha_value) {
  gsub("\\.", "p", format(alpha_value, trim = TRUE, scientific = FALSE))
}

plot_data_by_alpha <- lapply(alpha_values, extract_alpha_data)
names(plot_data_by_alpha) <- vapply(alpha_values, alpha_tag, character(1))

extracted_data <- do.call(rbind, plot_data_by_alpha)
extracted_data <- extracted_data[, c(
  "Plot_Alpha", "Source_Alpha", "Gap", "Display_Method", "CRM_r_model",
  "Metric", "Raw_Value", "Plot_Value", "Nmax_eff", "n_scenarios_found"
)]
write.csv(
  extracted_data,
  file.path(output_dir, "random_sce_gap_comparison_extracted_data.csv"),
  row.names = FALSE
)

all_alphas_pdf <- file.path(output_dir, "random_sce_gap_comparison_all_alphas.pdf")
pdf(all_alphas_pdf, width = 16, height = 10, onefile = TRUE, useDingbats = FALSE)
for (alpha_value in alpha_values) {
  draw_figure(alpha_value, plot_data_by_alpha[[alpha_tag(alpha_value)]])
}
dev.off()

for (alpha_value in alpha_values) {
  file_stub <- paste0("random_sce_gap_comparison_alpha", alpha_tag(alpha_value))
  plot_data <- plot_data_by_alpha[[alpha_tag(alpha_value)]]

  pdf(
    file.path(output_dir, paste0(file_stub, ".pdf")),
    width = 16,
    height = 10,
    onefile = FALSE,
    useDingbats = FALSE
  )
  draw_figure(alpha_value, plot_data)
  dev.off()

  png(
    file.path(output_dir, paste0(file_stub, ".png")),
    width = 3200,
    height = 2000,
    res = 200
  )
  draw_figure(alpha_value, plot_data)
  dev.off()
}

cat("Wrote extracted data and plots to:", normalizePath(output_dir, winslash = "/"), "\n")
