## ============================================================
## Random-scenario gap comparison plots for Presentation 7-20-2026
##
## Produces one nine-panel figure per alpha value.  The original six
## operational-characteristic metrics are supplemented with MCSE and the
## mean and standard deviation of paired MTD-selection differences.  CFO is
## included only for the 0.15 target-gap results, where it is available.
## ============================================================

input_dir <- "Presentation 7-20-2026"
output_dir <- file.path(input_dir, "Random Scenario")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  "0.05" = "randomsce_targetgap0p05_table_summary.csv",
  "0.10" = "randomsce_targetgap0p10_table_summary.csv",
  "0.15" = "randomsce_targetgap0p15_table_summary.csv"
)

gap_levels <- names(input_files)
alpha_values <- c(0, 0.3, 0.6, 0.9)

crm_method_labels <- c(
  r_fixed = "CRM",
  alpha_crm = "Alpha-CRM",
  cumu_crm = "IPCRM",
  random = "Discount CRM"
)
cfo_label <- "CFO (PRIDE; gap 0.15)"
oracle_label <- "Oracle (CRM, alpha = 0)"

method_colours <- c(
  "CRM" = "#F8766D",
  "Alpha-CRM" = "#A3A500",
  "IPCRM" = "#00BA38",
  "Discount CRM" = "#00B0F6",
  "CFO (PRIDE; gap 0.15)" = "#7B3294",
  "Oracle (CRM, alpha = 0)" = "#E76BF3"
)

metric_specs <- data.frame(
  Metric = c(
    "Average MTD selection %",
    "Average MTD allocation",
    "Average overdose selection %",
    "Average overdose allocation",
    "Average total unique patients",
    "Average trial duration",
    "MCSE of MTD selection %",
    "Mean paired MTD selection difference vs oracle %",
    "SD paired MTD selection difference vs oracle %"
  ),
  Panel = LETTERS[1:9],
  Title = c(
    "MTD selection",
    "MTD allocation",
    "Overdose selection",
    "Overdose allocation",
    "Sample size",
    "Average trial duration",
    "MTD-selection MCSE",
    "Paired MTD-selection mean difference",
    "Paired MTD-selection SD"
  ),
  Y_Label = c(
    "Percentage (%)",
    "Percentage (%)",
    "Percentage (%)",
    "Percentage (%)",
    "Patients",
    "Time (in weeks)",
    "Percentage points",
    "Percentage points",
    "Percentage points"
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
  "Model", "CRM_r_model", "CFO_Method", "Alpha_true", "Metric", "Total",
  "Duration", "Nmax_eff", "Gap", "n_scenarios_found",
  "n_scenarios_paired_with_oracle"
)
missing_columns <- setdiff(required_columns, names(summary_data))
if (length(missing_columns) > 0L) {
  stop("The summary data are missing columns: ", paste(missing_columns, collapse = ", "))
}

summary_data$Alpha_true <- as.numeric(summary_data$Alpha_true)
summary_data$Nmax_eff <- as.numeric(summary_data$Nmax_eff)
summary_data$Total <- as.numeric(summary_data$Total)
summary_data$Duration <- as.numeric(summary_data$Duration)

cfo_rows <- summary_data[summary_data$Model == "CFO", , drop = FALSE]
if (nrow(cfo_rows) == 0L) {
  stop("CFO results were expected in the 0.15-gap summary, but none were found.")
}
if (any(cfo_rows$Gap != "0.15")) {
  stop("CFO results must only appear in the 0.15-gap summary.")
}
if (any(cfo_rows$CFO_Method != "pride")) {
  stop("Only the CFO PRIDE results are supported by this plotting script.")
}

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

expected_combinations <- function(alpha_value) {
  expected <- expand.grid(
    Gap = gap_levels,
    Display_Method = unname(crm_method_labels),
    Metric = metric_specs$Metric,
    stringsAsFactors = FALSE
  )

  cfo_expected <- expand.grid(
    Gap = "0.15",
    Display_Method = cfo_label,
    Metric = metric_specs$Metric,
    stringsAsFactors = FALSE
  )
  expected <- rbind(expected, cfo_expected)

  if (alpha_value > 0) {
    oracle_expected <- expand.grid(
      Gap = gap_levels,
      Display_Method = oracle_label,
      Metric = metric_specs$Metric,
      stringsAsFactors = FALSE
    )
    expected <- rbind(expected, oracle_expected)
  }
  expected
}

validate_plot_data <- function(data, alpha_value) {
  expected <- expected_combinations(alpha_value)
  observed_key <- paste(
    as.character(data$Gap), as.character(data$Display_Method),
    as.character(data$Metric), sep = "__"
  )
  expected_key <- paste(
    expected$Gap, expected$Display_Method, expected$Metric, sep = "__"
  )

  unexpected_keys <- setdiff(observed_key, expected_key)
  if (length(unexpected_keys) > 0L) {
    stop("Unexpected plotting combinations at alpha = ", alpha_value, ": ",
         paste(unexpected_keys, collapse = "; "))
  }

  key_count <- table(factor(observed_key, levels = expected_key))
  if (any(key_count != 1L)) {
    bad_keys <- names(key_count)[key_count != 1L]
    stop(
      "Expected exactly one value for every required alpha/gap/method/metric ",
      "combination at alpha = ", alpha_value, ". Problem keys: ",
      paste(bad_keys, collapse = "; ")
    )
  }
  if (any(!is.finite(data$Raw_Value)) || any(!is.finite(data$Plot_Value))) {
    stop("Non-finite value found in the extracted plot data for alpha = ", alpha_value)
  }
}

extract_alpha_data <- function(alpha_value) {
  crm <- summary_data[
    summary_data$Model == "CRM" &
      summary_data$CRM_r_model %in% names(crm_method_labels) &
      abs(summary_data$Alpha_true - alpha_value) < 1e-12 &
      summary_data$Metric %in% metric_specs$Metric,
    , drop = FALSE
  ]
  crm$Display_Method <- unname(crm_method_labels[crm$CRM_r_model])
  crm$Source_Alpha <- crm$Alpha_true

  cfo <- summary_data[
    summary_data$Model == "CFO" &
      summary_data$CFO_Method == "pride" &
      summary_data$Gap == "0.15" &
      abs(summary_data$Alpha_true - alpha_value) < 1e-12 &
      summary_data$Metric %in% metric_specs$Metric,
    , drop = FALSE
  ]
  cfo$Display_Method <- cfo_label
  cfo$Source_Alpha <- cfo$Alpha_true

  core <- rbind(crm, cfo)
  if (alpha_value > 0) {
    oracle <- summary_data[
      summary_data$Model == "CRM" &
        summary_data$CRM_r_model == "r_fixed" &
        abs(summary_data$Alpha_true) < 1e-12 &
        summary_data$Metric %in% metric_specs$Metric,
      , drop = FALSE
    ]
    oracle$Display_Method <- oracle_label
    oracle$Source_Alpha <- oracle$Alpha_true
    core <- rbind(core, oracle)
  }

  display_methods <- c(unname(crm_method_labels), cfo_label)
  if (alpha_value > 0) {
    display_methods <- c(display_methods, oracle_label)
  }

  core$Plot_Alpha <- alpha_value
  core$Gap <- factor(core$Gap, levels = gap_levels)
  core <- make_plot_value(core)
  core$Display_Method <- factor(core$Display_Method, levels = display_methods)
  core$Metric <- factor(core$Metric, levels = metric_specs$Metric)
  core <- core[order(core$Metric, core$Gap, core$Display_Method), , drop = FALSE]

  validate_plot_data(core, alpha_value)
  core
}

methods_for_gap <- function(gap, display_methods) {
  if (gap == "0.15") display_methods else setdiff(display_methods, cfo_label)
}

pretty_y_limits <- function(metric, values) {
  max_value <- max(values, na.rm = TRUE)

  if (metric == "Mean paired MTD selection difference vs oracle %") {
    limit <- max(5, ceiling(max(abs(values), na.rm = TRUE) / 5) * 5)
    return(c(-limit, limit))
  }
  if (metric == "MCSE of MTD selection %") {
    return(c(0, max(0.06, ceiling(max_value / 0.01) * 0.01)))
  }
  if (metric == "SD paired MTD selection difference vs oracle %") {
    return(c(0, max(15, ceiling(max_value / 5) * 5)))
  }
  if (metric == "Average MTD selection %") return(c(0, max(70, ceiling(max_value / 10) * 10)))
  if (metric == "Average MTD allocation") return(c(0, max(50, ceiling(max_value / 10) * 10)))
  if (metric == "Average overdose selection %") return(c(0, max(35, ceiling(max_value / 5) * 5)))
  if (metric == "Average overdose allocation") return(c(0, max(30, ceiling(max_value / 5) * 5)))
  if (metric == "Average total unique patients") return(c(0, max(30, ceiling(max_value / 5) * 5)))
  if (metric == "Average trial duration") return(c(0, max(60, ceiling(max_value / 5) * 5)))

  c(0, ceiling(max_value))
}

format_y_ticks <- function(metric, ticks) {
  if (metric == "MCSE of MTD selection %") {
    return(formatC(ticks, format = "f", digits = 2))
  }
  formatC(ticks, format = "fg", flag = "#", digits = 3)
}

plot_one_panel <- function(data, metric_spec, display_methods) {
  panel_data <- data[data$Metric == metric_spec$Metric, , drop = FALSE]
  values <- panel_data$Plot_Value
  ylim <- pretty_y_limits(metric_spec$Metric, values)
  y_ticks <- pretty(ylim, n = 5)
  y_ticks <- y_ticks[y_ticks >= ylim[1] & y_ticks <= ylim[2]]

  plot(
    NA,
    xlim = c(0.4, length(gap_levels) + 0.6),
    ylim = ylim,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  abline(h = y_ticks, col = "white", lwd = 1.1)
  if (ylim[1] < 0) abline(h = 0, col = "#4D4D4D", lwd = 0.8)

  for (gap_index in seq_along(gap_levels)) {
    gap <- gap_levels[gap_index]
    gap_methods <- methods_for_gap(gap, display_methods)
    gap_rows <- panel_data[
      as.character(panel_data$Gap) == gap &
        as.character(panel_data$Display_Method) %in% gap_methods,
      , drop = FALSE
    ]
    gap_rows <- gap_rows[match(gap_methods, as.character(gap_rows$Display_Method)), , drop = FALSE]
    bar_width <- 0.78 / length(gap_methods)
    bar_midpoints <- gap_index +
      (seq_along(gap_methods) - (length(gap_methods) + 1) / 2) * bar_width

    rect(
      xleft = bar_midpoints - bar_width / 2,
      ybottom = pmin(0, gap_rows$Plot_Value),
      xright = bar_midpoints + bar_width / 2,
      ytop = pmax(0, gap_rows$Plot_Value),
      col = unname(method_colours[gap_methods]),
      border = NA
    )
  }

  axis(1, at = seq_along(gap_levels), labels = gap_levels, tck = -0.02)
  axis(
    2, at = y_ticks, labels = format_y_ticks(metric_spec$Metric, y_ticks),
    las = 1, tck = -0.02
  )
  box(col = "white")
  mtext(metric_spec$Y_Label, side = 2, line = 2.7, cex = 0.86)
  title(
    main = paste0("(", metric_spec$Panel, ") ", metric_spec$Title),
    cex.main = 0.95,
    font.main = 1
  )
}

draw_figure <- function(alpha_value, plot_data) {
  display_methods <- levels(plot_data$Display_Method)
  layout(
    matrix(c(1:9, rep(10, 3)), nrow = 4, byrow = TRUE),
    widths = c(1, 1, 1),
    heights = c(1, 1, 1, 0.18)
  )
  par(oma = c(0, 0, 2.1, 0), bg = "#EBEBEB", fg = "#4D4D4D")

  for (ii in seq_len(nrow(metric_specs))) {
    par(
      mar = c(3.4, 4.0, 2.5, 0.55),
      mgp = c(2.15, 0.55, 0),
      tcl = -0.25,
      col.axis = "#4D4D4D",
      col.lab = "#333333",
      col.main = "#111111"
    )
    plot_one_panel(plot_data, metric_specs[ii, ], display_methods)
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center",
    legend = display_methods,
    fill = unname(method_colours[display_methods]),
    border = NA,
    bty = "n",
    horiz = TRUE,
    cex = if (alpha_value == 0) 0.88 else 0.78,
    x.intersp = 0.8
  )

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
  "Plot_Alpha", "Source_Alpha", "Gap", "Display_Method", "Model",
  "CRM_r_model", "CFO_Method", "Metric", "Raw_Value", "Plot_Value",
  "Nmax_eff", "n_scenarios_found", "n_scenarios_paired_with_oracle"
)]
write.csv(
  extracted_data,
  file.path(output_dir, "random_sce_gap_comparison_extracted_data.csv"),
  row.names = FALSE
)

all_alphas_pdf <- file.path(output_dir, "random_sce_gap_comparison_all_alphas.pdf")
pdf(all_alphas_pdf, width = 18, height = 14, onefile = TRUE, useDingbats = FALSE)
for (alpha_value in alpha_values) {
  draw_figure(alpha_value, plot_data_by_alpha[[alpha_tag(alpha_value)]])
}
dev.off()

for (alpha_value in alpha_values) {
  file_stub <- paste0("random_sce_gap_comparison_alpha", alpha_tag(alpha_value))
  plot_data <- plot_data_by_alpha[[alpha_tag(alpha_value)]]

  pdf(
    file.path(output_dir, paste0(file_stub, ".pdf")),
    width = 18,
    height = 14,
    onefile = FALSE,
    useDingbats = FALSE
  )
  draw_figure(alpha_value, plot_data)
  dev.off()

  png(
    file.path(output_dir, paste0(file_stub, ".png")),
    width = 3600,
    height = 2800,
    res = 200
  )
  draw_figure(alpha_value, plot_data)
  dev.off()
}

cat("Wrote extracted data and plots to:", normalizePath(output_dir, winslash = "/"), "\n")
