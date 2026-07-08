## ============================================================
## Plot 5-dose high-gap scenarios from Set_5dose_adaptive_r_37
## ============================================================

target <- 0.30

out_dir <- file.path("Results", "Set_5dose_adaptive_r_37_gap_scenarios")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

scenario_meta <- data.frame(
  Scenario = 26:37,
  Scenario_Group = c(
    rep("higher_gap_above_MTD", 6L),
    rep("higher_gap_below_MTD", 6L)
  ),
  True_MTD = c(3L, 3L, 2L, 2L, 4L, 4L, 4L, 3L, 3L, 3L, 2L, 2L),
  Dose1 = c(0.15, 0.11, 0.20, 0.19, 0.05, 0.04, 0.02, 0.01, 0.05, 0.03, 0.13, 0.08),
  Dose2 = c(0.22, 0.17, 0.30, 0.30, 0.08, 0.06, 0.05, 0.02, 0.13, 0.09, 0.29, 0.30),
  Dose3 = c(0.30, 0.30, 0.50, 0.56, 0.17, 0.10, 0.11, 0.09, 0.30, 0.30, 0.35, 0.40),
  Dose4 = c(0.62, 0.50, 0.60, 0.74, 0.27, 0.20, 0.30, 0.32, 0.36, 0.40, 0.40, 0.45),
  Dose5 = c(0.80, 0.62, 0.77, 0.86, 0.54, 0.60, 0.36, 0.38, 0.42, 0.45, 0.43, 0.50),
  stringsAsFactors = FALSE
)

plot_gap_group <- function(rows, title, file_base, ylim) {
  dose_cols <- paste0("Dose", 1:5)
  curve_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  curve_pch <- c(16, 17, 15, 18, 3, 7)
  curve_lty <- c(1, 1, 1, 1, 2, 2)
  
  draw_plot <- function() {
    par(mar = c(4.6, 4.8, 3.5, 1.2), xpd = NA)
    plot(
      c(1, 5),
      ylim,
      type = "n",
      axes = FALSE,
      xlab = "Dose level",
      ylab = "True DLT probability",
      main = title
    )
    axis(1, at = 1:5)
    axis(2, las = 1)
    box()
    grid(nx = NA, ny = NULL, col = "#e6e6e6")
    abline(h = target, col = "#666666", lty = 2, lwd = 1.8)
    
    for (i in seq_len(nrow(rows))) {
      y <- as.numeric(rows[i, dose_cols])
      lines(
        1:5,
        y,
        type = "b",
        col = curve_cols[i],
        pch = curve_pch[i],
        lty = curve_lty[i],
        lwd = 2
      )
      mtd <- rows$True_MTD[i]
      points(
        mtd,
        y[mtd],
        pch = 21,
        bg = "white",
        col = curve_cols[i],
        lwd = 2.2,
        cex = 1.45
      )
    }
    
    legend(
      "topleft",
      legend = paste0("Scenario ", rows$Scenario, " (MTD ", rows$True_MTD, ")"),
      col = curve_cols[seq_len(nrow(rows))],
      pch = curve_pch[seq_len(nrow(rows))],
      lty = curve_lty[seq_len(nrow(rows))],
      lwd = 2,
      bty = "n",
      cex = 0.86
    )
    legend(
      "bottomright",
      legend = c("Target = 0.30", "Open marker = true MTD"),
      col = c("#666666", "black"),
      pch = c(NA, 21),
      lty = c(2, NA),
      lwd = c(1.8, NA),
      pt.bg = c(NA, "white"),
      bty = "n",
      cex = 0.86
    )
  }
  
  png_file <- file.path(out_dir, paste0(file_base, ".png"))
  pdf_file <- file.path(out_dir, paste0(file_base, ".pdf"))
  
  png(png_file, width = 1200, height = 800, res = 140)
  draw_plot()
  dev.off()
  
  pdf(pdf_file, width = 9.5, height = 6.5, onefile = FALSE, useDingbats = FALSE)
  draw_plot()
  dev.off()
  
  c(png = png_file, pdf = pdf_file)
}

above_rows <- scenario_meta[scenario_meta$Scenario %in% 26:31, ]
below_rows <- scenario_meta[scenario_meta$Scenario %in% 32:37, ]

out_files <- c(
  plot_gap_group(
    rows = above_rows,
    title = "5-Dose Scenarios 26-31: Larger Gaps Above MTD",
    file_base = "Set_5dose_scenarios_26_to_31_larger_gaps_above_MTD",
    ylim = c(0, 0.9)
  ),
  plot_gap_group(
    rows = below_rows,
    title = "5-Dose Scenarios 32-37: Larger Gaps Below MTD",
    file_base = "Set_5dose_scenarios_32_to_37_larger_gaps_below_MTD",
    ylim = c(0, 0.65)
  )
)

cat("Wrote:\n")
cat(paste0("  ", out_files, collapse = "\n"), "\n")
