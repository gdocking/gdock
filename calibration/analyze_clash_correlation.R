#!/usr/bin/env Rscript
#===========================================================================#
# Analyze correlation between energy weights and clash percentage
#
# Goal: Find weight ranges that minimize clashes in top-ranked models
#===========================================================================#

# Read grid search results
results_file <- "results/grid_search_results.tsv"
if (!file.exists(results_file)) {
  stop("Run grid_search.R first to generate results")
}

df <- read.delim(results_file, header = TRUE, sep = "\t")
cat(sprintf("Loaded %d weight combinations\n\n", nrow(df)))

#===========================================================================#
# Analyze correlation with each weight
#===========================================================================#

cat("=== CORRELATION ANALYSIS ===\n\n")

# Correlation between each weight and mean_top1_clash
cor_vdw <- cor(df$w_vdw, df$mean_top1_clash)
cor_elec <- cor(df$w_elec, df$mean_top1_clash)
cor_desolv <- cor(df$w_desolv, df$mean_top1_clash)

cat(sprintf("Correlation with mean_top1_clash:\n"))
cat(sprintf("  w_vdw:    %+.3f\n", cor_vdw))
cat(sprintf("  w_elec:   %+.3f\n", cor_elec))
cat(sprintf("  w_desolv: %+.3f\n", cor_desolv))
cat("\n")

#===========================================================================#
# Analyze clash % by w_desolv bins
#===========================================================================#

cat("=== CLASH % BY w_desolv ===\n\n")

# Group by w_desolv and calculate mean clash %
desolv_summary <- aggregate(
  cbind(mean_top1_clash, top50, top10, top1) ~ w_desolv,
  data = df,
  FUN = mean
)

# Sort by w_desolv
desolv_summary <- desolv_summary[order(desolv_summary$w_desolv), ]

cat("w_desolv\tmean_clash%\ttop50%\ttop10%\ttop1%\n")
for (i in 1:nrow(desolv_summary)) {
  cat(sprintf("%.2f\t\t%.2f\t\t%.1f\t%.1f\t%.1f\n",
              desolv_summary$w_desolv[i],
              desolv_summary$mean_top1_clash[i],
              desolv_summary$top50[i],
              desolv_summary$top10[i],
              desolv_summary$top1[i]))
}

#===========================================================================#
# Analyze clash % by w_vdw bins
#===========================================================================#

cat("\n=== CLASH % BY w_vdw ===\n\n")

vdw_summary <- aggregate(
  cbind(mean_top1_clash, top50, top10, top1) ~ w_vdw,
  data = df,
  FUN = mean
)
vdw_summary <- vdw_summary[order(vdw_summary$w_vdw), ]

cat("w_vdw\tmean_clash%\ttop50%\ttop10%\ttop1%\n")
for (i in 1:nrow(vdw_summary)) {
  cat(sprintf("%.1f\t%.2f\t\t%.1f\t%.1f\t%.1f\n",
              vdw_summary$w_vdw[i],
              vdw_summary$mean_top1_clash[i],
              vdw_summary$top50[i],
              vdw_summary$top10[i],
              vdw_summary$top1[i]))
}

#===========================================================================#
# Find best trade-off: low clash + high success rate
#===========================================================================#

cat("\n=== BEST TRADE-OFF ANALYSIS ===\n\n")

# Find combinations with lowest clash that still have reasonable success
# Sort by clash first, then by top50 (descending)
df_sorted <- df[order(df$mean_top1_clash, -df$top50), ]

cat("Top 10 combinations with lowest clash %:\n")
cat("w_vdw\tw_elec\tw_desolv\tclash%\ttop50%\ttop10%\ttop1%\n")
for (i in 1:min(10, nrow(df_sorted))) {
  row <- df_sorted[i, ]
  cat(sprintf("%.1f\t%.3f\t%.2f\t\t%.2f\t%.1f\t%.1f\t%.1f\n",
              row$w_vdw, row$w_elec, row$w_desolv,
              row$mean_top1_clash, row$top50, row$top10, row$top1))
}

cat("\n")

# Find the Pareto frontier: best top50 for each clash % bin
cat("Pareto frontier (best top50 at each clash level):\n")
clash_bins <- seq(floor(min(df$mean_top1_clash)), ceiling(max(df$mean_top1_clash)), by = 0.5)
cat("clash_bin\tbest_top50\tw_vdw\tw_elec\tw_desolv\n")

for (bin_start in clash_bins[1:(length(clash_bins)-1)]) {
  bin_end <- bin_start + 0.5
  subset <- df[df$mean_top1_clash >= bin_start & df$mean_top1_clash < bin_end, ]
  if (nrow(subset) > 0) {
    best_idx <- which.max(subset$top50)
    best <- subset[best_idx, ]
    cat(sprintf("%.1f-%.1f\t\t%.1f\t\t%.1f\t%.3f\t%.2f\n",
                bin_start, bin_end, best$top50,
                best$w_vdw, best$w_elec, best$w_desolv))
  }
}

#===========================================================================#
# Summary statistics
#===========================================================================#

cat("\n=== SUMMARY ===\n\n")
cat(sprintf("Clash %% range: %.2f - %.2f\n", min(df$mean_top1_clash), max(df$mean_top1_clash)))
cat(sprintf("Mean clash %%: %.2f\n", mean(df$mean_top1_clash)))
cat(sprintf("\n"))

# Best by top50 (ignoring clash)
best_top50 <- df[which.max(df$top50), ]
cat(sprintf("Best by top50 (ignoring clash):\n"))
cat(sprintf("  Weights: w_vdw=%.1f, w_elec=%.3f, w_desolv=%.2f\n",
            best_top50$w_vdw, best_top50$w_elec, best_top50$w_desolv))
cat(sprintf("  Performance: top50=%.1f%%, top1=%.1f%%, clash=%.2f%%\n",
            best_top50$top50, best_top50$top1, best_top50$mean_top1_clash))
