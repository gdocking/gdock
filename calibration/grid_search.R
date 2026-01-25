#!/usr/bin/env Rscript
#===========================================================================#
#
#  This grid search aims to calibrates gdock energy function weights
#   (w_vdw, w_elec, w_desolv) by finding the combination that best ranks
#   near-native structures (model #1) across a benchmark dataset of
#   protein-protein docking decoys.
#
#  IMPORTANT: Weight combinations are rejected if top-ranked models
#   exceed the clash threshold. This ensures the scoring function
#   produces physically reasonable structures.
#
#===========================================================================#

#===========================================================================#
# Configuration
#===========================================================================#
MAX_CLASH_PCT <- 1.0  # Maximum allowed clash % for top-ranked models

#===========================================================================#
# Grid search function
#===========================================================================#
evaluate_weights <- function(df, complexes, w_vdw, w_elec, w_desolv, max_clash_pct = MAX_CLASH_PCT) {

  # Re-score dataframe with new weights
  df$new_score <- w_vdw * df$vdw + w_elec * df$elec + w_desolv * df$desolv

  # Calculate metrics
  top1_count <- 0
  top5_count <- 0
  top10_count <- 0
  top20_count <- 0
  top50_count <- 0
  ranks <- c()
  clash_violations <- 0
  top1_clash_values <- c()

  for (comp in complexes) {
    # Get complex subset
    comp_data <- df[df$complex == comp, ]

    # Sort by new score
    comp_data <- comp_data[order(comp_data$new_score), ]
    comp_data$rank <- 1:nrow(comp_data)

    # Check clash % of top-ranked model
    top1_clash <- comp_data$clash_pct[1]
    top1_clash_values <- c(top1_clash_values, top1_clash)

    if (top1_clash > max_clash_pct) {
      clash_violations <- clash_violations + 1
    }

    # Find model #1 and where its ranked
    model1_rank <- comp_data$rank[comp_data$model == 1]
    ranks <- c(ranks, model1_rank)

    # Calculate success rate
    if (model1_rank == 1) top1_count <- top1_count + 1
    if (model1_rank <= 5) top5_count <- top5_count + 1
    if (model1_rank <= 10) top10_count <- top10_count + 1
    if (model1_rank <= 20) top20_count <- top20_count + 1
    if (model1_rank <= 50) top50_count <- top50_count + 1
  }

  # Calculate clash violation rate
  clash_violation_pct <- 100 * clash_violations / length(complexes)
  mean_top1_clash <- mean(top1_clash_values)

  return(c(
    top1 = 100 * top1_count / length(complexes),
    top5 = 100 * top5_count / length(complexes),
    top10 = 100 * top10_count / length(complexes),
    top20 = 100 * top20_count / length(complexes),
    top50 = 100 * top50_count / length(complexes),

    # average position of the native structure accross ALL complexes
    #  with a given combination of weights
    # ex: if mean of a combination is 42.3, it means that the
    #  native structure was (on average) on position 42 - ideally
    #  it should be on position 1
    mean_rank = mean(ranks),
    median_rank = median(ranks),

    # Clash metrics
    clash_violation_pct = clash_violation_pct,
    mean_top1_clash = mean_top1_clash
  ))
}

#===========================================================================#
# Read data
#===========================================================================#
data_file <- file.path("results/raw_energies.tsv")
df <- read.delim(data_file, header = TRUE, sep = "\t")

complexes <- unique(df$complex)
cat(sprintf("Loaded %d rows for %d complexes\n", nrow(df), length(complexes)))

#===========================================================================#
# Define grid
#===========================================================================#
w_vdw_values <- seq(0.1, 5.0, by = 0.1)
w_elec_values <- seq(0.01, 0.3, by = 0.005)
w_desolv_values <- seq(0.1, 3.0, by = 0.05)

weights_grid <- expand.grid(
  w_vdw = w_vdw_values,
  w_elec = w_elec_values,
  w_desolv = w_desolv_values
)

n_combinations <- nrow(weights_grid)

#===========================================================================#
# Run the search
#===========================================================================#

# Setup parallelism
library(parallel)
n_cores <- detectCores() - 1 # leave 1 core free
cl <- makeCluster(n_cores)

# Export necessary data and functions to all worker processes
clusterExport(cl, c("df", "complexes", "evaluate_weights", "weights_grid", "MAX_CLASH_PCT"))

# Evaluate combinations
cat(sprintf("Testing %d weight combinations with %d cores, this might take a while...\n\n", n_combinations, n_cores))
results_list <- parLapply(cl, 1:n_combinations, function(i) {

  w_vdw <- weights_grid$w_vdw[i]
  w_elec <- weights_grid$w_elec[i]
  w_desolv <- weights_grid$w_desolv[i]

  metrics <- evaluate_weights(df, complexes, w_vdw, w_elec, w_desolv)

  # return a row with the weights and the metrics
  c(w_vdw, w_elec, w_desolv, metrics)
})

stopCluster(cl)

#===========================================================================#
# Analyze
#===========================================================================#

results <- do.call(rbind, results_list)
colnames(results) <- c("w_vdw", "w_elec", "w_desolv", "top1", "top5",
                       "top10", "top20", "top50", "mean_rank", "median_rank",
                       "clash_violation_pct", "mean_top1_clash")
results_df <- as.data.frame(results)

# Report clash statistics
cat(sprintf("\n=== CLASH ANALYSIS ===\n"))
cat(sprintf("Max allowed clash %%: %.1f%%\n", MAX_CLASH_PCT))
n_valid <- sum(results_df$clash_violation_pct == 0)
cat(sprintf("Weight combinations with 0%% clash violations: %d / %d (%.1f%%)\n",
            n_valid, nrow(results_df), 100 * n_valid / nrow(results_df)))

# Filter to only include weight combinations with acceptable clash rates
valid_results <- results_df[results_df$clash_violation_pct == 0, ]

if (nrow(valid_results) == 0) {
  cat("\nWARNING: No weight combinations satisfy the clash constraint!\n")
  cat("Showing best results without clash filtering (for reference):\n\n")
  valid_results <- results_df
}

# Sort by best performance: top50 (desc), top20 (desc), top10 (desc), top5 (desc), mean_rank (asc)
#  --- here we are optimizing for the best model to be in the top50, if there is a tie, it select the one with
#  the higest top20, so on and so forth
valid_results <- valid_results[order(-valid_results$top50, -valid_results$top20,
                                     -valid_results$top10, -valid_results$top5,
                                     valid_results$mean_rank), ]

# Also sort full results for saving
results_df <- results_df[order(-results_df$top50, -results_df$top20,
                               -results_df$top10, -results_df$top5,
                               results_df$mean_rank, results_df$clash_violation_pct), ]

# Print best result details
cat("\n=== BEST RESULT (with clash constraint) ===\n")
cat(sprintf("Weight combination with highest top50 success AND 0%% clash violations\n\n"))
best <- valid_results[1, ]
cat(sprintf("Weights:\n"))
cat(sprintf("  w_vdw    = %.3f\n", best$w_vdw))
cat(sprintf("  w_elec   = %.3f\n", best$w_elec))
cat(sprintf("  w_desolv = %.3f\n", best$w_desolv))
cat(sprintf("\nPerformance:\n"))
cat(sprintf("  Top1:        %.2f%%\n", best$top1))
cat(sprintf("  Top5:        %.2f%%\n", best$top5))
cat(sprintf("  Top10:       %.2f%%\n", best$top10))
cat(sprintf("  Top20:       %.2f%%\n", best$top20))
cat(sprintf("  Top50:       %.2f%%\n", best$top50))
cat(sprintf("  Mean rank:   %.1f\n", best$mean_rank))
cat(sprintf("  Median rank: %.1f\n", best$median_rank))
cat(sprintf("\nClash metrics:\n"))
cat(sprintf("  Clash violations: %.1f%%\n", best$clash_violation_pct))
cat(sprintf("  Mean top1 clash:  %.2f%%\n", best$mean_top1_clash))

# Save results to file
output_file <- file.path("results/grid_search_results.tsv")
write.table(results_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nFull results saved to: %s\n", output_file))

# Save valid (clash-constrained) results separately
valid_output_file <- file.path("results/grid_search_valid_results.tsv")
write.table(valid_results, valid_output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Valid results (clash-constrained) saved to: %s\n", valid_output_file))

# Save best weights as JSON
best_weights_file <- file.path("results/optimized_weights.json")
best_weights_json <- sprintf('{
  "w_vdw": %.6f,
  "w_elec": %.6f,
  "w_desolv": %.6f,
  "w_air": 1.0,
  "performance": {
    "top1_percent": %.2f,
    "top5_percent": %.2f,
    "top10_percent": %.2f,
    "top20_percent": %.2f,
    "top50_percent": %.2f,
    "mean_rank": %.2f,
    "median_rank": %.2f
  },
  "clash_metrics": {
    "max_clash_pct_threshold": %.2f,
    "clash_violation_pct": %.2f,
    "mean_top1_clash_pct": %.2f
  },
  "method": "grid_search_with_clash_constraint",
  "n_complexes": %d,
  "n_combinations_tested": %d,
  "n_valid_combinations": %d
}', best$w_vdw, best$w_elec, best$w_desolv,
   best$top1, best$top5, best$top10, best$top20, best$top50,
   best$mean_rank, best$median_rank,
   MAX_CLASH_PCT, best$clash_violation_pct, best$mean_top1_clash,
   length(complexes), n_combinations, nrow(valid_results))

writeLines(best_weights_json, best_weights_file)
cat(sprintf("Best weights saved to: %s\n", best_weights_file))
