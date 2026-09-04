#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- arguments[grep("^--file=", arguments)][1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
root <- dirname(dirname(script_path))

system_names <- c(
  butyrate_assembly_2021 = "Butyrate",
  carlstrom_phyllosphere_2019 = "Carlström",
  schafer_phyllosphere_2022 = "Schäfer"
)

score_columns <- c(
  "PLNNetwork" = "pln_score_percentile",
  "Poisson GLMNet" = "glmnet_score_percentile",
  "SparCC" = "sparcc_score_percentile",
  "SPIEC-EASI" = "spieceasi_score_percentile",
  "OneNet" = "onenet_score"
)

canonicalize <- function(frame) {
  first <- as.character(frame$taxon_1)
  second <- as.character(frame$taxon_2)
  frame$taxon_1 <- pmin(first, second)
  frame$taxon_2 <- pmax(first, second)
  frame
}

average_precision <- function(labels, scores) {
  keep <- !is.na(labels) & !is.na(scores)
  labels <- as.integer(labels[keep])
  scores <- as.numeric(scores[keep])
  positives <- sum(labels == 1L)
  if (positives == 0L) {
    return(NA_real_)
  }

  ordering <- order(scores, decreasing = TRUE)
  labels <- labels[ordering]
  scores <- scores[ordering]
  threshold_end <- c(which(diff(scores) != 0), length(scores))
  cumulative_positive <- cumsum(labels == 1L)[threshold_end]
  cumulative_total <- threshold_end
  recall <- cumulative_positive / positives
  precision <- cumulative_positive / cumulative_total
  sum(c(recall[1], diff(recall)) * precision)
}

precision_recall_curve <- function(labels, scores) {
  keep <- !is.na(labels) & !is.na(scores)
  labels <- as.integer(labels[keep])
  scores <- as.numeric(scores[keep])
  positives <- sum(labels == 1L)
  ordering <- order(scores, decreasing = TRUE)
  labels <- labels[ordering]
  scores <- scores[ordering]
  threshold_end <- c(which(diff(scores) != 0), length(scores))
  cumulative_positive <- cumsum(labels == 1L)[threshold_end]
  data.frame(
    recall = c(0, cumulative_positive / positives),
    precision = c(1, cumulative_positive / threshold_end)
  )
}

load_system <- function(system) {
  truth <- canonicalize(read.csv(
    file.path(root, "cleaned_data", paste0(system, "_tested_pairs.csv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))
  features <- canonicalize(read.csv(
    file.path(root, "analysis_data", paste0(system, "_pair_features.csv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))
  onenet <- canonicalize(read.csv(
    file.path(root, "analysis_data", paste0(system, "_onenet_pair_scores.csv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))

  keys <- c("dataset", "taxon_1", "taxon_2")
  frame <- merge(truth, features, by = keys, all = FALSE)
  frame <- merge(frame, onenet, by = keys, all = FALSE)
  frame <- frame[!is.na(frame$interaction_label), ]
  frame$interaction_label <- as.integer(frame$interaction_label)
  frame
}

rows <- list()
curve_rows <- list()
for (system in names(system_names)) {
  frame <- load_system(system)
  pair_ids <- paste(frame$taxon_1, frame$taxon_2, sep = "||")
  for (method in names(score_columns)) {
    score_column <- unname(score_columns[method])
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = system,
      system = unname(system_names[system]),
      method = method,
      tested_pairs = length(unique(pair_ids)),
      evaluated_outcomes = nrow(frame),
      positives = sum(frame$interaction_label == 1L),
      prevalence = mean(frame$interaction_label == 1L),
      auprc = average_precision(frame$interaction_label, frame[[score_column]])
    )
    curve <- precision_recall_curve(
      frame$interaction_label,
      frame[[score_column]]
    )
    curve$analysis_set <- system
    curve$system <- unname(system_names[system])
    curve$method <- method
    curve_rows[[length(curve_rows) + 1L]] <- curve
  }
}
results <- do.call(rbind, rows)
curves <- do.call(rbind, curve_rows)

results$system <- factor(results$system, levels = unname(system_names))
results$method <- factor(results$method, levels = rev(names(score_columns)))
curves$system <- factor(curves$system, levels = unname(system_names))
curves$method <- factor(curves$method, levels = names(score_columns))
baselines <- unique(results[c("system", "prevalence")])
bootstrap_path <- file.path(root, "results", "abundance_bootstrap_5_auprc.csv")
bootstrap_results <- read.csv(bootstrap_path, stringsAsFactors = FALSE)
bootstrap_results <- bootstrap_results[
  bootstrap_results$fit_status == "complete" &
    !is.na(bootstrap_results$auprc) &
    bootstrap_results$analysis_set %in% names(system_names),
]
bootstrap_results$system <- factor(
  unname(system_names[bootstrap_results$analysis_set]),
  levels = unname(system_names)
)
bootstrap_results$method <- factor(
  bootstrap_results$method,
  levels = rev(names(score_columns))
)
bootstrap_summary <- aggregate(
  auprc ~ analysis_set + system + method,
  data = bootstrap_results,
  FUN = function(values) c(mean = mean(values), sd = sd(values))
)
bootstrap_summary <- data.frame(
  analysis_set = bootstrap_summary$analysis_set,
  system = bootstrap_summary$system,
  method = bootstrap_summary$method,
  mean_auprc = bootstrap_summary$auprc[, "mean"],
  sd_auprc = bootstrap_summary$auprc[, "sd"]
)
onenet_summary <- results[results$method == "OneNet", ]
onenet_summary <- data.frame(
  analysis_set = onenet_summary$analysis_set,
  system = onenet_summary$system,
  method = onenet_summary$method,
  mean_auprc = onenet_summary$auprc,
  sd_auprc = NA_real_
)
bar_data <- rbind(bootstrap_summary, onenet_summary)
bar_data$system <- factor(bar_data$system, levels = unname(system_names))
bar_data$method <- factor(bar_data$method, levels = rev(names(score_columns)))
bar_data$lower_auprc <- pmax(0, bar_data$mean_auprc - bar_data$sd_auprc)
bar_data$upper_auprc <- pmin(1, bar_data$mean_auprc + bar_data$sd_auprc)

bar_plot <- ggplot(bar_data, aes(x = mean_auprc, y = method)) +
  geom_col(width = 0.62, fill = "grey35") +
  geom_vline(
    data = baselines,
    aes(xintercept = prevalence),
    linetype = "11",
    linewidth = 1.0,
    colour = "#B2182B"
  ) +
  geom_errorbar(
    aes(xmin = lower_auprc, xmax = upper_auprc),
    width = 0.22,
    linewidth = 0.55,
    orientation = "y",
    na.rm = TRUE
  ) +
  facet_wrap(~system, nrow = 1, scales = "free_x") +
  scale_x_continuous(
    labels = function(values) sprintf("%.2f", values),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  labs(
    x = "AUPRC across all experimentally tested pairs",
    y = NULL
  ) +
  theme(panel.spacing.x = grid::unit(0.8, "lines"))

curve_plot <- ggplot(curves, aes(x = recall, y = precision, colour = method)) +
  geom_step(linewidth = 0.55, direction = "vh") +
  facet_wrap(~system, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Recall",
    y = "Precision",
    colour = "Method"
  ) +
  theme(legend.position = "right")

figure_dir <- file.path(root, "paper", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  results,
  file.path(root, "results", "unsupervised_full_tested_pair_auprc.csv"),
  row.names = FALSE
)
write.csv(
  bar_data,
  file.path(root, "results", "unsupervised_bootstrap_5_summary.csv"),
  row.names = FALSE
)
ggsave(
  file.path(figure_dir, "unsupervised_full_tested_pair_auprc.pdf"),
  bar_plot,
  width = 10.8,
  height = 2.20,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "unsupervised_full_tested_pair_auprc.png"),
  bar_plot,
  width = 10.8,
  height = 2.20,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(figure_dir, "unsupervised_full_tested_pair_pr_curves.pdf"),
  curve_plot,
  width = 10.8,
  height = 2.20,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "unsupervised_full_tested_pair_pr_curves.png"),
  curve_plot,
  width = 10.8,
  height = 2.20,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote full tested-pair AUPRC and precision--recall figures.")
