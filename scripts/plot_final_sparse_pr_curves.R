#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
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
method_columns <- c(
  "LENS" = "supervised_score",
  "OneNet" = "onenet_score",
  "SparCC" = "sparcc_score_percentile",
  "PLNNetwork" = "pln_score_percentile",
  "SPIEC-EASI" = "spieceasi_score_percentile"
)
recall_grid <- seq(0, 1, by = 0.01)

predictions <- fread(file.path(
  root, "results", "final_sparse_pr_predictions_80pct.csv"
))

feature_scores <- rbindlist(lapply(names(system_names), function(system_id) {
  frame <- fread(file.path(
    root, "analysis_data", paste0(system_id, "_pair_features.csv")
  ))
  first <- as.character(frame$taxon_1)
  second <- as.character(frame$taxon_2)
  frame[, `:=`(
    analysis_set = system_id,
    taxon_1 = pmin(first, second),
    taxon_2 = pmax(first, second)
  )]
  frame[, .(
    analysis_set,
    taxon_1,
    taxon_2,
    sparcc_score_percentile,
    pln_score_percentile,
    spieceasi_score_percentile
  )]
}))

first <- as.character(predictions$taxon_1)
second <- as.character(predictions$taxon_2)
predictions[, `:=`(
  taxon_1 = pmin(first, second),
  taxon_2 = pmax(first, second)
)]
predictions <- merge(
  predictions,
  feature_scores,
  by = c("analysis_set", "taxon_1", "taxon_2"),
  all.x = TRUE,
  sort = FALSE
)
stopifnot(
  !anyNA(predictions$sparcc_score_percentile),
  !anyNA(predictions$pln_score_percentile),
  !anyNA(predictions$spieceasi_score_percentile)
)

interpolated_curve <- function(labels, scores) {
  ordering <- order(scores, decreasing = TRUE)
  labels <- labels[ordering]
  true_positive <- cumsum(labels == 1L)
  false_positive <- cumsum(labels == 0L)
  recall <- true_positive / sum(labels == 1L)
  precision <- true_positive / (true_positive + false_positive)
  recall <- c(0, recall)
  precision <- c(1, precision)
  data.table(
    recall = recall_grid,
    precision = vapply(recall_grid, function(value) {
      eligible <- precision[recall >= value]
      if (length(eligible)) max(eligible) else NA_real_
    }, numeric(1))
  )
}

fold_curves <- rbindlist(lapply(names(system_names), function(system_id) {
  system_data <- predictions[analysis_set == system_id]
  rbindlist(lapply(sort(unique(system_data$fold)), function(fold_id) {
    fold_data <- system_data[fold == fold_id]
    rbindlist(lapply(names(method_columns), function(method_name) {
      curve <- interpolated_curve(
        fold_data$interaction_label,
        fold_data[[method_columns[[method_name]]]]
      )
      curve[, `:=`(
        analysis_set = system_id,
        fold = fold_id,
        method = method_name
      )]
      curve
    }))
  }))
}))

curve_summary <- fold_curves[, .(
  mean_precision = mean(precision),
  se_precision = sd(precision) / sqrt(.N),
  folds = .N
), by = .(analysis_set, method, recall)]
curve_summary[, `:=`(
  lower_se = pmax(0, mean_precision - se_precision),
  upper_se = pmin(1, mean_precision + se_precision),
  system = factor(unname(system_names[analysis_set]), levels = unname(system_names)),
  method = factor(method, levels = names(method_columns))
)]

plot <- ggplot(
  curve_summary,
  aes(
    x = recall,
    y = mean_precision,
    colour = method,
    fill = method,
    group = method
  )
) +
  geom_ribbon(
    aes(ymin = lower_se, ymax = upper_se),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~system,
    nrow = 1,
    scales = "fixed",
    axes = "margins",
    axis.labels = "margins"
  ) +
  scale_colour_manual(values = c(
    "LENS" = "#3333FF",
    "OneNet" = "#FF3333",
    "SparCC" = "#404040",
    "PLNNetwork" = "#777777",
    "SPIEC-EASI" = "#B0B0B0"
  )) +
  scale_fill_manual(values = c(
    "LENS" = "#3333FF",
    "OneNet" = "#FF3333",
    "SparCC" = "#404040",
    "PLNNetwork" = "#777777",
    "SPIEC-EASI" = "#B0B0B0"
  )) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Recall",
    y = "Precision",
    colour = "Method",
    fill = "Method"
  ) +
  theme(
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.position = "right"
  )

figure_dir <- file.path(root, "paper", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(
  curve_summary,
  file.path(root, "results", "final_sparse_pr_curve_summary.csv")
)
ggsave(
  file.path(figure_dir, "final_sparse_pr_curves.pdf"),
  plot,
  width = 10.8,
  height = 2.7,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "final_sparse_pr_curves.png"),
  plot,
  width = 10.8,
  height = 2.7,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote final LENS precision-recall curves.")
