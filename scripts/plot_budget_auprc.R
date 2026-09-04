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

results <- fread(file.path(root, "results", "pair_five_fold_cv.csv"))
results <- results[analysis_set %in% names(system_names)]

supervised <- results[method == "Supervised combined" & budget_fraction > 0]
summary <- supervised[, .(
  mean_auprc = mean(auprc),
  se_auprc = sd(auprc) / sqrt(.N),
  minimum_measured_pairs = min(budget_pairs),
  maximum_measured_pairs = max(budget_pairs),
  folds = .N
), by = .(analysis_set, budget_fraction)]

onenet <- results[budget_fraction == 0, .(
  mean_auprc = mean(auprc),
  se_auprc = sd(auprc) / sqrt(.N),
  folds = .N
), by = analysis_set]

summary[, `:=`(
  lower_se = pmax(0, mean_auprc - se_auprc),
  upper_se = pmin(1, mean_auprc + se_auprc)
)]
onenet[, `:=`(
  lower_se = pmax(0, mean_auprc - se_auprc),
  upper_se = pmin(1, mean_auprc + se_auprc)
)]

summary[, system := factor(
  unname(system_names[analysis_set]),
  levels = unname(system_names)
)]
summary[, representative_pairs := as.integer(round(
  (minimum_measured_pairs + maximum_measured_pairs) / 2
))]
summary[, pair_count_label := paste0(
  format(representative_pairs, big.mark = ",", scientific = FALSE, trim = TRUE)
)]
summary[, budget_label := pair_count_label]
label_levels <- summary[order(system, budget_fraction), unique(budget_label)]
summary[, budget_label := factor(budget_label, levels = label_levels)]
summary[, method_label := "Supervised\nModel"]

onenet <- merge(
  onenet,
  summary[, .(analysis_set, budget_fraction, budget_label, system)],
  by = "analysis_set",
  allow.cartesian = TRUE,
  all.x = TRUE
)
onenet[, method_label := "OneNet"]
onenet[, budget_label := factor(budget_label, levels = label_levels)]

plot_data <- rbindlist(list(
  summary[, .(
    analysis_set, system, budget_fraction, budget_label, method_label,
    mean_auprc, se_auprc, lower_se, upper_se, folds
  )],
  onenet[, .(
    analysis_set, system, budget_fraction, budget_label, method_label,
    mean_auprc, se_auprc, lower_se, upper_se, folds
  )]
), use.names = TRUE)

plot_data[, method_label := factor(
  method_label,
  levels = c("Supervised\nModel", "OneNet")
)]

plot <- ggplot(
  plot_data,
  aes(
    x = budget_label,
    y = mean_auprc,
    group = method_label,
    colour = method_label,
    fill = method_label
  )
) +
  geom_ribbon(
    aes(ymin = lower_se, ymax = upper_se),
    alpha = 0.22,
    colour = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~system, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(
    "Supervised\nModel" = "#3333FF",
    "OneNet" = "#FF3333"
  )) +
  scale_fill_manual(values = c(
    "Supervised\nModel" = "#3333FF",
    "OneNet" = "#FF3333"
  )) +
  scale_x_discrete(expand = expansion(add = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    x = "Number of measured pair identifiers (n)",
    y = "AUPRC",
    colour = "Method",
    fill = "Method"
  ) +
  theme(
    panel.spacing.x = grid::unit(0.8, "lines"),
    legend.position = "right"
  )

figure_dir <- file.path(root, "paper", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(root, "results", "budget_auprc_five_fold_figure_summary.csv"))
fwrite(onenet, file.path(root, "results", "budget_auprc_five_fold_references.csv"))
ggsave(
  file.path(figure_dir, "supervised_auprc_by_budget.pdf"),
  plot,
  width = 10.8,
  height = 2.55,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "supervised_auprc_by_budget.png"),
  plot,
  width = 10.8,
  height = 2.55,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote supervised AUPRC budget figure.")
