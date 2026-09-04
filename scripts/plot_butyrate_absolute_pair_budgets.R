#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- arguments[grep("^--file=", arguments)][1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
root <- dirname(dirname(script_path))

system_files <- c(
  Butyrate = "butyrate_absolute_pair_budgets_repeats.csv",
  `Carlström` = "carlstrom_absolute_pair_budgets_repeats.csv",
  `Schäfer` = "schafer_absolute_pair_budgets_repeats.csv"
)

repeats <- rbindlist(lapply(names(system_files), function(system_name) {
  result <- fread(file.path(root, "results", system_files[[system_name]]))
  result[, system := system_name]
  result
}))

maximum_labels <- c(
  Butyrate = "Maximum\n(83.2)",
  `Carlström` = "Maximum\n(791.2)",
  `Schäfer` = "Maximum\n(1219.2)"
)
repeats[, budget_display := fifelse(
  budget_label == "maximum", unname(maximum_labels[system]), budget_label
)]
display_levels <- c(
  "10", "20", "50", "Maximum\n(83.2)",
  "100", "200", "500", "Maximum\n(791.2)",
  "1000", "Maximum\n(1219.2)"
)
repeats[, budget_display := factor(budget_display, levels = display_levels)]

supervised <- repeats[, .(
  system,
  budget_display,
  method = "LENS",
  auprc = supervised_auprc
)]
onenet <- repeats[, .(
  system,
  budget_display,
  method = "OneNet",
  auprc = onenet_auprc
)]
plot_data <- rbindlist(list(supervised, onenet))
plot_data[, `:=`(
  system = factor(system, levels = names(system_files)),
  method = factor(method, levels = c("LENS", "OneNet"))
)]
summary_data <- plot_data[, .(
  mean_auprc = mean(auprc),
  se_auprc = sd(auprc) / sqrt(.N)
), by = .(system, budget_display, method)]
summary_data[, `:=`(
  lower_auprc = pmax(0, mean_auprc - se_auprc),
  upper_auprc = pmin(1, mean_auprc + se_auprc)
)]

plot <- ggplot(
  summary_data,
  aes(
    x = budget_display,
    y = mean_auprc,
    colour = method,
    fill = method,
    group = method
  )
) +
  geom_ribbon(
    aes(ymin = lower_auprc, ymax = upper_auprc),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~system, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(
    "LENS" = "#3333FF",
    "OneNet" = "#FF3333"
  )) +
  scale_fill_manual(values = c(
    "LENS" = "#3333FF",
    "OneNet" = "#FF3333"
  )) +
  scale_x_discrete(expand = expansion(add = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    x = "Number of measured pair identifiers (n)",
    y = "AUPRC",
    colour = "Method"
  ) +
  guides(fill = "none") +
  theme(
    panel.spacing.x = grid::unit(1.0, "lines"),
    legend.position = "right"
  )

figure_dir <- file.path(root, "paper", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(
  summary_data,
  file.path(root, "results", "absolute_pair_budgets_figure_data.csv")
)
ggsave(
  file.path(figure_dir, "butyrate_absolute_pair_budgets.pdf"),
  plot,
  width = 10.8,
  height = 3.0,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "butyrate_absolute_pair_budgets.png"),
  plot,
  width = 10.8,
  height = 3.0,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote absolute pair budget figure for Butyrate, Carlstrom, and Schafer.")
