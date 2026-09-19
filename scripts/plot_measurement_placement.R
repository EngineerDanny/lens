#!/usr/bin/env Rscript
# Three panels with identical test sets for both label placement rules.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]))
root <- dirname(dirname(script))
output <- file.path(root, "results", "measurement_placement")
metrics <- fread(file.path(output, "partition_metrics.csv"))
audit <- fread(file.path(output, "partition_audit.csv"))
names <- c(butyrate_assembly_2021 = "Butyrate", carlstrom_phyllosphere_2019 = "Carlström",
           schafer_phyllosphere_2022 = "Schäfer")
counts <- audit[, .(valid = sum(status == "ok"), total = .N), by = analysis_set]
labels <- setNames(sprintf("%s (%d/%d partitions)", names[counts$analysis_set], counts$valid, counts$total), counts$analysis_set)
methods <- c("LENS distributed", "LENS restricted", "OneNet")
plot_data <- metrics[method %in% methods & status == "ok"]
stopifnot(!anyDuplicated(plot_data[, .(analysis_set, `repeat`, budget_pairs, method)]),
          all(is.finite(plot_data$auprc)))
plot_data[, panel := factor(labels[analysis_set], levels = labels[names(names)])]
plot_data[, method := factor(method, levels = methods)]
plot_data[, budget := factor(budget_pairs, levels = sort(unique(budget_pairs)))]
colours <- c("LENS distributed" = "#3333FF", "LENS restricted" = "#D97706", "OneNet" = "#555555")
plot <- ggplot(plot_data, aes(budget, auprc, colour = method, fill = method)) +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"),
               width = 0.75, alpha = 0.18, linewidth = 0.45,
               outlier.shape = 1, outlier.size = 1.1, staplewidth = 0.5) +
  facet_wrap(~panel, nrow = 1, scales = "free") +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(x = "Measured training pairs (n)", y = "AUPRC", colour = NULL, fill = NULL) +
  theme(legend.position = "right", panel.spacing.x = grid::unit(0.8, "lines"))
fwrite(plot_data, file.path(output, "figure_data.csv"))
ggsave(file.path(output, "measurement_placement.png"), plot, width = 11.4, height = 3.6,
       units = "in", dpi = 300, bg = "white")
message("Wrote results/measurement_placement/measurement_placement.png")
