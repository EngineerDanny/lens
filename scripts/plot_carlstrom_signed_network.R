#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(igraph)
})

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), ".."), mustWork = FALSE)
if (!file.exists(file.path(root, "results", "final_sparse_pr_predictions_80pct.csv"))) {
  root <- normalizePath(".")
}

predictions <- read.csv(
  file.path(root, "results", "final_sparse_pr_predictions_80pct.csv"),
  stringsAsFactors = FALSE
)
features <- read.csv(
  file.path(root, "analysis_data", "carlstrom_phyllosphere_2019_pair_features.csv"),
  stringsAsFactors = FALSE
)
predictions <- predictions[predictions$analysis_set == "carlstrom_phyllosphere_2019", ]
truth <- read.csv(file.path(root, "cleaned_data", "carlstrom_phyllosphere_2019_tested_pairs.csv"))
pair_key <- function(d) paste(pmin(d$taxon_1, d$taxon_2), pmax(d$taxon_1, d$taxon_2), sep = "||")
stopifnot(nrow(predictions) == 989, !anyDuplicated(pair_key(predictions)), !anyDuplicated(pair_key(truth)))
ix <- match(pair_key(predictions), pair_key(truth))
stopifnot(!anyNA(ix), all(predictions$interaction_label == truth$interaction_label[ix]))
predictions$sign_label <- ifelse(predictions$interaction_label == 0, "neutral", truth$effect_sign[ix])
stopifnot(all(predictions$sign_label %in% c("positive", "negative", "neutral")))
fi <- match(pair_key(predictions), pair_key(features))
stopifnot(!anyNA(fi))
predictions$prevalence_1 <- features$prevalence_1[fi]
predictions$prevalence_2 <- features$prevalence_2[fi]
predictions$interaction_probability <- predictions$supervised_score
predictions$display_selected <- FALSE
predictions$onenet_display_selected <- FALSE
fold_audit <- list()
for (fold_id in sort(unique(predictions$fold))) {
  test <- which(predictions$fold == fold_id)
  train <- which(predictions$fold != fold_id)
  stopifnot(all(predictions$budget_pairs[test] == length(train)))
  k <- max(1L, as.integer(round(mean(predictions$interaction_label[train]) * length(test))))
  # Deterministic tie breaking by pair identifier; test outcomes never enter ranking.
  a <- order(-predictions$supervised_score[test], predictions$pair_id[test])
  b <- order(-predictions$onenet_score[test], predictions$pair_id[test])
  predictions$display_selected[test[a[seq_len(k)]]] <- TRUE
  predictions$onenet_display_selected[test[b[seq_len(k)]]] <- TRUE
  fold_audit[[as.character(fold_id)]] <- data.frame(fold = fold_id, train_pairs = length(train),
    test_pairs = length(test), train_interactions = sum(predictions$interaction_label[train]), selected_edges = k)
}
write.csv(predictions, file.path(root, "analysis_data", "carlstrom_binary_lens_network_predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, fold_audit), file.path(root, "results", "carlstrom_binary_network_fold_audit.csv"), row.names = FALSE)
metrics <- do.call(rbind, lapply(c("LENS", "OneNet"), function(method) {
  selected <- if (method == "LENS") predictions$display_selected else predictions$onenet_display_selected
  hit <- sum(predictions$interaction_label[selected])
  data.frame(method, selected_edges = sum(selected),
    selected_positive_effects = sum(predictions$sign_label[selected] == "positive"),
    selected_negative_effects = sum(predictions$sign_label[selected] == "negative"),
    selected_neutrals = sum(predictions$sign_label[selected] == "neutral"),
    precision = hit / sum(selected), recall = hit / sum(predictions$interaction_label))
}))
write.csv(metrics, file.path(root, "results", "carlstrom_binary_network_edge_summary.csv"), row.names = FALSE)
print(metrics)

canonical <- function(first, second) {
  data.frame(
    taxon_1 = pmin(as.character(first), as.character(second)),
    taxon_2 = pmax(as.character(first), as.character(second)),
    stringsAsFactors = FALSE
  )
}

ordered <- canonical(features$taxon_1, features$taxon_2)
features$taxon_1 <- ordered$taxon_1
features$taxon_2 <- ordered$taxon_2
features$layout_score <- rowMeans(
  features[, c("pln_score_percentile", "glmnet_score_percentile", "sparcc_score_percentile")],
  na.rm = TRUE
)
layout_edges <- features[
  rank(-features$layout_score, ties.method = "first") <= ceiling(0.10 * nrow(features)),
  c("taxon_1", "taxon_2")
]

taxa <- sort(unique(c(features$taxon_1, features$taxon_2)))
graph <- graph_from_data_frame(layout_edges, directed = FALSE, vertices = taxa)
set.seed(20260821)
coordinates <- layout_with_fr(graph, niter = 2000, grid = "nogrid")
nodes <- data.frame(
  taxon = V(graph)$name,
  x = coordinates[, 1],
  y = coordinates[, 2],
  stringsAsFactors = FALSE
)
code_index <- seq_along(taxa) - 1L
taxon_key <- data.frame(
  code = paste0(
    LETTERS[code_index %/% 26L + 1L],
    LETTERS[code_index %% 26L + 1L]
  ),
  taxon = taxa,
  stringsAsFactors = FALSE
)
nodes <- merge(nodes, taxon_key, by = "taxon", all.x = TRUE)

prevalence <- rbind(
  data.frame(taxon = predictions$taxon_1, prevalence = predictions$prevalence_1),
  data.frame(taxon = predictions$taxon_2, prevalence = predictions$prevalence_2)
)
prevalence <- aggregate(prevalence ~ taxon, prevalence, mean, na.rm = TRUE)
nodes <- merge(nodes, prevalence, by = "taxon", all.x = TRUE)
nodes$prevalence[is.na(nodes$prevalence)] <- 0

observed <- predictions[predictions$sign_label != "neutral", ]
observed$panel <- "A  Experimental effects"
observed$edge_sign <- observed$sign_label
observed$edge_alpha <- 0.78

predicted <- predictions[tolower(as.character(predictions$display_selected)) == "true", ]
predicted$panel <- "B  Held-out LENS selected edges"
predicted$edge_sign <- ifelse(
  predicted$sign_label == "neutral", "neutral", predicted$sign_label
)
predicted$edge_alpha <- 0.25 + 0.65 * predicted$interaction_probability

onenet <- predictions[
  tolower(as.character(predictions$onenet_display_selected)) == "true",
]
onenet$panel <- "C  OneNet selected edges"
onenet$edge_sign <- ifelse(
  onenet$sign_label == "neutral", "neutral", onenet$sign_label
)
onenet$edge_alpha <- 0.72

edges <- rbind(
  observed[, c("taxon_1", "taxon_2", "panel", "edge_sign", "edge_alpha")],
  predicted[, c("taxon_1", "taxon_2", "panel", "edge_sign", "edge_alpha")],
  onenet[, c("taxon_1", "taxon_2", "panel", "edge_sign", "edge_alpha")]
)
edges <- merge(edges, nodes[, c("taxon", "x", "y")], by.x = "taxon_1", by.y = "taxon")
names(edges)[names(edges) %in% c("x", "y")] <- c("x", "y")
edges <- merge(edges, nodes[, c("taxon", "x", "y")], by.x = "taxon_2", by.y = "taxon")
names(edges)[names(edges) %in% c("x.x", "y.x", "x.y", "y.y")] <- c("x", "y", "xend", "yend")

panel_levels <- c(
  "A  Experimental effects",
  "B  Held-out LENS selected edges",
  "C  OneNet selected edges"
)
edges$panel <- factor(edges$panel, levels = panel_levels)
node_panels <- merge(nodes, data.frame(panel = factor(panel_levels, levels = panel_levels)))

plot <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend,
        colour = edge_sign, linetype = edge_sign, alpha = edge_alpha),
    linewidth = 0.55,
    lineend = "round"
  ) +
  geom_point(
    data = node_panels,
    aes(x = x, y = y, size = prevalence),
    shape = 21,
    fill = "white",
    colour = "grey20",
    stroke = 0.45
  ) +
  geom_text(
    data = node_panels,
    aes(x = x, y = y, label = code),
    size = 1.65,
    colour = "grey15",
    fontface = "bold"
  ) +
  facet_wrap(~panel, nrow = 1) +
  coord_equal(clip = "off") +
  scale_colour_manual(
    values = c(negative = "#FF3B30", neutral = "grey65", positive = "#243BFF"),
    breaks = c("positive", "negative", "neutral"),
    labels = c("Positive effect", "Negative effect", "Tested neutral"),
    name = "Experimental outcome"
  ) +
  scale_linetype_manual(
    values = c(negative = "dashed", neutral = "dotted", positive = "solid"),
    breaks = c("positive", "negative", "neutral"),
    labels = c("Positive effect", "Negative effect", "Tested neutral"),
    name = "Experimental outcome"
  ) +
  scale_alpha_identity() +
  scale_size_continuous(range = c(3.2, 5.2), guide = "none") +
  theme_void(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing = grid::unit(1.3, "lines"),
    plot.margin = margin(6, 12, 6, 6),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA)
  )

output <- file.path(root, "paper", "figures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(
  taxon_key,
  file.path(root, "results", "carlstrom_taxon_code_key.csv"),
  row.names = FALSE
)

key_columns <- 5L
key_rows <- ceiling(nrow(taxon_key) / key_columns)
key_entries <- sprintf("%s = %s", taxon_key$code, taxon_key$taxon)
length(key_entries) <- key_rows * key_columns
key_matrix <- matrix(key_entries, nrow = key_rows, ncol = key_columns)
key_matrix[is.na(key_matrix)] <- ""
key_lines <- apply(
  key_matrix,
  1,
  function(values) paste(sprintf("%-18s", values), collapse = "")
)
key_title <- grid::textGrob(
  "Taxon codes",
  x = 0.015,
  hjust = 0,
  gp = grid::gpar(fontface = "bold", fontsize = 8)
)
key_grob <- grid::textGrob(
  paste(key_lines, collapse = "\n"),
  x = 0.015,
  hjust = 0,
  gp = grid::gpar(fontfamily = "mono", fontsize = 5.0, lineheight = 0.88)
)
combined <- gridExtra::arrangeGrob(
  plot,
  key_title,
  key_grob,
  ncol = 1,
  heights = c(1, 0.055, 0.30)
)
ggsave(
  file.path(output, "carlstrom_signed_network.png"),
  combined, width = 7.2, height = 4.85, units = "in", dpi = 320,
  bg = "white"
)

cat(sprintf(
  "Observed edges: %d; displayed LENS edges: %d; displayed OneNet edges: %d\n",
  nrow(observed), nrow(predicted), nrow(onenet)
))
