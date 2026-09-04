#!/usr/bin/env Rscript

.libPaths(c(normalizePath("r_library", mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

dataset <- "carlstrom_phyllosphere_2019"
input <- fread(file.path("external_data", dataset, "processed", "removal_counts.csv"), check.names = FALSE)
id_columns <- c("sample_id", "experiment", "treatment", "replicate")
taxa <- setdiff(names(input), id_columns)

directional <- list()
dropouts <- sort(unique(input[treatment != "ALL", treatment]))
for (dropout in dropouts) {
  omitted <- paste0("Leaf", sub("^-", "", dropout))
  subset <- input[treatment %in% c("ALL", dropout)]
  metadata <- data.frame(
    row.names = subset$sample_id,
    experiment = factor(subset$experiment),
    treatment = factor(ifelse(subset$treatment == "ALL", "control", "dropout"), levels = c("control", "dropout"))
  )
  counts <- t(as.matrix(subset[, ..taxa]))
  storage.mode(counts) <- "integer"
  colnames(counts) <- subset$sample_id
  keep <- rowSums(counts) > 0
  dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = metadata, design = ~ experiment + treatment)
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.table(results(dds, contrast = c("treatment", "dropout", "control")), keep.rownames = "target_taxon")
  result[is.na(padj), padj := 1]
  result[, `:=`(
    dataset = dataset,
    source_taxon = omitted,
    interaction_label = as.integer(padj < 0.05),
    tested_status = ifelse(padj < 0.05, "positive", "neutral"),
    effect_sign = fifelse(padj >= 0.05, "neutral", fifelse(log2FoldChange > 0, "positive", "negative")),
    n_evidence = nrow(subset)
  )]
  directional[[dropout]] <- result[target_taxon != omitted]
}

directional <- rbindlist(directional, fill = TRUE)
fwrite(directional, file.path("external_data", dataset, "processed", "truth_directional.csv"), na = "")

directional[, `:=`(
  taxon_1 = pmin(source_taxon, target_taxon),
  taxon_2 = pmax(source_taxon, target_taxon)
)]
undirected <- directional[, .(
  interaction_label = max(interaction_label),
  tested_status = ifelse(max(interaction_label) == 1L, "positive", "neutral"),
  effect_sign = if (all(effect_sign == "neutral")) "neutral" else if (length(unique(effect_sign[effect_sign != "neutral"])) == 1L) unique(effect_sign[effect_sign != "neutral"]) else "conflicting",
  n_evidence = .N,
  min_adjusted_p = min(padj),
  max_abs_log2_fold_change = max(abs(log2FoldChange), na.rm = TRUE),
  truth_type = "community_dropout",
  experimental_setting = "Arabidopsis_phyllosphere",
  label_rule = "DESeq2_BH_adjusted_p_below_0.05",
  mixed_evidence = length(unique(interaction_label)) > 1L
), by = .(dataset, taxon_1, taxon_2)]
setorder(undirected, taxon_1, taxon_2)
fwrite(undirected, file.path("cleaned_data", paste0(dataset, "_tested_pairs.csv")), na = "")
