#!/usr/bin/env Rscript
# Run from the repository root: Rscript scripts/plot_lens_feature_correlations.R
# Descriptive only: no fitting, feature selection, or outcome correlations.
library(ggplot2)
systems <- c("butyrate_assembly_2021", "carlstrom_phyllosphere_2019", "schafer_phyllosphere_2022")
names(systems) <- c("Butyrate", "Carlström", "Schäfer")
scores <- c("pln_score_percentile", "glmnet_score_percentile", "sparcc_score_percentile")
direct <- c("direct_prevalence_min", "direct_prevalence_max", "direct_joint_prevalence",
  "direct_presence_jaccard", "presence_phi_abs", "presence_log_odds_abs",
  "spearman_all_abs", "spearman_copresent_abs", "log_ratio_variance",
  "abundance_contrast_min", "abundance_contrast_max")
features <- c(scores, direct)
labels <- c("PLNNetwork", "GLMNet", "SparCC", "Min. prevalence", "Max. prevalence",
  "Joint prevalence", "Presence Jaccard", "Abs. presence phi", "Abs. log odds",
  "Abs. Spearman", "Abs. joint Spearman", "Log ratio variance", "Min. contrast", "Max. contrast")
canonical <- function(d) {
  a <- as.character(d$taxon_1); b <- as.character(d$taxon_2)
  d$taxon_1 <- pmin(a, b); d$taxon_2 <- pmax(a, b); d
}
safe_cor <- function(a, b) {
  if (length(a) < 3 || length(unique(a)) < 2 || length(unique(b)) < 2) return(NA_real_)
  cor(a, b, method = "spearman")
}
contrast <- function(v, g) {
  if (sum(g) < 2 || sum(!g) < 2) return(NA_real_)
  s <- sqrt((var(v[g]) + var(v[!g])) / 2)
  if (!is.finite(s) || s == 0) return(0)
  (mean(v[g]) - mean(v[!g])) / s
}
dir.create("results/lens_feature_correlations", recursive = TRUE, showWarnings = FALSE)
tiles <- list(); audits <- list()
for (title in names(systems)) {
  system <- systems[[title]]
  if (system == systems[[1]]) {
    d <- read.csv("training_data/all_labeled_pairs.csv", check.names = FALSE)
    d <- canonical(d[d$dataset == system & !is.na(d$interaction_label), ])
  } else {
    truth <- canonical(read.csv(paste0("cleaned_data/", system, "_tested_pairs.csv")))
    truth <- truth[!is.na(truth$interaction_label), ]
    f <- canonical(read.csv(paste0("analysis_data/", system, "_pair_features.csv")))
    stopifnot(!anyDuplicated(f[c("taxon_1", "taxon_2")]))
    d <- merge(truth, f, by = c("dataset", "taxon_1", "taxon_2"))
    message(system, ": ", nrow(d), " matched rows of ", nrow(truth),
      " truth rows (same inner join as the LENS loader).")
  }
  stopifnot(!anyDuplicated(d[c("taxon_1", "taxon_2")]))
  x <- read.csv(paste0("cleaned_data/", system, "_abundance.csv"), check.names = FALSE)
  x <- as.matrix(x[-1]); storage.mode(x) <- "double"
  stopifnot(all(is.finite(x)), all(x >= 0), all(c(d$taxon_1, d$taxon_2) %in% colnames(x)))
  pseudocount <- if (any(x > 0)) min(x[x > 0]) / 2 else 0.5
  lx <- log(x + pseudocount)
  z <- t(vapply(seq_len(nrow(d)), function(i) {
    a <- x[, d$taxon_1[i]]; b <- x[, d$taxon_2[i]]
    pa <- a > 0; pb <- b > 0; both <- sum(pa & pb)
    onlya <- sum(pa & !pb); onlyb <- sum(!pa & pb); neither <- sum(!pa & !pb)
    phi <- if (length(unique(pa)) > 1 && length(unique(pb)) > 1) cor(as.numeric(pa), as.numeric(pb)) else 0
    ct <- abs(c(contrast(lx[, d$taxon_1[i]], pb), contrast(lx[, d$taxon_2[i]], pa)))
    c(min(mean(pa), mean(pb)), max(mean(pa), mean(pb)), mean(pa & pb),
      if (sum(pa | pb)) both / sum(pa | pb) else 0, abs(phi),
      abs(log(((both + .5) * (neither + .5)) / ((onlya + .5) * (onlyb + .5)))),
      abs(safe_cor(a, b)), abs(safe_cor(a[pa & pb], b[pa & pb])),
      var(lx[, d$taxon_1[i]] - lx[, d$taxon_2[i]]),
      if (any(is.finite(ct))) min(ct, na.rm = TRUE) else NA_real_,
      if (any(is.finite(ct))) max(ct, na.rm = TRUE) else NA_real_)
  }, numeric(11)))
  colnames(z) <- direct
  inputs <- cbind(d[c("taxon_1", "taxon_2", scores)], z)
  write.csv(inputs, paste0("results/lens_feature_correlations/", system, "_inputs.csv"), row.names = FALSE)
  values <- as.matrix(inputs[features])
  # Pairwise complete observations; no imputation using the full system.
  rho <- suppressWarnings(cor(values, method = "spearman", use = "pairwise.complete.obs"))
  counts <- crossprod(1L * is.finite(values))
  stopifnot(isTRUE(all.equal(rho, t(rho))), all(abs(rho[is.finite(rho)]) <= 1 + 1e-12))
  g <- expand.grid(x = seq_along(features), y = seq_along(features))
  g$rho <- as.vector(rho); g$n_complete <- as.vector(counts)
  g$feature_x <- features[g$x]; g$feature_y <- features[g$y]
  g$system <- title; g$panel <- paste0(title, "\n", nrow(d), " tested pairs")
  tiles[[title]] <- g
  audits[[title]] <- data.frame(system, feature = features,
    missing = colSums(!is.finite(values)), distinct = apply(values, 2, function(v) length(unique(v[is.finite(v)]))))
}
all <- do.call(rbind, tiles)
all$panel <- factor(all$panel, levels = unique(all$panel))
write.csv(all, "results/lens_feature_correlations/correlations.csv", row.names = FALSE)
write.csv(do.call(rbind, audits), "results/lens_feature_correlations/feature_audit.csv", row.names = FALSE)
p <- ggplot(all, aes(x, y, fill = rho)) +
  geom_tile(colour = "white", linewidth = .15) +
  facet_wrap(~panel, nrow = 1) + coord_fixed() +
  scale_x_continuous(breaks = 1:14, labels = labels, expand = c(0, 0)) +
  scale_y_reverse(breaks = 1:14, labels = labels, expand = c(0, 0)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
    limits = c(-1, 1), midpoint = 0, na.value = "grey65", name = "Spearman\ncorrelation") +
  labs(x = NULL, y = NULL, title = "Correlations among the 14 LENS inputs",
    caption = "Tested pairs; no outcome labels used. Grey indicates an undefined correlation.\nCorrelations use available observations for each feature pair; no imputation or model fitting.") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10), strip.text = element_text(size = 12),
    plot.caption = element_text(hjust = 0), legend.position = "right")
ggsave("results/lens_feature_correlations/heatmap.png", p, width = 16, height = 6.8, dpi = 250, bg = "white")
print(do.call(rbind, audits))
