#!/usr/bin/env Rscript

base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth/host_fitness_2018", mustWork = TRUE)
proc_dir <- file.path(base_dir, "processed")
if (!dir.exists(proc_dir)) dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

species <- c("LP", "LB", "AP", "AT", "AO")
abund_cols <- paste0(species, " abund (CFU)")

cfu <- read.csv(file.path(base_dir, "FlygutCFUsData.csv"), check.names = FALSE, stringsAsFactors = FALSE)
trt <- read.csv(file.path(base_dir, "TreatmentSummary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
names(trt)[match(c("Lp", "Lb", "Ap", "At", "Ao"), names(trt))] <- species

cfu$sample_id <- sprintf("hostfit_t%02d_%04d", cfu$treatment, seq_len(nrow(cfu)))
abundance <- cfu[, c("sample_id", abund_cols), drop = FALSE]
colnames(abundance) <- c("sample_id", species)
write_tsv_gz(abundance, file.path(proc_dir, "abundance_matrix.tsv.gz"))

taxa_map <- data.frame(
  taxon_id = species,
  source_name = species,
  stringsAsFactors = FALSE
)
write_tsv_gz(taxa_map, file.path(proc_dir, "taxa_map.tsv.gz"))

single_mean <- setNames(rep(NA_real_, length(species)), species)
for (sp in species) {
  sid <- trt$treatment[trt[[sp]] == "Y" & rowSums(trt[, species] == "Y") == 1L]
  rows <- cfu$treatment %in% sid
  vals <- as.numeric(cfu[rows, paste0(sp, " abund (CFU)"), drop = TRUE])
  if (length(vals)) single_mean[sp] <- mean(vals, na.rm = TRUE)
}

eps <- 1
edge_rows <- list()
k <- 1L
pair_rows <- trt[rowSums(trt[, species] == "Y") == 2L, , drop = FALSE]
for (i in seq_len(nrow(pair_rows))) {
  members <- species[pair_rows[i, species] == "Y"]
  pair_id <- pair_rows$treatment[i]
  pair_data <- cfu[cfu$treatment == pair_id, , drop = FALSE]
  for (target in members) {
    source <- setdiff(members, target)
    base <- single_mean[target]
    vals <- as.numeric(pair_data[[paste0(target, " abund (CFU)")]])
    vals <- vals[is.finite(vals)]
    if (!is.finite(base) || !length(vals)) next
    for (v in vals) {
      log2fc <- log2((v + eps) / (base + eps))
      edge_rows[[k]] <- data.frame(
        partner_taxon = source,
        focal_taxon = target,
        effect_sign = ifelse(log2fc > 0, 1L, ifelse(log2fc < 0, -1L, 0L)),
        log2_ratio = log2fc,
        delta_abundance = v - base,
        context = "pair_vs_singleton",
        treatment = pair_id,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}

truth_raw <- do.call(rbind, edge_rows)
agg <- aggregate(
  cbind(log2_ratio, delta_abundance) ~ partner_taxon + focal_taxon,
  data = truth_raw,
  FUN = mean
)
counts <- aggregate(treatment ~ partner_taxon + focal_taxon, data = truth_raw, FUN = length)
colnames(counts)[3] <- "n_evidence"
truth_edges <- merge(agg, counts, by = c("partner_taxon", "focal_taxon"))
truth_edges$effect_sign <- ifelse(
  truth_edges$log2_ratio > 0, 1L,
  ifelse(truth_edges$log2_ratio < 0, -1L, 0L)
)
truth_edges$context <- "pair_vs_singleton_mean"
truth_edges <- truth_edges[abs(truth_edges$log2_ratio) >= 0.5, ]
truth_edges <- truth_edges[order(-abs(truth_edges$log2_ratio), truth_edges$partner_taxon, truth_edges$focal_taxon), ]

write_tsv_gz(truth_raw, file.path(proc_dir, "truth_edges_raw.tsv.gz"))
write_tsv_gz(truth_edges, file.path(proc_dir, "truth_edges.tsv.gz"))

message("Wrote:")
message(" - ", file.path(proc_dir, "abundance_matrix.tsv.gz"))
message(" - ", file.path(proc_dir, "taxa_map.tsv.gz"))
message(" - ", file.path(proc_dir, "truth_edges_raw.tsv.gz"))
message(" - ", file.path(proc_dir, "truth_edges.tsv.gz"))
