#!/usr/bin/env Rscript

base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth/butyrate_assembly_2021", mustWork = TRUE)
zip_path <- file.path(base_dir, "DesignSyntheticGutMicrobiomeAssemblyFunction-v1.0.zip")
proc_dir <- file.path(base_dir, "processed")
if (!dir.exists(proc_dir)) dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

species <- c(
  "ER","FP","AC","HB","CC","RI","DP","BH","CA","PC","DL","CG",
  "BF","EL","CH","BO","BT","BU","BV","BC","BY","PJ","DF","BL","BP","BA"
)
fraction_cols <- paste0(species, " Fraction")

tmpdir <- tempfile("butyrate_extract_")
dir.create(tmpdir)
on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

master_rel <- "RyanLincolnClark-DesignSyntheticGutMicrobiomeAssemblyFunction-8fbb777/commonfiles/2020_02_28_MasterDF.csv"
unzip(zip_path, files = master_rel, exdir = tmpdir)
master_path <- file.path(tmpdir, master_rel)
df <- read.csv(master_path, check.names = FALSE, stringsAsFactors = FALSE)

present <- as.matrix(df[, species, drop = FALSE])
present[is.na(present)] <- 0
richness <- rowSums(present > 0)

keep <- df[["Contamination?"]] == "No" & richness %in% c(1, 2)
df_small <- df[keep, c("Treatment", "Rep", "Experiment No.", "Sequenced", species, fraction_cols), drop = FALSE]
df_small$richness <- rowSums(as.matrix(df_small[, species, drop = FALSE]) > 0)
df_small$sample_id <- sprintf(
  "butyrate_exp%s_%s_rep%s_%04d",
  df_small[["Experiment No."]],
  gsub("[^A-Za-z0-9]+", "_", df_small$Treatment),
  df_small$Rep,
  seq_len(nrow(df_small))
)

abundance <- df_small[, c("sample_id", fraction_cols), drop = FALSE]
colnames(abundance) <- c("sample_id", species)
write_tsv_gz(abundance, file.path(proc_dir, "abundance_matrix.tsv.gz"))

taxa_map <- data.frame(
  taxon_id = species,
  source_name = species,
  stringsAsFactors = FALSE
)
write_tsv_gz(taxa_map, file.path(proc_dir, "taxa_map.tsv.gz"))

eps <- 1e-6
edge_rows <- list()
k <- 1L
for (exp_id in sort(unique(df_small[["Experiment No."]]))) {
  de <- df_small[df_small[["Experiment No."]] == exp_id, , drop = FALSE]
  singles <- de[de$richness == 1, , drop = FALSE]
  pairs <- de[de$richness == 2, , drop = FALSE]
  if (!nrow(singles) || !nrow(pairs)) next

  single_mean <- setNames(rep(NA_real_, length(species)), species)
  for (sp in species) {
    frac_col <- paste0(sp, " Fraction")
    x <- singles[singles[[sp]] > 0, frac_col, drop = TRUE]
    if (length(x)) single_mean[sp] <- mean(as.numeric(x), na.rm = TRUE)
  }

  for (i in seq_len(nrow(pairs))) {
    members <- species[as.numeric(pairs[i, species, drop = TRUE]) > 0]
    if (length(members) != 2L) next
    for (target in members) {
      source <- setdiff(members, target)
      base <- single_mean[target]
      pair_value <- as.numeric(pairs[i, paste0(target, " Fraction"), drop = TRUE])
      if (!is.finite(base) || !is.finite(pair_value)) next
      log2fc <- log2((pair_value + eps) / (base + eps))
      delta <- pair_value - base
      edge_rows[[k]] <- data.frame(
        partner_taxon = source,
        focal_taxon = target,
        effect_sign = ifelse(log2fc > 0, 1L, ifelse(log2fc < 0, -1L, 0L)),
        log2_ratio = log2fc,
        delta_fraction = delta,
        context = sprintf("Experiment_%s", exp_id),
        treatment = pairs$Treatment[i],
        sample_id = pairs$sample_id[i],
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}

truth_raw <- do.call(rbind, edge_rows)
agg <- aggregate(
  cbind(log2_ratio, delta_fraction) ~ partner_taxon + focal_taxon,
  data = truth_raw,
  FUN = mean
)
counts <- aggregate(sample_id ~ partner_taxon + focal_taxon, data = truth_raw, FUN = length)
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
