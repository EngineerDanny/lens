#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
})

base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth/cdiff_inhibition_2024", mustWork = TRUE)
out_dir <- file.path(base_dir, "processed")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

parse_group <- function(x) {
  parts <- unlist(strsplit(as.character(x), "[-+]"))
  strain <- intersect(parts, c("DSM", "MS001", "MS008", "MS014"))
  species <- setdiff(parts, c("DSM", "MS001", "MS008", "MS014", "a", "b", "c"))
  species <- species[species != ""]
  list(
    strain = if (length(strain)) strain[[1]] else NA_character_,
    species = species,
    community_id = paste(sort(species), collapse = "+")
  )
}

allowed_species <- c("BU", "CA", "CS", "CH", "DP")

# Direct abundance table used for network fitting.
s15c <- read.csv(file.path(base_dir, "Source Data", "Figure S15C", "gLV_validation_processed.csv"), check.names = FALSE, stringsAsFactors = FALSE)
s15c_meta <- lapply(s15c$Group, parse_group)
s15c$strain <- vapply(s15c_meta, `[[`, "", "strain")
s15c$community_id <- vapply(s15c_meta, `[[`, "", "community_id")
s15c$species_csv <- vapply(s15c_meta, function(x) paste(x$species, collapse = ","), "")

abundance <- s15c[, c("Sample", "Group", "strain", "community_id", "species_csv", "BU", "CA", "CS", "CH", "DP", "CD")]
names(abundance)[1:2] <- c("sample_id", "raw_group")
write_tsv_gz(abundance, file.path(out_dir, "validation_abundance.tsv.gz"))

# Combined abundance table for focused incoming-edge recovery.
fig3 <- as.data.frame(read_excel(file.path(base_dir, "Source Data", "Figure 3A", "EXP0028_gLV_validation.xlsx")))
fig3_meta <- lapply(fig3$Group, parse_group)
fig3$strain <- vapply(fig3_meta, `[[`, "", "strain")
fig3$community_id <- vapply(fig3_meta, `[[`, "", "community_id")
fig3$species_csv <- vapply(fig3_meta, function(x) paste(x$species, collapse = ","), "")
for (nm in c("BU", "CA", "CS", "CH", "DP", "CD")) {
  if (!nm %in% names(fig3)) fig3[[nm]] <- 0
}
fig3$source_dataset <- "Figure_3A"
s15c$source_dataset <- "Figure_S15C"
incoming_abundance <- rbind(
  s15c[, c("Sample", "Group", "strain", "community_id", "species_csv", "source_dataset", "BU", "CA", "CS", "CH", "DP", "CD")],
  fig3[, c("Sample", "Group", "strain", "community_id", "species_csv", "source_dataset", "BU", "CA", "CS", "CH", "DP", "CD")]
)
names(incoming_abundance)[1:2] <- c("sample_id", "raw_group")
write_tsv_gz(incoming_abundance, file.path(out_dir, "cdiff_incoming_abundance.tsv.gz"))

# Direct truth table from the same validation experiment.
truth_s15c <- aggregate(`CD abs` ~ strain + community_id + species_csv, data = s15c[!is.na(s15c$strain), ], FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
truth_s15c <- data.frame(
  strain = truth_s15c$strain,
  community_id = truth_s15c$community_id,
  species_csv = truth_s15c$species_csv,
  observed_cd_abs_mean = truth_s15c$`CD abs`[, "mean"],
  observed_cd_abs_sd = truth_s15c$`CD abs`[, "sd"],
  n_rep = truth_s15c$`CD abs`[, "n"],
  source = "Figure_S15C",
  stringsAsFactors = FALSE
)

# Additional experimentally validated stronger-inhibition communities.
s13a <- read.csv(file.path(base_dir, "Source Data", "Figure S13A", "Combined_higherinhibitorythanCDCH.csv"), check.names = FALSE, stringsAsFactors = FALSE)
s13a <- s13a[!(s13a$Community %in% c("a", "b", "c")), , drop = FALSE]
s13a_meta <- lapply(s13a$Community, parse_group)
s13a$strain <- vapply(s13a_meta, `[[`, "", "strain")
s13a$community_id <- vapply(s13a_meta, `[[`, "", "community_id")
s13a$species_csv <- vapply(s13a_meta, function(x) paste(x$species, collapse = ","), "")
s13a$all_allowed <- vapply(s13a_meta, function(x) all(x$species %in% allowed_species), FALSE)
s13a <- s13a[s13a$all_allowed, , drop = FALSE]

truth_s13a <- aggregate(`CD abs at 24h` ~ strain + community_id + species_csv, data = s13a, FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
truth_s13a <- data.frame(
  strain = truth_s13a$strain,
  community_id = truth_s13a$community_id,
  species_csv = truth_s13a$species_csv,
  observed_cd_abs_mean = truth_s13a$`CD abs at 24h`[, "mean"],
  observed_cd_abs_sd = truth_s13a$`CD abs at 24h`[, "sd"],
  n_rep = truth_s13a$`CD abs at 24h`[, "n"],
  source = "Figure_S13A",
  stringsAsFactors = FALSE
)

truth <- rbind(truth_s15c, truth_s13a)
truth <- truth[order(truth$strain, truth$observed_cd_abs_mean, truth$community_id), ]
write_tsv_gz(truth, file.path(out_dir, "cdiff_inhibition_truth.tsv.gz"))

# Focused incoming-edge truth for CD.
candidate_taxa <- c("BU", "CA", "CS", "CH", "DP")
strains <- c("DSM", "MS001", "MS008", "MS014")
incoming_truth <- expand.grid(
  strain = strains,
  source_taxon = candidate_taxa,
  stringsAsFactors = FALSE
)
incoming_truth$target_taxon <- "CD"
incoming_truth$label <- as.integer(incoming_truth$source_taxon %in% c("CS", "CH"))
incoming_truth$label_type <- ifelse(
  incoming_truth$label == 1L,
  "validated_positive",
  "candidate_background"
)
incoming_truth$evidence_source <- ifelse(
  incoming_truth$source_taxon == "CS",
  "Figure_3A_direct_pair",
  ifelse(
    incoming_truth$source_taxon == "CH",
    "Figure_S13A_direct_pair",
    "candidate_background_only"
  )
)
incoming_truth <- incoming_truth[, c(
  "strain", "source_taxon", "target_taxon",
  "label", "label_type", "evidence_source"
)]
write_tsv_gz(incoming_truth, file.path(out_dir, "cdiff_incoming_truth.tsv.gz"))

message("Wrote:")
message(" - ", file.path(out_dir, "validation_abundance.tsv.gz"))
message(" - ", file.path(out_dir, "cdiff_incoming_abundance.tsv.gz"))
message(" - ", file.path(out_dir, "cdiff_inhibition_truth.tsv.gz"))
message(" - ", file.path(out_dir, "cdiff_incoming_truth.tsv.gz"))
