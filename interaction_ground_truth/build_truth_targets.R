#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth", mustWork = TRUE)

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

sort_pair <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

build_undirected <- function(directed) {
  directed$key <- sort_pair(directed$source_taxon, directed$target_taxon)
  split_rows <- split(directed, directed$key)
  rows <- lapply(split_rows, function(df) {
    nodes <- strsplit(df$key[1], "\\|\\|")[[1]]
    signs <- unique(df$effect_sign[df$effect_sign != 0])
    sign_consensus <- if (length(signs) == 1L) signs else NA_integer_
    strength <- suppressWarnings(max(abs(df$strength), na.rm = TRUE))
    if (!is.finite(strength)) strength <- NA_real_
    data.frame(
      taxon_1 = nodes[1],
      taxon_2 = nodes[2],
      n_directions = nrow(unique(df[c("source_taxon", "target_taxon")])),
      sign_consensus = sign_consensus,
      max_abs_strength = strength,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_directed <- function(directed) {
  directed$key <- paste(directed$source_taxon, directed$target_taxon, sep = "||")
  split_rows <- split(directed, directed$key)
  rows <- lapply(split_rows, function(df) {
    signs <- unique(df$effect_sign[df$effect_sign != 0])
    sign_consensus <- if (length(signs) == 1L) signs else NA_integer_
    strength <- suppressWarnings(max(abs(df$strength), na.rm = TRUE))
    if (!is.finite(strength)) strength <- NA_real_
    data.frame(
      source_taxon = df$source_taxon[1],
      target_taxon = df$target_taxon[1],
      effect_sign = sign_consensus,
      max_abs_strength = strength,
      n_evidence = nrow(df),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_pairinterax <- function() {
  proc_dir <- file.path(base_dir, "pairinterax", "processed")
  truth <- read_tsv_gz(file.path(proc_dir, "truth_edges.tsv.gz"))

  dir1 <- truth[truth$phenotype_a != 0, c("canonical_b", "canonical_a", "phenotype_a", "interaction_label")]
  names(dir1) <- c("source_taxon", "target_taxon", "effect_sign", "evidence")
  dir2 <- truth[truth$phenotype_b != 0, c("canonical_a", "canonical_b", "phenotype_b", "interaction_label")]
  names(dir2) <- c("source_taxon", "target_taxon", "effect_sign", "evidence")
  directed <- rbind(dir1, dir2)
  directed$effect_sign <- as.integer(directed$effect_sign)
  directed$strength <- abs(directed$effect_sign)
  directed <- unique(directed)
  directed <- directed[directed$source_taxon != directed$target_taxon, , drop = FALSE]

  directed_std <- build_directed(directed)
  write_tsv_gz(directed_std, file.path(proc_dir, "truth_directed.tsv.gz"))
  write_tsv_gz(build_undirected(directed), file.path(proc_dir, "truth_undirected.tsv.gz"))
}

build_omm12 <- function() {
  proc_dir <- file.path(base_dir, "omm12", "processed")
  truth <- read_tsv_gz(file.path(proc_dir, "truth_edges.tsv.gz"))

  directed <- truth[, c("partner_taxon", "focal_taxon", "effect_sign", "log2_rbm", "coku")]
  names(directed) <- c("source_taxon", "target_taxon", "effect_sign", "strength", "evidence")
  directed$effect_sign <- as.integer(directed$effect_sign)
  directed$strength <- as.numeric(directed$strength)
  directed <- directed[!is.na(directed$effect_sign) & directed$effect_sign != 0, , drop = FALSE]
  directed <- directed[directed$source_taxon != directed$target_taxon, , drop = FALSE]

  directed_std <- build_directed(directed)
  write_tsv_gz(directed_std, file.path(proc_dir, "truth_directed.tsv.gz"))
  write_tsv_gz(build_undirected(directed), file.path(proc_dir, "truth_undirected.tsv.gz"))
}

build_omm12_keystone_2023 <- function() {
  proc_dir <- file.path(base_dir, "omm12_keystone_2023", "processed")

  truth <- read_tsv_gz(file.path(proc_dir, "truth_edges_in_vitro.tsv.gz"))
  directed <- truth[, c("partner_taxon", "focal_taxon", "effect_sign", "log2_ratio", "context")]
  names(directed) <- c("source_taxon", "target_taxon", "effect_sign", "strength", "evidence")
  directed$effect_sign <- as.integer(directed$effect_sign)
  directed$strength <- as.numeric(directed$strength)
  directed <- directed[!is.na(directed$effect_sign) & directed$effect_sign != 0, , drop = FALSE]
  directed <- directed[directed$source_taxon != directed$target_taxon, , drop = FALSE]

  directed_std <- build_directed(directed)
  write_tsv_gz(directed_std, file.path(proc_dir, "truth_directed.tsv.gz"))
  write_tsv_gz(build_undirected(directed), file.path(proc_dir, "truth_undirected.tsv.gz"))
}

build_venturelli <- function() {
  proc_dir <- file.path(base_dir, "venturelli_2018", "processed")
  truth <- read_tsv_gz(file.path(proc_dir, "truth_edges.tsv.gz"))

  directed <- truth[, c("predictor_taxon", "response_taxon", "effect_sign", "coefficient", "model_id")]
  names(directed) <- c("source_taxon", "target_taxon", "effect_sign", "strength", "evidence")
  directed$effect_sign <- as.integer(directed$effect_sign)
  directed$strength <- as.numeric(directed$strength)
  directed <- directed[!is.na(directed$effect_sign) & directed$effect_sign != 0, , drop = FALSE]
  directed <- directed[directed$source_taxon != directed$target_taxon, , drop = FALSE]

  directed_std <- build_directed(directed)
  write_tsv_gz(directed_std, file.path(proc_dir, "truth_directed.tsv.gz"))
  write_tsv_gz(build_undirected(directed), file.path(proc_dir, "truth_undirected.tsv.gz"))
}

build_butyrate_assembly_2021 <- function() {
  proc_dir <- file.path(base_dir, "butyrate_assembly_2021", "processed")
  truth <- read_tsv_gz(file.path(proc_dir, "truth_edges.tsv.gz"))

  directed <- truth[, c("partner_taxon", "focal_taxon", "effect_sign", "log2_ratio", "context")]
  names(directed) <- c("source_taxon", "target_taxon", "effect_sign", "strength", "evidence")
  directed$effect_sign <- as.integer(directed$effect_sign)
  directed$strength <- as.numeric(directed$strength)
  directed <- directed[!is.na(directed$effect_sign) & directed$effect_sign != 0, , drop = FALSE]
  directed <- directed[directed$source_taxon != directed$target_taxon, , drop = FALSE]

  directed_std <- build_directed(directed)
  write_tsv_gz(directed_std, file.path(proc_dir, "truth_directed.tsv.gz"))
  write_tsv_gz(build_undirected(directed), file.path(proc_dir, "truth_undirected.tsv.gz"))
}

build_host_fitness_2018 <- function() {
  proc_dir <- file.path(base_dir, "host_fitness_2018", "processed")
  truth <- read_tsv_gz(file.path(proc_dir, "truth_edges.tsv.gz"))

  directed <- truth[, c("partner_taxon", "focal_taxon", "effect_sign", "log2_ratio", "context")]
  names(directed) <- c("source_taxon", "target_taxon", "effect_sign", "strength", "evidence")
  directed$effect_sign <- as.integer(directed$effect_sign)
  directed$strength <- as.numeric(directed$strength)
  directed <- directed[!is.na(directed$effect_sign) & directed$effect_sign != 0, , drop = FALSE]
  directed <- directed[directed$source_taxon != directed$target_taxon, , drop = FALSE]

  directed_std <- build_directed(directed)
  write_tsv_gz(directed_std, file.path(proc_dir, "truth_directed.tsv.gz"))
  write_tsv_gz(build_undirected(directed), file.path(proc_dir, "truth_undirected.tsv.gz"))
}

targets <- if (length(args)) args else c("pairinterax", "omm12", "omm12_keystone_2023")
if ("all" %in% targets) targets <- c("pairinterax", "omm12", "omm12_keystone_2023", "butyrate_assembly_2021", "host_fitness_2018")

for (target in targets) {
  message("Building truth targets for ", target, " ...")
  switch(
    target,
    pairinterax = build_pairinterax(),
    omm12 = build_omm12(),
    omm12_keystone_2023 = build_omm12_keystone_2023(),
    butyrate_assembly_2021 = build_butyrate_assembly_2021(),
    host_fitness_2018 = build_host_fitness_2018(),
    venturelli_2018 = build_venturelli(),
    stop("Unknown target: ", target)
  )
}

message("Done.")
