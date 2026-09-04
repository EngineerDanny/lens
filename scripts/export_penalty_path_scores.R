#!/usr/bin/env Rscript

local_lib <- normalizePath("r_library", mustWork = TRUE)
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

datasets <- commandArgs(trailingOnly = TRUE)
if (!length(datasets)) {
  datasets <- c(
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
    "friedman_microcosm_2017"
  )
}

partial_from_precision <- function(Omega) {
  Omega <- as.matrix(Omega)
  denom <- sqrt(outer(diag(Omega), diag(Omega)))
  out <- -Omega / denom
  out[!is.finite(out)] <- 0
  diag(out) <- 0
  out
}

path_summaries <- function(score_path, pair_index, prefix, stability_path = NULL) {
  n_pairs <- length(pair_index)
  values <- vapply(score_path, function(M) {
    M <- as.matrix(M)
    abs(M[pair_index])
  }, numeric(n_pairs))
  if (is.null(dim(values))) values <- matrix(values, ncol = 1)
  present <- values > 1e-10
  first <- apply(present, 1, function(x) {
    hit <- which(x)
    if (length(hit)) hit[[1]] else NA_integer_
  })
  n_path <- ncol(values)
  entry <- ifelse(is.na(first), 0, (n_path - first + 1) / n_path)
  out <- data.table(
    entry = entry,
    presence = rowMeans(present),
    mean_abs = rowMeans(values),
    max_abs = apply(values, 1, max),
    terminal_abs = values[, n_path]
  )
  setnames(out, names(out), paste0(prefix, "_", names(out)))

  if (!is.null(stability_path)) {
    stability <- vapply(stability_path, function(M) {
      M <- as.matrix(M)
      abs(M[pair_index])
    }, numeric(n_pairs))
    if (is.null(dim(stability))) stability <- matrix(stability, ncol = 1)
    out[, paste0(prefix, "_stability_mean") := rowMeans(stability)]
    out[, paste0(prefix, "_stability_max") := apply(stability, 1, max)]
  }
  out
}

dir.create("analysis_data", recursive = TRUE, showWarnings = FALSE)

for (dataset in datasets) {
  message("Exporting penalty path scores for ", dataset)
  abundance <- fread(file.path("cleaned_data", paste0(dataset, "_abundance.csv")), check.names = FALSE)
  taxa <- names(abundance)[-1]
  p <- length(taxa)
  pair_index <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  pair_linear <- which(upper.tri(matrix(0, p, p)))
  out <- data.table(
    dataset = dataset,
    taxon_1 = taxa[pair_index[, 1]],
    taxon_2 = taxa[pair_index[, 2]]
  )

  pln <- readRDS(file.path("analysis_cache", dataset, "pln.rds"))
  pln_path <- lapply(pln$models, function(model) partial_from_precision(model$model_par$Omega))
  out <- cbind(out, path_summaries(pln_path, pair_linear, "pln_path"))

  spiec <- readRDS(file.path("analysis_cache", dataset, "spieceasi.rds"))
  spiec_path <- lapply(spiec$est$icov, partial_from_precision)
  out <- cbind(out, path_summaries(
    spiec_path, pair_linear, "spieceasi_path", spiec$select$stars$merge
  ))

  spring <- readRDS(file.path("analysis_cache", dataset, "spring.rds"))
  spring_path <- lapply(spring$fit$est$beta, function(beta) {
    beta <- as.matrix(beta)
    score <- (beta + t(beta)) / 2
    diag(score) <- 0
    score
  })
  out <- cbind(out, path_summaries(
    spring_path, pair_linear, "spring_path", spring$output$stars$merge
  ))

  fwrite(out, file.path("analysis_data", paste0(dataset, "_penalty_path_scores.csv")))
}
