#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- arguments[grep("^--file=", arguments)][1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
root <- dirname(dirname(script_path))
build_script <- file.path(root, "scripts", "build_pair_features.R")
local_library <- file.path(root, "r_library")

systems <- c(
  "butyrate_assembly_2021",
  "carlstrom_phyllosphere_2019",
  "schafer_phyllosphere_2022"
)
methods <- c(
  PLNNetwork = "pln",
  `Poisson GLMNet` = "glmnet",
  SparCC = "sparcc",
  `SPIEC-EASI` = "spieceasi"
)
replicates <- 5L
seed <- 20260825L

canonicalize <- function(frame) {
  first <- as.character(frame$taxon_1)
  second <- as.character(frame$taxon_2)
  frame$taxon_1 <- pmin(first, second)
  frame$taxon_2 <- pmax(first, second)
  frame
}

average_precision <- function(labels, scores) {
  ordering <- order(scores, decreasing = TRUE)
  labels <- as.integer(labels[ordering])
  scores <- as.numeric(scores[ordering])
  positives <- sum(labels == 1L)
  threshold_end <- c(which(diff(scores) != 0), length(scores))
  cumulative_positive <- cumsum(labels == 1L)[threshold_end]
  recall <- cumulative_positive / positives
  precision <- cumulative_positive / threshold_end
  sum(c(recall[1], diff(recall)) * precision)
}

load_truth <- function(system) {
  suffix <- "_tested_pairs.csv"
  truth <- fread(file.path(root, "cleaned_data", paste0(system, suffix)))
  truth <- canonicalize(truth)
  truth[!is.na(interaction_label)]
}

run_one <- function(task) {
  system <- task$system
  bootstrap_id <- task$bootstrap
  set.seed(seed + match(system, systems) * 100L + bootstrap_id)
  abundance <- fread(file.path(root, "cleaned_data", paste0(system, "_abundance.csv")))
  sampled_rows <- sample.int(nrow(abundance), nrow(abundance), replace = TRUE)
  sampled <- copy(abundance[sampled_rows])
  setnames(sampled, 1L, "sample_id")
  sampled[[1L]] <- sprintf("bootstrap_%02d_sample_%05d", bootstrap_id, seq_len(nrow(sampled)))

  variable <- vapply(sampled[, -1L, with = FALSE], function(values) {
    length(unique(values)) > 1L
  }, logical(1))
  omitted_taxa <- names(variable)[!variable]
  sampled <- sampled[, c(TRUE, variable), with = FALSE]

  work <- tempfile(pattern = paste0("network_bootstrap_", system, "_"))
  dir.create(file.path(work, "cleaned_data"), recursive = TRUE)
  dir.create(file.path(work, "analysis_data"), recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
  file.symlink(local_library, file.path(work, "r_library"))
  fwrite(sampled, file.path(work, "cleaned_data", paste0(system, "_abundance.csv")))

  previous_directory <- getwd()
  setwd(work)
  command_output <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = c(build_script, system),
    stdout = TRUE,
    stderr = TRUE,
    env = "PAIR_METHODS=pln,glmnet,spieceasi,sparcc",
    wait = TRUE
  )
  setwd(previous_directory)
  exit_status <- attr(command_output, "status")
  if (is.null(exit_status)) exit_status <- 0L

  status_path <- file.path(work, "analysis_data", paste0(system, "_method_status.csv"))
  feature_path <- file.path(work, "analysis_data", paste0(system, "_pair_features.csv"))
  if (exit_status != 0L || !file.exists(feature_path)) {
    return(data.table(
      analysis_set = system,
      bootstrap = bootstrap_id,
      method = names(methods),
      auprc = NA_real_,
      omitted_taxa = length(omitted_taxa),
      fit_status = "failed",
      detail = paste(tail(command_output, 1L), collapse = " ")
    ))
  }

  fit_status <- fread(status_path)
  features <- canonicalize(fread(feature_path))
  truth <- load_truth(system)
  all_taxa <- names(abundance)[-1L]
  all_pairs <- as.data.table(t(combn(all_taxa, 2L)))
  setnames(all_pairs, c("taxon_1", "taxon_2"))
  all_pairs <- canonicalize(all_pairs)

  output <- list()
  for (method_name in names(methods)) {
    method_key <- unname(methods[method_name])
    raw_column <- paste0(method_key, "_abs_score")
    status <- fit_status[method == method_key]
    if (!nrow(status) || status$status != "complete" || !raw_column %in% names(features)) {
      output[[length(output) + 1L]] <- data.table(
        analysis_set = system, bootstrap = bootstrap_id, method = method_name,
        auprc = NA_real_, omitted_taxa = length(omitted_taxa),
        fit_status = "failed", detail = if (nrow(status)) status$detail else "missing fit"
      )
      next
    }
    scores <- merge(
      all_pairs,
      features[, c("taxon_1", "taxon_2", raw_column), with = FALSE],
      by = c("taxon_1", "taxon_2"),
      all.x = TRUE
    )
    set(scores, which(is.na(scores[[raw_column]])), raw_column, 0)
    scores[, score_percentile := frank(get(raw_column), ties.method = "average") / .N]
    evaluated <- merge(
      truth[, .(taxon_1, taxon_2, interaction_label)],
      scores[, .(taxon_1, taxon_2, score_percentile)],
      by = c("taxon_1", "taxon_2"),
      all = FALSE,
      allow.cartesian = TRUE
    )
    output[[length(output) + 1L]] <- data.table(
      analysis_set = system, bootstrap = bootstrap_id, method = method_name,
      auprc = average_precision(evaluated$interaction_label, evaluated$score_percentile),
      omitted_taxa = length(omitted_taxa), fit_status = "complete", detail = ""
    )
  }
  rbindlist(output, fill = TRUE)
}

tasks <- CJ(system = systems, bootstrap = seq_len(replicates), sorted = FALSE)
task_list <- lapply(seq_len(nrow(tasks)), function(index) as.list(tasks[index]))
cores <- min(4L, detectCores(), length(task_list))
message("Running ", nrow(tasks), " bootstrap fits with ", cores, " workers")
results <- rbindlist(mclapply(task_list, run_one, mc.cores = cores), fill = TRUE)
setorder(results, analysis_set, method, bootstrap)
fwrite(results, file.path(root, "results", "abundance_bootstrap_5_auprc.csv"))
message("Wrote abundance bootstrap results")
