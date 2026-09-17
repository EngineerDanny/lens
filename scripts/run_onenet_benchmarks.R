#!/usr/bin/env Rscript

.libPaths(c(file.path(getwd(), "r_library"), .libPaths()))
suppressPackageStartupMessages(library(OneNet))
options(future.globals.maxSize = 4 * 1024^3)

root <- normalizePath(getwd())
systems <- c(
  "butyrate_assembly_2021",
  "carlstrom_phyllosphere_2019",
  "schafer_phyllosphere_2022"
)
rep_num <- as.integer(Sys.getenv("ONENET_REPS", "30"))
cores <- as.integer(Sys.getenv("ONENET_CORES", "4"))
seed <- 20260821L

prepare_count_like <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x)) || any(x < 0)) stop("Abundance values must be finite and nonnegative")
  totals <- rowSums(x)
  keep <- totals > 0
  if (!all(keep)) x <- x[keep, , drop = FALSE]
  round(x / rowSums(x) * 10000)
}

edge_table <- function(taxa, score) {
  shell <- matrix(0, length(taxa), length(taxa), dimnames = list(taxa, taxa))
  stopifnot(length(score) == sum(upper.tri(shell)))
  shell[upper.tri(shell)] <- score
  indices <- which(upper.tri(shell), arr.ind = TRUE)
  data.frame(
    taxon_1 = rownames(shell)[indices[, 1]],
    taxon_2 = colnames(shell)[indices[, 2]],
    onenet_score = score,
    stringsAsFactors = FALSE
  )
}

adapt_mean_stability_safe <- function(inference_collection, mean_stability = 0.8) {
  ne <- inference_collection[[1]]$nedges
  methods <- names(inference_collection)
  stability_path <- do.call(rbind, lapply(seq_along(inference_collection), function(i) {
    cbind(inference_collection[[i]]$lambda_stab, method = methods[i])
  }))
  candidates <- if (ne < 40) seq_len(max(1, floor(ne / 4))) else seq(10, ne / 4, 10)
  curve <- do.call(rbind, lapply(candidates, function(target) {
    chosen <- do.call(rbind, lapply(methods, function(method) {
      part <- stability_path[stability_path$method == method, , drop = FALSE]
      part[which.min(abs(part$nedges - target)), , drop = FALSE]
    }))
    data.frame(mean_stability = mean(chosen$stability), edges = target)
  }))
  minimum <- curve$edges[which.min(curve$mean_stability)]
  eligible <- curve[curve$edges < minimum, , drop = FALSE]
  if (!nrow(eligible)) eligible <- curve
  final_edges <- eligible$edges[which.min(abs(eligible$mean_stability - mean_stability))]
  selected <- do.call(rbind, lapply(methods, function(method) {
    part <- stability_path[stability_path$method == method, , drop = FALSE]
    part[which.min(abs(part$nedges - final_edges)), , drop = FALSE]
  }))
  frequencies <- do.call(cbind, lapply(seq_along(inference_collection), function(i) {
    method <- methods[i]
    lambda <- selected$lambda[selected$method == method][1]
    OneNet:::get_vec(lambda, inference_collection[[i]], method)
  }))
  colnames(frequencies) <- methods
  list(
    freqs = as.data.frame(frequencies),
    stab_data = selected[, c("stability", "nedges", "method"), drop = FALSE]
  )
}

dir.create(file.path(root, "analysis_data", "onenet_cache"), recursive = TRUE, showWarnings = FALSE)

for (index in seq_along(systems)) {
  system <- systems[index]
  message("OneNet: ", system, " with ", rep_num, " resamples")
  abundance <- read.csv(
    file.path(root, "cleaned_data", paste0(system, "_abundance.csv")),
    check.names = FALSE
  )
  taxa <- names(abundance)[-1]
  counts <- prepare_count_like(abundance[, -1, drop = FALSE])
  colnames(counts) <- taxa

  cache <- file.path(root, "analysis_data", "onenet_cache",
                     paste0(system, "_six_methods_reps", rep_num, "_inferences.rds"))
  if (file.exists(cache)) {
    inference <- readRDS(cache)
  } else {
    set.seed(seed + index - 1L)
    inference <- all_inferences_new(
      data = counts,
      rep.num = rep_num,
      methods = c("PLNnetwork", "SpiecEasi", "gCoda", "EMtree", "Magma", "SPRING"),
      parallel = FALSE,
      cores = cores,
      seed = seed + index - 1L
    )
    saveRDS(inference, cache)
  }

  adapted <- adapt_mean_stability_safe(inference, mean_stability = 0.8)
  aggregate <- compute_aggreg_measures(adapted$freqs)
  scores <- edge_table(taxa, aggregate$mean)
  scores$dataset <- system
  scores <- scores[, c("dataset", "taxon_1", "taxon_2", "onenet_score")]
  write.csv(
    scores,
    file.path(root, "analysis_data", paste0(system, "_onenet_pair_scores.csv")),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(adapted$stab_data),
    file.path(root, "analysis_data", paste0(system, "_onenet_stability.csv")),
    row.names = FALSE
  )
}

versions <- data.frame(
  package = c("OneNet", "PLNmodels", "SpiecEasi", "SPRING", "EMtree", "rMAGMA"),
  version = vapply(c("OneNet", "PLNmodels", "SpiecEasi", "SPRING", "EMtree", "rMAGMA"),
                   function(x) as.character(packageVersion(x)), character(1))
)
write.csv(versions, file.path(root, "results", "onenet_package_versions.csv"), row.names = FALSE)
