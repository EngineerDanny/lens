#!/usr/bin/env Rscript

local_lib <- normalizePath("r_library", mustWork = TRUE)
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(PLNmodels)
  library(glmnet)
  library(SpiecEasi)
  library(SPRING)
})

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) args[[1]] else stop("Provide a dataset name")

root <- normalizePath(".", mustWork = TRUE)
input_path <- file.path(root, "cleaned_data", paste0(dataset, "_abundance.csv"))
output_dir <- file.path(root, "analysis_data")
cache_dir <- file.path(root, "analysis_cache", dataset)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 426001L + match(dataset, c(
  "omm12", "omm12_keystone_2023", "pairinterax",
  "butyrate_assembly_2021", "host_fitness_2018", "wortel_syncom_2026",
  "schafer_phyllosphere_2022", "carlstrom_phyllosphere_2019",
  "friedman_microcosm_2017"
))
set.seed(seed)

abundance <- fread(input_path, check.names = FALSE)
sample_id <- abundance[[1]]
X_native <- as.matrix(abundance[, -1, with = FALSE])
storage.mode(X_native) <- "double"
rownames(X_native) <- sample_id
taxa <- colnames(X_native)
n <- nrow(X_native)
p <- ncol(X_native)

if (any(!is.finite(X_native)) || any(X_native < 0)) stop("Invalid abundance values")
if (any(rowSums(X_native) <= 0)) stop("Zero-total rows remain")
if (any(apply(X_native, 2, function(x) length(unique(x)) <= 1))) stop("Constant taxa remain")

measurement <- switch(
  dataset,
  omm12 = "absolute_abundance",
  omm12_keystone_2023 = "absolute_abundance",
  pairinterax = "relative_abundance",
  butyrate_assembly_2021 = "proportion",
  host_fitness_2018 = "cfu_abundance",
  wortel_syncom_2026 = "qpcr_abundance",
  schafer_phyllosphere_2022 = "amplicon_count",
  carlstrom_phyllosphere_2019 = "amplicon_count",
  friedman_microcosm_2017 = "absolute_abundance"
)

prepare_count_like <- function(X, target_depth = 10000L) {
  is_integerish <- all(abs(X - round(X)) < 1e-8)
  if (is_integerish) return(round(X))
  rs <- rowSums(X)
  out <- round(X / rs * target_depth)
  dimnames(out) <- dimnames(X)
  out
}

X_count <- prepare_count_like(X_native)

pair_index <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
pair_table <- data.table(
  dataset = dataset,
  taxon_1 = taxa[pair_index[, 1]],
  taxon_2 = taxa[pair_index[, 2]]
)

extract_upper <- function(M) as.numeric(M[pair_index])

partial_from_precision <- function(Omega) {
  Omega <- as.matrix(Omega)
  denom <- sqrt(outer(diag(Omega), diag(Omega)))
  out <- -Omega / denom
  out[!is.finite(out)] <- 0
  diag(out) <- 0
  out
}

score_percentile <- function(values) {
  values <- abs(values)
  rank(values, ties.method = "average") / length(values)
}

record_method <- function(method, status, seconds, detail = "") {
  data.table(
    dataset = dataset,
    method = method,
    status = status,
    elapsed_seconds = round(seconds, 3),
    detail = detail
  )
}

method_rows <- list()

add_method <- function(method, score_matrix, selected_matrix, stability_matrix = NULL) {
  raw <- extract_upper(score_matrix)
  selected <- as.integer(extract_upper(selected_matrix) != 0)
  pair_table[, paste0(method, "_score") := raw]
  pair_table[, paste0(method, "_abs_score") := abs(raw)]
  pair_table[, paste0(method, "_score_percentile") := score_percentile(raw)]
  pair_table[, paste0(method, "_selected") := selected]
  if (!is.null(stability_matrix)) {
    pair_table[, paste0(method, "_stability") := extract_upper(stability_matrix)]
  }
}

pln_prepare <- function(counts) {
  prepare_data(
    counts = counts,
    covariates = data.frame(Intercept = rep(1, nrow(counts)), row.names = rownames(counts)),
    offset = "none"
  )
}

fit_pln <- function() {
  init <- PLN(
    Abundance ~ 1,
    data = pln_prepare(X_count),
    control = PLN_param(backend = "nlopt", trace = 0)
  )
  sigma <- init$model_par$Sigma
  max_rho <- max(abs(sigma[upper.tri(sigma)]))
  penalties <- exp(seq(log(max_rho), log(max_rho * 0.005), length.out = 30))
  fit <- PLNnetwork(
    Abundance ~ 1,
    data = pln_prepare(X_count),
    penalties = penalties,
    control = PLNnetwork_param(
      backend = "nlopt", trace = 0, penalize_diagonal = FALSE, inception = init
    )
  )
  Omega <- fit$getBestModel(crit = "EBIC")$model_par$Omega
  list(score = partial_from_precision(Omega), selected = abs(Omega) > 1e-8, fit = fit)
}

fit_glmnet <- function() {
  X_log <- log1p(X_count)
  coefficients <- matrix(0, p, p, dimnames = list(taxa, taxa))
  selected_directed <- matrix(FALSE, p, p, dimnames = list(taxa, taxa))
  for (j in seq_len(p)) {
    y <- X_count[, j]
    predictors <- X_log[, -j, drop = FALSE]
    if (length(unique(y)) < 2 || !ncol(predictors)) next
    set.seed(seed + j)
    fit <- tryCatch(
      cv.glmnet(
        predictors, y, family = "poisson", alpha = 1,
        nfolds = max(3L, min(5L, nrow(predictors)))
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    beta <- as.numeric(coef(fit, s = "lambda.1se")[-1])
    coefficients[j, -j] <- beta
    selected_directed[j, -j] <- abs(beta) > 1e-10
  }
  score <- (coefficients + t(coefficients)) / 2
  selected <- selected_directed & t(selected_directed)
  diag(score) <- 0
  diag(selected) <- FALSE
  list(score = score, selected = selected, fit = coefficients)
}

fit_spiec <- function() {
  set.seed(seed)
  fit <- spiec.easi(
    X_count,
    method = "glasso",
    nlambda = 30,
    lambda.min.ratio = 0.05,
    verbose = FALSE,
    pulsar.params = list(rep.num = 20, ncores = 1)
  )
  opt <- getOptInd(fit)
  stability <- as.matrix(fit$select$stars$merge[[opt]])
  list(
    score = partial_from_precision(getOptiCov(fit)),
    selected = as.matrix(getRefit(fit)) != 0,
    stability = stability,
    fit = fit
  )
}

fit_spring <- function() {
  quantitative <- measurement %in% c("absolute_abundance", "cfu_abundance") || !any(X_native == 0)
  set.seed(seed)
  fit <- SPRING(
    X_native,
    quantitative = quantitative,
    nlambda = 20,
    lambdaseq = "data-specific",
    seed = seed,
    ncores = 1,
    rep.num = 20,
    verbose = FALSE,
    verboseR = FALSE,
    Rmethod = "approx"
  )
  opt <- fit$output$stars$opt.index
  if (is.null(opt) || opt < 1) opt <- 1L
  beta <- as.matrix(fit$fit$est$beta[[opt]])
  score <- (beta + t(beta)) / 2
  selected <- as.matrix(fit$fit$refit$stars) != 0
  stability <- as.matrix(fit$output$stars$merge[[opt]])
  list(score = score, selected = selected, stability = stability, fit = fit)
}

fit_sparcc <- function() {
  set.seed(seed)
  fit <- SpiecEasi::sparcc(X_count, iter = 20, inner_iter = 10, th = 0.1)
  score <- as.matrix(fit$Cor)
  diag(score) <- 0
  selected <- abs(score) >= 0.3
  diag(selected) <- FALSE
  list(score = score, selected = selected, fit = fit)
}

methods <- list(
  pln = fit_pln,
  glmnet = fit_glmnet,
  spieceasi = fit_spiec,
  spring = fit_spring,
  sparcc = fit_sparcc
)

requested_methods <- Sys.getenv("PAIR_METHODS", "")
if (nzchar(requested_methods)) {
  requested_methods <- strsplit(requested_methods, ",", fixed = TRUE)[[1]]
  unknown_methods <- setdiff(requested_methods, names(methods))
  if (length(unknown_methods)) {
    stop("Unknown methods in PAIR_METHODS: ", paste(unknown_methods, collapse = ", "))
  }
  methods <- methods[requested_methods]
}

for (method in names(methods)) {
  message("Fitting ", method, " for ", dataset)
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(methods[[method]](), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(result, "error")) {
    method_rows[[method]] <- record_method(method, "failed", elapsed, conditionMessage(result))
    message("  failed: ", conditionMessage(result))
    next
  }
  add_method(method, result$score, result$selected, result$stability)
  saveRDS(result$fit, file.path(cache_dir, paste0(method, ".rds")))
  method_rows[[method]] <- record_method(method, "complete", elapsed)
  message("  complete in ", round(elapsed, 1), " seconds")
}

presence <- X_native > 0
prevalence <- colMeans(presence)
positive_mean <- vapply(seq_len(p), function(j) mean(X_native[presence[, j], j]), numeric(1))
dispersion <- vapply(seq_len(p), function(j) {
  values <- X_native[, j]
  mean_value <- mean(values)
  if (mean_value > 0) var(values) / mean_value else NA_real_
}, numeric(1))

i <- pair_index[, 1]
j <- pair_index[, 2]
joint_count <- vapply(seq_len(nrow(pair_index)), function(k) sum(presence[, i[k]] & presence[, j[k]]), integer(1))
only_1 <- vapply(seq_len(nrow(pair_index)), function(k) sum(presence[, i[k]] & !presence[, j[k]]), integer(1))
only_2 <- vapply(seq_len(nrow(pair_index)), function(k) sum(!presence[, i[k]] & presence[, j[k]]), integer(1))
neither <- n - joint_count - only_1 - only_2
union_count <- joint_count + only_1 + only_2

pair_table[, `:=`(
  prevalence_1 = prevalence[i],
  prevalence_2 = prevalence[j],
  prevalence_min = pmin(prevalence[i], prevalence[j]),
  prevalence_max = pmax(prevalence[i], prevalence[j]),
  joint_prevalence = joint_count / n,
  n_joint_positive = joint_count,
  n_taxon_1_only = only_1,
  n_taxon_2_only = only_2,
  n_neither = neither,
  presence_jaccard = ifelse(union_count > 0, joint_count / union_count, 0),
  mean_positive_abundance_1 = positive_mean[i],
  mean_positive_abundance_2 = positive_mean[j],
  dispersion_1 = dispersion[i],
  dispersion_2 = dispersion[j],
  n_samples = n,
  n_taxa = p,
  overall_zero_frequency = mean(X_native == 0),
  measurement_scale = measurement,
  count_conversion = ifelse(all(abs(X_native - round(X_native)) < 1e-8), "native_integer", "fixed_total_10000")
)]

setorder(pair_table, taxon_1, taxon_2)
fwrite(pair_table, file.path(output_dir, paste0(dataset, "_pair_features.csv")), na = "")
fwrite(rbindlist(method_rows, fill = TRUE), file.path(output_dir, paste0(dataset, "_method_status.csv")), na = "")

message("Wrote ", nrow(pair_table), " possible pairs for ", dataset)
