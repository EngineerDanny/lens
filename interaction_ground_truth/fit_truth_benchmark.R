#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(PLNmodels)
  library(glassoFast)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[1] else "omm12"

base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth", mustWork = TRUE)
proc_dir <- file.path(base_dir, dataset, "processed")
out_dir <- file.path(proc_dir, "benchmark_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

sort_pair <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

prepare_count_like <- function(X, target_depth = 10000L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0

  is_integerish <- all(abs(X - round(X)) < 1e-8)
  if (is_integerish) {
    return(round(X))
  }

  rs <- rowSums(X)
  scaled <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  keep <- rs > 0
  scaled[keep, ] <- round(X[keep, , drop = FALSE] / rs[keep] * target_depth)
  scaled
}

load_benchmark_matrix <- function(dataset) {
  if (dataset == "pairinterax") {
    df <- read_tsv_gz(file.path(proc_dir, "abundance_matrix.tsv.gz"))
    X <- as.matrix(df[, -1, drop = FALSE])
    rownames(X) <- df$sample_id
    return(X)
  }

  if (dataset == "omm12") {
    df <- read_tsv_gz(file.path(proc_dir, "community_absabundance_in_vitro.tsv.gz"))
    X <- as.matrix(df[, -1, drop = FALSE])
    rownames(X) <- df$sample_id
    return(X)
  }

  if (dataset == "omm12_keystone_2023") {
    df <- read_tsv_gz(file.path(proc_dir, "community_absabundance_in_vitro.tsv.gz"))
    X <- as.matrix(df[, -1, drop = FALSE])
    rownames(X) <- df$sample_id
    return(X)
  }

  if (dataset == "butyrate_assembly_2021") {
    df <- read_tsv_gz(file.path(proc_dir, "abundance_matrix.tsv.gz"))
    X <- as.matrix(df[, -1, drop = FALSE])
    rownames(X) <- df$sample_id
    return(X)
  }

  if (dataset == "host_fitness_2018") {
    df <- read_tsv_gz(file.path(proc_dir, "abundance_matrix.tsv.gz"))
    X <- as.matrix(df[, -1, drop = FALSE])
    rownames(X) <- df$sample_id
    return(X)
  }

  if (dataset == "venturelli_2018") {
    ev_paths <- c(
      EV1 = file.path(proc_dir, "pairwise_time_series_EV1.tsv.gz"),
      EV2 = file.path(proc_dir, "pairwise_time_series_EV2.tsv.gz")
    )
    taxon_ids <- read_tsv_gz(file.path(proc_dir, "taxa_map.tsv.gz"))$taxon_id
    mats <- lapply(names(ev_paths), function(tag) {
      long <- read_tsv_gz(ev_paths[[tag]])
      long$sample_id <- paste(tag, long$community_code, sprintf("t%02d", long$time_index), sep = "_")
      samples <- unique(long$sample_id)
      X <- matrix(0, nrow = length(samples), ncol = length(taxon_ids), dimnames = list(samples, taxon_ids))
      for (i in seq_len(nrow(long))) {
        X[long$sample_id[i], long$focal_taxon[i]] <- long$abundance[i]
      }
      X
    })
    return(do.call(rbind, mats))
  }

  stop("Unsupported dataset: ", dataset)
}

pln_prepare_counts <- function(counts) {
  prepare_data(
    counts = counts,
    covariates = data.frame(Intercept = rep(1, nrow(counts)), row.names = rownames(counts)),
    offset = "none"
  )
}

support_from_precision <- function(Omega, tol = 1e-8) {
  support <- abs(Omega) > tol
  diag(support) <- FALSE
  support | t(support)
}

support_from_coef_mat <- function(coef_mat, require_mutual = TRUE, tol = 1e-8) {
  nz <- abs(coef_mat) > tol
  diag(nz) <- FALSE
  support <- if (require_mutual) nz & t(nz) else nz | t(nz)
  diag(support) <- FALSE
  support
}

loto_to_omega <- function(coef_matrix, sigma_diag) {
  p <- nrow(coef_matrix)
  Omega <- matrix(0, p, p)
  for (j in seq_len(p)) {
    Omega[j, j] <- 1 / sigma_diag[j]
    for (k in which(coef_matrix[j, ] != 0)) {
      Omega[j, k] <- -coef_matrix[j, k] * Omega[j, j]
    }
  }
  Omega <- (Omega + t(Omega)) / 2
  eig <- eigen(Omega, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-6)
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

fit_plnnetwork <- function(X_train_raw) {
  pln_init <- PLN(
    Abundance ~ 1,
    data = pln_prepare_counts(X_train_raw),
    control = PLN_param(backend = "nlopt", trace = 0)
  )
  sigma_pln <- pln_init$model_par$Sigma
  max_rho <- max(abs(sigma_pln[upper.tri(sigma_pln)]))
  rho_seq <- exp(seq(log(max_rho), log(max_rho * 0.005), length.out = 30))

  fit <- PLNnetwork(
    Abundance ~ 1,
    data = pln_prepare_counts(X_train_raw),
    penalties = rho_seq,
    control = PLNnetwork_param(
      backend = "nlopt",
      trace = 0,
      penalize_diagonal = FALSE,
      inception = pln_init
    )
  )
  fit$getBestModel(crit = "EBIC")$model_par$Omega
}

fit_loto_pln_glasso <- function(X_train_raw, X_train_log, rho_grid_len = 12L) {
  n <- nrow(X_train_raw)
  p <- ncol(X_train_raw)
  coef_mat <- matrix(0, p, p)
  sigma_diag <- pmax(apply(X_train_log, 2, var), 1e-6)
  for (target_idx in seq_len(p)) {
    order_idx <- c(target_idx, setdiff(seq_len(p), target_idx))
    counts_target <- X_train_raw[, order_idx, drop = FALSE]
    log_target <- X_train_log[, order_idx, drop = FALSE]
    fit <- PLN(
      Abundance ~ 1,
      data = pln_prepare_counts(counts_target),
      control = PLN_param(backend = "nlopt", trace = 0)
    )
    sigma_pln <- fit$model_par$Sigma
    max_rho <- max(abs(sigma_pln[upper.tri(sigma_pln)]))
    if (!is.finite(max_rho) || max_rho <= 0) next
    rho_seq <- exp(seq(log(max_rho), log(max_rho * 0.01), length.out = rho_grid_len))
    centered <- sweep(log_target, 2, colMeans(log_target))
    best_bic <- Inf
    best_omega <- NULL
    for (rho in rho_seq) {
      tryCatch({
        gl <- glassoFast(sigma_pln, rho = rho)
        omega <- gl$wi
        omega_jj <- omega[1, 1]
        if (!is.finite(omega_jj) || omega_jj <= 0) return(NULL)
        cond_mean <- -centered[, -1, drop = FALSE] %*% omega[-1, 1] / omega_jj
        resid <- centered[, 1] - as.vector(cond_mean)
        cond_ll_total <- sum(0.5 * log(omega_jj) - 0.5 * omega_jj * resid^2 - 0.5 * log(2 * pi))
        df <- sum(omega[1, -1] != 0)
        bic <- -2 * cond_ll_total + df * log(n)
        if (bic < best_bic) {
          best_bic <- bic
          best_omega <- omega
        }
      }, error = function(e) NULL)
    }
    if (is.null(best_omega)) next
    coef_mat[target_idx, setdiff(seq_len(p), target_idx)] <- -best_omega[1, -1] / best_omega[1, 1]
    sigma_diag[target_idx] <- 1 / best_omega[1, 1]
  }
  list(omega = loto_to_omega(coef_mat, sigma_diag), coef_mat = coef_mat)
}

safe_glmnet_coefs <- function(x_train, y_train) {
  p_minus_1 <- ncol(x_train)
  zero_coefs <- rep(0, p_minus_1)
  if (p_minus_1 == 0L || length(unique(y_train)) < 2L) return(zero_coefs)
  x_train <- as.matrix(x_train)
  cv_nfolds <- max(3L, min(5L, nrow(x_train)))
  fit <- tryCatch(cv.glmnet(x_train, y_train, family = "poisson", alpha = 1, nfolds = cv_nfolds), error = function(e) NULL)
  if (is.null(fit)) return(zero_coefs)
  coefs <- tryCatch(as.numeric(coef(fit, s = "lambda.1se")[-1]), error = function(e) zero_coefs)
  if (length(coefs) != p_minus_1) zero_coefs else coefs
}

fit_loto_glmnet <- function(X_train_raw, X_train_log) {
  p <- ncol(X_train_raw)
  coef_mat <- matrix(0, p, p)
  for (j in seq_len(p)) {
    y_train <- X_train_raw[, j]
    x_train <- X_train_log[, -j, drop = FALSE]
    coef_mat[j, -j] <- safe_glmnet_coefs(x_train, y_train)
  }
  sigma_diag <- pmax(apply(X_train_log, 2, var), 1e-6)
  list(omega = loto_to_omega(coef_mat, sigma_diag), coef_mat = coef_mat)
}

omega_to_prediction <- function(Omega, taxa) {
  idx <- which(upper.tri(Omega) & Omega != 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.frame(taxon_1 = character(), taxon_2 = character(), sign = integer(), weight = numeric(), stringsAsFactors = FALSE))
  }
  w <- Omega[idx]
  data.frame(
    taxon_1 = taxa[idx[, 1]],
    taxon_2 = taxa[idx[, 2]],
    sign = ifelse(w > 0, 1L, -1L),
    weight = as.numeric(w),
    stringsAsFactors = FALSE
  )
}

support_to_prediction <- function(support, weight_mat, taxa) {
  idx <- which(upper.tri(support) & support, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.frame(taxon_1 = character(), taxon_2 = character(), sign = integer(), weight = numeric(), stringsAsFactors = FALSE))
  }
  w <- weight_mat[idx]
  data.frame(
    taxon_1 = taxa[idx[, 1]],
    taxon_2 = taxa[idx[, 2]],
    sign = ifelse(w > 0, 1L, -1L),
    weight = as.numeric(w),
    stringsAsFactors = FALSE
  )
}

compute_metrics <- function(truth, pred) {
  truth$key <- sort_pair(truth$taxon_1, truth$taxon_2)
  pred$key <- sort_pair(pred$taxon_1, pred$taxon_2)
  truth_keys <- unique(truth$key)
  pred_keys <- unique(pred$key)
  tp_keys <- intersect(truth_keys, pred_keys)
  fp_keys <- setdiff(pred_keys, truth_keys)
  fn_keys <- setdiff(truth_keys, pred_keys)
  tp <- length(tp_keys); fp <- length(fp_keys); fn <- length(fn_keys)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1 <- if (is.finite(precision + recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_
  sign_accuracy <- NA_real_
  if ("sign_consensus" %in% names(truth) && "sign" %in% names(pred)) {
    tsub <- truth[!is.na(truth$sign_consensus), c("key", "sign_consensus"), drop = FALSE]
    psub <- pred[!is.na(pred$sign), c("key", "sign"), drop = FALSE]
    merged <- merge(tsub, psub, by = "key")
    if (nrow(merged) > 0) sign_accuracy <- mean(merged$sign_consensus == merged$sign)
  }
  data.frame(
    dataset = dataset,
    n_truth = length(truth_keys),
    n_pred = length(pred_keys),
    tp = tp, fp = fp, fn = fn,
    precision = precision,
    recall = recall,
    f1 = f1,
    sign_accuracy = sign_accuracy,
    stringsAsFactors = FALSE
  )
}

message("Loading dataset: ", dataset)
X_input <- load_benchmark_matrix(dataset)
X_counts <- prepare_count_like(X_input)
X_counts <- X_counts[rowSums(X_counts) > 0, colSums(X_counts) > 0, drop = FALSE]
X_log <- log1p(X_counts)
taxa <- colnames(X_counts)
truth <- read_tsv_gz(file.path(proc_dir, "truth_undirected.tsv.gz"))

methods <- c("Baseline_diagonal", "PLNnetwork", "LOTO_glmnet_CV1se")
metric_rows <- list()

for (method in methods) {
  message("Fitting ", method, " ...")
  t0 <- proc.time()[["elapsed"]]
  if (method == "Baseline_diagonal") {
    pred <- data.frame(taxon_1 = character(), taxon_2 = character(), sign = integer(), weight = numeric(), stringsAsFactors = FALSE)
  } else if (method == "PLNnetwork") {
    omega <- fit_plnnetwork(X_counts)
    pred <- omega_to_prediction(omega, taxa)
  } else if (method == "LOTO_glmnet_CV1se") {
    fit <- fit_loto_glmnet(X_counts, X_log)
    support <- support_from_coef_mat(fit$coef_mat, require_mutual = TRUE)
    weight_mat <- (fit$coef_mat + t(fit$coef_mat)) / 2
    pred <- support_to_prediction(support, weight_mat, taxa)
  } else {
    stop("Unknown method: ", method)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  pred_path <- file.path(out_dir, paste0(method, "_predictions.tsv.gz"))
  write_tsv_gz(pred, pred_path)
  metrics <- compute_metrics(truth, pred)
  metrics$method <- method
  metrics$elapsed_sec <- elapsed
  metric_rows[[method]] <- metrics
}

summary_df <- rbindlist(metric_rows, fill = TRUE)
fwrite(summary_df[, c("dataset", "method", "n_truth", "n_pred", "tp", "fp", "fn", "precision", "recall", "f1", "sign_accuracy", "elapsed_sec")], file.path(out_dir, "metrics_summary.tsv"), sep = "\t")
print(summary_df[, c("dataset", "method", "precision", "recall", "f1", "sign_accuracy", "elapsed_sec")])
