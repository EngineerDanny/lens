#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(PLNmodels)
  library(glassoFast)
  library(glmnet)
})

base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth/cdiff_inhibition_2024", mustWork = TRUE)
proc_dir <- file.path(base_dir, "processed")
out_dir <- file.path(proc_dir, "incoming_edge_benchmark")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

prepare_count_like <- function(X, target_depth = 10000L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  is_integerish <- all(abs(X - round(X)) < 1e-8)
  if (is_integerish) return(round(X))
  rs <- rowSums(X)
  scaled <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  keep <- rs > 0
  scaled[keep, ] <- round(X[keep, , drop = FALSE] / rs[keep] * target_depth)
  scaled
}

pln_prepare_counts <- function(counts) {
  prepare_data(
    counts = counts,
    covariates = data.frame(Intercept = rep(1, nrow(counts)), row.names = rownames(counts)),
    offset = "none"
  )
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

fit_pln_glasso <- function(X_train_raw) {
  n <- nrow(X_train_raw)
  fit <- PLN(
    Abundance ~ 1,
    data = pln_prepare_counts(X_train_raw),
    control = PLN_param(backend = "nlopt", trace = 0)
  )
  sigma_pln <- fit$model_par$Sigma
  p <- ncol(sigma_pln)
  max_rho <- max(abs(sigma_pln[upper.tri(sigma_pln)]))
  rho_seq <- exp(seq(log(max_rho), log(max(max_rho * 0.01, 1e-6)), length.out = 20))
  best_rho <- rho_seq[1]
  best_ebic <- Inf
  for (rho in rho_seq) {
    tryCatch({
      gl <- glassoFast(sigma_pln, rho = rho)
      logdet <- determinant(gl$wi, logarithm = TRUE)$modulus
      loglik_train <- (n / 2) * (logdet - sum(diag(gl$wi %*% sigma_pln)))
      n_edges <- sum(gl$wi[upper.tri(gl$wi)] != 0)
      ebic <- -2 * loglik_train + n_edges * log(n) + 4 * 0.5 * n_edges * log(p)
      if (ebic < best_ebic) {
        best_ebic <- ebic
        best_rho <- rho
      }
    }, error = function(e) NULL)
  }
  glassoFast(sigma_pln, rho = best_rho)$wi
}

fit_loto_pln_glasso <- function(X_train_raw, X_train_log, rho_grid_len = 10L) {
  n <- nrow(X_train_raw)
  p <- ncol(X_train_raw)
  coef_mat <- matrix(0, p, p, dimnames = list(colnames(X_train_raw), colnames(X_train_raw)))
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
    rho_seq <- exp(seq(log(max_rho), log(max(max_rho * 0.01, 1e-6)), length.out = rho_grid_len))
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
  coef_mat <- matrix(0, p, p, dimnames = list(colnames(X_train_raw), colnames(X_train_raw)))
  for (j in seq_len(p)) {
    y_train <- X_train_raw[, j]
    x_train <- X_train_log[, -j, drop = FALSE]
    coef_mat[j, -j] <- safe_glmnet_coefs(x_train, y_train)
  }
  sigma_diag <- pmax(apply(X_train_log, 2, var), 1e-6)
  list(omega = loto_to_omega(coef_mat, sigma_diag), coef_mat = coef_mat)
}

pln_inhibition_vector <- function(omega, cd_taxon, candidate_taxa) {
  idx <- match(cd_taxon, colnames(omega))
  out <- setNames(rep(0, length(candidate_taxa)), candidate_taxa)
  off <- seq_len(ncol(omega))[-idx]
  pc <- -omega[idx, off] / sqrt(omega[idx, idx] * diag(omega)[off])
  names(pc) <- colnames(omega)[off]
  out[names(pc)] <- pmax(-pc, 0)
  out
}

loto_inhibition_vector <- function(coef_mat, cd_taxon, candidate_taxa) {
  idx <- match(cd_taxon, rownames(coef_mat))
  out <- setNames(rep(0, length(candidate_taxa)), candidate_taxa)
  eff <- coef_mat[idx, setdiff(colnames(coef_mat), cd_taxon)]
  out[names(eff)] <- pmax(-eff, 0)
  out
}

eval_incoming_edges <- function(scores_df, k = 2L) {
  scores <- scores_df$predicted_inhibition_score
  labels <- scores_df$label
  n_pos <- sum(labels == 1L)
  if (!n_pos) stop("No positives in incoming-edge truth.")

  if (length(unique(scores)) <= 1L) {
    top_hits <- 0L
    precision_at_k <- 0
    recall_at_k <- 0
    f1_at_k <- 0
    average_precision <- 0
    positive_ranks <- rep(length(scores) + 1, n_pos)
  } else {
    ord <- order(-scores, scores_df$source_taxon)
    labels_ord <- labels[ord]
    top_idx <- ord[seq_len(min(k, length(ord)))]
    top_hits <- sum(labels[top_idx] == 1L)
    precision_at_k <- top_hits / min(k, length(ord))
    recall_at_k <- top_hits / n_pos
    f1_at_k <- if ((precision_at_k + recall_at_k) == 0) 0 else {
      2 * precision_at_k * recall_at_k / (precision_at_k + recall_at_k)
    }
    cum_tp <- cumsum(labels_ord == 1L)
    prec_seq <- cum_tp / seq_along(labels_ord)
    average_precision <- sum(prec_seq[labels_ord == 1L]) / n_pos
    positive_ranks <- match(which(labels == 1L), ord)
  }

  data.frame(
    top_hits = top_hits,
    precision_at_2 = precision_at_k,
    recall_at_2 = recall_at_k,
    f1_at_2 = f1_at_k,
    average_precision = average_precision,
    mean_positive_rank = mean(positive_ranks),
    stringsAsFactors = FALSE
  )
}

abund <- read_tsv_gz(file.path(proc_dir, "cdiff_incoming_abundance.tsv.gz"))
truth <- read_tsv_gz(file.path(proc_dir, "cdiff_incoming_truth.tsv.gz"))

candidate_taxa <- c("BU", "CA", "CS", "CH", "DP")
strains <- c("DSM", "MS001", "MS008", "MS014")
methods <- c("Baseline_diagonal", "PLN_glasso", "LOTO_PLN_glasso", "LOTO_glmnet_CV1se")

pred_rows <- list()
metric_rows <- list()

for (strain in strains) {
  fit_df <- abund[is.na(abund$strain) | abund$strain == "" | abund$strain == strain, , drop = FALSE]
  X_counts <- prepare_count_like(as.matrix(fit_df[, c(candidate_taxa, "CD"), drop = FALSE]))
  X_counts <- X_counts[rowSums(X_counts) > 0, colSums(X_counts) > 0, drop = FALSE]
  X_log <- log1p(X_counts)

  truth_strain <- truth[truth$strain == strain, , drop = FALSE]
  if (!nrow(truth_strain)) next

  for (method in methods) {
    message("Fitting ", strain, " / ", method, " ...")
    t0 <- proc.time()[["elapsed"]]
    inhib_vec <- setNames(rep(0, length(candidate_taxa)), candidate_taxa)

    if (method == "PLN_glasso") {
      omega <- fit_pln_glasso(X_counts)
      inhib_vec <- pln_inhibition_vector(omega, "CD", candidate_taxa)
    } else if (method == "LOTO_PLN_glasso") {
      fit <- fit_loto_pln_glasso(X_counts, X_log)
      inhib_vec <- loto_inhibition_vector(fit$coef_mat, "CD", candidate_taxa)
    } else if (method == "LOTO_glmnet_CV1se") {
      fit <- fit_loto_glmnet(X_counts, X_log)
      inhib_vec <- loto_inhibition_vector(fit$coef_mat, "CD", candidate_taxa)
    }
    elapsed <- proc.time()[["elapsed"]] - t0

    pred_df <- truth_strain
    pred_df$method <- method
    pred_df$predicted_inhibition_score <- inhib_vec[pred_df$source_taxon]
    pred_df$predicted_positive <- as.integer(pred_df$predicted_inhibition_score > 0)
    pred_rows[[paste(strain, method, sep = "_")]] <- pred_df

    met <- eval_incoming_edges(pred_df, k = 2L)
    met$strain <- strain
    met$method <- method
    met$n_candidate_edges <- nrow(pred_df)
    met$elapsed_sec <- elapsed
    metric_rows[[paste(strain, method, sep = "_")]] <- met
  }
}

pred_all <- rbindlist(pred_rows, fill = TRUE)
metrics <- rbindlist(metric_rows, fill = TRUE)
summary_df <- metrics[, .(
  mean_precision_at_2 = mean(precision_at_2, na.rm = TRUE),
  mean_recall_at_2 = mean(recall_at_2, na.rm = TRUE),
  mean_f1_at_2 = mean(f1_at_2, na.rm = TRUE),
  mean_average_precision = mean(average_precision, na.rm = TRUE),
  mean_positive_rank = mean(mean_positive_rank, na.rm = TRUE),
  perfect_top2 = sum(top_hits == 2L, na.rm = TRUE),
  strains_scored = .N,
  mean_elapsed_sec = mean(elapsed_sec, na.rm = TRUE)
), by = method]

write_tsv_gz(pred_all, file.path(out_dir, "incoming_edge_predictions.tsv.gz"))
fwrite(metrics, file.path(out_dir, "per_strain_metrics.tsv"), sep = "\t")
fwrite(summary_df, file.path(out_dir, "summary_metrics.tsv"), sep = "\t")

print(summary_df)
