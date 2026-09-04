#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: evaluate_network_predictions.R <dataset> <prediction_file>")
}

dataset <- args[1]
prediction_file <- normalizePath(args[2], mustWork = TRUE)
base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth", mustWork = TRUE)

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

sort_pair <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

infer_prediction <- function(df) {
  nm <- names(df)
  if (all(c("source_taxon", "target_taxon") %in% nm)) {
    df$mode <- "directed"
    df$source <- df$source_taxon
    df$target <- df$target_taxon
  } else if (all(c("predictor_taxon", "response_taxon") %in% nm)) {
    df$mode <- "directed"
    df$source <- df$predictor_taxon
    df$target <- df$response_taxon
  } else if (all(c("focal_taxon", "partner_taxon") %in% nm)) {
    df$mode <- "undirected"
    df$a <- df$focal_taxon
    df$b <- df$partner_taxon
  } else if (all(c("taxon_i", "taxon_j") %in% nm)) {
    df$mode <- "undirected"
    df$a <- df$taxon_i
    df$b <- df$taxon_j
  } else if (all(c("taxon_1", "taxon_2") %in% nm)) {
    df$mode <- "undirected"
    df$a <- df$taxon_1
    df$b <- df$taxon_2
  } else {
    stop("Cannot infer prediction format from columns: ", paste(nm, collapse = ", "))
  }

  sign_col <- intersect(c("effect_sign", "sign", "sign_consensus", "edge_sign"), names(df))
  if (length(sign_col)) {
    df$pred_sign <- suppressWarnings(as.integer(df[[sign_col[1]]]))
  } else {
    df$pred_sign <- NA_integer_
  }

  if (unique(df$mode) == "directed") {
    df$key <- paste(df$source, df$target, sep = "||")
    unique(df[, c("key", "source", "target", "pred_sign")])
  } else {
    df$key <- sort_pair(df$a, df$b)
    nodes <- do.call(rbind, strsplit(df$key, "\\|\\|"))
    out <- unique(df[, c("key", "pred_sign")])
    parts <- do.call(rbind, strsplit(out$key, "\\|\\|"))
    out$taxon_1 <- parts[, 1]
    out$taxon_2 <- parts[, 2]
    out
  }
}

compute_metrics <- function(truth_df, pred_df, mode) {
  truth_keys <- unique(truth_df$key)
  pred_keys <- unique(pred_df$key)
  tp_keys <- intersect(truth_keys, pred_keys)
  fp_keys <- setdiff(pred_keys, truth_keys)
  fn_keys <- setdiff(truth_keys, pred_keys)

  tp <- length(tp_keys)
  fp <- length(fp_keys)
  fn <- length(fn_keys)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1 <- if (is.finite(precision + recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_

  sign_accuracy <- NA_real_
  signed_tp <- NA_integer_
  truth_signable <- truth_df[!is.na(truth_df$truth_sign), c("key", "truth_sign"), drop = FALSE]
  pred_signable <- pred_df[!is.na(pred_df$pred_sign), c("key", "pred_sign"), drop = FALSE]
  if (nrow(truth_signable) > 0 && nrow(pred_signable) > 0) {
    merged <- merge(truth_signable, pred_signable, by = "key")
    if (nrow(merged) > 0) {
      sign_accuracy <- mean(merged$truth_sign == merged$pred_sign)
      signed_tp <- sum(merged$truth_sign == merged$pred_sign)
    }
  }

  data.frame(
    dataset = dataset,
    mode = mode,
    n_truth = length(truth_keys),
    n_pred = length(pred_keys),
    tp = tp,
    fp = fp,
    fn = fn,
    precision = precision,
    recall = recall,
    f1 = f1,
    signed_tp = signed_tp,
    sign_accuracy = sign_accuracy,
    stringsAsFactors = FALSE
  )
}

truth_dir <- file.path(base_dir, dataset, "processed")
truth_directed <- read_tsv_gz(file.path(truth_dir, "truth_directed.tsv.gz"))
truth_undirected <- read_tsv_gz(file.path(truth_dir, "truth_undirected.tsv.gz"))

truth_directed$key <- paste(truth_directed$source_taxon, truth_directed$target_taxon, sep = "||")
truth_directed$truth_sign <- suppressWarnings(as.integer(truth_directed$effect_sign))
truth_undirected$key <- sort_pair(truth_undirected$taxon_1, truth_undirected$taxon_2)
truth_undirected$truth_sign <- suppressWarnings(as.integer(truth_undirected$sign_consensus))

pred_path_lower <- tolower(prediction_file)
pred <- if (grepl("\\.gz$", pred_path_lower)) read_tsv_gz(prediction_file) else {
  read.delim(prediction_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}
pred_std <- infer_prediction(pred)

mode <- if ("source" %in% names(pred_std)) "directed" else "undirected"
metrics <- if (mode == "directed") {
  compute_metrics(truth_directed, pred_std, mode)
} else {
  compute_metrics(truth_undirected, pred_std, mode)
}

out_dir <- file.path(truth_dir, "evaluation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(prediction_file)), "_metrics.tsv"))
write_tsv(metrics, out_file)
print(metrics)
