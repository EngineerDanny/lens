#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

base_dir <- normalizePath("/projects/genomic-ml/da2343/PLN/pln_eval/data/interaction_ground_truth", mustWork = TRUE)

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
}

clean_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s*\\([^\\)]*\\)$", "", x)
  x <- gsub("^'+|'+$", "", x)
  trimws(x)
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

extract_pairinterax <- function() {
  src_dir <- file.path(base_dir, "pairinterax")
  out_dir <- file.path(src_dir, "processed")
  dir_create(out_dir)

  interactions <- suppressMessages(read_excel(file.path(src_dir, "pairinterax_interactions.xlsx"), skip = 2))
  prevalence <- suppressMessages(read_excel(file.path(src_dir, "pairinterax_mapping_prevalence.xlsx"), skip = 2))

  interactions <- as.data.frame(interactions, stringsAsFactors = FALSE)
  prevalence <- as.data.frame(prevalence, stringsAsFactors = FALSE)

  interactions$taxon_a <- interactions$A
  interactions$taxon_b <- interactions$B
  interactions$canonical_a <- clean_name(interactions$taxon_a)
  interactions$canonical_b <- clean_name(interactions$taxon_b)
  interactions$phenotype_a <- as.numeric(interactions$`Phenotype A`)
  interactions$phenotype_b <- as.numeric(interactions$`Phenotype B`)
  interactions$interaction_label <- interactions$Interaction
  interactions_out <- interactions[, c(
    "taxon_a", "taxon_b", "canonical_a", "canonical_b",
    "phenotype_a", "phenotype_b", "interaction_label"
  )]

  prevalence <- prevalence[!is.na(prevalence$Taxa) & prevalence$Taxa != "Enterotypes", , drop = FALSE]
  prevalence$canonical_taxon <- clean_name(prevalence$Taxa)
  sample_cols <- names(prevalence)[seq.int(5L, ncol(prevalence))]
  abundance_mat <- prevalence[, sample_cols, drop = FALSE]
  abundance_mat[] <- lapply(abundance_mat, function(x) suppressWarnings(as.numeric(x)))
  abundance_matrix <- data.frame(
    sample_id = sample_cols,
    t(as.matrix(abundance_mat)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(abundance_matrix)[-1] <- prevalence$canonical_taxon

  taxa_sources <- unique(c(prevalence$Taxa, interactions$taxon_a, interactions$taxon_b))
  taxa_map <- data.frame(
    original_taxon = taxa_sources,
    canonical_taxon = clean_name(taxa_sources),
    stringsAsFactors = FALSE
  )
  taxa_map$in_abundance_table <- taxa_map$canonical_taxon %in% prevalence$canonical_taxon
  taxa_map$in_interaction_table <- taxa_map$canonical_taxon %in% c(interactions$canonical_a, interactions$canonical_b)

  write_tsv_gz(interactions_out, file.path(out_dir, "truth_edges.tsv.gz"))
  write_tsv_gz(abundance_matrix, file.path(out_dir, "abundance_matrix.tsv.gz"))
  write_tsv_gz(taxa_map, file.path(out_dir, "taxa_map.tsv.gz"))
}

extract_omm12 <- function() {
  src_dir <- file.path(base_dir, "omm12")
  out_dir <- file.path(src_dir, "processed")
  dir_create(out_dir)

  wb <- file.path(src_dir, "41396_2021_1153_moesm2_esm.xlsx")

  taxa_map <- data.frame(
    taxon_id = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58"),
    full_name = c(
      "Enterococcus faecalis",
      "Bifidobacterium animalis",
      "Acutalibacter muris",
      "Muribaculum intestinale",
      "Flavonifractor plautii",
      "Enterocloster clostridioformis",
      "Akkermansia muciniphila",
      "Turicimonas muris",
      "Clostridium innocuum",
      "Bacteroides caecimuris",
      "Limosilactobacillus reuteri",
      "Blautia coccoides"
    ),
    mapping_note = c(
      rep("from model archive filenames", 2),
      "species inferred from missing code in model archive plus workbook code KB18",
      rep("from model archive filenames", 9)
    ),
    stringsAsFactors = FALSE
  )

  abundance_in_vitro <- as.data.frame(suppressMessages(read_excel(wb, sheet = "community absabundance in vitro")), stringsAsFactors = FALSE)
  abundance_in_vitro <- abundance_in_vitro[, c("sample", taxa_map$taxon_id), drop = FALSE]
  names(abundance_in_vitro)[1] <- "sample_id"

  abundance_adult <- as.data.frame(suppressMessages(read_excel(wb, sheet = "in vivo adult mice")), stringsAsFactors = FALSE)
  names(abundance_adult)[1] <- "sample_id"
  abundance_adult <- abundance_adult[, c("sample_id", taxa_map$taxon_id), drop = FALSE]

  abundance_infant <- as.data.frame(suppressMessages(read_excel(wb, sheet = "in vivo infant mice")), stringsAsFactors = FALSE)
  names(abundance_infant)[1] <- "sample_id"
  abundance_infant <- abundance_infant[, c("sample_id", taxa_map$taxon_id), drop = FALSE]

  rbm <- as.data.frame(suppressMessages(read_excel(wb, sheet = "rbm 72h mean")), stringsAsFactors = FALSE)
  rbm <- rbm[rbm$group == "co" & !grepl(" Mono$", rbm$coku), c("coku", "probe", "mean_rbm", "sd"), drop = FALSE]

  decode_coku <- function(x, codes) {
    for (code in codes[order(nchar(codes), decreasing = TRUE)]) {
      if (startsWith(x, code)) {
        rest <- substring(x, nchar(code) + 1)
        if (rest %in% codes) return(c(code, rest))
      }
    }
    c(NA_character_, NA_character_)
  }

  pairs <- t(vapply(rbm$coku, decode_coku, character(2), codes = taxa_map$taxon_id))
  rbm$taxon_1 <- pairs[, 1]
  rbm$taxon_2 <- pairs[, 2]
  rbm$partner_taxon <- ifelse(rbm$probe == rbm$taxon_1, rbm$taxon_2,
                              ifelse(rbm$probe == rbm$taxon_2, rbm$taxon_1, NA_character_))
  rbm$effect_sign <- ifelse(rbm$mean_rbm > 1, 1L, ifelse(rbm$mean_rbm < 1, -1L, 0L))
  rbm$log2_rbm <- log2(rbm$mean_rbm)
  truth_edges <- rbm[, c("coku", "probe", "partner_taxon", "mean_rbm", "sd", "log2_rbm", "effect_sign")]
  names(truth_edges)[2] <- "focal_taxon"

  write_tsv_gz(abundance_in_vitro, file.path(out_dir, "community_absabundance_in_vitro.tsv.gz"))
  write_tsv_gz(abundance_adult, file.path(out_dir, "in_vivo_adult_mice.tsv.gz"))
  write_tsv_gz(abundance_infant, file.path(out_dir, "in_vivo_infant_mice.tsv.gz"))
  write_tsv_gz(truth_edges, file.path(out_dir, "truth_edges.tsv.gz"))
  write_tsv_gz(taxa_map, file.path(out_dir, "taxa_map.tsv.gz"))
}

extract_omm12_keystone_2023 <- function() {
  src_dir <- file.path(base_dir, "omm12_keystone_2023")
  out_dir <- file.path(src_dir, "processed")
  dir_create(out_dir)

  wb <- file.path(src_dir, "source-data_table_1_revision.xlsx")

  taxa_map <- data.frame(
    taxon_id = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58"),
    full_name = c(
      "Enterococcus faecalis",
      "Bifidobacterium animalis",
      "Acutalibacter muris",
      "Muribaculum intestinale",
      "Flavonifractor plautii",
      "Enterocloster clostridioformis",
      "Akkermansia muciniphila",
      "Turicimonas muris",
      "Clostridium innocuum",
      "Bacteroides caecimuris",
      "Limosilactobacillus reuteri",
      "Blautia coccoides"
    ),
    mapping_note = "Reused OMM12 code-to-species mapping",
    stringsAsFactors = FALSE
  )

  abundance_in_vitro <- as.data.frame(suppressMessages(read_excel(wb, sheet = "Tab1")), stringsAsFactors = FALSE)
  abundance_in_vitro <- abundance_in_vitro[, c("16Scorr", taxa_map$taxon_id), drop = FALSE]
  names(abundance_in_vitro)[1] <- "sample_id"

  abundance_in_vivo <- as.data.frame(suppressMessages(read_excel(wb, sheet = "Tab13")), stringsAsFactors = FALSE)
  abundance_in_vivo <- abundance_in_vivo[, c("absab_16Scorr", taxa_map$taxon_id), drop = FALSE]
  names(abundance_in_vivo)[1] <- "sample_id"

  build_truth_from_ratio <- function(df, context_col, deletion_prefix) {
    ratio_cols <- setdiff(names(df), c(context_col, "probe"))
    rows <- lapply(ratio_cols, function(col) {
      out <- data.frame(
        context = df[[context_col]],
        focal_taxon = df$probe,
        partner_taxon = sub(paste0("^", deletion_prefix), "", col),
        ratio = as.numeric(df[[col]]),
        stringsAsFactors = FALSE
      )
      out$effect_sign <- ifelse(out$ratio > 1, 1L, ifelse(out$ratio < 1, -1L, 0L))
      out$log2_ratio <- log2(out$ratio)
      out
    })
    out <- do.call(rbind, rows)
    out[!is.na(out$effect_sign), , drop = FALSE]
  }

  truth_in_vitro <- build_truth_from_ratio(
    as.data.frame(suppressMessages(read_excel(wb, sheet = "Tab5")), stringsAsFactors = FALSE),
    context_col = "medium",
    deletion_prefix = "OMM11-"
  )

  truth_in_vivo <- build_truth_from_ratio(
    as.data.frame(suppressMessages(read_excel(wb, sheet = "Tab15")), stringsAsFactors = FALSE),
    context_col = "region",
    deletion_prefix = "OMM-"
  )

  write_tsv_gz(abundance_in_vitro, file.path(out_dir, "community_absabundance_in_vitro.tsv.gz"))
  write_tsv_gz(abundance_in_vivo, file.path(out_dir, "community_absabundance_in_vivo.tsv.gz"))
  write_tsv_gz(truth_in_vitro, file.path(out_dir, "truth_edges_in_vitro.tsv.gz"))
  write_tsv_gz(truth_in_vivo, file.path(out_dir, "truth_edges_in_vivo.tsv.gz"))
  write_tsv_gz(taxa_map, file.path(out_dir, "taxa_map.tsv.gz"))
}

extract_venturelli <- function() {
  src_dir <- file.path(base_dir, "venturelli_2018")
  out_dir <- file.path(src_dir, "processed")
  dir_create(out_dir)

  parse_pairwise_sheet <- function(path, sheet) {
    df <- as.data.frame(suppressMessages(read_excel(path, sheet = sheet, col_names = FALSE)), stringsAsFactors = FALSE)
    names(df) <- c("community_code", paste0("t", seq_len(ncol(df) - 1L)))
    df$community_code <- gsub("[^A-Z]", "", trimws(df$community_code))
    df <- df[nzchar(df$community_code), , drop = FALSE]
    df$focal_taxon <- substr(df$community_code, 1, 2)
    df$partner_taxon <- substr(df$community_code, 3, 4)

    long <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
      vals <- unlist(df[i, grep("^t", names(df)), drop = TRUE], use.names = FALSE)
      data.frame(
        community_code = df$community_code[i],
        focal_taxon = df$focal_taxon[i],
        partner_taxon = df$partner_taxon[i],
        time_index = seq_along(vals),
        abundance = suppressWarnings(as.numeric(vals)),
        stringsAsFactors = FALSE
      )
    }))
    long
  }

  ev1 <- parse_pairwise_sheet(file.path(src_dir, "dataset_EV1.xlsx"), "PW1")
  ev2 <- parse_pairwise_sheet(file.path(src_dir, "dataset_EV2.xlsx"), "Dataset EV2")

  zip_path <- file.path(src_dir, "code_EV1.zip")
  zip_list <- unzip(zip_path, list = TRUE)
  xml_files <- zip_list$Name[grepl("\\.xml$", zip_list$Name)]

  parse_xml_truth <- function(xml_file) {
    lines <- unzip(zip_path, files = xml_file, exdir = tempdir(), overwrite = TRUE)
    txt <- paste(readLines(lines, warn = FALSE), collapse = "")
    matches <- gregexpr('<parameter[^>]*id="a([A-Z]{2})([A-Z]{2})"[^>]*value="([^"]+)"', txt, perl = TRUE)
    m <- regmatches(txt, matches)[[1]]
    if (!length(m)) return(data.frame())
    parts <- do.call(rbind, lapply(m, function(x) {
      parsed <- regmatches(x, regexec('id="a([A-Z]{2})([A-Z]{2})"[^>]*value="([^"]+)"', x, perl = TRUE))[[1]]
      data.frame(
        response_taxon = parsed[2],
        predictor_taxon = parsed[3],
        coefficient = as.numeric(parsed[4]),
        stringsAsFactors = FALSE
      )
    }))
    parts$model_id <- tools::file_path_sans_ext(basename(xml_file))
    parts$effect_sign <- ifelse(parts$coefficient > 0, 1L, ifelse(parts$coefficient < 0, -1L, 0L))
    parts
  }

  truth_all <- do.call(rbind, lapply(xml_files, parse_xml_truth))
  truth_nonzero <- truth_all[truth_all$coefficient != 0, , drop = FALSE]
  taxon_ids <- sort(unique(c(truth_all$response_taxon, truth_all$predictor_taxon, ev1$focal_taxon, ev1$partner_taxon, ev2$focal_taxon, ev2$partner_taxon)))
  taxa_map <- data.frame(
    taxon_id = taxon_ids,
    taxon_label = taxon_ids,
    full_name = NA_character_,
    mapping_note = "full species name not recovered from supplementary files; article uses abbreviations",
    stringsAsFactors = FALSE
  )

  write_tsv_gz(ev1, file.path(out_dir, "pairwise_time_series_EV1.tsv.gz"))
  write_tsv_gz(ev2, file.path(out_dir, "pairwise_time_series_EV2.tsv.gz"))
  write_tsv_gz(truth_nonzero, file.path(out_dir, "truth_edges.tsv.gz"))
  write_tsv_gz(taxa_map, file.path(out_dir, "taxa_map.tsv.gz"))
}

targets <- if (length(args)) args else c("pairinterax", "omm12", "omm12_keystone_2023")
if ("all" %in% targets) targets <- c("pairinterax", "omm12", "omm12_keystone_2023")

for (target in targets) {
  message("Extracting ", target, " ...")
  switch(
    target,
    pairinterax = extract_pairinterax(),
    omm12 = extract_omm12(),
    omm12_keystone_2023 = extract_omm12_keystone_2023(),
    venturelli_2018 = extract_venturelli(),
    stop("Unknown target: ", target)
  )
}

message("Done.")
