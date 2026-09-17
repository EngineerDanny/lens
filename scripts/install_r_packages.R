#!/usr/bin/env Rscript

dir.create("r_library", recursive = TRUE, showWarnings = FALSE)
local_lib <- normalizePath("r_library", mustWork = TRUE)
.libPaths(c(local_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c("remotes", "pulsar", "huge", "PLNmodels", "glmnet", "data.table", "ggplot2", "igraph")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, lib = local_lib, Ncpus = 4)
}

if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
  remotes::install_github(
    "zdk123/SpiecEasi@03b96da",
    lib = local_lib,
    dependencies = TRUE,
    upgrade = "never",
    build_vignettes = FALSE
  )
}

if (!requireNamespace("SPRING", quietly = TRUE)) {
  remotes::install_github(
    "GraceYoon/SPRING@3d641a4b",
    lib = local_lib,
    dependencies = NA,
    upgrade = "never",
    build_vignettes = FALSE
  )
}

# OneNet brings additional estimators and is optional when using saved scores.
if (Sys.getenv("INSTALL_ONENET", "false") == "true" &&
    !requireNamespace("OneNet", quietly = TRUE)) {
  remotes::install_github(
    "metagenopolis/OneNet@d969e9df41139199619e0fe43be368724292403d",
    lib = local_lib,
    dependencies = TRUE,
    upgrade = "never",
    build_vignettes = FALSE
  )
}

used <- c("PLNmodels", "glmnet", "SpiecEasi", "SPRING", "pulsar", "huge")
versions <- data.frame(
  package = used,
  version = vapply(used, function(package) as.character(packageVersion(package)), character(1)),
  stringsAsFactors = FALSE
)
dir.create("analysis_data", recursive = TRUE, showWarnings = FALSE)
write.csv(versions, "analysis_data/package_versions.csv", row.names = FALSE, quote = TRUE)
