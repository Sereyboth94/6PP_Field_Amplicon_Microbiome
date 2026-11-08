# Run once in RStudio: source("tools/Install_Dependencies.R")
# BiocManager chooses the Bioconductor release appropriate for installed R.
install_amplicon_dependencies <- function(include_legacy = FALSE, include_optional = TRUE) {
  if (getRversion() < "4.3.0") stop("Use R 4.3 or newer; the user's reported R version is 4.4.0.")
  cran <- c("ggplot2", "dplyr", "tidyr", "tibble", "vegan", "patchwork")
  bioc <- "phyloseq"
  if (include_optional) cran <- c(cran, "ragg", "pdftools")
  if (include_legacy) {
    cran <- c(cran, "ape", "RColorBrewer", "forcats"); bioc <- c(bioc, "DESeq2")
  }
  missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) utils::install.packages(missing, repos = "https://cloud.r-project.org")
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
  expected <- c(cran, bioc)
  missing_after <- expected[!vapply(expected, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_after)) stop("Packages still missing: ", paste(missing_after, collapse = ", "))
  message("Dependencies available. Cairo support: ", capabilities("cairo"))
  invisible(data.frame(Package = expected, Version = vapply(expected, function(p) as.character(utils::packageVersion(p)), character(1))))
}
if (!isTRUE(getOption("amplicon.no_autorun", FALSE))) {
  install_amplicon_dependencies(include_legacy = "--legacy" %in% commandArgs(trailingOnly = TRUE))
}
