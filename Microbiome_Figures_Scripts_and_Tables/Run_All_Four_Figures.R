# Source this saved file in RStudio, or use Rscript from any working directory.
# Add --check-only on the command line to validate inputs without making figures.
.amp_source_paths <- vapply(sys.frames(), function(fr) {
  if (!is.null(fr$ofile)) as.character(fr$ofile)[1L] else ""
}, character(1))
.amp_source_paths <- .amp_source_paths[nzchar(.amp_source_paths)]
.amp_cli_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.amp_this_file <- if (length(.amp_source_paths)) tail(.amp_source_paths, 1L) else {
  if (length(.amp_cli_file)) .amp_cli_file[1L] else ""
}
.amp_runner_root <- if (nzchar(.amp_this_file)) dirname(dirname(normalizePath(.amp_this_file, winslash = "/", mustWork = TRUE))) else getwd()
source(file.path(.amp_runner_root, "R", "project_utils.R"), local = TRUE, encoding = "UTF-8")

run_all_four_figures <- function(check_only = FALSE, repo_root = .amp_runner_root, data_root = NULL, output_root = NULL) {
  repo_root <- amplicon_find_root(repo_root)
  old_options <- options(amplicon.repo_root = repo_root)
  on.exit(options(old_options), add = TRUE)
  if (!is.null(data_root)) {
    previous <- options(amplicon.data_root = normalizePath(data_root, winslash = "/", mustWork = TRUE))
    on.exit(options(previous), add = TRUE)
  }
  if (!is.null(output_root)) {
    previous_output <- options(amplicon.output_root = output_root)
    on.exit(options(previous_output), add = TRUE)
  }
  datasets <- amplicon_datasets()
  scripts <- file.path(repo_root, "Microbiome_Figures_Scripts_and_Tables", datasets$Script)
  if (any(!file.exists(scripts))) stop("Missing analysis scripts:\n", paste(scripts[!file.exists(scripts)], collapse = "\n"))
  checks <- amplicon_validate_inputs(repo_root)
  print(checks, row.names = FALSE)
  if (isTRUE(check_only)) {
    message("Input checks passed for all four datasets. No figures were generated.")
    return(invisible(checks))
  }
  required <- c("phyloseq", "ggplot2", "dplyr", "tidyr", "tibble", "vegan", "patchwork")
  if (isTRUE(getOption("amplicon.run_legacy", FALSE))) required <- c(required, "ape", "DESeq2", "RColorBrewer", "forcats")
  absent <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing packages: ", paste(absent, collapse = ", "), ".\nSource tools/Install_Dependencies.R once, then run this file again.", call. = FALSE)
  if (!isTRUE(capabilities("cairo"))) stop("Cairo graphics support is required for the requested PDF exports.")
  out <- amplicon_paths(repo_root, "Pasture", "Bacteria")$output_root
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = out)
  if (!dir.create(log_dir)) stop("Cannot create run log folder: ", log_dir)
  utils::write.csv(checks, file.path(log_dir, "Input_checks.csv"), row.names = FALSE)
  inputs <- unlist(lapply(seq_len(nrow(datasets)), function(i) {
    p <- amplicon_paths(repo_root, datasets$Plant[i], datasets$Marker[i]); c(p$asv, p$taxonomy, p$metadata)
  }), use.names = FALSE)
  utils::write.csv(data.frame(File = inputs, MD5 = unname(tools::md5sum(inputs))), file.path(log_dir, "Input_MD5.csv"), row.names = FALSE)
  completed <- character()
  for (i in seq_len(nrow(datasets))) {
    message("Running ", datasets$Dataset[i], " (", i, "/4)")
    analysis_env <- new.env(parent = globalenv())
    tryCatch(source(scripts[i], local = analysis_env, chdir = FALSE, encoding = "UTF-8"), error = function(e) {
      writeLines(c(paste("Dataset:", datasets$Dataset[i]), conditionMessage(e)), file.path(log_dir, "ERROR.txt"))
      stop("Analysis stopped at ", datasets$Dataset[i], ": ", conditionMessage(e), "\nEarlier completed datasets: ",
           paste(completed, collapse = ", "), "\nRun log: ", log_dir, call. = FALSE)
    })
    completed <- c(completed, datasets$Dataset[i])
    writeLines(completed, file.path(log_dir, "Completed_datasets.txt"))
    rm(analysis_env); invisible(gc())
  }
  writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
  message("All four datasets finished. Outputs: ", out)
  invisible(list(checks = checks, completed = completed, output_root = out, log_dir = log_dir))
}
if (!isTRUE(getOption("amplicon.no_autorun", FALSE))) {
  run_all_four_figures(check_only = "--check-only" %in% commandArgs(trailingOnly = TRUE))
}
