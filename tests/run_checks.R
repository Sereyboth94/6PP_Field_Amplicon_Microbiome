# Base-R checks; temporary test tables never become analysis inputs or figures.
.paths <- vapply(sys.frames(), function(fr) if (!is.null(fr$ofile)) as.character(fr$ofile)[1L] else "", character(1))
.paths <- .paths[nzchar(.paths)]
.cli <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.this <- if (length(.paths)) tail(.paths, 1L) else if (length(.cli)) .cli[1L] else ""
.repo <- if (nzchar(.this)) dirname(dirname(normalizePath(.this, winslash = "/", mustWork = TRUE))) else getwd()

run_repository_checks <- function(root) {
  files <- list.files(root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  for (path in files) parse(path, encoding = "UTF-8")
  message("PASS: all R files parse.")
  source(file.path(root, "R", "project_utils.R"), local = TRUE, encoding = "UTF-8")
  datasets <- amplicon_datasets()
  fixture <- tempfile("amplicon repo test with spaces ")
  dir.create(fixture)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  prior <- options(amplicon.repo_root = NULL, amplicon.data_root = NULL,
                   amplicon.output_root = NULL, amplicon.no_autorun = TRUE)
  on.exit(options(prior), add = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  copy_relative <- function(relative) {
    target <- file.path(fixture, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(root, relative), target))
  }
  copy_relative("R/project_utils.R")
  copy_relative("Microbiome_Figures_Scripts_and_Tables/Run_All_Four_Figures.R")
  for (i in seq_len(nrow(datasets))) {
    copy_relative(file.path("Microbiome_Figures_Scripts_and_Tables", datasets$Script[i]))
    relative <- file.path(datasets$Plant[i], datasets$Marker[i], datasets$Metadata[i])
    copy_relative(relative)
    md <- utils::read.delim(file.path(fixture, relative), check.names = FALSE, stringsAsFactors = FALSE)
    ids <- md[[1L]]
    counts <- as.data.frame(matrix(1000L, nrow = 3L, ncol = length(ids)), check.names = FALSE)
    colnames(counts) <- ids
    asv <- data.frame(ASV_ID = paste0("TEST_ASV_", 1:3), counts, check.names = FALSE)
    tax <- data.frame(ASV_ID = asv$ASV_ID, Family = paste0("TEST_FAMILY_", 1:3))
    folder <- dirname(file.path(fixture, relative))
    utils::write.table(asv, file.path(folder, "ASV_filtered.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(tax, file.path(folder, "Taxonomy_filtered.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  validated <- amplicon_validate_inputs(fixture)
  stopifnot(nrow(validated) == 4L, all(validated$Samples == 32L),
            all(validated$Rhizosphere == 16L), all(validated$Root == 16L), all(validated$Status == "PASS"))
  message("PASS: four supplied metadata tables and temporary ASV/taxonomy fixtures match the expected design.")
  setwd(tempdir())
  runner <- new.env(parent = globalenv())
  source(file.path(fixture, "Microbiome_Figures_Scripts_and_Tables", "Run_All_Four_Figures.R"), local = runner)
  result <- runner$run_all_four_figures(check_only = TRUE)
  stopifnot(nrow(result) == 4L, !dir.exists(file.path(fixture, "results")))
  message("PASS: runner finds a repository with spaces from an unrelated working directory; check-only creates no figures.")
  for (i in seq_len(nrow(datasets))) {
    original_script <- file.path(fixture, "Microbiome_Figures_Scripts_and_Tables", datasets$Script[i])
    lines <- readLines(original_script, warn = FALSE)
    begin <- grep("^# Locate this saved script", lines)
    end <- grep("^base_dir <- \\.amp_paths\\$base_dir$", lines)
    stopifnot(length(begin) == 1L, length(end) == 1L)
    probe <- file.path(dirname(original_script), "temporary_path_probe.R")
    writeLines(c(lines[seq.int(begin, end)], "main_out_dir <- .amp_paths$main_out_dir"), probe)
    env <- new.env(parent = globalenv())
    source(probe, local = env)
    expected <- amplicon_paths(fixture, datasets$Plant[i], datasets$Marker[i])
    stopifnot(identical(env$base_dir, expected$base_dir), identical(env$main_out_dir, expected$main_out_dir))
    unlink(probe)
  }
  message("PASS: all four standalone path blocks locate the correct input and output folders.")
  expect_error <- function(expression, text) {
    error <- tryCatch({force(expression); NULL}, error = function(e) conditionMessage(e))
    if (is.null(error) || !grepl(text, error, fixed = TRUE)) stop("Expected validation failure was not observed: ", text)
  }
  paths <- amplicon_paths(fixture, "Pasture", "Bacteria")
  abundance <- readLines(paths$asv)
  file.rename(paths$asv, paste0(paths$asv, ".temporarily_missing"))
  expect_error(amplicon_validate_inputs(fixture), "missing input files")
  file.rename(paste0(paths$asv, ".temporarily_missing"), paths$asv)
  message("PASS: missing inputs stop the run before figure generation.")
  asv <- utils::read.delim(paths$asv, check.names = FALSE)
  bad <- asv; bad[1L, 2L] <- -1
  utils::write.table(bad, paths$asv, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(amplicon_validate_inputs(fixture), "non-negative whole numbers")
  bad <- asv; names(bad)[2L] <- "MISSING_SAMPLE"
  utils::write.table(bad, paths$asv, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(amplicon_validate_inputs(fixture), "metadata samples absent")
  writeLines(abundance, paths$asv)
  message("PASS: invalid counts and missing metadata samples are rejected.")
  md <- utils::read.delim(paths$metadata, check.names = FALSE, stringsAsFactors = FALSE)
  metadata_lines <- readLines(paths$metadata)
  bad_meta <- rbind(md, md[1L, , drop = FALSE])
  utils::write.table(bad_meta, paths$metadata, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(amplicon_validate_inputs(fixture), "Metadata sample IDs are missing or duplicated")
  bad_meta <- md; bad_meta$Treatment[1L] <- "UNKNOWN"
  utils::write.table(bad_meta, paths$metadata, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(amplicon_validate_inputs(fixture), "unknown or missing treatment")
  writeLines(metadata_lines, paths$metadata)
  message("PASS: duplicated IDs and invalid experimental labels are rejected.")
  taxonomy_lines <- readLines(paths$taxonomy)
  writeLines(sub("Family", "NotFamily", taxonomy_lines, fixed = TRUE), paths$taxonomy)
  expect_error(amplicon_validate_inputs(fixture), "Family column")
  writeLines(taxonomy_lines, paths$taxonomy)
  stopifnot(all(amplicon_validate_inputs(fixture)$Status == "PASS"))
  message("PASS: missing family taxonomy is rejected; valid inputs pass again after restoration.")
  message("All repository checks passed. These checks do not recompute manuscript results.")
  invisible(TRUE)
}
run_repository_checks(.repo)
