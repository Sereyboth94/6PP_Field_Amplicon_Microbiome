# Base-R validation only. No package installation or figure generation.
.source_paths <- vapply(sys.frames(), function(fr) if (!is.null(fr$ofile)) as.character(fr$ofile)[1L] else "", character(1))
.source_paths <- .source_paths[nzchar(.source_paths)]
.cli <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.file <- if (length(.source_paths)) tail(.source_paths, 1L) else if (length(.cli)) .cli[1L] else ""
.root <- if (nzchar(.file)) dirname(dirname(normalizePath(.file, winslash = "/", mustWork = TRUE))) else getwd()
source(file.path(.root, "R", "project_utils.R"), local = TRUE, encoding = "UTF-8")
print(amplicon_validate_inputs(.root), row.names = FALSE)
message("All processed inputs passed validation.")
