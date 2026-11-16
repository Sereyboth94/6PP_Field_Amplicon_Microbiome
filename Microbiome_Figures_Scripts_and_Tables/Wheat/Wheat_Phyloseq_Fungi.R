# Generate the approved tall figures. Run with source(..., encoding = "UTF-8").
local({
  paths <- vapply(sys.frames(), function(fr) if (!is.null(fr$ofile)) as.character(fr$ofile)[1] else "", character(1))
  paths <- paths[nzchar(paths)]
  cli <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  self <- if (length(paths)) tail(paths, 1) else if (length(cli)) cli[1] else ""
  if (!nzchar(self)) stop("Source this complete script or use Rscript; do not run selected lines.")
  root <- dirname(normalizePath(self, winslash = "/", mustWork = TRUE))
  root <- dirname(root)
  runtime <- new.env(parent = globalenv())
  sys.source(file.path(root, "R", "Publication_Workflow.R"), envir = runtime)
  runtime$run_publication_workflow(root, datasets = "Wheat_Fungi",
    figures = c("main", "supplementary"), replot = FALSE,
    cli_args = commandArgs(trailingOnly = TRUE))
})
