# Shared base-R project paths and processed-input checks. No automatic execution.
amplicon_datasets <- function() {
  data.frame(
    Dataset = c("Pasture_Bacteria", "Pasture_Fungi", "Wheat_Bacteria", "Wheat_Fungi"),
    Plant = c("Pasture", "Pasture", "Wheat", "Wheat"),
    Marker = c("Bacteria", "Fungi", "Bacteria", "Fungi"),
    Script = c("Pasture/Phyloseq_Bacteria.R", "Pasture/Phyloseq_Fungi.R",
               "Wheat/Wheat_Phyloseq_Bacteria.R", "Wheat/Wheat_Phyloseq_Fungi.R"),
    Metadata = c("Metadat_Bacteria_Field.txt", "Metadat_Fungi_Pasture.txt",
                 "Metadat_Bacteria_Field.txt", "Metadat_Fungal_Wheat.txt"),
    Figure = 1:4, Min_depth = c(500L, 500L, 100L, 500L), stringsAsFactors = FALSE)
}

amplicon_find_root <- function(start = getwd()) {
  override <- getOption("amplicon.repo_root", "")
  if (nzchar(override)) start <- override
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  if (!dir.exists(start)) start <- dirname(start)
  repeat {
    if (file.exists(file.path(start, "R", "project_utils.R")) &&
        file.exists(file.path(start, "Microbiome_Figures_Scripts_and_Tables", "Run_All_Four_Figures.R"))) return(start)
    parent <- dirname(start)
    if (identical(parent, start)) break
    start <- parent
  }
  stop("Repository root not found. Open the .Rproj file or Source the saved runner.", call. = FALSE)
}

amplicon_paths <- function(repo_root, plant, marker) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  datasets <- amplicon_datasets()
  row <- datasets[datasets$Plant == plant & datasets$Marker == marker, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown dataset: ", plant, "/", marker, call. = FALSE)
  data_root <- getOption("amplicon.data_root", Sys.getenv("AMPLICON_DATA_ROOT", repo_root))
  output_root <- getOption("amplicon.output_root", Sys.getenv("AMPLICON_OUTPUT_ROOT", file.path(repo_root, "results")))
  if (!nzchar(data_root)) data_root <- repo_root
  if (!nzchar(output_root)) output_root <- file.path(repo_root, "results")
  base_dir <- file.path(data_root, plant, marker)
  list(repo_root = repo_root, data_root = data_root, output_root = output_root,
       project_dir = file.path(data_root, plant), base_dir = base_dir,
       main_out_dir = file.path(output_root, plant, marker, paste0("Phyloseq_Field_", marker)),
       asv = file.path(base_dir, "ASV_filtered.txt"), taxonomy = file.path(base_dir, "Taxonomy_filtered.txt"),
       metadata = file.path(base_dir, row$Metadata), dataset = row)
}

amplicon_expected_ids <- function(plant) {
  columns <- if (identical(plant, "Pasture")) 5:8 else 1:4
  as.vector(outer(LETTERS[1:8], sprintf("%02d", columns), paste0))
}

amplicon_validate_dataset <- function(paths, check_design = TRUE) {
  label <- paths$dataset$Dataset
  fail <- function(...) stop(label, ": ", ..., call. = FALSE)
  files <- c(paths$asv, paths$taxonomy, paths$metadata)
  absent <- !file.exists(files) | dir.exists(files)
  if (any(absent)) fail("missing input files:\n", paste(files[absent], collapse = "\n"),
                       "\nRun tools/Import_Original_Data.ps1 to copy the original processed tables.")
  if (any(file.info(files)$size <= 0)) fail("an input file is empty.")
  read_input <- function(path) tryCatch(
    utils::read.delim(path, header = TRUE, sep = "\t", check.names = FALSE,
                      stringsAsFactors = FALSE, strip.white = TRUE),
    error = function(e) fail("cannot read ", basename(path), ": ", conditionMessage(e)))
  asv <- read_input(paths$asv)
  tax <- read_input(paths$taxonomy)
  meta <- read_input(paths$metadata)
  if (ncol(asv) < 2L || nrow(asv) < 1L) fail("the ASV table must have an ID column and sample counts.")
  if (ncol(tax) < 2L || nrow(tax) < 1L) fail("the taxonomy table must have an ID column and ranks.")
  if (ncol(meta) < 5L || nrow(meta) < 1L) fail("the metadata table must have sample IDs and four attributes.")
  check_ids <- function(x, description) {
    x <- as.character(x)
    if (anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x)) fail(description, " are missing or duplicated.")
    if (any(x != trimws(x))) fail(description, " contain leading or trailing spaces.")
    x
  }
  asv_ids <- check_ids(asv[[1L]], "ASV IDs")
  tax_ids <- check_ids(tax[[1L]], "Taxonomy ASV IDs")
  sample_ids <- check_ids(meta[[1L]], "Metadata sample IDs")
  asv_samples <- check_ids(colnames(asv)[-1L], "ASV sample columns")
  if (!"Family" %in% colnames(tax)) fail("the taxonomy table must include a Family column.")
  required <- c("Treatment", "Moisture", "Source", "Plant")
  if (!all(required %in% names(meta))) fail("metadata columns required: ", paste(required, collapse = ", "))
  if (!all(sample_ids %in% asv_samples)) fail("metadata samples absent from ASV table: ", paste(setdiff(sample_ids, asv_samples), collapse = ", "))
  counts <- asv[, sample_ids, drop = FALSE]
  if (!all(vapply(counts, is.numeric, logical(1)))) fail("ASV sample counts must be numeric.")
  counts <- as.matrix(counts)
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) fail("ASV counts must be finite, non-negative whole numbers.")
  common <- intersect(asv_ids, tax_ids)
  if (!length(common)) fail("no ASV IDs are shared by abundance and taxonomy tables.")
  treatment <- trimws(as.character(meta$Treatment))
  treatment[treatment %in% c("FVOC", "6-PP")] <- "6PP"
  moisture <- trimws(as.character(meta$Moisture)); source <- trimws(as.character(meta$Source))
  plant <- trimws(as.character(meta$Plant))
  if (anyNA(c(treatment, moisture, source, plant)) || !all(treatment %in% c("Control", "6PP")) ||
      !all(moisture %in% c("Irrigated", "Drought")) || !all(source %in% c("Root", "Rhizosphere")) ||
      !all(plant == paths$dataset$Plant)) fail("unknown or missing treatment, moisture, compartment, or plant values.")
  if (check_design) {
    if (!setequal(sample_ids, amplicon_expected_ids(paths$dataset$Plant))) fail("the expected 32 sample IDs for this plant are not all present.")
    groups <- table(factor(source, c("Root", "Rhizosphere")), factor(moisture, c("Irrigated", "Drought")), factor(treatment, c("Control", "6PP")))
    if (any(groups != 4L)) fail("expected four input samples in every compartment x moisture x treatment group.")
  }
  data.frame(Dataset = label, Samples = length(sample_ids), Rhizosphere = sum(source == "Rhizosphere"),
             Root = sum(source == "Root"), ASVs_in_table = length(asv_ids), ASVs_with_taxonomy = length(common),
             Extra_ASV_sample_columns = length(setdiff(asv_samples, sample_ids)), Status = "PASS", stringsAsFactors = FALSE)
}

amplicon_validate_inputs <- function(repo_root, check_design = TRUE) {
  datasets <- amplicon_datasets(); results <- list(); errors <- character()
  for (i in seq_len(nrow(datasets))) {
    paths <- amplicon_paths(repo_root, datasets$Plant[i], datasets$Marker[i])
    result <- tryCatch(amplicon_validate_dataset(paths, check_design), error = function(e) conditionMessage(e))
    if (is.character(result)) errors <- c(errors, result) else results[[length(results) + 1L]] <- result
  }
  if (length(errors)) stop("Input validation failed before figure generation:\n\n", paste(errors, collapse = "\n\n"), call. = FALSE)
  do.call(rbind, results)
}
