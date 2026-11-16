# Shared publication workflow; source an entry script, not this module.
run_publication_workflow <- function(script_root, datasets = NULL,
                                     figures = c("main", "supplementary"),
                                     replot = FALSE, cli_args = character()) {
  required <- c("ggplot2", "dplyr", "tidyr", "patchwork", "vegan", "ape")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
    ". Install them with install.packages(c(", paste(sprintf('"%s"', missing), collapse = ", "), ")).")
  `%>%` <- getExportedValue("dplyr", "%>%")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  script_root <- normalizePath(script_root, winslash = "/", mustWork = TRUE)
  source(file.path(script_root, "Figure_Settings.R"), local = TRUE, encoding = "UTF-8")
  extra_settings <- getOption("amplicon.settings_file", NULL)
  if (!is.null(extra_settings)) source(extra_settings, local = TRUE, encoding = "UTF-8")
  positive <- c(publication_dpi, publication_main_width, publication_main_height,
    publication_pie_width, publication_pie_height, publication_base_size,
    publication_axis_size, publication_main_family_size, publication_pie_family_size)
  if (any(!is.finite(positive)) || any(positive <= 0)) stop("Figure sizes and fonts must be positive finite numbers.")
  if (!publication_raster_backend %in% c("cairo", "ragg")) stop("Choose cairo or ragg for the raster backend.")
  if (!isTRUE(capabilities("cairo"))) stop("This R installation needs Cairo support for embedded-font PDFs.")
  if (publication_raster_backend == "ragg" && !requireNamespace("ragg", quietly = TRUE)) stop("Install ragg or use the default cairo backend.")
  if (isTRUE(publication_export_svg) && !requireNamespace("svglite", quietly = TRUE)) stop("Install svglite or set publication_export_svg <- FALSE.")
  if (!publication_pie_order %in% c("within_panel", "compartment")) stop("Invalid pie ordering setting.")
  if (any(c(pie_top_n, stacked_top_n) < 1) || any(c(pie_top_n, stacked_top_n) %% 1 != 0)) stop("Family numbers must be positive integers.")
  if (publication_pie_label_min < 0 || publication_pie_label_min > 100) stop("Pie label threshold must be between 0 and 100.")
  if (length(publication_main_column_widths) != 3 || length(publication_main_row_heights) != 4 ||
      length(publication_pie_panel_heights) != 2 || any(c(publication_main_column_widths,
      publication_main_row_heights, publication_pie_panel_heights) <= 0)) stop("Invalid panel layout weights.")
  source(file.path(script_root, "R", "figure_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(script_root, "R", "Figure_Layouts.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(script_root, "R", "Compute_Study.R"), local = TRUE, encoding = "UTF-8")
  slugs <- c("Pasture_Bacteria", "Pasture_Fungi", "Wheat_Bacteria", "Wheat_Fungi")
  if (is.null(datasets)) datasets <- slugs
  if (!length(datasets) || any(!datasets %in% slugs)) stop("Unknown dataset selection.")
  if (!length(figures) || any(!figures %in% c("main", "supplementary"))) stop("Unknown figure selection.")
  input_root <- if (length(cli_args) >= 1) cli_args[1] else getOption("amplicon.input_root", dirname(script_root))
  output_parent <- if (length(cli_args) >= 2) cli_args[2] else getOption("amplicon.output_root", file.path(script_root, "Tall_Figure_Outputs"))
  input_root <- normalizePath(input_root, winslash = "/", mustWork = TRUE)
  metadata_files <- c("Metadat_Bacteria_Field.txt", "Metadat_Fungi_Pasture.txt",
    "Metadat_Bacteria_Field.txt", "Metadat_Fungal_Wheat.txt")
  relative_inputs <- unlist(lapply(seq_along(slugs), function(i) file.path(
    sub("_", "/", slugs[i]), c("ASV_filtered.txt", "Taxonomy_filtered.txt", metadata_files[i]))))
  input_files <- file.path(input_root, relative_inputs)
  if (any(!file.exists(input_files))) stop("Missing processed input tables:\n",
    paste(input_files[!file.exists(input_files)], collapse = "\n"),
    "\nExpected the Pasture and Wheat data folders beside Microbiome_Figures_Scripts_and_Tables, or set options(amplicon.input_root = 'path/to/github').")
  input_md5 <- unname(tools::md5sum(input_files));names(input_md5) <- relative_inputs
  if (anyNA(input_md5)) stop("An input table is unreadable; ensure OneDrive files are available locally.")
  dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
  output_parent <- normalizePath(output_parent, winslash = "/", mustWork = TRUE)
  cache_file <- getOption("amplicon.cache_file", NULL)
  if (isTRUE(replot) && is.null(cache_file)) {
    latest <- file.path(output_parent, "LATEST_ANALYSIS_RUN.txt")
    if (!file.exists(latest)) stop("Run Run_All_Four_Figures.R once before Replot_Figures.R.")
    previous <- readLines(latest, warn = FALSE, n = 1)
    cache_file <- file.path(output_parent, previous, "Cache", "study_analyses.rds")
  }
  if (!is.null(cache_file)) {
    cache_file <- normalizePath(cache_file, winslash = "/", mustWork = TRUE)
    cache <- readRDS(cache_file)
    if (!is.list(cache) || !identical(cache$format, 1L) || !identical(names(cache$objects), slugs)) stop("Incompatible analysis cache; run the full workflow again.")
    if (!identical(cache$input_md5, input_md5)) stop("Input tables differ from the cached analysis. Run Run_All_Four_Figures.R again without amplicon.cache_file.")
    analysis_mode <- "Reused existing analyses; statistical models were not rerun."
  } else {
    computed <- compute_study(input_root)
    cache <- c(list(format = 1L, input_md5 = input_md5), computed)
    analysis_mode <- "Recomputed all eight datasets with the recorded fixed seeds."
  }
  run_name <- paste0("Run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  if (dir.exists(file.path(output_parent, run_name))) {
    suffix <- 1L
    while (dir.exists(file.path(output_parent, paste0(run_name, "_", suffix)))) suffix <- suffix + 1L
    run_name <- paste0(run_name, "_", suffix)
  }
  run_dir <- file.path(output_parent, run_name)
  for (folder in c("Figures", "Supporting_Data", "Tables", "Cache")) dir.create(file.path(run_dir, folder), recursive = TRUE)
  for (nm in names(cache$tables)) write.csv(cache$tables[[nm]], file.path(run_dir, "Supporting_Data", paste0(nm, ".csv")), row.names = FALSE)
  write.csv(data.frame(File = relative_inputs, Bytes = file.info(input_files)$size, MD5 = unname(input_md5)),
    file.path(run_dir, "Supporting_Data", "Input_Checksums.csv"), row.names = FALSE)
  saveRDS(cache, file.path(run_dir, "Cache", "study_analyses.rds"))
  captions_main <- character();captions_supp <- character();profiles <- list()
  for (i in which(slugs %in% datasets)) {
    slug <- slugs[i];object <- cache$objects[[slug]]
    results <- object$results;families <- object$families
    publication_number <- i;publication_slug <- slug
    publication_title <- paste(sub("_.*", "", slug), if (grepl("Bacteria", slug)) "bacterial" else "fungal", "communities")
    publication_family_title <- if (grepl("Bacteria", slug)) "Bacterial families" else "Fungal families"
    publication_out_dir <- file.path(run_dir, "Figures")
    publication_data_dir <- file.path(run_dir, "Supporting_Data", slug)
    dir.create(publication_data_dir, recursive = TRUE)
    for (comp in names(families)) write.csv(families[[comp]], file.path(publication_data_dir,
      paste0(comp, "_Per_Sample_Family_Percent.csv")), row.names = FALSE)
    rank <- dplyr::bind_rows(families) %>% dplyr::group_by(Source, Family) %>%
      dplyr::summarise(Mean = mean(Abundance), .groups = "drop") %>%
      tidyr::complete(Source = c("Rhizosphere", "Root"), Family, fill = list(Mean = 0)) %>%
      dplyr::group_by(Family) %>% dplyr::summarise(Mean = mean(Mean), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(Mean), Family)
    selected <- head(rank$Family, stacked_top_n)
    pie_selected <- lapply(families, pub_select_pie_families)
    current_colours <- pub_palette(unique(c(selected, unlist(pie_selected, use.names = FALSE))))
    if (length(publication_family_colours)) {
      if (is.null(names(publication_family_colours)) || anyDuplicated(names(publication_family_colours))) stop("Family colour overrides need unique names.")
      grDevices::col2rgb(publication_family_colours)
      present <- intersect(names(publication_family_colours), names(current_colours))
      current_colours[present] <- publication_family_colours[present]
    }
    write.csv(data.frame(Family = names(current_colours), Colour = unname(current_colours)),
      file.path(publication_data_dir, "Family_Colour_Key.csv"), row.names = FALSE)
    if ("main" %in% figures) {
      bars <- setNames(lapply(seq_along(families), function(k)
        pub_bar_panel(families[[k]], selected, current_colours, names(families)[k], c("C", "F")[k])), names(families))
      main <- journal_main(results, bars)
      pub_save_three(main, file.path(publication_out_dir, paste0("Figure", i, "_", slug, "_Main_A-F")),
        publication_main_width, publication_main_height)
      captions_main <- c(captions_main, paste0("Figure ", i, ". ", publication_title,
        ". (A,D) Shannon diversity in rhizosphere (A) and roots (D). Boxes show medians and interquartile ranges; whiskers extend to observations within 1.5 times the interquartile range. Points show individual retained samples. ",
        if (i == 1) "Rhizosphere Control-Drought has n = 3; other groups have n = 4. " else "Each treatment-by-moisture group has n = 4. ",
        "(B,E) PCoA of Bray-Curtis dissimilarities of square-root relative ASV abundances. Axis labels show variation percentages. Orange and blue denote Control and 6-PP; circles and triangles denote irrigated and rainfall-excluded samples. ",
        "(C,F) Mean family proportions by treatment and moisture. The top ", stacked_top_n,
        " groups are selected using an equal-weighted average of compartment means; remaining groups form Others. Keys are below their respective rows. The 6-PP label denotes the seed formulation; controls received water. Drought denotes rainfall exclusion. Family comparisons are descriptive."))
    }
    if ("supplementary" %in% figures) {
      pies <- pub_build_supplementary_pies(families, current_colours, pie_selected, write_outputs = TRUE)
      profiles[[slug]] <- pies$data;captions_supp <- c(captions_supp, pies$caption)
    }
    a <- cache$tables$alpha_tests;p <- cache$tables$permanova_tests
    rows <- do.call(rbind, lapply(c("Rhizosphere", "Root"), function(comp) do.call(rbind,
      lapply(c("Treatment", "Moisture", "Treatment:Moisture"), function(term) {
        x <- a[a$Dataset == slug & a$Source == comp & a$Term == term, ]
        y <- p[p$Dataset == slug & p$Source == comp & p$Term == term, ]
        data.frame(Compartment = comp, Effect = term, Shannon_F = x$`F value`,
          Shannon_P = x$`Pr(>F)`, PERMANOVA_R2 = y$R2, PERMANOVA_P = y$`Pr(>F)`,
          Shannon_BH_q = x$BH_q, PERMANOVA_BH_q = y$BH_q)
      }))))
    write.csv(rows, file.path(run_dir, "Tables", paste0("Table", i, "_", slug, ".csv")), row.names = FALSE)
    message("Finished: ", slug)
  }
  writeLines(paste(c(captions_main, captions_supp), collapse = "\n\n"),
    file.path(run_dir, "Figures", "Figure_captions.txt"), useBytes = TRUE)
  if (length(profiles)) write.csv(dplyr::bind_rows(profiles), file.path(run_dir, "Supporting_Data", "six_panel_pie_data.csv"), row.names = FALSE)
  expected_figures <- length(datasets) * length(figures)
  for (extension in c("pdf", "png", "tiff")) {
    written <- list.files(file.path(run_dir, "Figures"), pattern = paste0("\\.", extension, "$"), full.names = TRUE)
    if (length(written) != expected_figures || any(file.info(written)$size <= 0)) stop("Incomplete ", extension, " export.")
  }
  if (!identical(unname(tools::md5sum(input_files)), unname(input_md5))) stop("Input files changed during processing; review the outputs.")
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "Supporting_Data", "sessionInfo.txt"))
  file.copy(file.path(script_root, "Figure_Settings.R"), file.path(run_dir, "Supporting_Data", "Figure_Settings_Used.R"))
  if (!is.null(extra_settings)) file.copy(extra_settings, file.path(run_dir, "Supporting_Data", "Additional_Settings_Used.R"))
  writeLines(c(analysis_mode, paste("Datasets:", paste(datasets, collapse = ", ")),
    paste("Figure types:", paste(figures, collapse = ", ")),
    paste("Raster dpi:", publication_dpi), paste("Raster backend:", publication_raster_backend),
    "Default main figures: 190 mm x 9 in; default supplementary figures: 190 mm x 9.5 in.",
    "The original unrestricted analysis is reproduced; this does not resolve field block/shelter allocation or pooled-marker provenance.",
    "Input files were verified unchanged."), file.path(run_dir, "RUN_COMPLETE.txt"))
  writeLines(run_name, file.path(output_parent, "LATEST_ANALYSIS_RUN.txt"))
  message("\nCompleted: ", expected_figures, " figures in PDF, PNG and ", publication_dpi, "-dpi TIFF.")
  message("Figures: ", normalizePath(file.path(run_dir, "Figures"), winslash = "/"))
  invisible(list(run_dir = run_dir, figures = file.path(run_dir, "Figures"), analysis_mode = analysis_mode))
}
