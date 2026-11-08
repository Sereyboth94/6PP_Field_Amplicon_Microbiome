###############################################################################
# Unified Phyloseq Analysis – Wheat Field Fungi
# Main A-F figure, supplemental A-B pies and one main statistical table
# Vector PDF with embedded text; PNG and LZW TIFF at 600 dpi
# Original Combined/Root/Rhizosphere workflow retained as an optional section
#
# Includes:
#   - Alpha diversity
#   - Beta diversity + PERMANOVA
#   - DESeq2 differential abundance
#   - Portrait family stacked plots
#   - Overall family pie charts
#   - Enrichment vs depletion tables
#   - Top-5 family summaries
###############################################################################

###############################################################################
run_legacy_outputs <- isTRUE(getOption("amplicon.run_legacy", FALSE))  # Explicit opt-in.

# 0) Load packages
###############################################################################
required_pkgs <- c("phyloseq", "ggplot2", "dplyr", "tidyr", "tibble", "vegan", "patchwork")
if (run_legacy_outputs) {
  required_pkgs <- unique(c(required_pkgs, "ape", "DESeq2", "RColorBrewer", "forcats"))
}

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste(
      "Missing packages:", paste(missing_pkgs, collapse = ", "),
      "\nPlease install them manually before running this script."
    )
  )
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

###############################################################################
# 1) User settings
###############################################################################
# Locate this saved script when sourced in RStudio or invoked with Rscript.
.amp_sources <- vapply(sys.frames(), function(fr) {
  if (!is.null(fr$ofile)) as.character(fr$ofile)[1L] else ""
}, character(1))
.amp_sources <- .amp_sources[nzchar(.amp_sources)]
.amp_cli <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.amp_file <- if (length(.amp_sources)) tail(.amp_sources, 1L) else if (length(.amp_cli)) .amp_cli[1L] else ""
.amp_root <- getOption("amplicon.repo_root", if (nzchar(.amp_file)) {
  dirname(dirname(dirname(normalizePath(.amp_file, winslash = "/", mustWork = TRUE))))
} else getwd())
source(file.path(.amp_root, "R", "project_utils.R"), local = TRUE, encoding = "UTF-8")
.amp_paths <- amplicon_paths(.amp_root, "Wheat", "Fungi")
project_dir <- .amp_paths$project_dir
base_dir <- .amp_paths$base_dir

# Change these if your real files include suffixes like (3)
asv_file  <- file.path(base_dir, "ASV_filtered.txt")
tax_file  <- file.path(base_dir, "Taxonomy_filtered.txt")
meta_file <- file.path(base_dir, "Metadat_Fungal_Wheat.txt")

invisible(amplicon_validate_dataset(.amp_paths))
main_out_dir <- .amp_paths$main_out_dir
if (!dir.exists(main_out_dir)) dir.create(main_out_dir, recursive = TRUE)

# filtering thresholds
min_total <- 10
min_prev  <- 0.05
min_depth <- 500

# plotting / summaries
stacked_top_n <- 15
pie_top_n     <- 15
top_n_summary <- 5

# Publication settings for the requested main figure and supplemental pies.
publication_number <- 4
publication_slug <- "Wheat_Fungi"
publication_title <- "Wheat fungal communities"
publication_family_title <- "Fungal families"
publication_out_dir <- file.path(main_out_dir, "Manuscript_Figures_Tables")
publication_font <- getOption("amplicon.font", "Arial")
publication_base_size <- 20
publication_axis_size <- 18
publication_legend_size <- 17.5
publication_family_size <- 20
publication_main_family_size <- 18
publication_bar_axis_size <- 17
publication_bar_title_size <- 20
publication_main_legend_columns <- 3
publication_tag_size <- 28
publication_dpi <- 600
publication_main_width <- 21
publication_main_height <- 14.5
publication_pie_width <- 18
publication_pie_height <- 12
publication_pie_label_min <- 5  # Show percentages on slices of at least 5%.
publication_pie_label_size <- 6.4  # ggplot geom_text size in millimetres, about 18 pt.
publication_permutations <- 999
publication_seed <- 20261105
publication_treatment_colours <- c("Control" = "#D55E00", "6PP" = "#0072B2")


###############################################################################
# 2) Helper functions
###############################################################################
save_text_output <- function(object, filename) {
  capture.output(object, file = filename)
}

percent_transform <- function(x) {
  100 * x / sum(x)
}

sqrt_rel_transform <- function(x) {
  sqrt(x / sum(x))
}

clean_taxonomy_labels <- function(x, rank = NULL) {
  
  x <- as.character(x)
  
  # Handle missing values first
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  
  # Remove common taxonomy prefixes
  x <- gsub("^[a-zA-Z]__", "", x)
  x <- gsub("^[a-zA-Z]_", "", x)
  x <- gsub("^[a-zA-Z]-", "", x)
  
  # Replace underscores with spaces
  x <- gsub("_", " ", x)
  
  # Remove incertae sedis annotations for cleaner plotting
  x <- gsub("\\s*fam\\s*\\(incertae sedis\\)", "", x, ignore.case = TRUE)
  x <- gsub("\\s*\\(incertae sedis\\)", "", x, ignore.case = TRUE)
  x <- gsub("\\s*incertae sedis", "", x, ignore.case = TRUE)
  
  # Remove extra whitespace
  x <- trimws(x)
  
  # Final missing-value handling
  if (!is.null(rank)) {
    x[is.na(x) | x == ""] <- paste0("Unclassified ", rank)
  } else {
    x[is.na(x) | x == ""] <- "Unclassified"
  }
  
  return(x)
}


safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

make_nature_palette <- function(n) {
  base_cols <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#B2DF8A",
    "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99", "#8DD3C7",
    "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5",
    "#FFED6F", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F"
  )
  
  if (n <= length(base_cols)) {
    pal <- base_cols[seq_len(n)]
  } else {
    pal <- grDevices::colorRampPalette(base_cols)(n)
  }
  pal
}

treatment_labels <- c(
  "Control" = "Control",
  "6PP" = "6-PP"
)

###############################################################################
# 3) Load data
###############################################################################
# Check exact input filenames before trying to load data.
input_files <- c(asv_file, tax_file, meta_file)
if (any(!file.exists(input_files))) {
  stop("Missing input files. Check base_dir and the filenames in section 1:\n",
       paste(input_files[!file.exists(input_files)], collapse = "\n"))
}

# ASV table
asv <- read.delim(asv_file, sep = "\t", header = TRUE, check.names = FALSE)
colnames(asv)[1] <- "ASV_ID"
rownames(asv) <- asv$ASV_ID
asv <- asv[, -1, drop = FALSE]

# Taxonomy table
tax <- read.delim(tax_file, sep = "\t", header = TRUE, check.names = FALSE)
colnames(tax)[1] <- "ASV_ID"
rownames(tax) <- tax$ASV_ID
tax <- tax[, -1, drop = FALSE]

# Metadata
meta <- read.delim(
  meta_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  strip.white = TRUE
)

colnames(meta)[1] <- "SampleID"
meta$SampleID  <- trimws(meta$SampleID)
meta$Treatment <- trimws(meta$Treatment)

meta$Treatment <- dplyr::recode(
  meta$Treatment,
  "FVOC" = "6PP",
  "6-PP" = "6PP",
  "6PP" = "6PP",
  "Control" = "Control"
)

meta$Treatment <- factor(meta$Treatment, levels = c("Control", "6PP"))
meta$Moisture  <- factor(trimws(meta$Moisture), levels = c("Irrigated", "Drought"))
meta$Source    <- factor(trimws(meta$Source), levels = c("Root", "Rhizosphere"))
meta$Plant     <- factor(trimws(meta$Plant))
rownames(meta) <- meta$SampleID

###############################################################################
# 4) Match samples and ASVs
###############################################################################
common_samples <- intersect(colnames(asv), rownames(meta))
if (length(common_samples) == 0) {
  stop("No matching sample IDs between ASV table and metadata.")
}

asv  <- asv[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

common_asvs <- intersect(rownames(asv), rownames(tax))
asv <- asv[common_asvs, , drop = FALSE]
tax <- tax[common_asvs, , drop = FALSE]

###############################################################################
# 5) Build full phyloseq object
###############################################################################
OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
SAM <- sample_data(meta)

ps_full <- phyloseq(OTU, TAX, SAM)

clean_family_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- "Unclassified"
  x <- gsub("^[a-z]__+", "", x)   # removes f__, g__, p__, etc.
  x <- gsub("^[a-z]_", "", x)     # removes f_, g_, p_, etc.
  x <- gsub("_fam_Incertae_sedis", " incertae_sedis", x)
  x
}

###############################################################################
# 6A) Pie chart function
###############################################################################
plot_family_pie <- function(ps, analysis_name, out_dir, top_families) {
  ps_family <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
  ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))
  
  df <- psmelt(ps_family_rel)
  df$Family <- clean_taxonomy_labels(df$Family, rank = "Family")
  
  df_sum <- df %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(perc = mean_ab * 100) %>%
    dplyr::arrange(dplyr::desc(perc))
  
  top_fams <- top_families
  
  df_pie <- df_sum %>%
    dplyr::mutate(Family = ifelse(Family %in% top_fams, Family, "Others")) %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(perc = sum(perc), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(perc))
  
  family_levels <- df_pie$Family
  family_levels <- c(setdiff(family_levels, "Others"), "Others")
  df_pie$Family <- factor(df_pie$Family, levels = family_levels)
  
  pal <- make_nature_palette(length(family_levels))
  names(pal) <- family_levels
  if ("Others" %in% names(pal)) pal["Others"] <- "#BDBDBD"
  
  df_pie$label <- ifelse(df_pie$perc >= 3, paste0(round(df_pie$perc, 1), "%"), "")
  
  p <- ggplot(df_pie, aes(x = "", y = perc, fill = Family)) +
    geom_col(width = 1, color = "white", linewidth = 0.2) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_fill_manual(values = pal, breaks = family_levels) +
    theme_void(base_size = 13) +
    labs(
      title = paste("Overall family composition -", analysis_name),
      fill = "Family"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    file.path(out_dir, paste0("Family_Pie_", analysis_name, ".png")),
    p, width = 10, height = 7, dpi = 300
  )
  
  write.csv(
    df_pie,
    file.path(out_dir, paste0("Family_Pie_Table_", analysis_name, ".csv")),
    row.names = FALSE
  )
}

###############################################################################
# 6B) Overall family composition bar plot
###############################################################################
plot_family_overall_bar <- function(ps, analysis_name, out_dir, top_families) {
  ps_family <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
  ps_family_rel <- transform_sample_counts(ps_family, function(x) 100 * x / sum(x))
  
  df <- psmelt(ps_family_rel)
  df$Family <- clean_taxonomy_labels(df$Family, rank = "Family")
  
  df_sum <- df %>%
    dplyr::group_by(Moisture, Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  top_fams <- top_families
  
  df_bar <- df_sum %>%
    dplyr::mutate(Family = ifelse(Family %in% top_fams, Family, "Others")) %>%
    dplyr::group_by(Moisture, Family) %>%
    dplyr::summarise(mean_ab = sum(mean_ab), .groups = "drop")
  
  family_levels <- df_bar %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(overall_mean = mean(mean_ab), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(overall_mean)) %>%
    dplyr::pull(Family)
  
  family_levels <- c(setdiff(family_levels, "Others"), "Others")
  df_bar$Family <- factor(df_bar$Family, levels = family_levels)
  df_bar$Moisture <- factor(df_bar$Moisture, levels = c("Irrigated", "Drought"))
  
  df_bar <- df_bar %>%
    dplyr::group_by(Moisture) %>%
    dplyr::mutate(
      total = sum(mean_ab, na.rm = TRUE),
      mean_ab = ifelse(total > 0, mean_ab / total * 100, 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-total)
  
  df_bar <- df_bar %>%
    dplyr::mutate(
      label = ifelse(mean_ab >= 4, paste0(round(mean_ab, 1), "%"), "")
    )
  
  pal <- make_nature_palette(length(family_levels))
  names(pal) <- family_levels
  if ("Others" %in% names(pal)) pal["Others"] <- "#BDBDBD"
  
  p <- ggplot(df_bar, aes(x = Moisture, y = mean_ab, fill = Family)) +
    geom_col(width = 0.7, color = NA) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_fill_manual(values = pal, breaks = family_levels) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      title = paste("Overall family composition -", analysis_name),
      x = NULL,
      y = "Mean relative abundance (%)",
      fill = "Family"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0)
    )
  
  ggsave(
    file.path(out_dir, paste0("Family_OverallBar_", analysis_name, ".png")),
    p, width = 7, height = 7, dpi = 300
  )
  
  write.csv(
    df_bar,
    file.path(out_dir, paste0("Family_OverallBar_Table_", analysis_name, ".csv")),
    row.names = FALSE
  )
  return(p)
}

###############################################################################
# 7) Portrait stacked family plot
###############################################################################
plot_family_stacked_portrait <- function(ps, analysis_name, out_dir, top_families) {
  ps_family <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
  ps_family_rel <- transform_sample_counts(ps_family, function(x) 100 * x / sum(x))
  
  df <- psmelt(ps_family_rel)
  df$Family <- clean_taxonomy_labels(df$Family, rank = "Family")
  
  meta_df <- data.frame(sample_data(ps))
  meta_df$Treatment <- as.character(meta_df$Treatment)
  meta_df$Source    <- as.character(meta_df$Source)
  meta_df$Moisture  <- as.character(meta_df$Moisture)
  
  df_sum <- df %>%
    dplyr::mutate(
      Treatment = as.character(Treatment),
      Source    = as.character(Source),
      Moisture  = as.character(Moisture)
    ) %>%
    dplyr::group_by(Treatment, Source, Moisture, Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  top_fams <- top_families
  
  df_sum2 <- df_sum %>%
    dplyr::mutate(Family = ifelse(Family %in% top_fams, Family, "Others")) %>%
    dplyr::group_by(Treatment, Source, Moisture, Family) %>%
    dplyr::summarise(mean_ab = sum(mean_ab, na.rm = TRUE), .groups = "drop")
  
  group_template <- meta_df %>%
    dplyr::distinct(Treatment, Source, Moisture)
  
  family_template <- data.frame(
    Family = sort(unique(df_sum2$Family)),
    stringsAsFactors = FALSE
  )
  
  df_sum2 <- tidyr::crossing(group_template, family_template) %>%
    dplyr::left_join(
      df_sum2,
      by = c("Treatment", "Source", "Moisture", "Family")
    ) %>%
    dplyr::mutate(
      mean_ab = ifelse(is.na(mean_ab) | !is.finite(mean_ab), 0, mean_ab)
    ) %>%
    dplyr::group_by(Treatment, Source, Moisture) %>%
    dplyr::mutate(
      total_ab = sum(mean_ab, na.rm = TRUE),
      mean_ab  = ifelse(total_ab > 0, mean_ab / total_ab * 100, 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-total_ab)
  
  df_sum2 <- df_sum2 %>%
    dplyr::group_by(Treatment, Source, Moisture) %>%
    dplyr::mutate(
      total_bar = sum(mean_ab, na.rm = TRUE),
      deficit   = 100 - total_bar
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      mean_ab = ifelse(Family == "Others", mean_ab + deficit, mean_ab),
      mean_ab = ifelse(is.na(mean_ab) | !is.finite(mean_ab), 0, mean_ab),
      mean_ab = pmax(mean_ab, 0),
      mean_ab = round(mean_ab, 8)
    ) %>%
    dplyr::select(-total_bar, -deficit)
  
  fam_order2 <- df_sum2 %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(overall_mean = mean(mean_ab, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(overall_mean))
  
  family_levels <- fam_order2$Family
  family_levels <- c(setdiff(family_levels, "Others"), "Others")
  df_sum2$Family <- factor(df_sum2$Family, levels = family_levels)
  
  df_sum2$Treatment <- factor(df_sum2$Treatment, levels = c("Control", "6PP"))
  df_sum2$Source    <- factor(df_sum2$Source, levels = c("Root", "Rhizosphere"))
  df_sum2$Moisture  <- factor(df_sum2$Moisture, levels = c("Irrigated", "Drought"))
  
  pal <- make_nature_palette(length(family_levels))
  names(pal) <- family_levels
  if ("Others" %in% names(pal)) pal["Others"] <- "#BDBDBD"
  
  present_sources <- unique(as.character(df_sum2$Source))
  single_source <- length(present_sources) == 1
  
  if (single_source) {
    p <- ggplot(df_sum2, aes(x = Treatment, y = mean_ab, fill = Family)) +
      geom_col(width = 0.8, color = NA) +
      facet_wrap(~ Moisture, nrow = 1) +
      scale_x_discrete(labels = treatment_labels) +
      scale_fill_manual(values = pal, breaks = family_levels) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(
        title = paste("Family-level composition -", analysis_name),
        x = NULL,
        y = "Mean relative abundance (%)",
        fill = "Family"
      ) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0)
      )
  } else {
    df_sum2 <- df_sum2 %>%
      dplyr::mutate(
        Group = factor(
          paste(Source, ifelse(Moisture == "Irrigated", "I", "D"), sep = "\n"),
          levels = c("Root\nI", "Root\nD", "Rhizosphere\nI", "Rhizosphere\nD")
        )
      )
    
    p <- ggplot(df_sum2, aes(x = Group, y = mean_ab, fill = Family)) +
      geom_col(width = 0.8, color = NA) +
      facet_wrap(
        ~ Treatment,
        nrow = 1,
        labeller = labeller(Treatment = treatment_labels)
      ) +
      scale_fill_manual(values = pal, breaks = family_levels) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(
        title = paste("Family-level composition -", analysis_name),
        x = NULL,
        y = "Mean relative abundance (%)",
        fill = "Family"
      ) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0)
      )
  }
  
  ggsave(
    file.path(out_dir, paste0("Family_StackedBar_Portrait_", analysis_name, ".png")),
    p, width = 10, height = 7, dpi = 300
  )
  
  write.csv(
    df_sum2,
    file.path(out_dir, paste0("Family_Composition_Table_Portrait_", analysis_name, ".csv")),
    row.names = FALSE
  )
  
  if (single_source) {
    check_totals <- df_sum2 %>%
      dplyr::group_by(Treatment, Moisture) %>%
      dplyr::summarise(total = sum(mean_ab), .groups = "drop")
  } else {
    check_totals <- df_sum2 %>%
      dplyr::group_by(Treatment, Source, Moisture) %>%
      dplyr::summarise(total = sum(mean_ab), .groups = "drop")
  }
  
  write.csv(
    check_totals,
    file.path(out_dir, paste0("Family_Composition_Check_", analysis_name, ".csv")),
    row.names = FALSE
  )
  
  print(check_totals)
  return(p)
}

###############################################################################
# 8) Enrichment / depletion summary
###############################################################################
make_enrichment_table <- function(ps, analysis_name, out_dir, top_n = 5, include_source = TRUE) {
  ps_family <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
  ps_family_rel <- transform_sample_counts(ps_family, percent_transform)
  
  df <- psmelt(ps_family_rel)
  df$Family <- clean_taxonomy_labels(df$Family, rank = "Family")
  
  top_fams <- df %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_ab)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(Family)
  
  df_top <- df %>% dplyr::filter(Family %in% top_fams)
  
  if (include_source) {
    summary_tbl <- df_top %>%
      dplyr::group_by(Source, Moisture, Treatment, Family) %>%
      dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Treatment, values_from = mean_ab) %>%
      dplyr::mutate(
        PercentChange_6PP_vs_Control = ifelse(
          .data[["Control"]] > 0,
          ((.data[["6PP"]] - .data[["Control"]]) / .data[["Control"]]) * 100,
          NA_real_
        ),
        Direction = dplyr::case_when(
          is.na(PercentChange_6PP_vs_Control) ~ "Undefined",
          PercentChange_6PP_vs_Control > 0 ~ "Enriched",
          PercentChange_6PP_vs_Control < 0 ~ "Depleted",
          TRUE ~ "No change"
        )
      ) %>%
      dplyr::arrange(Source, Moisture, dplyr::desc(abs(PercentChange_6PP_vs_Control)))
  } else {
    summary_tbl <- df_top %>%
      dplyr::group_by(Moisture, Treatment, Family) %>%
      dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Treatment, values_from = mean_ab) %>%
      dplyr::mutate(
        PercentChange_6PP_vs_Control = ifelse(
          .data[["Control"]] > 0,
          ((.data[["6PP"]] - .data[["Control"]]) / .data[["Control"]]) * 100,
          NA_real_
        ),
        Direction = dplyr::case_when(
          is.na(PercentChange_6PP_vs_Control) ~ "Undefined",
          PercentChange_6PP_vs_Control > 0 ~ "Enriched",
          PercentChange_6PP_vs_Control < 0 ~ "Depleted",
          TRUE ~ "No change"
        )
      ) %>%
      dplyr::arrange(Moisture, dplyr::desc(abs(PercentChange_6PP_vs_Control)))
  }
  
  write.csv(
    summary_tbl,
    file.path(out_dir, paste0("Family_Enrichment_Depletion_", analysis_name, ".csv")),
    row.names = FALSE
  )
}

###############################################################################
# 9) Main analysis function
###############################################################################
run_phyloseq_analysis <- function(ps_in, analysis_name, out_dir, include_source_in_model = TRUE) {
  safe_dir(out_dir)
  
  ps <- prune_samples(sample_sums(ps_in) >= min_depth, ps_in)
  
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  
  keep_taxa <- rowSums(otu) >= min_total &
    rowMeans(otu > 0) >= min_prev
  
  ps <- prune_taxa(keep_taxa, ps)
  
  saveRDS(ps, file.path(out_dir, paste0("phyloseq_", analysis_name, "_filtered.rds")))
  
  meta_df <- data.frame(sample_data(ps))
  if (!"SampleID" %in% colnames(meta_df)) {
    meta_df <- meta_df %>% tibble::rownames_to_column("SampleID")
  }
  
  ###########################################################################
  # Alpha diversity
  ###########################################################################
  alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon", "Chao1")) %>%
    tibble::rownames_to_column("SampleID") %>%
    left_join(meta_df, by = "SampleID")
  
  write.csv(
    alpha_df,
    file.path(out_dir, paste0("AlphaDiversity_", analysis_name, ".csv")),
    row.names = FALSE
  )
  
  treatment_labels <- c(
    "Control" = "Control",
    "6PP" = "6-PP"
  )
  
  treatment_cols <- c(
    "Control" = "#F8766D",
    "6PP"     = "#00BA38"
  )
  
  alpha_theme <- theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )
  
  if (include_source_in_model) {
    alpha_formula <- as.formula("Shannon ~ Treatment * Moisture * Source")
    
    p_alpha <- ggplot(alpha_df, aes(x = Treatment, y = Shannon, fill = Treatment)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
      facet_grid(Source ~ Moisture) +
      scale_x_discrete(labels = treatment_labels) +
      scale_fill_manual(
        values = treatment_cols,
        labels = treatment_labels,
        name = "Treatment"
      ) +
      theme_bw(base_size = 14) +
      alpha_theme +
      labs(
        title = paste("Shannon diversity -", analysis_name),
        x = NULL,
        y = "Shannon index"
      )
    
  } else {
    alpha_formula <- as.formula("Shannon ~ Treatment * Moisture")
    
    p_alpha <- ggplot(alpha_df, aes(x = Treatment, y = Shannon, fill = Treatment)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
      facet_wrap(~ Moisture) +
      scale_x_discrete(labels = treatment_labels) +
      scale_fill_manual(
        values = treatment_cols,
        labels = treatment_labels,
        name = "Treatment"
      ) +
      theme_bw(base_size = 14) +
      alpha_theme +
      labs(
        title = paste("Shannon diversity -", analysis_name),
        x = NULL,
        y = "Shannon index"
      )
  }
  
  alpha_aov <- aov(alpha_formula, data = alpha_df)
  
  save_text_output(
    summary(alpha_aov),
    file.path(out_dir, paste0("Alpha_Shannon_ANOVA_", analysis_name, ".txt"))
  )
  
  ggsave(
    file.path(out_dir, paste0("Alpha_Shannon_Boxplot_", analysis_name, ".png")),
    p_alpha,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ###########################################################################
  # Beta diversity
  ###########################################################################
  ps_rel_sqrt <- transform_sample_counts(ps, sqrt_rel_transform)
  dist_bc <- phyloseq::distance(ps_rel_sqrt, method = "bray")
  ord <- ordinate(ps_rel_sqrt, method = "PCoA", distance = dist_bc)
  
  beta_theme <- theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )
  
  if (include_source_in_model) {
    p_pcoa <- plot_ordination(
      ps_rel_sqrt,
      ord,
      color = "Treatment",
      shape = "Moisture"
    ) +
      geom_point(size = 3, alpha = 0.9) +
      facet_wrap(~ Source) +
      scale_color_manual(
        values = treatment_cols,
        labels = treatment_labels,
        name = "Treatment"
      ) +
      theme_bw(base_size = 14) +
      beta_theme +
      labs(title = paste("PCoA (Bray–Curtis) -", analysis_name))
    
  } else {
    p_pcoa <- plot_ordination(
      ps_rel_sqrt,
      ord,
      color = "Treatment",
      shape = "Moisture"
    ) +
      geom_point(size = 3, alpha = 0.9) +
      scale_color_manual(
        values = treatment_cols,
        labels = treatment_labels,
        name = "Treatment"
      ) +
      theme_bw(base_size = 14) +
      beta_theme +
      labs(title = paste("PCoA (Bray–Curtis) -", analysis_name))
  }
  
  ggsave(
    file.path(out_dir, paste0("PCoA_BrayCurtis_", analysis_name, ".png")),
    p_pcoa,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  beta_meta <- data.frame(sample_data(ps))
  
  if (include_source_in_model) {
    adon <- adonis2(
      dist_bc ~ Treatment * Moisture * Source,
      data = beta_meta,
      permutations = 999,
      by = "terms"
    )
  } else {
    adon <- adonis2(
      dist_bc ~ Treatment * Moisture,
      data = beta_meta,
      permutations = 999,
      by = "terms"
    )
  }
  
  save_text_output(
    adon,
    file.path(out_dir, paste0("PERMANOVA_BrayCurtis_", analysis_name, ".txt"))
  )
  
  ###########################################################################
  # Combined Alpha + Beta figure
  ###########################################################################
  combined_alpha_beta <- (p_alpha + p_pcoa) +
    patchwork::plot_layout(widths = c(1, 1)) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(
        size = 20,
        face = "bold"
      )
    )
  
  ggsave(
    file.path(out_dir, paste0("Alpha_Beta_Combined_", analysis_name, ".png")),
    combined_alpha_beta,
    width = 14,
    height = 6,
    dpi = 300
  )
  
  ###########################################################################
  # DESeq2
  ###########################################################################
  if (include_source_in_model) {
    dds <- phyloseq_to_deseq2(ps, ~ Source + Moisture + Treatment)
  } else {
    dds <- phyloseq_to_deseq2(ps, ~ Moisture + Treatment)
  }
  
  dds <- DESeq(dds, fitType = "local", sfType = "poscounts")
  res_6PP <- results(dds, contrast = c("Treatment", "6PP", "Control"))
  
  res_6PP_df <- as.data.frame(res_6PP) %>%
    tibble::rownames_to_column("ASV_ID") %>%
    arrange(padj)
  
  write.csv(
    res_6PP_df,
    file.path(out_dir, paste0("DESeq2_6PP_vs_Control_", analysis_name, ".csv")),
    row.names = FALSE
  )
  
  tax_df <- as.data.frame(tax_table(ps)) %>%
    tibble::rownames_to_column("ASV_ID")
  
  res_6PP_tax <- left_join(res_6PP_df, tax_df, by = "ASV_ID")
  write.csv(
    res_6PP_tax,
    file.path(out_dir, paste0("DESeq2_6PP_vs_Control_withTaxonomy_", analysis_name, ".csv")),
    row.names = FALSE
  )
  
  ###########################################################################
  # Shared top families for consistent family plots
  ###########################################################################
  ps_family_shared <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
  ps_family_rel_shared <- transform_sample_counts(ps_family_shared, function(x) 100 * x / sum(x))
  
  df_shared <- psmelt(ps_family_rel_shared)
  df_shared$Family <- clean_taxonomy_labels(df_shared$Family, rank = "Family")
  
  top_families <- df_shared %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_ab)) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::pull(Family)
  
  ###########################################################################
  # Family portrait plot
  ###########################################################################
  plot_family_stacked_portrait(
    ps = ps,
    analysis_name = analysis_name,
    out_dir = out_dir,
    top_families = top_families
  )
  
  ###########################################################################
  # Family pie chart
  ###########################################################################
  plot_family_pie(
    ps = ps,
    analysis_name = analysis_name,
    out_dir = out_dir,
    top_families = top_families
  )
  
  ###########################################################################
  # Overall stacked bar
  ###########################################################################
  plot_family_overall_bar(
    ps = ps,
    analysis_name = analysis_name,
    out_dir = out_dir,
    top_families = top_families
  )
  
  ###########################################################################
  # Enrichment vs depletion summary
  ###########################################################################
  make_enrichment_table(
    ps = ps,
    analysis_name = analysis_name,
    out_dir = out_dir,
    top_n = top_n_summary,
    include_source = include_source_in_model
  )
  
  ###########################################################################
  # Top-5 families
  ###########################################################################
  ps_family <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
  ps_family_rel <- transform_sample_counts(ps_family, percent_transform)
  df_family <- psmelt(ps_family_rel)
  df_family$Family <- clean_taxonomy_labels(df_family$Family, rank = "Family")
  
  top5 <- df_family %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_ab)) %>%
    dplyr::slice_head(n = top_n_summary)
  
  write.csv(
    top5,
    file.path(out_dir, paste0("Top5_Families_Overall_", analysis_name, ".csv")),
    row.names = FALSE
  )
  
  top5_names <- top5$Family
  
  if (include_source_in_model) {
    top5_by_group <- df_family %>%
      dplyr::filter(Family %in% top5_names) %>%
      dplyr::group_by(Source, Moisture, Treatment, Family) %>%
      dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Treatment, values_from = mean_ab)
  } else {
    top5_by_group <- df_family %>%
      dplyr::filter(Family %in% top5_names) %>%
      dplyr::group_by(Moisture, Treatment, Family) %>%
      dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Treatment, values_from = mean_ab)
  }
  
  write.csv(
    top5_by_group,
    file.path(out_dir, paste0("Top5_Families_ByGroup_", analysis_name, ".csv")),
    row.names = FALSE
  )
  
  ###########################################################################
  # Session info
  ###########################################################################
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, paste0("sessionInfo_", analysis_name, ".txt"))
  )
  
  message("Finished: ", analysis_name)
}


###############################################################################
# 10) Publication figures and table
###############################################################################
# Publication outputs shared in content, but embedded in each standalone script.
# Statistical models and filtering thresholds match the original scripts.

pub_theme <- function() {
  ggplot2::theme_bw(base_size = publication_base_size, base_family = publication_font) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(size = publication_axis_size, colour = "black"),
      axis.title = ggplot2::element_text(size = publication_axis_size + 1),
      strip.text = ggplot2::element_text(size = publication_axis_size, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F0F0F0", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#EAEAEA", linewidth = 0.25),
      plot.title = ggplot2::element_text(size = publication_base_size + 1, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = publication_axis_size, margin = ggplot2::margin(b = 10)),
      plot.tag = ggplot2::element_text(size = publication_tag_size, face = "bold"),
      plot.margin = ggplot2::margin(14, 12, 14, 12),
      legend.title = ggplot2::element_text(size = publication_legend_size, face = "bold"),
      legend.text = ggplot2::element_text(size = publication_legend_size, lineheight = 1.05),
      legend.key.height = grid::unit(0.50, "cm"),
      legend.key.width = grid::unit(0.52, "cm"),
      legend.spacing.y = grid::unit(0.20, "cm"),
      legend.box = "vertical"
    )
}

pub_clean_family <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("^[A-Za-z]__?", "", x)
  x <- trimws(gsub("[[:space:]]+", " ", gsub("_", " ", x)))
  x[x == "" | toupper(x) == "NA"] <- "Unclassified family"
  x
}

pub_wrap <- function(x, width = 28) {
  vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

pub_save_three <- function(plot, stem, width, height) {
  if (!isTRUE(capabilities("cairo"))) {
    stop("Cairo graphics support is required for PDFs with embedded, selectable text. ",
         "Check capabilities('cairo') in your R installation.")
  }
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  draw_plot <- function() {
    if (inherits(plot, c("ggplot", "patchwork"))) print(plot) else {
      grid::grid.newpage()
      grid::grid.draw(plot)
    }
  }
  write_device <- function(ext) {
    filename <- paste0(stem, ".", ext)
    if (ext == "pdf") {
      grDevices::cairo_pdf(filename, width = width, height = height,
        family = publication_font, bg = "white", onefile = TRUE,
        fallback_resolution = publication_dpi)
    } else if (requireNamespace("ragg", quietly = TRUE)) {
      if (ext == "png") ragg::agg_png(filename, width = width, height = height,
        units = "in", res = publication_dpi, background = "white")
      else ragg::agg_tiff(filename, width = width, height = height,
        units = "in", res = publication_dpi, background = "white", compression = "lzw")
    } else {
      if (ext == "png") grDevices::png(filename, width = width, height = height,
        units = "in", res = publication_dpi, bg = "white", type = "cairo")
      else grDevices::tiff(filename, width = width, height = height,
        units = "in", res = publication_dpi, bg = "white", type = "cairo", compression = "lzw")
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_plot()
    invisible(filename)
  }
  files <- vapply(c("pdf", "png", "tiff"), write_device, character(1))
  if (any(!file.exists(files)) || any(file.info(files)$size == 0)) {
    stop("An export is missing or empty: ", stem)
  }
  # Optional embedding audit; pdftools is not needed to create embedded fonts.
  if (requireNamespace("pdftools", quietly = TRUE)) {
    font_audit <- pdftools::pdf_fonts(paste0(stem, ".pdf"))
    utils::write.csv(font_audit, paste0(stem, "_PDF_Fonts.csv"), row.names = FALSE)
    if (nrow(font_audit) > 0 && any(!font_audit$embedded)) {
      warning("A PDF font was not embedded; inspect the PDF_Fonts.csv audit: ", stem)
    }
  }
  invisible(files)
}

pub_filter_compartment <- function(ps_in, compartment) {
  meta0 <- data.frame(phyloseq::sample_data(ps_in), check.names = FALSE)
  ids <- rownames(meta0)[as.character(meta0$Source) == compartment & !is.na(meta0$Source)]
  if (length(ids) == 0) stop("No samples for ", compartment)
  ps <- phyloseq::prune_samples(ids, ps_in)
  ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) >= min_depth, ps)
  if (phyloseq::nsamples(ps) < 5) stop("Too few samples after depth filtering: ", compartment)
  counts <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)
  keep <- rowSums(counts) >= min_total & rowMeans(counts > 0) >= min_prev
  if (!any(keep)) stop("No ASVs remain after filtering: ", compartment)
  ps <- phyloseq::prune_taxa(keep, ps)
  if (any(phyloseq::sample_sums(ps) == 0)) {
    stop("Zero-count samples after ASV filtering in ", compartment, ": ",
         paste(phyloseq::sample_names(ps)[phyloseq::sample_sums(ps) == 0], collapse = ", "))
  }
  ps
}

pub_analyse_compartment <- function(ps, compartment, seed) {
  md <- data.frame(phyloseq::sample_data(ps), check.names = FALSE)
  md$SampleID <- rownames(md)
  md$Treatment <- factor(as.character(md$Treatment), levels = c("Control", "6PP"))
  md$Moisture <- factor(as.character(md$Moisture), levels = c("Irrigated", "Drought"))
  if (anyNA(md[, c("Treatment", "Moisture")])) stop("Unknown/missing treatment or moisture: ", compartment)
  if (any(table(md$Treatment, md$Moisture) == 0)) {
    stop("All four treatment-by-moisture groups must be represented: ", compartment)
  }
  alpha <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon"))
  alpha$SampleID <- rownames(alpha)
  alpha <- dplyr::left_join(alpha, md, by = "SampleID")
  alpha_fit <- stats::aov(Shannon ~ Treatment * Moisture, data = alpha)
  alpha_anova <- as.data.frame(summary(alpha_fit)[[1]])
  alpha_anova$Term <- trimws(rownames(alpha_anova))

  ps_sqrt <- phyloseq::transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
  dist_bc <- phyloseq::distance(ps_sqrt, method = "bray")
  ord <- phyloseq::ordinate(ps_sqrt, method = "PCoA", distance = dist_bc)
  if (ncol(ord$vectors) < 2) stop("Fewer than two PCoA axes: ", compartment)
  ord_data <- data.frame(SampleID = rownames(ord$vectors),
    PCoA1 = ord$vectors[, 1], PCoA2 = ord$vectors[, 2], check.names = FALSE)
  ord_data <- dplyr::left_join(ord_data, md, by = "SampleID")
  # Explicit ordering aligns the distance matrix and model metadata.
  md <- md[attr(dist_bc, "Labels"), , drop = FALSE]
  set.seed(seed)
  permanova <- vegan::adonis2(dist_bc ~ Treatment * Moisture,
    data = md, permutations = publication_permutations, by = "terms")
  permanova <- as.data.frame(permanova)
  permanova$Term <- rownames(permanova)
  list(ps = ps, metadata = md, alpha = alpha, alpha_anova = alpha_anova,
       permanova = permanova, ordination = ord_data,
       eigenvalues = ord$values, eigen_percent = 100 * ord$values$Relative_eig[1:2],
       compartment = compartment, seed = seed)
}

pub_family_data <- function(ps, compartment) {
  counts <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)
  taxonomy <- as(phyloseq::tax_table(ps), "matrix")
  if (!"Family" %in% colnames(taxonomy)) stop("Taxonomy must include a Family column.")
  family_counts <- rowsum(counts, group = pub_clean_family(taxonomy[rownames(counts), "Family"]), reorder = FALSE)
  family_pct <- 100 * sweep(family_counts, 2, colSums(family_counts), "/")
  df <- as.data.frame(as.table(family_pct), responseName = "Abundance", stringsAsFactors = FALSE)
  names(df)[1:2] <- c("Family", "SampleID")
  md <- data.frame(phyloseq::sample_data(ps), check.names = FALSE)
  md$SampleID <- rownames(md)
  df <- dplyr::left_join(df, md[, c("SampleID", "Treatment", "Moisture"), drop = FALSE], by = "SampleID")
  df$Source <- compartment
  check <- df %>% dplyr::group_by(SampleID) %>% dplyr::summarise(Total = sum(Abundance), .groups = "drop")
  if (any(abs(check$Total - 100) > 1e-7)) stop("Family percentages do not total 100 in ", compartment)
  df
}

pub_palette <- function(families) {
  # A named mapping guarantees that a family keeps its colour in every panel.
  base <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#B2DF8A", "#FB9A99",
    "#FDBF6F", "#CAB2D6", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD",
    "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
    "#DC0000", "#7E6148", "#B09C85", "#8DD3C7", "#B15928", "#6A3D9A",
    "#33A02C", "#E31A1C", "#1F78B4", "#FF7F00", "#B79F00", "#00BFC4")
  if (length(families) > length(base)) base <- grDevices::hcl.colors(length(families), "Dynamic")
  out <- stats::setNames(base[seq_along(families)], families)
  out["Others"] <- "#BDBDBD"
  out
}

pub_collapse_families <- function(df, selected) {
  df %>%
    dplyr::mutate(Family = ifelse(Family %in% selected, Family, "Others")) %>%
    dplyr::group_by(SampleID, Treatment, Moisture, Source, Family) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
}

pub_alpha_panel <- function(result, tag) {
  ggplot2::ggplot(result$alpha, ggplot2::aes(x = Treatment, y = Shannon, fill = Treatment)) +
    ggplot2::geom_boxplot(width = 0.68, outlier.shape = NA, colour = "#222222", linewidth = 0.6, show.legend = FALSE) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.13, height = 0, seed = result$seed),
      colour = "#222222", size = 2.4, show.legend = FALSE) +
    ggplot2::facet_wrap(~ Moisture, nrow = 1) +
    ggplot2::scale_x_discrete(labels = treatment_labels) +
    ggplot2::scale_fill_manual(values = publication_treatment_colours, guide = "none") +
    ggplot2::labs(title = result$compartment, subtitle = "Alpha diversity", x = NULL, y = "Shannon index", tag = tag) +
    pub_theme() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, vjust = 1))
}

pub_beta_panel <- function(result, tag) {
  axis_label <- function(k) sprintf("PCoA %d (%.1f%%)", k, result$eigen_percent[k])
  ggplot2::ggplot(result$ordination, ggplot2::aes(x = PCoA1, y = PCoA2, colour = Treatment, shape = Moisture)) +
    ggplot2::geom_point(size = 3.8) +
    ggplot2::scale_colour_manual(name = "Seed treatment", values = publication_treatment_colours,
      limits = c("Control", "6PP"), labels = treatment_labels, drop = FALSE) +
    ggplot2::scale_shape_manual(name = "Moisture", values = c("Irrigated" = 16, "Drought" = 17),
      limits = c("Irrigated", "Drought"), drop = FALSE) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, nrow = 1, title.position = "top"),
      shape = ggplot2::guide_legend(order = 2, nrow = 1, title.position = "top")) +
    ggplot2::labs(title = result$compartment, subtitle = "Beta diversity", x = axis_label(1), y = axis_label(2), tag = tag) +
    pub_theme()
}

pub_family_legend <- function(families, colours, title, ncol = 2, font_size = NULL) {
  # Draw the full named colour key explicitly. Its entries do not depend on
  # whether a family is present in the current compartment or ggplot2 version.
  labels <- pub_wrap(families, width = 22)
  row <- ceiling(seq_along(families) / ncol)
  column <- (seq_along(families) - 1) %% ncol
  line_count <- lengths(strsplit(labels, "\n", fixed = TRUE))
  row_height <- vapply(seq_len(max(row)), function(i) max(line_count[row == i]) * 1.10 + 0.28, numeric(1))
  boundaries <- c(0, cumsum(row_height))
  total <- sum(row_height)
  dat <- data.frame(Family = families, Label = labels,
    Colour = unname(colours[families]), Column = column,
    Y = total - (head(boundaries, -1)[row] + tail(boundaries, -1)[row]) / 2)
  if (anyNA(dat$Colour)) stop("A family legend is missing a named colour mapping.")
  family_size <- if (!is.null(font_size)) font_size else if (exists("publication_family_size", inherits = TRUE)) publication_family_size else publication_legend_size
  text_size <- family_size / (72.27 / 25.4)
  ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(x = Column * 1.25 + 0.035, y = Y, fill = Colour),
      shape = 22, size = 4.5, stroke = 0, show.legend = FALSE) +
    ggplot2::scale_fill_identity(guide = "none") +
    ggplot2::geom_text(ggplot2::aes(x = Column * 1.25 + 0.13, y = Y, label = Label),
      hjust = 0, size = text_size, family = publication_font, lineheight = 1.05) +
    ggplot2::annotate("text", x = 0.015, y = total + 0.62, label = title,
      hjust = 0, family = publication_font, fontface = "bold", size = text_size) +
    ggplot2::coord_cartesian(xlim = c(0, ncol * 1.25), ylim = c(0, total + 1.25), expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_family = publication_font) +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 8, 7, 8))
}

pub_bar_panel <- function(df, selected, colours, compartment, tag) {
  legend_columns <- if (exists("publication_main_legend_columns", inherits = TRUE)) publication_main_legend_columns else 3
  family_size <- if (exists("publication_main_family_size", inherits = TRUE)) publication_main_family_size else 18
  axis_size <- if (exists("publication_bar_axis_size", inherits = TRUE)) publication_bar_axis_size else 17
  title_size <- if (exists("publication_bar_title_size", inherits = TRUE)) publication_bar_title_size else 20
  dat <- pub_collapse_families(df, selected) %>%
    dplyr::group_by(Moisture, Treatment, Family) %>%
    dplyr::summarise(Mean_relative_abundance = mean(Abundance), .groups = "drop")
  # Means include zero-abundance samples: family aggregation occurs per sample first.
  levels <- c(selected, "Others")
  dat$Family <- factor(dat$Family, levels = levels)
  dat$Treatment <- factor(dat$Treatment, levels = c("Control", "6PP"))
  dat$Moisture <- factor(dat$Moisture, levels = c("Irrigated", "Drought"))
  totals <- dat %>% dplyr::group_by(Moisture, Treatment) %>%
    dplyr::summarise(Total = sum(Mean_relative_abundance), .groups = "drop")
  if (any(abs(totals$Total - 100) > 1e-7)) stop("Stacked bars do not total 100: ", compartment)
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Treatment, y = Mean_relative_abundance, fill = Family)) +
    ggplot2::geom_col(width = 0.76, show.legend = TRUE) +
    ggplot2::facet_wrap(~ Moisture, nrow = 1) +
    ggplot2::scale_x_discrete(labels = treatment_labels) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::scale_fill_manual(name = publication_family_title, values = colours, limits = levels,
      breaks = levels, labels = pub_wrap(levels, width = 22), drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 3, ncol = legend_columns, byrow = TRUE,
      title.position = "top", override.aes = list(colour = NA))) +
    ggplot2::labs(title = paste(compartment, "community composition"), subtitle = NULL, x = NULL,
      y = "Mean relative\nabundance (%)", tag = tag) +
    pub_theme() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = axis_size),
      axis.title = ggplot2::element_text(size = axis_size + 1),
      strip.text = ggplot2::element_text(size = axis_size),
      plot.title = ggplot2::element_text(size = title_size),
      legend.text = ggplot2::element_text(size = family_size),
      legend.title = ggplot2::element_text(size = family_size),
      legend.position = "bottom", legend.box = "vertical", legend.justification = "left")
  list(plot = p, data = dat, totals = totals,
       legend = pub_family_legend(levels, colours, publication_family_title, ncol = legend_columns, font_size = family_size))
}

pub_pie_panel <- function(df, selected, colours, compartment, tag) {
  dat <- pub_collapse_families(df, selected) %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(Percent = mean(Abundance), .groups = "drop") %>%
    dplyr::arrange(Family == "Others", dplyr::desc(Percent), Family)
  if (abs(sum(dat$Percent) - 100) > 1e-7) stop("Pie does not total 100: ", compartment)
  # Rank slices independently within each compartment. Colour mappings remain
  # shared, but a shared factor order would scramble the other pie's ranks.
  legend_families <- c(setdiff(as.character(dat$Family), "Others"), "Others")
  dat$Family <- factor(dat$Family, levels = legend_families)
  dat$Display_order <- seq_len(nrow(dat))
  dat$Label <- ifelse(dat$Percent >= publication_pie_label_min,
    sprintf("%.1f%%", dat$Percent), "")
  cumulative <- c(0, cumsum(dat$Percent)) / 100
  dat$Mid_angle <- pi / 2 - 2 * pi * (head(cumulative, -1) + tail(cumulative, -1)) / 2
  wedges <- dplyr::bind_rows(lapply(seq_len(nrow(dat)), function(i) {
    theta <- seq(pi / 2 - 2 * pi * cumulative[i], pi / 2 - 2 * pi * cumulative[i + 1],
                 length.out = max(12, ceiling(dat$Percent[i] / 100 * 720)))
    data.frame(x = c(0, cos(theta), 0), y = c(0, sin(theta), 0), Family = dat$Family[i])
  }))
  wedges$Family <- factor(wedges$Family, levels = legend_families)
  labels <- dat[nzchar(dat$Label), , drop = FALSE]
  labels$Side <- ifelse(cos(labels$Mid_angle) >= 0, 1, -1)
  labels$Anchor_x <- 0.98 * cos(labels$Mid_angle)
  labels$Anchor_y <- 0.98 * sin(labels$Mid_angle)
  labels$Label_y <- sin(labels$Mid_angle)
  # Preserve order on each side and enforce enough vertical separation for 18 pt text.
  spread_labels <- function(y, gap = 0.19, lower = -1.10, upper = 1.10) {
    if (length(y) < 2) return(pmax(lower, pmin(upper, y)))
    o <- order(y)
    z <- pmax(lower, pmin(upper, y[o]))
    for (i in 2:length(z)) z[i] <- max(z[i], z[i - 1] + gap)
    if (tail(z, 1) > upper) z <- z - (tail(z, 1) - upper)
    if (z[1] < lower) z <- seq(lower, upper, length.out = length(z))
    out <- numeric(length(y)); out[o] <- z; out
  }
  for (side in c(-1, 1)) {
    idx <- which(labels$Side == side)
    labels$Label_y[idx] <- spread_labels(labels$Label_y[idx])
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = wedges, ggplot2::aes(x = x, y = y, fill = Family, group = Family),
      colour = "white", linewidth = 0.45, show.legend = TRUE) +
    ggplot2::geom_segment(data = labels,
      ggplot2::aes(x = Anchor_x, y = Anchor_y, xend = Side * 1.13, yend = Label_y),
      colour = "#666666", linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_segment(data = labels,
      ggplot2::aes(x = Side * 1.13, xend = Side * 1.19, y = Label_y, yend = Label_y),
      colour = "#666666", linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_text(data = labels,
      ggplot2::aes(x = Side * 1.23, y = Label_y, label = Label, hjust = ifelse(Side > 0, 0, 1)),
      size = publication_pie_label_size, family = publication_font, colour = "black", show.legend = FALSE) +
    ggplot2::coord_fixed(xlim = c(-1.72, 1.72), ylim = c(-1.22, 1.22), expand = FALSE, clip = "off") +
    ggplot2::scale_fill_manual(name = publication_family_title, values = colours, limits = legend_families,
      breaks = legend_families, labels = pub_wrap(legend_families, width = 22), drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4, byrow = TRUE,
      title.position = "left", override.aes = list(colour = NA))) +
    ggplot2::labs(title = compartment, tag = tag) +
    ggplot2::theme_void(base_size = publication_base_size, base_family = publication_font) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = publication_base_size + 1, face = "bold", hjust = 0.5),
      plot.tag = ggplot2::element_text(size = publication_tag_size, face = "bold"),
      plot.margin = ggplot2::margin(15, 12, 15, 12),
      legend.title = ggplot2::element_text(size = publication_legend_size, face = "bold"),
      legend.text = ggplot2::element_text(size = publication_legend_size, lineheight = 1.05),
      legend.key.height = grid::unit(0.5, "cm"),
      legend.key.width = grid::unit(0.5, "cm"),
      legend.spacing.x = grid::unit(0.3, "cm")
    )
  list(plot = p, data = dat, family_order = legend_families)
}

pub_main_table <- function(results) {
  terms <- c("Treatment", "Moisture", "Treatment:Moisture")
  labels <- c("Seed treatment", "Moisture", "Treatment x moisture")
  dplyr::bind_rows(lapply(results, function(r) {
    a <- r$alpha_anova[match(terms, r$alpha_anova$Term), , drop = FALSE]
    b <- r$permanova[match(terms, r$permanova$Term), , drop = FALSE]
    data.frame(Compartment = r$compartment, Effect = labels,
      Alpha_F = a[["F value"]], Alpha_P = a[["Pr(>F)"]],
      PERMANOVA_R2 = b[["R2"]], PERMANOVA_P = b[["Pr(>F)"]], check.names = FALSE)
  }))
}

pub_display_table <- function(tbl) {
  p_format <- function(x) format.pval(x, digits = 3, eps = 0.001)
  out <- data.frame(tbl$Compartment, tbl$Effect, sprintf("%.2f", tbl$Alpha_F),
    p_format(tbl$Alpha_P), sprintf("%.3f", tbl$PERMANOVA_R2), p_format(tbl$PERMANOVA_P),
    stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- c("Compartment", "Effect", "Alpha F", "Alpha P", "PERMANOVA R2", "PERMANOVA P")
  out
}

pub_table_plot <- function(tbl, title, note) {
  display <- pub_display_table(tbl)
  headers <- c("Compartment", "Effect", "Alpha\ndiversity F", "Alpha\ndiversity P", "PERMANOVA\nR2", "PERMANOVA\nP")
  widths <- c(0.17, 0.29, 0.13, 0.13, 0.14, 0.14)
  edges <- c(0, cumsum(widths))
  centres <- (head(edges, -1) + tail(edges, -1)) / 2
  n <- nrow(display)
  cells <- expand.grid(row = seq_len(n), column = seq_len(ncol(display)))
  cells$label <- as.matrix(display)[cbind(cells$row, cells$column)]
  cells$x <- centres[cells$column]
  cells$y <- n - cells$row + 0.5
  cells$face <- "plain"
  cells$face[cells$column == 4 & !is.na(tbl$Alpha_P[cells$row]) & tbl$Alpha_P[cells$row] < 0.05] <- "bold"
  cells$face[cells$column == 6 & !is.na(tbl$PERMANOVA_P[cells$row]) & tbl$PERMANOVA_P[cells$row] < 0.05] <- "bold"
  bg <- data.frame(ymin = 0:(n - 1), ymax = 1:n, fill = rep(c("white", "#F4F6F8"), length.out = n))
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = bg, ggplot2::aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = fill), colour = NA) +
    ggplot2::scale_fill_identity() +
    ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = n, ymax = n + 1.35, fill = "#E7EDF2", colour = NA) +
    ggplot2::geom_hline(yintercept = c(0:n, n + 1.35), colour = "#D9D9D9", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = edges, colour = "#D9D9D9", linewidth = 0.4) +
    ggplot2::annotate("text", x = centres, y = n + 0.65, label = headers,
      family = publication_font, fontface = "bold", size = 5.0, lineheight = 1.05) +
    ggplot2::geom_text(data = cells, ggplot2::aes(x = x, y = y, label = label, fontface = face),
      family = publication_font, size = 4.9) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, n + 1.35), clip = "off", expand = FALSE) +
    ggplot2::labs(title = title, caption = paste(strwrap(note, width = 140), collapse = "\n")) +
    ggplot2::theme_void(base_family = publication_font) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 19, face = "bold", margin = ggplot2::margin(b = 18)),
      plot.caption = ggplot2::element_text(size = 13, hjust = 0, lineheight = 1.2, margin = ggplot2::margin(t = 16)),
      plot.margin = ggplot2::margin(20, 22, 20, 22))
}

pub_write_html_table <- function(tbl, path, title, note) {
  esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
  display <- pub_display_table(tbl)
  rows <- vapply(seq_len(nrow(display)), function(i) {
    cells <- vapply(seq_len(ncol(display)), function(j) {
      value <- esc(display[i, j])
      significant <- (j == 4 && !is.na(tbl$Alpha_P[i]) && tbl$Alpha_P[i] < 0.05) ||
        (j == 6 && !is.na(tbl$PERMANOVA_P[i]) && tbl$PERMANOVA_P[i] < 0.05)
      if (significant) value <- paste0("<strong>", value, "</strong>")
      paste0("<td>", value, "</td>")
    }, character(1))
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  }, character(1))
  html <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>', esc(title),
    '</title><style>body{font:12pt Arial,sans-serif;margin:2cm;color:#000}table{border-collapse:collapse;width:100%}',
    'th,td{border:1px solid #d9d9d9;padding:9px;text-align:center}th{background:#e7edf2}',
    'td:nth-child(1),td:nth-child(2){text-align:left}p{line-height:1.4}</style></head><body><h2>', esc(title),
    '</h2><table><thead><tr>', paste0('<th>', esc(names(display)), '</th>', collapse = ''),
    '</tr></thead><tbody>', paste(rows, collapse = '\n'), '</tbody></table><p>', esc(note), '</p></body></html>')
  writeLines(html, path, useBytes = TRUE)
}

pub_main_figure <- function(results, bars) {
  # Treat each complete panel independently so a tall family legend cannot
  # squash the alpha-diversity and ordination plotting areas in the same row.
  standalone <- function(p) patchwork::wrap_elements(full = ggplot2::ggplotGrob(p))
  composition <- function(bar) {
    patchwork::wrap_plots(
      standalone(bar$plot + ggplot2::theme(legend.position = "none")), standalone(bar$legend),
      ncol = 1, heights = c(1, 0.70), guides = "keep")
  }
  patchwork::wrap_plots(
    A = standalone(pub_alpha_panel(results$Rhizosphere, "A")),
    B = standalone(pub_beta_panel(results$Rhizosphere, "B") + ggplot2::theme(legend.position = "bottom")),
    C = composition(bars$Rhizosphere),
    D = standalone(pub_alpha_panel(results$Root, "D")),
    E = standalone(pub_beta_panel(results$Root, "E") + ggplot2::theme(legend.position = "bottom")),
    F = composition(bars$Root),
    ncol = 3, widths = c(1, 1.02, 1.75), guides = "keep")
}

pub_supplement_figure <- function(pies, colours) {
  panel <- function(compartment) {
    key <- pub_family_legend(pies[[compartment]]$family_order, colours,
      publication_family_title, ncol = 2)
    patchwork::wrap_plots(
      pies[[compartment]]$plot + ggplot2::theme(legend.position = "none"), key,
      ncol = 1, heights = c(1, 0.70), guides = "keep")
  }
  patchwork::wrap_plots(panel("Rhizosphere"), panel("Root"), ncol = 2, guides = "keep")
}

pub_build_outputs <- function(ps_input) {
  dir.create(publication_out_dir, recursive = TRUE, showWarnings = FALSE)
  data_dir <- file.path(publication_out_dir, "Supporting_Data")
  dir.create(data_dir, showWarnings = FALSE)
  compartments <- c("Rhizosphere", "Root")
  results <- stats::setNames(lapply(seq_along(compartments), function(i) {
    ps <- pub_filter_compartment(ps_input, compartments[i])
    pub_analyse_compartment(ps, compartments[i], publication_seed + i)
  }), compartments)
  families <- stats::setNames(lapply(compartments, function(compartment) {
    pub_family_data(results[[compartment]]$ps, compartment)
  }), compartments)
  all_family <- dplyr::bind_rows(families)
  rank_table <- all_family %>% dplyr::group_by(Source, Family) %>%
    dplyr::summarise(Mean_percent = mean(Abundance), .groups = "drop") %>%
    tidyr::complete(Source = compartments, Family, fill = list(Mean_percent = 0)) %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(Mean_across_compartments = mean(Mean_percent), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Mean_across_compartments), Family)
  selected <- head(rank_table$Family, stacked_top_n)
  pie_selected <- lapply(families, function(df) {
    rank <- df %>% dplyr::group_by(Family) %>%
      dplyr::summarise(Percent = mean(Abundance), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(Percent), Family)
    head(rank$Family, pie_top_n)
  })
  all_names <- unique(c(selected, unlist(pie_selected, use.names = FALSE)))
  all_names <- setdiff(all_names, "Others")
  colours <- pub_palette(all_names)
  bars <- stats::setNames(lapply(seq_along(compartments), function(i) {
    compartment <- compartments[i]
    pub_bar_panel(families[[compartment]], selected, colours, compartment, c("C", "F")[i])
  }), compartments)
  main <- pub_main_figure(results, bars)
  main_stem <- file.path(publication_out_dir, paste0("Figure", publication_number, "_", publication_slug, "_Main_A-F"))
  pub_save_three(main, main_stem, publication_main_width, publication_main_height)

  pies <- stats::setNames(lapply(seq_along(compartments), function(i) {
    compartment <- compartments[i]
    pub_pie_panel(families[[compartment]], pie_selected[[compartment]], colours, compartment, LETTERS[i])
  }), compartments)
  supp <- pub_supplement_figure(pies, colours)
  supp_stem <- file.path(publication_out_dir, paste0("FigureS", publication_number, "_", publication_slug, "_Pies_A-B"))
  pub_save_three(supp, supp_stem, publication_pie_width, publication_pie_height)

  tbl <- pub_main_table(results)
  table_stem <- file.path(publication_out_dir, paste0("Table", publication_number, "_", publication_slug, "_Main"))
  note <- paste0("Shannon alpha diversity: two-way ANOVA. PERMANOVA: Bray-Curtis dissimilarity of square-root relative abundances; ",
    publication_permutations, " unrestricted permutations; sequential term tests. P < 0.05 is bold. Models retain the original scripts' Treatment * Moisture design.")
  utils::write.csv(tbl, paste0(table_stem, ".csv"), row.names = FALSE)
  pub_write_html_table(tbl, paste0(table_stem, ".html"), paste0("Table ", publication_number, ". ", publication_title), note)
  pub_save_three(pub_table_plot(tbl, paste0("Table ", publication_number, ". ", publication_title), note), table_stem, 14, 5.2)

  for (compartment in compartments) {
    r <- results[[compartment]]
    saveRDS(r$ps, file.path(data_dir, paste0(compartment, "_filtered_phyloseq.rds")))
    utils::write.csv(r$alpha, file.path(data_dir, paste0(compartment, "_Alpha_Data.csv")), row.names = FALSE)
    utils::write.csv(r$ordination, file.path(data_dir, paste0(compartment, "_PCoA_Coordinates.csv")), row.names = FALSE)
    utils::write.csv(r$eigenvalues, file.path(data_dir, paste0(compartment, "_PCoA_Eigenvalues.csv")), row.names = FALSE)
    utils::write.csv(r$alpha_anova, file.path(data_dir, paste0(compartment, "_Alpha_ANOVA_Full.csv")), row.names = FALSE)
    utils::write.csv(r$permanova, file.path(data_dir, paste0(compartment, "_PERMANOVA_Full.csv")), row.names = FALSE)
    utils::write.csv(bars[[compartment]]$data, file.path(data_dir, paste0(compartment, "_Stacked_Bar_Data.csv")), row.names = FALSE)
    utils::write.csv(bars[[compartment]]$totals, file.path(data_dir, paste0(compartment, "_Bar_Percent_Check.csv")), row.names = FALSE)
    utils::write.csv(pies[[compartment]]$data, file.path(data_dir, paste0(compartment, "_Pie_Data.csv")), row.names = FALSE)
    utils::write.csv(families[[compartment]], file.path(data_dir, paste0(compartment, "_Per_Sample_Family_Percent.csv")), row.names = FALSE)
  }
  counts <- dplyr::bind_rows(lapply(results, function(r) {
    r$metadata %>% dplyr::count(Source, Moisture, Treatment, name = "N")
  }))
  utils::write.csv(counts, file.path(data_dir, "Sample_Counts_After_Filtering.csv"), row.names = FALSE)
  rank_table$Shown_in_main <- rank_table$Family %in% selected
  utils::write.csv(rank_table, file.path(data_dir, "Family_Ranking.csv"), row.names = FALSE)
  utils::write.csv(data.frame(Family = names(colours), Colour = unname(colours)), file.path(data_dir, "Family_Colour_Key.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(publication_out_dir, "sessionInfo.txt"))
  saveRDS(list(Main = main, Supplemental = supp, Table = tbl), file.path(data_dir, "Publication_Plots_and_Table.rds"))
  message("Main A-F figure, supplemental A-B pies and one main table exported to: ", publication_out_dir)
  invisible(list(results = results, main = main, supplemental = supp, table = tbl))
}

publication_results <- pub_build_outputs(ps_full)

# Optional original outputs, including DESeq2 and combined-compartment analyses.
if (run_legacy_outputs) {
  ###############################################################################
  # 10) Run analyses
  ###############################################################################
  # Combined
  run_phyloseq_analysis(
    ps_in = ps_full,
    analysis_name = "Combined",
    out_dir = file.path(main_out_dir, "Combined"),
    include_source_in_model = TRUE
  )

  # Root-only
  ps_root <- subset_samples(ps_full, Source == "Root")
  run_phyloseq_analysis(
    ps_in = ps_root,
    analysis_name = "Root",
    out_dir = file.path(main_out_dir, "Root"),
    include_source_in_model = FALSE
  )

  # Rhizosphere-only
  ps_rhizo <- subset_samples(ps_full, Source == "Rhizosphere")
  run_phyloseq_analysis(
    ps_in = ps_rhizo,
    analysis_name = "Rhizosphere",
    out_dir = file.path(main_out_dir, "Rhizosphere"),
    include_source_in_model = FALSE
  )


  # =========================
  # RHIZOSPHERE COMBINED PLOT
  # =========================
  ps_rhizo <- subset_samples(ps_full, Source == "Rhizosphere")

  ps_rhizo_family <- tax_glom(ps_rhizo, taxrank = "Family", NArm = FALSE)
  ps_rhizo_rel <- transform_sample_counts(ps_rhizo_family, function(x) 100 * x / sum(x))

  df_rhizo <- psmelt(ps_rhizo_rel)
  df_rhizo$Family <- clean_taxonomy_labels(df_rhizo$Family, rank = "Family")

  top_families_rhizo <- df_rhizo %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_ab)) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::pull(Family)

  pA_rhizo <- plot_family_overall_bar(
    ps = ps_rhizo,
    analysis_name = "Rhizosphere",
    out_dir = file.path(main_out_dir, "Rhizosphere"),
    top_families = top_families_rhizo
  )

  pB_rhizo <- plot_family_stacked_portrait(
    ps = ps_rhizo,
    analysis_name = "Rhizosphere",
    out_dir = file.path(main_out_dir, "Rhizosphere"),
    top_families = top_families_rhizo
  )

  combined_rhizo <- patchwork::wrap_plots(
    pA_rhizo + theme(legend.position = "none"),
    pB_rhizo,
    ncol = 2,
    widths = c(1, 1.2)
  ) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        vjust = 1
      )
    )

  ggsave(
    file.path(main_out_dir, "Rhizosphere", "Family_Combined_Rhizosphere.png"),
    combined_rhizo,
    width = 14,
    height = 8,
    dpi = 300
  )

  # ==================
  # ROOT COMBINED PLOT
  # ==================
  ps_root <- subset_samples(ps_full, Source == "Root")

  ps_root_family <- tax_glom(ps_root, taxrank = "Family", NArm = FALSE)
  ps_root_rel <- transform_sample_counts(ps_root_family, function(x) 100 * x / sum(x))

  df_root <- psmelt(ps_root_rel)
  df_root$Family <- clean_taxonomy_labels(df_root$Family, rank = "Family")

  top_families_root <- df_root %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_ab)) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::pull(Family)

  pA_root <- plot_family_overall_bar(
    ps = ps_root,
    analysis_name = "Root",
    out_dir = file.path(main_out_dir, "Root"),
    top_families = top_families_root
  )

  pB_root <- plot_family_stacked_portrait(
    ps = ps_root,
    analysis_name = "Root",
    out_dir = file.path(main_out_dir, "Root"),
    top_families = top_families_root
  )

  combined_root <- patchwork::wrap_plots(
    pA_root + theme(legend.position = "none"),
    pB_root,
    ncol = 2,
    widths = c(1, 1.2)
  ) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        vjust = 1
      )
    )

  ggsave(
    file.path(main_out_dir, "Root", "Family_Combined_Root.png"),
    combined_root,
    width = 14,
    height = 8,
    dpi = 300
  )

  message("All analyses completed successfully.")
  message("Output folder: ", main_out_dir)
}
