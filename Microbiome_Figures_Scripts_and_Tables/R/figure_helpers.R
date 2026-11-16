# Shared plotting functions for journal figures.
# Derived from the supplied field microbiome scripts; no analysis is called by loading this file.

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
    } else if (identical(publication_raster_backend, "ragg")) {
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
  if (isTRUE(publication_export_svg)) {
    filename <- paste0(stem, ".svg")
    svglite::svglite(filename, width = width, height = height, bg = "white")
    tryCatch(draw_plot(), finally = grDevices::dev.off())
    files <- c(files, filename)
  }
  if (any(!file.exists(files)) || any(file.info(files)$size == 0)) {
    stop("An export is missing or empty: ", stem)
  }
  # Optional embedding audit; pdftools is not needed to create embedded fonts.
  if (requireNamespace("pdftools", quietly = TRUE)) {
    font_audit <- pdftools::pdf_fonts(paste0(stem, ".pdf"))
    utils::write.csv(font_audit, file.path(publication_data_dir, paste0(basename(stem), "_PDF_Fonts.csv")), row.names = FALSE)
    if (nrow(font_audit) > 0 && any(!font_audit$embedded)) {
      warning("A PDF font was not embedded; inspect the PDF_Fonts.csv audit: ", stem)
    }
  }
  invisible(files)
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
  row_height <- vapply(seq_len(max(row)), function(i)
    max(line_count[row == i]) * 1.75 - 0.25, numeric(1))
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
      shape = 22, size = 2.5, stroke = 0, show.legend = FALSE) +
    ggplot2::scale_fill_identity(guide = "none") +
    ggplot2::geom_text(ggplot2::aes(x = Column * 1.25 + 0.13, y = Y, label = Label),
      hjust = 0, size = text_size, family = publication_font, lineheight = 1.05) +
    ggplot2::annotate("text", x = 0.015, y = total + 1.15, label = title,
      hjust = 0, family = publication_font, fontface = "bold", size = text_size) +
    ggplot2::coord_cartesian(xlim = c(0, ncol * 1.25), ylim = c(0, total + 2.20), expand = FALSE, clip = "off") +
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

pub_complete_family_percent <- function(df, compartment) {
  required <- c("SampleID", "Treatment", "Moisture", "Source", "Family", "Abundance")
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("Missing family-data columns for ", compartment, ": ", paste(missing, collapse = ", "))
  if (!nrow(df)) stop("No family data for ", compartment)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in c("SampleID", "Treatment", "Moisture", "Source")) df[[nm]] <- as.character(df[[nm]])
  if (anyNA(df[, c("SampleID", "Treatment", "Moisture", "Source")])) stop("Missing sample metadata in ", compartment)
  if (any(!nzchar(df$SampleID))) stop("Empty sample identifier in ", compartment)
  if (any(df$Source != compartment)) stop("Mixed compartment labels in ", compartment, " family data")
  if (any(!df$Treatment %in% c("Control", "6PP"))) stop("Unknown treatment in ", compartment, "; expected Control or 6PP")
  if (any(!df$Moisture %in% c("Irrigated", "Drought"))) stop("Unknown moisture label in ", compartment)
  if (!is.numeric(df$Abundance) || any(!is.finite(df$Abundance)) || any(df$Abundance < 0)) {
    stop("Family percentages must be finite, non-negative numbers in ", compartment)
  }
  df$Family <- pub_clean_family(df$Family)
  if (any(df$Family == "Others")) stop("Use the uncategorised per-sample family table, not previously collapsed pie data, for ", compartment)
  metadata <- unique(df[, c("SampleID", "Treatment", "Moisture", "Source"), drop = FALSE])
  if (anyDuplicated(metadata$SampleID)) stop("Conflicting metadata for a sample in ", compartment)
  collapsed <- df %>% dplyr::group_by(SampleID, Family) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
  totals <- collapsed %>% dplyr::group_by(SampleID) %>%
    dplyr::summarise(Total = sum(Abundance), .groups = "drop")
  if (any(abs(totals$Total - 100) > 1e-6)) {
    stop("Per-sample family percentages do not total 100 in ", compartment,
      ". Use *_Per_Sample_Family_Percent.csv from the publication workflow.")
  }
  # A missing sample-family row is zero abundance, not a reason to reduce
  # the denominator used for a family's mean.
  tidyr::expand_grid(SampleID = metadata$SampleID, Family = sort(unique(collapsed$Family))) %>%
    dplyr::left_join(collapsed, by = c("SampleID", "Family")) %>%
    dplyr::mutate(Abundance = dplyr::coalesce(Abundance, 0)) %>%
    dplyr::left_join(metadata, by = "SampleID")
}

pub_select_pie_families <- function(df) {
  rank <- df %>% dplyr::group_by(Family) %>%
    dplyr::summarise(Percent = mean(Abundance), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Percent), Family)
  head(as.character(rank$Family), pie_top_n)
}

pub_pie_panel <- function(df, selected, colours, compartment, tag,
                          treatment = NULL, family_order = NULL) {
  df <- pub_complete_family_percent(df, compartment)
  selected <- unique(setdiff(as.character(selected), "Others"))
  if (!length(selected)) stop("No named families selected for ", compartment)
  if (!is.null(treatment)) {
    if (length(treatment) != 1L || !treatment %in% c("Control", "6PP")) stop("Invalid pie treatment subset")
    df <- df[df$Treatment == treatment, , drop = FALSE]
  }
  metadata <- unique(df[, c("SampleID", "Treatment", "Moisture", "Source"), drop = FALSE])
  if (!nrow(metadata)) stop("No retained samples for ", compartment, " / ", treatment)
  n_samples <- nrow(metadata)
  counts <- c(Irrigated = sum(metadata$Moisture == "Irrigated"), Drought = sum(metadata$Moisture == "Drought"))
  families <- c(selected, "Others")
  collapsed <- pub_collapse_families(df, selected)
  dat <- tidyr::expand_grid(SampleID = metadata$SampleID, Family = families) %>%
    dplyr::left_join(collapsed[, c("SampleID", "Family", "Abundance")], by = c("SampleID", "Family")) %>%
    dplyr::mutate(Abundance = dplyr::coalesce(Abundance, 0)) %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(Percent = mean(Abundance), .groups = "drop")
  if (abs(sum(dat$Percent) - 100) > 1e-7) stop("Pie does not total 100: ", compartment, " / ", treatment)
  sort_mode <- if (exists("publication_pie_order", inherits = TRUE)) publication_pie_order else "within_panel"
  if (!sort_mode %in% c("within_panel", "compartment")) stop("publication_pie_order must be within_panel or compartment")
  if (sort_mode == "within_panel") {
    dat <- dat %>% dplyr::arrange(Family == "Others", dplyr::desc(Percent), Family)
    legend_families <- c(setdiff(as.character(dat$Family), "Others"), "Others")
  } else {
    if (is.null(family_order)) family_order <- families
    legend_families <- c(setdiff(as.character(family_order), "Others"), "Others")
    if (!setequal(legend_families, families)) stop("Inconsistent common family order")
    dat <- dat[match(legend_families, dat$Family), , drop = FALSE]
  }
  if (anyNA(colours[legend_families])) stop("Missing named family colour for ", compartment)
  dat$Family <- factor(dat$Family, levels = legend_families)
  dat$Display_order <- seq_len(nrow(dat))
  dat$Label <- ifelse(dat$Percent >= publication_pie_label_min, sprintf("%.1f%%", dat$Percent), "")
  cumulative <- c(0, cumsum(dat$Percent)) / 100
  dat$Mid_angle <- pi / 2 - 2 * pi * (head(cumulative, -1) + tail(cumulative, -1)) / 2
  wedges <- dplyr::bind_rows(lapply(which(dat$Percent > 0), function(i) {
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
  spread_labels <- function(y, gap = 0.36, lower = -1.30, upper = 1.30) {
    if (length(y) < 2) return(pmax(lower, pmin(upper, y)))
    o <- order(y); z <- pmax(lower, pmin(upper, y[o]))
    for (i in 2:length(z)) z[i] <- max(z[i], z[i - 1] + gap)
    if (tail(z, 1) > upper) z <- z - (tail(z, 1) - upper)
    if (z[1] < lower) z <- seq(lower, upper, length.out = length(z))
    out <- numeric(length(y)); out[o] <- z; out
  }
  for (side in c(-1, 1)) {
    idx <- which(labels$Side == side)
    labels$Label_y[idx] <- spread_labels(labels$Label_y[idx])
  }
  group_label <- if (is.null(treatment)) "All treatments" else if (treatment == "6PP") "6-PP" else "Control"
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
    ggplot2::coord_fixed(xlim = c(-1.72, 1.72), ylim = c(-1.40, 1.40), expand = FALSE, clip = "off") +
    ggplot2::scale_fill_manual(name = publication_family_title, values = colours, limits = legend_families,
      breaks = legend_families, labels = pub_wrap(legend_families, width = 22), drop = FALSE) +
    ggplot2::labs(title = compartment, subtitle = sprintf("%s (n = %d)", group_label, n_samples), tag = tag) +
    ggplot2::theme_void(base_size = publication_base_size, base_family = publication_font) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = publication_base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = publication_axis_size, hjust = 0.5,
        margin = ggplot2::margin(t = 4, b = 8)),
      plot.tag = ggplot2::element_text(size = publication_tag_size, face = "bold"),
      plot.margin = ggplot2::margin(15, 12, 15, 12),
      legend.position = "none"
    )
  subset_label <- if (is.null(treatment)) "All" else treatment
  dat$Panel <- tag; dat$Source <- compartment; dat$Treatment_subset <- subset_label
  dat$N_samples <- n_samples; dat$N_Irrigated <- unname(counts["Irrigated"]); dat$N_Drought <- unname(counts["Drought"])
  metadata$Panel <- tag
  list(plot = p, data = dat, family_order = legend_families, metadata = metadata,
    n = n_samples, group = subset_label, compartment = compartment, tag = tag)
}

pub_build_supplementary_pies <- function(families, colours = NULL,
                                         selected_by_compartment = NULL,
                                         write_outputs = TRUE) {
  compartments <- c("Rhizosphere", "Root")
  if (!all(compartments %in% names(families))) stop("Both Rhizosphere and Root family tables are required")
  families <- stats::setNames(lapply(compartments, function(nm) {
    pub_complete_family_percent(families[[nm]], nm)
  }), compartments)
  if (is.null(selected_by_compartment)) selected_by_compartment <- lapply(families, pub_select_pie_families)
  if (!all(compartments %in% names(selected_by_compartment))) stop("Missing compartment-specific family selection")
  selected_by_compartment <- lapply(selected_by_compartment[compartments], function(x) unique(setdiff(as.character(x), "Others")))
  required_colours <- unique(c(unlist(selected_by_compartment, use.names = FALSE), "Others"))
  if (is.null(colours)) colours <- pub_palette(setdiff(required_colours, "Others"))
  if (is.null(names(colours)) || anyDuplicated(names(colours)) || anyNA(colours[required_colours])) {
    stop("The named colour key does not cover all selected pie families")
  }
  specification <- data.frame(Panel = LETTERS[1:6],
    Source = rep(compartments, each = 3),
    Treatment_subset = rep(c("All", "Control", "6PP"), 2), stringsAsFactors = FALSE)
  pies <- stats::setNames(lapply(seq_len(nrow(specification)), function(i) {
    s <- specification[i, ]
    selected <- selected_by_compartment[[s$Source]]
    pub_pie_panel(families[[s$Source]], selected, colours, s$Source, s$Panel,
      treatment = if (s$Treatment_subset == "All") NULL else s$Treatment_subset,
      family_order = c(selected, "Others"))
  }), specification$Panel)
  combined <- dplyr::bind_rows(lapply(pies, function(x) {
    d <- x$data; d$Family <- as.character(d$Family); d
  }))
  combined$Figure <- paste0("Figure S", publication_number)
  combined$Dataset <- publication_slug
  combined <- combined[, c("Figure", "Dataset", "Panel", "Source", "Treatment_subset",
    "N_samples", "N_Irrigated", "N_Drought", "Display_order", "Family", "Percent", "Label", "Mid_angle")]
  panel_counts <- unique(combined[, c("Figure", "Dataset", "Panel", "Source", "Treatment_subset",
    "N_samples", "N_Irrigated", "N_Drought")])
  selection <- dplyr::bind_rows(lapply(compartments, function(nm) {
    data.frame(Source = nm, Overview_rank = seq_along(selected_by_compartment[[nm]]),
      Family = selected_by_compartment[[nm]], stringsAsFactors = FALSE)
  }))
  family_word <- if (grepl("Fungi", publication_slug)) "families or family-level taxonomic groups" else "families"
  order_text <- if (publication_pie_order == "within_panel") {
    "Named groups are ordered clockwise from largest to smallest mean contribution within each panel, with Others last."
  } else {
    "Named groups follow a common clockwise order within each compartment, ranked by the overall profile in A or D, with Others last."
  }
  caption <- paste0("Figure S", publication_number, ". ", publication_title,
    ": overall and treatment-specific family composition. ",
    "(A) Overall rhizosphere community, (B) rhizosphere Control, (C) rhizosphere 6-PP, ",
    "(D) overall root community, (E) root Control, and (F) root 6-PP. ",
    "Each pie shows the arithmetic mean of per-sample relative abundances, giving every retained sample equal weight. ",
    "Overall panels pool both seed treatments and moisture regimes; treatment-specific panels pool irrigated and drought samples within the indicated seed treatment. ",
    "Up to ", pie_top_n, " of the most abundant ", family_word,
    " are selected independently for each compartment from its overall profile and retained in all three panels for that compartment. Remaining groups are pooled as Others. ",
    order_text, " Family colours are shared across panels and match the corresponding main figure. ",
    "Percentage labels are shown for components contributing at least ", publication_pie_label_min,
    "%. Retained sample numbers are shown in each panel and moisture-specific counts are provided with the supporting data. ",
    if (grepl("Fungi", publication_slug)) "Labels containing fam. incertae sedis represent unresolved family-level assignments within the named order. " else "",
    "These relative-abundance summaries are descriptive and do not establish family-specific differential abundance, absolute population changes, or moisture-specific treatment effects. ",
    "6-PP, 6-pentyl-alpha-pyrone.")
  figure <- pub_supplement_figure(pies, colours)
  if (isTRUE(write_outputs)) {
    dir.create(publication_out_dir, recursive = TRUE, showWarnings = FALSE)
    data_dir <- publication_data_dir
    dir.create(data_dir, showWarnings = FALSE)
    stem <- file.path(publication_out_dir, paste0("FigureS", publication_number, "_", publication_slug, "_Pies_A-F"))
    pub_save_three(figure, stem, publication_pie_width, publication_pie_height)
    writeLines(caption, file.path(data_dir, paste0(basename(stem), "_Legend.txt")), useBytes = TRUE)
    utils::write.csv(combined, file.path(data_dir, "Supplementary_Pies_A-F_Data.csv"), row.names = FALSE)
    utils::write.csv(panel_counts, file.path(data_dir, "Supplementary_Pies_A-F_Sample_Counts.csv"), row.names = FALSE)
    utils::write.csv(selection, file.path(data_dir, "Supplementary_Pies_A-F_Family_Selection.csv"), row.names = FALSE)
    for (i in seq_len(nrow(specification))) {
      s <- specification[i, ]
      suffix <- if (s$Treatment_subset == "All") "Overall" else s$Treatment_subset
      utils::write.csv(pies[[s$Panel]]$data,
        file.path(data_dir, paste0(s$Source, "_", suffix, "_Pie_Data.csv")), row.names = FALSE)
    }
    message("Supplementary A-F pies exported to: ", publication_out_dir)
  }
  invisible(list(figure = figure, panels = pies, data = combined,
    sample_counts = panel_counts, selection = selection, colours = colours, caption = caption))
}
