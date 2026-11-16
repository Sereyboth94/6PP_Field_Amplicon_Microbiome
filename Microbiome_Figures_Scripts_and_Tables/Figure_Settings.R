# Journal artwork is drawn at 190 mm width; text sizes are physical points.
publication_font <- "Arial"
publication_base_size <- 9
publication_axis_size <- 8
publication_legend_size <- 8
publication_family_size <- 8.5
publication_main_family_size <- 8.5
publication_bar_axis_size <- 8
publication_bar_title_size <- 9
publication_main_legend_columns <- 3
publication_tag_size <- 12
publication_dpi <- 600
publication_main_width <- 190 / 25.4
publication_main_height <- 9
publication_pie_width <- 190 / 25.4
publication_pie_height <- 9.5
publication_pie_family_size <- 8.5
publication_pie_order <- "within_panel"
publication_pie_label_min <- 5
publication_pie_label_size <- 2.8
publication_treatment_colours <- c(Control = "#D55E00", "6PP" = "#0072B2")
treatment_labels <- c(Control = "Control", "6PP" = "6-PP")
pie_top_n <- 15
stacked_top_n <- 15


# Cairo matches the approved reference exports. ragg is an optional alternative.
publication_raster_backend <- "cairo"  # "cairo" or "ragg"
publication_export_svg <- FALSE        # TRUE adds editable SVG; requires svglite
publication_main_column_widths <- c(1, 1.02, 1.04)
publication_main_row_heights <- c(3.0, 1.8, 3.0, 1.8)
publication_pie_panel_heights <- c(1, 1.9)

# Optional named family-colour overrides, shared by bars and pies within a dataset.
# Example: publication_family_colours <- c(Rhizobiaceae = "#D95F02")
publication_family_colours <- character(0)

# These controls change presentation, not Shannon/PERMANOVA model fitting.
# Keep the supplied defaults to reproduce the approved tall figures.
