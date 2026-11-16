# Use one full-width family key beneath each compartment row in main figures.
journal_main <- function(results, bars) {
  wrap <- function(p) patchwork::wrap_elements(full = p)
  compact <- function(p) p + ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5),
    legend.key.size = grid::unit(3, "mm"), legend.spacing.y = grid::unit(1.5, "mm"))
  a <- compact(pub_alpha_panel(results$Rhizosphere, "A"))
  b <- compact(pub_beta_panel(results$Rhizosphere, "B") + ggplot2::theme(legend.position = "bottom"))
  c <- compact(bars$Rhizosphere$plot + ggplot2::labs(title = "Rhizosphere\ncommunity composition") + ggplot2::theme(legend.position = "none"))
  d <- compact(pub_alpha_panel(results$Root, "D"))
  e <- compact(pub_beta_panel(results$Root, "E") + ggplot2::theme(legend.position = "bottom"))
  f <- compact(bars$Root$plot + ggplot2::labs(title = "Root\ncommunity composition") + ggplot2::theme(legend.position = "none"))
  compact_key <- ggplot2::theme(legend.spacing.y = grid::unit(0, "pt"),
    legend.margin = ggplot2::margin(0), legend.box.margin = ggplot2::margin(0),
    legend.box.spacing = grid::unit(2, "pt"))
  compact_guides <- ggplot2::guides(colour = ggplot2::guide_legend(title = NULL, nrow = 1, order = 1),
    shape = ggplot2::guide_legend(title = NULL, nrow = 1, order = 2))
  b <- b + compact_key + compact_guides
  e <- e + compact_key + compact_guides
  # Keep point sizes proportional to the journal-sized drawing area.
  a$layers[[2]]$aes_params$size <- 1.2; d$layers[[2]]$aes_params$size <- 1.2
  b$layers[[1]]$aes_params$size <- 1.8; e$layers[[1]]$aes_params$size <- 1.8
  row1 <- patchwork::wrap_plots(wrap(a), wrap(b), wrap(c), ncol = 3,
    widths = publication_main_column_widths, guides = "keep")
  row2 <- patchwork::wrap_plots(wrap(d), wrap(e), wrap(f), ncol = 3,
    widths = publication_main_column_widths, guides = "keep")
  key1 <- pub_family_legend(levels(bars$Rhizosphere$data$Family), current_colours,
    publication_family_title, ncol = publication_main_legend_columns, font_size = publication_main_family_size)
  key2 <- pub_family_legend(levels(bars$Root$data$Family), current_colours,
    publication_family_title, ncol = publication_main_legend_columns, font_size = publication_main_family_size)
  patchwork::wrap_plots(row1, wrap(key1), row2, wrap(key2), ncol = 1,
    heights = publication_main_row_heights, guides = "keep")
}

# A single-column legend under each pie keeps family text readable at 190 mm.
pub_supplement_figure <- function(pies, colours) {
  stopifnot(identical(names(pies), LETTERS[1:6]))
  wrap <- function(p) patchwork::wrap_elements(full = p)
  panels <- lapply(pies, function(pie) {
    key <- pub_family_legend(pie$family_order, colours, publication_family_title,
      ncol = 1, font_size = publication_pie_family_size)
    p <- pie$plot + ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4))
    patchwork::wrap_plots(wrap(p), wrap(key), ncol = 1, heights = publication_pie_panel_heights, guides = "keep")
  })
  patchwork::wrap_plots(panels, ncol = 3, nrow = 2, guides = "keep")
}

