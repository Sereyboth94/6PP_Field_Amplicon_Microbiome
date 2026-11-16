# The filtering, ANOVA, distance, permutation and diagnostic calculations
# below are retained from the validated journal-figure pipeline.
compute_study <- function(input_root) {
  output_root <- ""  # No input or output files are changed by this function.
  publication_permutations <- 999
  publication_seed <- 20261105
spec <- data.frame(Host = c("Pasture", "Pasture", "Wheat", "Wheat"),
  Marker = c("Bacteria", "Fungi", "Bacteria", "Fungi"),
  Metadata = c("Metadat_Bacteria_Field.txt", "Metadat_Fungi_Pasture.txt",
    "Metadat_Bacteria_Field.txt", "Metadat_Fungal_Wheat.txt"), stringsAsFactors = FALSE)
all_counts <- list(); all_alpha <- list(); all_permanova <- list(); all_diagnostics <- list()
all_libraries <- list(); all_excluded <- list(); all_profiles <- list(); all_objects <- list()

for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]; slug <- paste(s$Host, s$Marker, sep = "_")
  publication_number <- i; publication_slug <- slug
  publication_title <- paste(s$Host, if (s$Marker == "Bacteria") "bacterial" else "fungal", "communities")
  publication_family_title <- if (s$Marker == "Bacteria") "Bacterial families" else "Fungal families"
  publication_out_dir <- file.path(output_root, slug)
  data_dir <- file.path(publication_out_dir, "Supporting_Data")
  input_dir <- file.path(input_root, s$Host, s$Marker)
  counts <- as.matrix(read.delim(file.path(input_dir, "ASV_filtered.txt"), row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE))
  storage.mode(counts) <- "numeric"
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == round(counts)))
  taxonomy <- read.delim(file.path(input_dir, "Taxonomy_filtered.txt"), row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE)
  md <- read.delim(file.path(input_dir, s$Metadata), check.names = FALSE, stringsAsFactors = FALSE)
  names(md)[1] <- "SampleID"; rownames(md) <- md$SampleID
  md$Treatment[md$Treatment %in% c("FVOC", "6-PP", "6PP")] <- "6PP"
  stopifnot(all(md$Treatment %in% c("Control", "6PP")), all(colnames(counts) %in% md$SampleID),
    all(rownames(counts) %in% rownames(taxonomy)))
  min_depth <- if (slug == "Wheat_Bacteria") 100 else 500
  results <- list(); families <- list()
  for (j in seq_along(c("Rhizosphere", "Root"))) {
    compartment <- c("Rhizosphere", "Root")[j]
    ids <- md$SampleID[md$Source == compartment]
    x <- counts[, ids, drop = FALSE]
    before <- colSums(x); excluded <- names(before)[before < min_depth]
    if (length(excluded)) all_excluded[[paste(slug, compartment)]] <- data.frame(Dataset = slug,
      Source = compartment, SampleID = excluded, Input_ASV_count = before[excluded], Threshold = min_depth)
    x <- x[, before >= min_depth, drop = FALSE]
    x <- x[rowSums(x) >= 10 & rowMeans(x > 0) >= 0.05, , drop = FALSE]
    stopifnot(all(colSums(x) > 0))
    metadata <- md[colnames(x), , drop = FALSE]
    metadata$Treatment <- factor(metadata$Treatment, levels = c("Control", "6PP"))
    metadata$Moisture <- factor(metadata$Moisture, levels = c("Irrigated", "Drought"))
    stopifnot(!anyNA(metadata[, c("Treatment", "Moisture")]))
    alpha <- metadata; alpha$Shannon <- vegan::diversity(t(x), index = "shannon")
    fit <- aov(Shannon ~ Treatment * Moisture, data = alpha)
    a <- as.data.frame(summary(fit)[[1]]); a$Term <- trimws(rownames(a))
    a$Dataset <- slug; a$Source <- compartment; a$Residual_df <- df.residual(fit)
    transformed <- sqrt(sweep(x, 2, colSums(x), "/"))
    distance <- vegan::vegdist(t(transformed), method = "bray")
    ord <- ape::pcoa(distance)
    coords <- data.frame(SampleID = rownames(ord$vectors), PCoA1 = ord$vectors[, 1], PCoA2 = ord$vectors[, 2])
    coords <- dplyr::left_join(coords, metadata, by = "SampleID")
    set.seed(publication_seed + j)
    p <- as.data.frame(vegan::adonis2(distance ~ Treatment * Moisture,
      data = metadata, permutations = publication_permutations, by = "terms"))
    p$Term <- rownames(p); p$Dataset <- slug; p$Source <- compartment
    results[[compartment]] <- list(alpha = alpha, ordination = coords, alpha_anova = a,
      permanova = p, eigen_percent = 100 * ord$values$Relative_eig[1:2],
      compartment = compartment, seed = publication_seed + j)
    family_counts <- rowsum(x, group = pub_clean_family(taxonomy[rownames(x), "Family"]), reorder = FALSE)
    family_percent <- 100 * sweep(family_counts, 2, colSums(family_counts), "/")
    df <- as.data.frame(as.table(family_percent), responseName = "Abundance", stringsAsFactors = FALSE)
    names(df)[1:2] <- c("Family", "SampleID")
    df <- dplyr::left_join(df, metadata[, c("SampleID", "Treatment", "Moisture"), drop = FALSE], by = "SampleID")
    df$Source <- compartment; families[[compartment]] <- df
    groups <- metadata %>% dplyr::count(Treatment, Moisture, name = "N")
    groups$Dataset <- slug; groups$Source <- compartment; all_counts[[paste(slug, compartment)]] <- groups
    libraries <- metadata[, c("SampleID", "Treatment", "Moisture", "Source")]
    libraries$Dataset <- slug; libraries$Input_ASV_count <- unname(before[metadata$SampleID])
    libraries$Retained_ASV_count <- unname(colSums(x)); libraries$Retained_ASVs_in_compartment <- nrow(x)
    all_libraries[[paste(slug, compartment)]] <- libraries
    all_alpha[[paste(slug, compartment)]] <- a; all_permanova[[paste(slug, compartment)]] <- p
    # Independent-group dispersion diagnostic; block restrictions remain unresolved.
    grouping <- interaction(metadata$Treatment, metadata$Moisture, drop = TRUE)
    bd <- vegan::betadisper(distance, grouping, type = "median", bias.adjust = TRUE)
    set.seed(publication_seed + 100 + j)
    bd_test <- vegan::permutest(bd, permutations = 999)$tab
    deviations <- abs(alpha$Shannon - ave(alpha$Shannon, grouping, FUN = median))
    lev <- summary(aov(deviations ~ grouping))[[1]]
    all_diagnostics[[paste(slug, compartment)]] <- data.frame(Dataset = slug, Source = compartment,
      Levene_F = lev[1, "F value"], Levene_P = lev[1, "Pr(>F)"],
      PERMDISP_F = bd_test[1, "F"], PERMDISP_P = bd_test[1, "Pr(>F)"],
      Shapiro_W = unname(shapiro.test(residuals(fit))$statistic), Shapiro_P = shapiro.test(residuals(fit))$p.value)
  }
  all_objects[[slug]] <- list(results = results, families = families)
}
sample_counts <- dplyr::bind_rows(all_counts)
alpha_tests <- dplyr::bind_rows(all_alpha)
permanova_tests <- dplyr::bind_rows(all_permanova)
library_counts <- dplyr::bind_rows(all_libraries)
diagnostics <- dplyr::bind_rows(all_diagnostics)
exclusions <- dplyr::bind_rows(all_excluded)
# Multiplicity is a sensitivity analysis within each model term and test type.
alpha_tests$BH_q <- NA_real_; permanova_tests$BH_q <- NA_real_
for (term in c("Treatment", "Moisture", "Treatment:Moisture")) {
 a <- which(alpha_tests$Term == term); alpha_tests$BH_q[a] <- p.adjust(alpha_tests$`Pr(>F)`[a], "BH")
 p <- which(permanova_tests$Term == term); permanova_tests$BH_q[p] <- p.adjust(permanova_tests$`Pr(>F)`[p], "BH")
}

  list(objects = all_objects, tables = list(sample_counts = sample_counts,
    alpha_tests = alpha_tests, permanova_tests = permanova_tests,
    library_counts = library_counts, diagnostics = diagnostics, exclusions = exclusions))
}
