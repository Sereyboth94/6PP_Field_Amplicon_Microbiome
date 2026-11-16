# Final tall microbiome figure scripts

This update reproduces the four main figures and four six-panel supplementary figures approved on 12 November 2026. It retains the four dataset script names and provides one master runner. All four entry scripts share one scientific calculation module and one set of drawing functions, so their defaults remain consistent.

## Install into your existing folder

Download the ZIP into Downloads, then extract and install it with the PowerShell commands below. Keep its Microbiome_Figures_Scripts_and_Tables folder together. The installer backs up only the script/documentation files it replaces and verifies the backups and every copied file. It does not change data tables or existing figures, and does not require administrator privileges.

```powershell
$updateFolder = Join-Path $env:USERPROFILE 'Downloads\Tall_Figure_Scripts_Update'
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\Microbiome_Figures_Scripts_and_Tables_Tall_Final.zip" -DestinationPath $updateFolder -Force
powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $updateFolder 'Microbiome_Figures_Scripts_and_Tables\Install_Updated_Scripts.ps1')
```

If your browser renames the extracted folder, adjust that source path. The default target is your existing OneDrive project. You may instead copy the new folder contents into the existing Microbiome_Figures_Scripts_and_Tables folder manually, retaining a backup first. The shared R directory and Figure_Settings.R are required; do not copy only the four short entry scripts.

## Run all figures in R or RStudio

```r
source("C:/Users/soths/OneDrive - Lincoln University/Writing/Amplicon_Field_Experiment/github/Microbiome_Figures_Scripts_and_Tables/Run_All_Four_Figures.R", encoding = "UTF-8")
```

No setwd() is required. Use source() on the complete file. R code must be entered in R/RStudio, not a PowerShell prompt. Install missing dependencies once in R with:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "vegan", "ape"))
```

The tested versions are recorded in VALIDATION.txt. New package versions and different system fonts can slightly change rendering. The approved Linux reference PDFs embed Liberation Sans, the Arial-compatible fallback used there; Windows usually resolves the default Arial request to Arial. Scientific values and configured physical layout do not depend on changing the font. Cairo is the default export backend. Optional ragg and svglite are used only when deliberately enabled in Figure_Settings.R. pdftools is optional and adds a PDF font audit when installed.

## Data paths

By default, the input root is the parent github folder. These existing data folders must be beside the scripts folder:

| Folder | Required files |
| --- | --- |
| Pasture/Bacteria | ASV_filtered.txt, Taxonomy_filtered.txt, Metadat_Bacteria_Field.txt |
| Pasture/Fungi | ASV_filtered.txt, Taxonomy_filtered.txt, Metadat_Fungi_Pasture.txt |
| Wheat/Bacteria | ASV_filtered.txt, Taxonomy_filtered.txt, Metadat_Bacteria_Field.txt |
| Wheat/Fungi | ASV_filtered.txt, Taxonomy_filtered.txt, Metadat_Fungal_Wheat.txt |

Input tables are read-only. For another location, set options(amplicon.input_root = "C:/path/to/data-parent") before sourcing. All four datasets are analysed even when a single-dataset entry script is run, because BH adjustment is defined across all eight compartment datasets.

## Output locations and formats

Each run creates a new directory under Microbiome_Figures_Scripts_and_Tables/Tall_Figure_Outputs/Run_DATE_TIME_PROCESS. The completion message prints its full path. Existing figures are retained. To use another output parent, set options(amplicon.output_root = "C:/path/to/output-parent").

Figures contains the 24 PDF/PNG/TIFF files with the approved filenames and Figure_captions.txt. Supporting_Data contains full numerical tests, sample counts, family proportions, plotting settings and session information. Tables contains four editable CSV statistical tables. Cache contains the numerical analysis used for fast replotting. No RAW FASTQ processing is performed.

| Figure type | Width | Height | Raster resolution |
| --- | --- | --- | --- |
| Main figures 1–4 | 190 mm | 9 in / 228.6 mm | 600 dpi |
| Supplementary figures S1–S4 | 190 mm | 9.5 in / 241.3 mm | 600 dpi |

TIFF uses LZW compression. PDF uses vector drawing and Cairo font embedding. At the default device rounding the TIFFs have 4488 × 5400 pixels (main) and 4488 × 5700 pixels (supplementary), tagged at 600 × 600 dpi. PNG is also exported at 600 dpi. Optional SVG export adds one editable vector file per figure.

## Panel arrangement

Main figures: A alpha diversity, B beta diversity and C composition in rhizosphere; D alpha diversity, E beta diversity and F composition in roots. A full-width three-column family key is under each corresponding row.

Supplementary figures: A rhizosphere overall, B rhizosphere Control, C rhizosphere 6-PP, D root overall, E root Control and F root 6-PP. Every pie has its own single-column key below it. The top 15 families are retained within each compartment. Named groups descend in abundance within each pie; Others is last. Colours stay consistent within a dataset. Means include zeros and give every retained sample equal weight.

## Supplementary pies in one script

To generate only the four six-panel supplementary figures, source Run_Supplementary_Pies_Only.R. This writes 12 image files. Run_All_Four_Figures.R generates both main and supplementary figures. Individual Pasture and Wheat entry scripts generate their matching main and supplementary pair without renumbering them.

## Edit the figures reproducibly

Edit Figure_Settings.R for font sizes, height/width, legend space, colours, family number or percentage-label threshold. The supplied defaults reproduce the approved tall layout. Increasing font size or family number may require increasing legend space or height.

After a completed run, source Replot_Figures.R to reuse the latest numerical analysis and redraw the figures without refitting models. The input checksums must still match. To replot a specific earlier cache, set options(amplicon.cache_file = "C:/path/to/Run_.../Cache/study_analyses.rds"). Any entry script can use that explicit cache. To return to a fresh analysis, clear it with options(amplicon.cache_file = NULL).

For changes to panel structure, edit R/Figure_Layouts.R. For drawing details, edit R/figure_helpers.R. Data and statistical changes belong in the input tables or scientific calculation module, and require a fresh analysis. For cosmetic edits outside R, the vector PDF or optional SVG preserves drawing/text objects; TIFF and PNG are raster exports. The most reliable final figures come from editing R and regenerating all formats together.

## Scientific specification

The scientific calculation reproduces the deposited workflow: compartment-specific depth thresholds of 500 counts (100 for wheat bacteria), ASV total counts ≥10 and prevalence ≥5%, Shannon from filtered counts without rarefaction, sequential Type I two-way ANOVA, square-root relative abundance and Bray–Curtis PCoA/PERMANOVA with 999 unrestricted permutations. Seeds are 20261106 for rhizosphere and 20261107 for root PERMANOVA. Diagnostics use seeds 20261206/20261207. BH sensitivity uses separate eight-dataset families for each term and test type.

The reference input tables are from Zenodo v1.0.0, https://doi.org/10.5281/zenodo.24147075. Sample D07 is excluded only from pasture rhizosphere bacteria. These figure scripts do not resolve the field block/shelter allocation, carrier-control details or pooled-marker processing questions already recorded in the manuscript review. No block identities are inferred from well names. Family comparisons remain descriptive.

These scripts are a subsequent update and should be archived in a reviewed new code release. They are not claimed to be part of Zenodo v1.0.0. The installer and scripts do not commit, push or publish anything to GitHub.
