# 6-PP field microbiome analysis in pasture and wheat

R workflows for bacterial and fungal communities in root-associated and
rhizosphere samples from a field experiment comparing control and
6-pentyl-alpha-pyrone (6-PP) seed treatments under irrigated and drought
conditions. Repository owner: [Sereyboth Soth](https://github.com/Sereyboth94).

The four current analysis scripts generate four main figures, four main
statistical tables, and four supplementary pie-chart figures. Each analysis
uses its own processed ASV, taxonomy, and sample-metadata tables.

## Quick start

Open `6PP_Field_Amplicon_Microbiome.Rproj` in RStudio. From the repository root:

```r
# Install missing dependencies once.
source("tools/Install_Dependencies.R")

# Check all twelve processed input tables without generating figures.
source("tools/Validate_Inputs.R")

# Generate all four datasets' manuscript outputs.
source("Microbiome_Figures_Scripts_and_Tables/Run_All_Four_Figures.R")
```

Alternatively, from a terminal with Rscript on PATH:

```bash
Rscript Microbiome_Figures_Scripts_and_Tables/Run_All_Four_Figures.R --check-only
Rscript Microbiome_Figures_Scripts_and_Tables/Run_All_Four_Figures.R
```

The saved runner can be sourced or invoked by its full path from another
working directory. It locates this repository automatically. No username,
OneDrive location, drive letter, or `setwd()` is embedded in the analysis paths.

## Initial setup from the original Windows folders

The prepared download includes the current code and four supplied metadata
tables. The eight ASV and taxonomy tables are on the original analysis computer
and must be copied before the first full run or GitHub publication.

If the new repository is inside the original `github` directory, run this
saved file in Windows PowerShell:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File ".\tools\Import_Original_Data.ps1"
```

This copies exactly the twelve files listed below from the original sibling
`Pasture` and `Wheat` folders. It verifies SHA256 checksums, writes
`data_checksums.csv`, and leaves the originals unchanged. Existing repository
inputs that differ are backed up locally before replacement. Nothing is
uploaded by the import script. Use `-SourceRoot` to select another original
`github` directory, or `-CheckOnly` to check readability without copying.

The separate `Setup_Amplicon_GitHub_Repo.ps1` download can extract this package,
perform that import, and initialise a new local repository in one command.
Do not initialise Git in the original parent `github` directory: that directory
also contains older scripts and unrelated results.

## Repository layout

| Path | Contents |
| --- | --- |
| `Microbiome_Figures_Scripts_and_Tables/Run_All_Four_Figures.R` | Main runner |
| `Microbiome_Figures_Scripts_and_Tables/Pasture/` | Current pasture bacterial and fungal scripts |
| `Microbiome_Figures_Scripts_and_Tables/Wheat/` | Current wheat bacterial and fungal scripts |
| `Pasture/Bacteria/`, `Pasture/Fungi/` | Processed pasture inputs |
| `Wheat/Bacteria/`, `Wheat/Fungi/` | Processed wheat inputs |
| `R/project_utils.R` | Shared paths, dataset mapping, and input checks |
| `tools/` | Dependency installation, input import/validation, and GitHub publication |
| `tests/run_checks.R` | Base-R validation and portability tests |
| `.github/workflows/checks.yml` | Automated syntax, validation, and input checks |
| `docs/` | Data dictionary, methods, and validation notes |
| `results/` | Generated figures, tables, supporting data, and logs; excluded from Git |

## Required input files

Every dataset folder contains `ASV_filtered.txt`, `Taxonomy_filtered.txt`, and
its metadata file. All are tab-delimited, despite the `.txt` extension.

| Dataset folder | Metadata filename | Sample IDs |
| --- | --- | --- |
| `Pasture/Bacteria/` | `Metadat_Bacteria_Field.txt` | A05–H08, columns 05–08 |
| `Pasture/Fungi/` | `Metadat_Fungi_Pasture.txt` | A05–H08, columns 05–08 |
| `Wheat/Bacteria/` | `Metadat_Bacteria_Field.txt` | A01–H04, columns 01–04 |
| `Wheat/Fungi/` | `Metadat_Fungal_Wheat.txt` | A01–H04, columns 01–04 |

There are 32 biological samples per plant system and 64 across the experiment;
bacterial and fungal analyses refer to the corresponding biological specimens.
Input checks require all expected sample IDs and four samples in each
compartment × moisture × treatment group. Sample retention after filtering is
reported separately by each analysis. See [the data dictionary](docs/DATA_DICTIONARY.md).

## Figures and tables

| Number | Dataset | Script under `Microbiome_Figures_Scripts_and_Tables/` |
| --- | --- | --- |
| Figure 1, Table 1, Figure S1 | Pasture bacteria | `Pasture/Phyloseq_Bacteria.R` |
| Figure 2, Table 2, Figure S2 | Pasture fungi | `Pasture/Phyloseq_Fungi.R` |
| Figure 3, Table 3, Figure S3 | Wheat bacteria | `Wheat/Wheat_Phyloseq_Bacteria.R` |
| Figure 4, Table 4, Figure S4 | Wheat fungi | `Wheat/Wheat_Phyloseq_Fungi.R` |

Each main figure contains six panels:

| Panel | Compartment | Display |
| --- | --- | --- |
| A | Rhizosphere | Shannon alpha diversity |
| B | Rhizosphere | Bray–Curtis PCoA |
| C | Rhizosphere | Mean family-level relative abundance |
| D | Root | Shannon alpha diversity |
| E | Root | Bray–Curtis PCoA |
| F | Root | Mean family-level relative abundance |

Each supplementary figure has **A, Rhizosphere** and **B, Root** pie charts.
Family legends are directly below their own composition or pie panels.
Legend titles are **Bacterial families** or **Fungal families**. Panels C and F
use three-column family legends and 18-point family labels. Main figures are
21 × 14.5 inches; supplementary figures are 18 × 12 inches.

Main bars display the top 15 families ranked jointly across the two
compartments, with equal compartment weighting, followed by Others. Pie charts
select and rank the top 15 families independently within each compartment,
with Others last. Colours remain consistent within each dataset. Families are
summed per sample before means are calculated, including zero abundances.

## Output files

Outputs are saved under:

```text
results/<Plant>/<Marker>/Phyloseq_Field_<Marker>/Manuscript_Figures_Tables/
```

For example, pasture bacterial outputs are in
`results/Pasture/Bacteria/Phyloseq_Field_Bacteria/Manuscript_Figures_Tables/`.

| Filename pattern | Formats |
| --- | --- |
| `FigureN_DATASET_Main_A-F` | PDF, PNG, TIFF |
| `FigureSN_DATASET_Pies_A-B` | PDF, PNG, TIFF |
| `TableN_DATASET_Main` | CSV, HTML, PDF, PNG, TIFF |
| `Supporting_Data/` | Plotted observations, model results, PCoA coordinates/eigenvalues, family percentages, sample counts, and RDS objects |
| `sessionInfo.txt` | Actual R and package versions for that dataset run |

PDFs use `cairo_pdf()` for vector graphics and embedded, selectable text. They
are not raster screenshots or outlined text. PNG and TIFF exports are 600 dpi;
TIFF uses lossless LZW compression. Cairo support is required. Arial is the
default font; install it for consistent rendering or set, for example,
`options(amplicon.font = "sans")` before running. PDF appearance may involve
system font substitution; optional `pdftools` font-audit CSVs report the fonts
actually embedded. Raster export uses `ragg` when available.

Run-level logs under `results/run_*/` contain input checks, input MD5 hashes,
completed datasets, and session information. A failed dataset is recorded in
`ERROR.txt`. The runner stops at the first analysis failure.

## Analysis settings

The portability update retains the four scripts' statistical methods and
publication defaults. Compartment analyses apply a minimum total ASV count of
10 and prevalence of 5%. Minimum sample depth is 500, except wheat bacteria,
which use 100. Shannon diversity is modelled as
`Shannon ~ Treatment * Moisture`. Community composition uses square-root
relative abundances, Bray–Curtis dissimilarity, PCoA, and PERMANOVA with
999 unrestricted permutations and sequential term tests (`by = "terms"`).
The base seed is 20261105, with a compartment-specific offset.

These scripts preserve the existing two-factor models. They do not add a
block term, restricted permutations, or new differential-abundance tests to the
main figures. Relative-abundance differences are descriptive. See
[methods and output details](docs/METHODS_AND_OUTPUTS.md).

Legacy analyses are disabled by default. To additionally run the original
legacy outputs, install their dependencies and opt in:

```r
options(amplicon.no_autorun = TRUE)
source("tools/Install_Dependencies.R")
install_amplicon_dependencies(include_legacy = TRUE)
options(amplicon.no_autorun = FALSE, amplicon.run_legacy = TRUE)
source("Microbiome_Figures_Scripts_and_Tables/Run_All_Four_Figures.R")
```

Core packages are phyloseq, ggplot2, dplyr, tidyr, tibble, vegan, and patchwork.
Legacy outputs additionally require DESeq2, ape, RColorBrewer, and forcats.
The installer uses CRAN and BiocManager; it does not force a Bioconductor release
incompatible with the installed R. R 4.3 or newer is expected; the user's
reported Windows R version is 4.4.0. No unverified package-version lockfile is included.

## Publish to GitHub

After importing the twelve real input files, run in Windows PowerShell:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File ".\tools\Publish_GitHub.ps1"
```

The default target is `Sereyboth94/6PP_Field_Amplicon_Microbiome`, with **private**
visibility. The script requires Git and [GitHub CLI](https://cli.github.com/),
checks the authenticated account, commits only the repository contents and
twelve processed inputs, and creates/pushes the new repository. If needed,
GitHub CLI opens its normal browser sign-in flow. An existing remote repository
is left unchanged. Use `-Visibility Public` only when ready to publish publicly,
or `-RepositoryName` to choose another name.

Raw FASTQ files, generated outputs, R histories, and credentials are excluded.
The processed tables are tracked with exact input bytes preserved. Raw-read
SRA accession information can be added once the archive record is finalised;
no accession is claimed in this repository.

## Checks and citation

```bash
Rscript tests/run_checks.R
Rscript tools/Validate_Inputs.R
```

The first command tests syntax and input/path checks using temporary test
tables; those tables are never used for manuscript figures. The second checks
the real processed inputs. CI runs both after the input import. Full figure
generation is a separate analysis step requiring the listed R packages.

See `CITATION.cff` for software citation metadata. Add the final manuscript
citation and a versioned archive DOI when available. Neither is invented here.

References: [Bioconductor installation](https://bioconductor.org/install/),
[GitHub CLI repository creation](https://cli.github.com/manual/gh_repo_create),
and [R Cairo graphics devices](https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/cairo.html).
