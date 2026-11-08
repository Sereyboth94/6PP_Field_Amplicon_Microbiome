# Methods and output interpretation

This repository packages the four supplied figure workflows without changing
their statistical model expressions or publication plotting functions.

Filtering is applied separately within Root and Rhizosphere. Samples below the
dataset's minimum depth are removed, then ASVs require at least 10 total reads
and presence in at least 5% of retained samples. A zero-count sample remaining
after ASV filtering stops the analysis. The original post-filter sample counts
and group representation checks are retained.

Shannon diversity is estimated from the filtered phyloseq object and analysed
using two-way ANOVA with treatment, moisture, and their interaction. Observed
richness is also written in the supporting data but is not the main alpha panel.

Bray–Curtis distances are calculated from square-root relative abundances.
PCoA supplies the beta-diversity panels. PERMANOVA uses 999 unrestricted
permutations with sequential term tests in the order Treatment, Moisture,
Treatment:Moisture. A fixed seed is set separately for each compartment.
Sequential effects can depend on term order; these scripts do not substitute
marginal tests or a blocked design.

Main tables combine the three model terms for each of the two compartments.
Full model outputs and actual sample counts are exported alongside the plotted
observations. The pipeline recalculates tables from the processed inputs;
manually transcribed manuscript statistics are not used as analysis inputs.

Family relative abundances are calculated by summing ASVs within each family
and sample before expressing those sums as percentages. Means include zeros.
Main-bar family selection ranks the mean abundance across compartments with
equal compartment weighting. Supplemental pies select families within each
compartment. Others pools the unselected families. Displayed relative-abundance
patterns do not establish absolute population changes or statistical
differential abundance.

PDFs use Cairo. Text remains selectable, although PDF editors may need the
corresponding local font to edit it faithfully. Optional `pdftools` audit files
record actual embedded fonts. PNG/TIFF files use 600 dpi; TIFF compression is
lossless LZW. The supplied plotting functions and dimensions are retained.

`Supporting_Data/Publication_Plots_and_Table.rds` stores the composed R plots
and main table. Other RDS files store filtered phyloseq objects. CSV outputs
contain alpha observations, ordination coordinates and eigenvalues, model
tables, per-sample family percentages, stacked-bar/pie values, sample counts,
rankings, and colour mappings.

The legacy branch is optional and disabled by default. Its DESeq2 analyses and
older plots are not required to generate the requested main figures or pies.
