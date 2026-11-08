# Processed input data

The repository consumes processed tables. It does not reproduce demultiplexing,
read filtering, ASV inference, or taxonomic classification from raw FASTQ files.

## ASV_filtered.txt

Tab-delimited table with an ASV identifier in the first column and sample IDs
in the remaining columns. Each cell is a non-negative integer read count.
Sample names must match the corresponding metadata exactly. Extra sample
columns can remain in the table; the analysis selects the metadata samples.
Missing metadata samples, duplicated identifiers, nonnumeric counts, negative
counts, and non-finite values fail the input check.

## Taxonomy_filtered.txt

Tab-delimited table with ASV identifiers in the first column and taxonomy ranks
in the remaining columns. A `Family` column is required. The analysis uses
ASVs shared by abundance and taxonomy tables. The preflight reports the number
of shared ASVs. Blank/unassigned families are displayed as Unclassified family
by the plotting code; family-level colour keys are saved with the outputs.

## Sample metadata

| Column | Accepted values and meaning |
| --- | --- |
| First column (`SampleID` or `Sample`) | Unique well identifier, e.g. A05; standardised to SampleID during reading |
| Treatment | Control, 6PP, 6-PP, or FVOC; treated labels are standardised to 6PP and displayed as 6-PP |
| Moisture | Irrigated or Drought |
| Source | Root or Rhizosphere |
| Plant | Pasture or Wheat, matching the input folder |

The expected design contains four biological replicates in each
plant × compartment × treatment × moisture combination. Letters A–D denote
Control and E–H denote treated samples in the supplied metadata.

| Plant | Column | Compartment | Moisture |
| --- | --- | --- | --- |
| Wheat | 01 | Root | Drought |
| Wheat | 02 | Root | Irrigated |
| Wheat | 03 | Rhizosphere | Drought |
| Wheat | 04 | Rhizosphere | Irrigated |
| Pasture | 05 | Root | Drought |
| Pasture | 06 | Root | Irrigated |
| Pasture | 07 | Rhizosphere | Drought |
| Pasture | 08 | Rhizosphere | Irrigated |

Sample letters are not assigned block identities by this repository. The
`Root` category refers to the root-associated sample compartment described in
the study; it should not be relabelled as an endosphere-only dataset without
supporting sampling information.

## Data import and provenance

`tools/Import_Original_Data.ps1` copies exactly three files per dataset from the
original directories. `data_checksums.csv` records relative filenames, file
sizes, and SHA256 hashes. The four supplied metadata files are included in the
prepared package; the original computer supplies the eight ASV/taxonomy files
and refreshes all four metadata files during import.

No artificial ASV or taxonomy table is placed in the analysis input folders.
Automated tests generate small temporary tables outside the repository and
delete them after testing. They test software behaviour, not study results.
