# Validation of the repository package

Package checks completed on 8 November 2026 using R 4.3.3 and PowerShell 7.4.6
in a temporary Linux test environment. The Windows helper scripts use
PowerShell 5.1-compatible syntax and standard cmdlets; they have not been run
on the user's Windows computer by the package authoring environment.

## Source preservation

The supplied scripts were matched to the current source bundle and the file
sizes in the user's directory listings. The analysis and plotting function
blocks are unchanged. Filtering values, model expressions, permutation counts,
seed, family selections, figure dimensions, label sizes, and export resolution
retain the supplied defaults. Changes are limited to portable input/output
paths, shared validation, configurable font/legacy options, and the runner.

## Checks completed

- Parsed every R and PowerShell file.
- Validated all four supplied metadata tables against the expected 32 samples
  and eight groups of four per plant/marker dataset, using temporary count and
  taxonomy tables for the software checks.
- Tested the runner and all four standalone path blocks from an unrelated
  working directory and a repository path containing spaces.
- Checked that missing inputs, missing samples, duplicate identifiers, invalid
  counts, missing Family taxonomy, and unknown experimental labels stop validation.
- Ran the Windows setup/import logic under PowerShell with temporary source
  directories. All twelve files were copied and SHA256-verified without changing
  source bytes. Repeated setup preserved matching files; changed destinations
  were refused by setup and backed up by an explicit import.
- Tested GitHub publication using a local stand-in for the GitHub CLI, with no
  remote account writes. Wrong accounts, existing repositories, and changed
  input hashes were refused. The successful simulation selected private
  visibility and the specified Sereyboth94 owner.
- Created a real local test Git commit and verified that all twelve processed
  files retained their exact byte hashes. Test raw-read and credential files
  were excluded.
- Parsed the citation and workflow YAML and checked read-only CI permissions.

## What still runs on the analysis computer

The eight original ASV and taxonomy tables were not supplied to the authoring
environment. They are not replaced with demonstration data. The setup script
copies them from the user's existing folders before publication. Full analysis
and 600-dpi exports must be run against those real tables in the user's R
environment. `sessionInfo.txt` and run logs record the actual package versions
and retained samples. No new manuscript results or figures are claimed from
the temporary software tests.

The prepared package is not evidence that a hosted GitHub repository has been
created. That operation occurs only when the publication helper succeeds using
the authenticated account.
