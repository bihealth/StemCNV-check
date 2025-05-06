# Changelog

All notable changes to this project will be documented in this file.

Note: change log is only added with version 0.5.1

## [0.5.1] - 2025-05-06

### Added

- Added automated detection and usage of available cores and - on WSL only - memory, added `--memory` option

### Fixed

- Fixed the annotated number of call in the precision estimation data
- Fixed crashes of make_summary_table when parsing PennCNV logs generated form wrongly annotated sex
- Fixed selection, log output and Maximum number restriction of sample_id collection for SNP clustering
- Fixed crashes caused by usage of snakemake pipe in WSl to/from /mnt/ drives (by not using pipe in those cases)
- Fixed counting of CNVs/LOhs for summary table to properly use the "call_count_excl_labels" config setting
- Fixed splitting of combined CNVs into individual callers for report and summary stats

### Changed

- Updated mehari to version 0.35.1 and its transcript db to 0.10.3
- Updated descriptions and annotations in config file and report, renamed CNV and SNV labels to be more consistent
- Updated SNP clustering on WSL only to not parallelize vcf io, since memory overhead was considerable
- Improved report summary table by splitting it into Data QC (incl. the reference) and sample QC
- Updated usage of global array definition to ignore arrays not used in the current sampletable

## [0.5.0] and before

TBD