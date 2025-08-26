# Changelog

All notable changes to this project will be documented in this file.

Note: change log is only added with version 0.5.1

## [1.0.0] - TBD

### Added

- Documentation
- CI tests for building docs
- Hotspot tables hosted via github-pages

### Changed

- SNV labels: ROI-overlap is now reportable (from critical)
- Updated format and references for stem-cell hotspot lists (CNVs and SNVs)

## [0.5.4] - 2025-08-11

### Changed

- Updated stem-cell hotspot lists (CNVs and SNVs)


## [0.5.3] - 2025-07-29

### Changed

- Updated precision estimation reference data
- Update CNV labels


## [0.5.2] - 2025-07-08

### Fixed

- Adapted to and pinned vembrane version to >=2.1,<3 (fixes penncnv input errors)
- Fixed sampletable parsing if xlsx file with comment lines on top was used

### Changed

- Updated defaults dendrogram plotting in report

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