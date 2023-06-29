# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2023-06-29

### Added
- Added average read length calculation
- Added filtering reads by sequencing time

## [1.0.0] - 2023-06-15

### Added
- Accept input as a folder path containing subfolders or files with prefix name barcode, if subfolders found, multiple fastq files within each barcode folder will be concatenated to a compressed barcodeXX.fastq.gz
- Filter reads by length
- Remove human host using minimap2
- Map the dehosted reads against a target ref either fna or mmi
- Screen AMR genes with KMA and AMRFinderPlus DB
- Generate ready view report for visualisation with Pavian and a summary in excel.
