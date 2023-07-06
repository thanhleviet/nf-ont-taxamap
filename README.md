# ONT-Taxonomic-Profiling

Taxonomic profiling ONT metagenomic reads using minimap2
## Usage
```
nextflow run thanhleviet/nf-ont-taxamap --help
```

```bash
Typical pipeline command:

  nextflow run thanhleviet/nf-ont-taxamap --input /path/to/a/folder/of/barcodes

Input/Output options
  --email            [string]  Email to send notification to [default: @quadram.ac.uk]
  --input            [string]  Must be an absolute path to a folder of either sub-folders or files starting having prefix name barcode [default: /bart]
  --outdir           [string]  Name of the output folder [default: results]

Filtering reads by sequencing time
  --skip_filter_time [boolean] Skip filtering by sequencing time ? [default: true]
  --ontime_from      [string]  Extract reads generated from... [default: 0h]
  --ontime_to        [string]  Extract reads generated until... [default: 1h]

Pipeline Options
  --skip_dehuman     [boolean] Skip human seq removal ?
  --human_ref        [string]  null
  --ref              [string]  Path to Taxa database [default: None]
  --score            [integer] Score to include the mapped reads in [default: 30]
  --min_read_length  [integer] Minimum length of reads to be included [default: 200]

AMR
  --skip_amr         [boolean] Skip screening AMR
  --amr_db           [string]  Path to AMR database [default: None]
```
## Input
See [Input](docs/input.md)
## Database
See [database](docs/database.md)
## Output
See [Output](docs/output.md)

## Authors

Thanh.Le-Viet <thanh.le-viet@quadram.ac.uk>
