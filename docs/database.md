# Building database

As minimap2 is used for mapping ONT reads, the database can be a multiple sequences fasta file, but if it's large, it's best to build an index for it.
## Requirements
To build and format the database for this pipeline, we will need to have `kraken2` and `seqkit`

```
micromamba install -c conda-forge -c bioconda kraken2 seqkit
```

## NCBI RefSeq

1. Download genome sequences with kraken2

```bash
mkdir -p refseq
export DB=$PWD/refseq
kraken2-build --download-library bacteria --no-masking --db $DB
kraken2-build --download-library viral --no-masking --db $DB
kraken2-build --download-library fungi --no-masking --db $DB
kraken2-build --download-library protozoa --no-masking --db $DB
kraken2-build --download-library archaea --no-masking --db $DB

#Once they are downloaded successfully, combine them all
cat {bacteria,viral,fungi,protozoa}/library.fna | seqkit replace -p 'kraken::taxid\|' -r > refseq.fna
```

Sequence headers in the `refseq.fasta` should look like then

```
455631|NZ_CM000441.1 Clostridioides difficile QCD-66c26 chromosome, whole genome shotgun sequence
455631|NZ_ABFD02000018.1 Clostridioides difficile QCD-66c26 contig00018_2, whole genome shotgun sequence
```
`455631` is taxid

2. As this is a large sequence file, it's best to build index for it

```bash
minimap2 -I8G -d refseq.mmi refseq.fna
```

The parameter I8G directs minimap2 to partition the sequences by every 8Gb, so the memory required for building the index should be at least 8*2.5 GB.
