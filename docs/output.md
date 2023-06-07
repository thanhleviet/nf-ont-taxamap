# Output
## Folders
- `concatenate`:Concatenated fastq file by barcode if the input is a folder of barcode sub folders
- `dehuman2`: reads after having human reads removed
- `filter`: filtered reads by length
- `minimap2`:mapped paf format
- `clean`: the mapped paf file after cleaning some bits
- `pipeline_info`: run metrics of the pipeline

## Files
###  **\*_reads.csv**: List of mapped reads
```
read,taxid,read_length,closet_reference,read_coverage,alignment_coverage,alignment_identity,score
e8eee79e-2c91-4305-913b-8d79679693c8,537011,3741,NZ_CP085932.1,94.97,91.05,94.48,92.73
```
Columns:
  - read: read id
  - taxid
  - read_length
  - closet_reference: the closest one that the read mapped to.
  - read_coverage
  - alignment_coverage
  - alignment_identity
  - score: is a harmonic function `2 x alignment_identity*alignment_coverage/(alignment_identity+alignment_coverage)`

###  **\*_report.csv**: Summary of mapped reads as per taxid

```
taxid,count,kingdom,phylum,class,order,family,genus,species,strain
562,8502,k__Bacteria,p__Proteobacteria,c__Gammaproteobacteria,o__Enterobacterales,f__Enterobacteriaceae,g__Escherichia,s__Escherichia coli,
287,3808,k__Bacteria,p__Proteobacteria,c__Gammaproteobacteria,o__Pseudomonadales,f__Pseudomonadaceae,g__Pseudomonas,s__Pseudomonas aeruginosa,
```
Columns:

  - taxid
  - count
  - kingdom
  - phylum
  - class
  - order
  - family
  - genus
  - species
  - strain

### **\*_metaphlan_report.csv**: Summary of mapped reads as per taxid by MPA format.

This file is generated for visualisation with the [Pavian Shinny app](https://fbreitwieser.shinyapps.io/pavian/)

```
#mpa_v3
#clade_name	NCBI_tax_id	relative_abundance	additional_species
k__Bacteria	2	33609
k__Eukaryota	2759	2
k__Viruses	10239	9
k__Bacteria|p__Actinobacteria	201174	660
k__Bacteria|p__Bacteroidetes	976	14405
```
