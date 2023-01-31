# /qib/platforms/Informatics/transfer/incoming/Lluis_Moragues/04Mar22_03Feb22_batch_adaptive/no_sample/20220304_1335_MN35050_FAP78209_528a17ec/fastq_HAC/barcode*.fastq.gz
# REF="/dev/shm/top_blood_infectious_species/blood_species.mmi"
# export NXF_EXECUTOR="local";
REF="/share/minimap_refseq/refseq.mmi"
REF_SAM="/share/minimap_refseq/refseq_ids.mmi"
HUMAN_REF="/qib/platforms/Informatics/transfer/outgoing/databases/CHM13/chm13.mmi"
# HUMAN_REF="/share/hg38_bowtie2/hg38.mmi"
# OUTDIR="04Mar22_03Feb22_batch_adaptive_top_species_taxonomy_score_30"
nextflow run minimap2.nf \
--reads "/qib/platforms/Informatics/transfer/incoming/Lluis_Moragues/1mL Standard results compilation/S. aureus/Replicate 2 (100 CFU missing)/basecalling/pass/barcode*.fastq.gz" \
--ref $REF \
--ref_sam $REF_SAM \
--human_ref $HUMAN_REF \
--score 85 \
--amr_db "/share/amrfinderplus/AMR" \
--min_read_length 250 \
--sam2lca_db "/share1/conda/lca/.sam2lca" \
--outdir "/share/Lluis/1mL_Standard_results_compilation/S_aureus/Rep2" \
-resume