nextflow run main.nf \
--samples 's3://czbiohub-mosquito/contigs/**.m9' \
--outdir "s3://czbiohub-mosquito/contig_quality_test" \
-profile aws
