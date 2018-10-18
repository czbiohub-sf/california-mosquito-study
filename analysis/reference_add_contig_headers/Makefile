.PHONY: upload

upload: scratch/mosquito_genomes2_withContigHeaders.fa.gz
	aws s3 cp $< s3://czbiohub-mosquito/references/mosquito_genomes2/mosquito_genomes2_withContigHeaders.fa.gz

scratch/mosquito_genomes2_withContigHeaders.fa.gz: scratch/ncbi_download.log
	zcat scratch/ncbi_download/* | bgzip -c > $@

scratch/ncbi_download.log:
	python download_ncbi_refs.py | tee $@

scratch/mosquito_genomes2_sequences.txt: scratch/mosquito_genomes2.fa.gz
	zcat $< | grep '^>' > $@

scratch/mosquito_genomes2.fa.gz:
	aws s3 cp s3://czbiohub-mosquito/references/mosquito_genomes2/mosquito_genomes2.fa.gz $@
