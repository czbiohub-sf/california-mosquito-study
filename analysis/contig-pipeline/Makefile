#scratch/CMS_040_RNA_A_S21.bowtie2.sam.picard_metrics:
#scratch/%.sam.picard_metrics:
#	picard CollectAlignmentSummaryMetrics R=scratch/mosquito_genomes2/mosquito_genomes2.fa.gz I=scratch/$*.sam O=$@

# add headers back into minimap2 output
scratch/CMS_040_RNA_A_S21.minimap2_sr.bam:
scratch/%.minimap2_sr.bam: scratch/%.sam
	samtools view -b -t scratch/mosquito_genomes2/mosquito_genomes2_withContigHeaders.fa.gz.fai -o $@ $<

scratch/CMS_040_RNA_A_S21.star.Aligned.out.sam: scratch/CMS_040_RNA_A_S21_R1_001.fastq.gz scratch/CMS_040_RNA_A_S21_R2_001.fastq.gz
	STAR --runThreadN 45 \
	 --genomeDir scratch/mosquito_genomes2/star_index \
	 --readFilesIn $^ \
	 --outFileNamePrefix scratch/CMS_040_RNA_A_S21.star. \
	 --readFilesCommand zcat \
	 --outFilterMultimapNmax 99999 \
	 --outFilterScoreMinOverLread .5 \
	 --outFilterMatchNminOverLread .5 \
	 --outReadsUnmapped Fastx \
	 --outFilterMismatchNmax 999

scratch/mosquito_genomes2/star_index: scratch/mosquito_genomes2/mosquito_genomes2.fa
	mkdir -p $@
	STAR --runThreadN 45 --runMode genomeGenerate --genomeDir $@ --genomeFastaFiles $<

scratch/mosquito_genomes2/mosquito_genomes2.fa: | scratch/mosquito_genomes2
	gunzip -c $@.gz > $@

scratch/CMS_040_RNA_A_S21.bowtie2.sam:
scratch/%.bowtie2.sam:
	bowtie2 -p 6 -x scratch/mosquito_genomes2/mosquito_genomes2 --very-sensitive-local -S $@ -1 scratch/$*_R1_001.fastq.gz -2 scratch/$*_R2_001.fastq.gz

scratch/CMS_040_RNA_A_S21.minimap2_sr.sam:
scratch/%.minimap2_sr.sam:
	minimap2 -ax sr scratch/mosquito_genomes2/mosquito_genomes2_withContigHeaders.fa.gz scratch/$*_R1_001.fastq.gz scratch/$*_R2_001.fastq.gz > $@

scratch/mosquito_genomes2/mosquito_genomes2_withContigHeaders.fa.gz.fai:
	samtools faidx scratch/mosquito_genomes2/mosquito_genomes2_withContigHeaders.fa.gz

scratch/mosquito_genomes2:
	aws s3 sync s3://czbiohub-mosquito/references/mosquito_genomes2 $@

scratch/CMS_040_RNA_A_S21_R1_001_fastqc.html:
scratch/CMS_040_RNA_A_S21_R2_001_fastqc.html:
%_fastqc.html: %.fastq.gz
	fastqc $<

scratch/CMS_040_RNA_A_S21_R2_001.fastq.gz:
scratch/CMS_040_RNA_A_S21_R1_001.fastq.gz:
scratch/CMS_040_RNA_A_S21_R%_001.fastq.gz:
	aws s3 cp s3://czbiohub-mosquito/sequences/CMS001_fastq.gz/$(@F) $@
