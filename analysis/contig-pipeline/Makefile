.PHONY: scratch/CMS_040_RNA_A_S21.star.Unmapped.contig_miner
scratch/CMS_040_RNA_A_S21.star.Unmapped.contig_miner: scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_1.fa scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_2.fa scratch/mosquito_genomes2/mosquito_genomes2_truncatedNames.fa
	mkdir -p $@
	cd $@ && ../../../../scripts/contig_miner_v3.py -i ../CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_1.fa ../CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_2.fa -r ../mosquito_genomes2/mosquito_genomes2_truncatedNames.fa -t 20 &> log.cm.out


scratch/mosquito_genomes2/mosquito_genomes2_truncatedNames.fa: scratch/mosquito_genomes2/mosquito_genomes2.fa
	python truncate_fasta_names.py < $< > $@

scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_1.fa: 
scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_2.fa: 
scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_%.fa: scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_%.fq
#	sed -n '1~4s/^@/>/p;2~4p' $< | sed -r 's/^(>.*[^\s])\s+00\s*$$/\1/' > $@
	sed -n '1~4s/^@/>/p;2~4p' $< > $@

scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_2.fq: | scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_1.fq

scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_1.fq:
	PriceSeqFilter -fp scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate1.fq scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate2.fq -op scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_1.fq scratch/CMS_040_RNA_A_S21.star.Unmapped.priceseqfilter_2.fq -a 20 -rnf 90 -log c -rqf 85 0.98

.PHONY: scratch/CMS_040_RNA_A_S21.star.Unmapped.spades
scratch/CMS_040_RNA_A_S21.star.Unmapped.spades: scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate1.fq scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate2.fq
	mkdir -p $@
	spades.py --meta -t 20 -1 scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate1.fq -2 scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate2.fq -o $@

scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate%.fq:
	cp scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate$* $@

scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate2.wc_l:
scratch/CMS_040_RNA_A_S21.star.Unmapped.out.mate1.wc_l:
scratch/CMS_040_RNA_A_S21.star.Aligned.out.sam.wc_l:
scratch/CMS_040_RNA_A_S21.minimap2_sr.sam.wc_l:
scratch/%.wc_l: scratch/%
	grep -v '^@' $< | grep -v '^\s*$$' | wc -l > $@

scratch/CMS_040_RNA_A_S21.minimap2_sr.bam.picard_metrics:
#scratch/CMS_040_RNA_A_S21.bowtie2.sam.picard_metrics:
scratch/%.picard_metrics: scratch/%
	picard CollectAlignmentSummaryMetrics R=scratch/mosquito_genomes2/mosquito_genomes2.fa.gz I=$< O=$@

scratch/CMS_040_RNA_A_S21.bowtie2.sam.samtools_flagstat:
scratch/CMS_040_RNA_A_S21.minimap2_sr.bam.samtools_flagstat:
# Fails with: [W::sam_parse1] urecognized reference name; treated as unmapped
#scratch/CMS_040_RNA_A_S21.star.Aligned.out.sam.samtools_flagstat:
scratch/%.samtools_flagstat: scratch/%
	samtools flagstat $< > $@

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

scratch/mosquito_genomes2/star_index_withContigHeaders: scratch/mosquito_genomes2/mosquito_genomes2_withContigHeaders.fa
	mkdir -p $@
	STAR --runThreadN 45 --runMode genomeGenerate --genomeDir $@ --genomeFastaFiles $<

scratch/mosquito_genomes2/star_index: scratch/mosquito_genomes2/mosquito_genomes2.fa
	mkdir -p $@
	STAR --runThreadN 45 --runMode genomeGenerate --genomeDir $@ --genomeFastaFiles $<

scratch/mosquito_genomes2/mosquito_genomes2_withContigHeaders.fa: | scratch/mosquito_genomes2
	gunzip -c $@.gz > $@

scratch/mosquito_genomes2/mosquito_genomes2.fa: | scratch/mosquito_genomes2
	gunzip -c $@.gz > $@

scratch/CMS_040_RNA_A_S21.bowtie2.sam:
scratch/%.bowtie2.sam:
	bowtie2 -p 20 -x scratch/mosquito_genomes2/mosquito_genomes2_withContigHeaders --very-sensitive-local -S $@ -1 scratch/$*_R1_001.fastq.gz -2 scratch/$*_R2_001.fastq.gz

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
