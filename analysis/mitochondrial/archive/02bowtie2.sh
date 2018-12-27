bowtie2-build aaegypti_EU352212.fasta aedes
bowtie2-build cpipiens_HQ724614.fasta culex

# bowtie2-build ../../s3/references/mosquito_mito/Culiseta_mitochondrion_genome_curated.fasta culiseta
# bowtie2-build ../../s3/references/mosquito_mito/Anopheles_mitochondrion_genome_curated.fasta anopheles

#for f in `cat mito_cms1_Culex.txt`
mkdir cpipiens_HQ724614
for fastq in `ls /mnt/data/s3/sequences/*/CMS_0*R1*`
do
  fn=$(basename ${fastq})
  prefix="${fn%.*.*}"
  prefix=${prefix/_R1_001/}
  dir='cpipiens_HQ724614/'
  prefix=${dir}${prefix}
  if [ ! -f ${prefix}.sam ]; then
    bowtie2 --threads 8 --qupto 1000000 -x culex -1 $fastq -2 ${fastq/R1/R2} -S ${prefix}.sam
  fi
  if [ ! -f ${prefix}.bam ]; then
    samtools view -Sb ${prefix}.sam | samtools sort -o ${prefix}.bam
  fi
  if [ ! -f ${prefix}.vcf.gz ]; then
    if [ -f ${prefix}.bam ]; then
      bcftools mpileup -Ou -f cpipiens_HQ724614.fasta ${prefix}.bam | bcftools call -mv -Oz -o ${prefix}.vcf.gz
      tabix ${prefix}.vcf.gz
      cat cpipiens_HQ724614.fasta | bcftools consensus ${prefix}.vcf.gz > ${prefix}_consensus.fasta
    fi 
  fi
done

for f in `cat mito_cms2_Culex.txt`
do
  if [ ! -f ${f}.sam ]; then
    bowtie2 --threads 8  --qupto 1000000 -x culex -1 /mnt/data/s3/sequences/CMS002_fastq.gz/${f}_R1_001.fastq.gz -2 /mnt/data/s3/sequences/CMS002_fastq.gz/${f}_R2_001.fastq.gz -S ${f}.sam
  fi
  if [ ! -f ${f}.bam ]; then
    samtools view -Sb ${f}.sam | samtools view -F 4 | samtools sort -o ${f}.bam
  fi
  if [ ! -f ${f}.vcf.gz ]; then
    if [ -f ${f}.bam ]; then
      bcftools mpileup -Ou -f cpipiens_HQ724614.fasta ${f}.bam | bcftools call -mv -Oz -o ${f}.vcf.gz
      tabix ${f}.vcf.gz
      cat cpipiens_HQ724614.fasta | bcftools consensus ${f}.vcf.gz > ${f}_consensus.fasta
    fi
  fi
done


for f in `cat mito_cms2_Aedes.txt`
do
  if [ ! -f ${f}.sam ]; then
    bowtie2 --threads 8  --qupto 1000000 -x aedes -1 /mnt/data/s3/sequences/CMS002_fastq.gz/${f}_R1_001.fastq.gz -2 /mnt/data/s3/sequences/CMS002_fastq.gz/${f}_R2_001.fastq.gz -S ${f}.sam
  fi
  fn=$(basename $f)
  if [ ! -f ${f}.bam ]; then
    samtools view -Sb ${f}.sam | samtools view -F 4 | samtools sort -o ${f}.bam
  fi
  if [ ! -f ${f}.vcf.gz ]; then
    if [ -f ${f}.bam ]; then
      bcftools mpileup -Ou -f aaegypti_EU352212.fasta ${f}.bam | bcftools call -mv -Oz -o ${f}.vcf.gz
      tabix ${f}.vcf.gz
      cat aaegypti_EU352212.fasta | bcftools consensus ${f}.vcf.gz > ${f}_consensus.fasta
    fi
  fi
done





samtools faidx ../cpipiens_HQ724614_COI.fasta
for i in `ls *.vcf.gz`
do
  prefix="${i%.*.*}"
  samtools faidx ../cpipiens_HQ724614_COI.fasta HQ724614.1:1446-2982 | bcftools consensus $i -o ${prefix}_COI.fasta
done

for i in `ls *COI.fasta`
do
  prefix="${i%.*.*}"
  A=">"
  prefix=$A$prefix
  sed -i "1s/.*/$prefix/" $i
done