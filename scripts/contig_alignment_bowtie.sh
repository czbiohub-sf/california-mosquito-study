# need to create /mnt/data/foo.py file and /mnt/data/names.sh file, files are below

mkdir -p /mnt/data/names /mnt/data/index; cd /mnt/data/names; mnt/data/names.sh; cd /mnt/data
for sample in `ls /mnt/data/names`
  do
    mkdir /mnt/data/reads /mnt/data/contigs /mnt/data/outputs
    #download material
    aws s3 cp s3://czbiohub-mosquito/sequences/CMS_cleaned_reads/${sample}_R1_cdh_lzw_trim30_PF.fastq.gz /mnt/data/reads
    aws s3 cp s3://czbiohub-mosquito/sequences/CMS_cleaned_reads/${sample}_R2_cdh_lzw_trim30_PF.fastq.gz /mnt/data/reads
    aws s3 cp s3://czbiohub-mosquito/contigs/${sample}/contigs.fasta /mnt/data/contigs
    cd /mnt/data/index
    #run bowtie at default frag length = 500bp
    bowtie2-build /mnt/data/contigs/contigs.fasta ${sample}
    echo "${sample} ALIGNMENT" >> /mnt/data/log.out
    bowtie2 -p4 -x ${sample} -q -1 /mnt/data/reads/${sample}_R1_cdh_lzw_trim30_PF.fastq.gz \
    -2 /mnt/data/reads/${sample}_R2_cdh_lzw_trim30_PF.fastq.gz  \
    --very-sensitive-local -S /mnt/data/outputs/${sample}.sam --no-unal --al-conc /mnt/data/outputs/${sample}_${base}.fq \
    2>> /mnt/data/log.out
    #sam to bam
    samtools view -S -b /mnt/data/outputs/${sample}.sam > /mnt/data/outputs/${sample}.bam
    #upload bam and sam to aws
    aws s3 cp /mnt/data/outputs/${sample}.sam s3://czbiohub-mosquito/alignment/sam/${sample}.sam
    aws s3 cp /mnt/data/outputs/${sample}.bam s3://czbiohub-mosquito/alignment/bam/${sample}.bam
    echo "${sample} done" >> /mnt/data/log2.out
    rm -r /mnt/data/reads /mnt/data/contigs /mnt/data/outputs
  done

#counting reads. need to create foo.py
#information on count code: https://www.biostars.org/p/145053/ & https://www.biostars.org/p/95929/
for sample in `ls /mnt/data/names_2`
  do
    mkdir /mnt/data/bam /mnt/data/outputs
    aws s3 cp s3://czbiohub-mosquito/alignment/bam/${sample}.bam /mnt/data/bam/
    samtools view -f 0x2 /mnt/data/bam/${sample}.bam | grep -v "XS:i:" | ./foo.py | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/outputs/bowtie_csp_counts.txt
    samtools view -F 260 /mnt/data/bam/${sample}.bam | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}'  > /mnt/data/outputs/bowtie_all_counts.txt
    samtools view -f 0x2 /mnt/data/bam/${sample}.bam | grep -v "XS:i:" | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/outputs/bowtie_cs_counts.txt
    aws s3 cp /mnt/data/outputs/bowtie_csp_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_csp_counts.txt
    aws s3 cp /mnt/data/outputs/bowtie_all_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_all_counts.txt
    aws s3 cp /mnt/data/outputs/bowtie_cs_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_cs_counts.txt
    echo "${sample} done" >> log3.out
    rm -r /mnt/data/bam  /mnt/data/outputs
  done


#bowtie and counting for fragment lengths of 700 and 1,000
mkdir -p /mnt/data/index2
for sample in `ls /mnt/data/names`
  do
    mkdir /mnt/data/reads /mnt/data/contigs /mnt/data/bam /mnt/data/sam
    aws s3 cp s3://czbiohub-mosquito/sequences/CMS_cleaned_reads/${sample}_R1_cdh_lzw_trim30_PF.fastq.gz /mnt/data/reads
    aws s3 cp s3://czbiohub-mosquito/sequences/CMS_cleaned_reads/${sample}_R2_cdh_lzw_trim30_PF.fastq.gz /mnt/data/reads
    aws s3 cp s3://czbiohub-mosquito/contigs/${sample}/contigs.fasta /mnt/data/contigs
    cd /mnt/data/index2

    bowtie2-build /mnt/data/contigs/contigs.fasta ${sample}
    echo "${sample} ALIGNMENT 700" >> /mnt/data/log.out
    bowtie2 -p4 -t -x ${sample} -q -1 /mnt/data/reads/${sample}_R1_cdh_lzw_trim30_PF.fastq.gz \
    -2 /mnt/data/reads/${sample}_R2_cdh_lzw_trim30_PF.fastq.gz  \
    --very-sensitive-local -X 700 -S /mnt/data/sam/${sample}_700.sam --no-unal --al-conc /mnt/data/sam/${sample}_conc_700.fq \
    2>> /mnt/data/log_may.out

    echo "${sample} ALIGNMENT 1000" >> /mnt/data/log.out
    bowtie2 -p4 -t -x ${sample} -q -1 /mnt/data/reads/${sample}_R1_cdh_lzw_trim30_PF.fastq.gz \
    -2 /mnt/data/reads/${sample}_R2_cdh_lzw_trim30_PF.fastq.gz  \
    --very-sensitive-local -X 1000 -S /mnt/data/sam/${sample}_1000.sam --no-unal --al-conc /mnt/data/sam/${sample}_conc_1000.fq \
    2>> /mnt/data/log_may.out

    samtools view -S -b /mnt/data/sam/${sample}_700.sam > /mnt/data/bam/${sample}_700.bam
    samtools view -S -b /mnt/data/sam/${sample}_1000.sam > /mnt/data/bam/${sample}_1000.bam
    rm -r /mnt/data/reads /mnt/data/contigs
    samtools view -f 0x2 /mnt/data/bam/${sample}_700.bam | grep -v "XS:i:" | /mnt/data/foo.py | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/counts/${sample}_bowtie_csp_700_counts.txt
    samtools view -F 260 /mnt/data/bam/${sample}_700.bam | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}'  > /mnt/data/counts/${sample}_bowtie_all_700_counts.txt
    samtools view -f 0x2 /mnt/data/bam/${sample}_700.bam | grep -v "XS:i:" | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/counts/${sample}_bowtie_cs_700_counts.txt

    samtools view -f 0x2 /mnt/data/bam/${sample}_1000.bam | grep -v "XS:i:" | /mnt/data/foo.py | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/counts/${sample}_bowtie_csp_1000_counts.txt
    samtools view -F 260 /mnt/data/bam/${sample}_1000.bam | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}'  > /mnt/data/counts/${sample}_bowtie_all_1000_counts.txt
    samtools view -f 0x2 /mnt/data/bam/${sample}_1000.bam | grep -v "XS:i:" | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/counts/${sample}_bowtie_cs_1000_counts.txt

    aws s3 cp /mnt/data/bam/${sample}_700.bam s3://czbiohub-mosquito/alignment/bam_700/${sample}.bam
    aws s3 cp /mnt/data/bam/${sample}_1000.bam s3://czbiohub-mosquito/alignment/bam_1000/${sample}.bam
    echo "${sample} done" >> /mnt/data/log2.out

    rm -r /mnt/data/bam /mnt/data/sam
  done

#upload all count files to s3
for sample in `ls /mnt/data/names`
  do
    aws s3 cp /mnt/data/counts/${sample}_bowtie_csp_700_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_csp_counts_700.txt
    aws s3 cp /mnt/data/counts/${sample}_bowtie_all_700_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_all_counts_700.txt
    aws s3 cp /mnt/data/counts/${sample}_bowtie_cs_700_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_cs_counts_700.txt
    aws s3 cp /mnt/data/counts/${sample}_bowtie_csp_1000_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_csp_counts_1000.txt
    aws s3 cp /mnt/data/counts/${sample}_bowtie_all_1000_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_all_counts_1000.txt
    aws s3 cp /mnt/data/counts/${sample}_bowtie_cs_1000_counts.txt s3://czbiohub-mosquito/contig_quality/${sample}/bowtie_cs_counts_1000.txt
    echo "${sample} done" >> /mnt/data/log3.out
  done



#foo.py script
#!/usr/bin/env python
import csv
import sys

f = csv.reader(sys.stdin, dialect="excel-tab")
of = csv.writer(sys.stdout, dialect="excel-tab")
last_read = None
for line in f :
    #take care of the header
    if(line[0][0] == "@") :
        of.writerow(line)
        continue

    if(last_read == None) :
        last_read = line
    else :
        if(last_read[0] == line[0]) :
            of.writerow(last_read)
            of.writerow(line)
            last_read = None
        else :
            last_read = line



#/mnt/data/names.sh
touch CMS001_001_Ra_S1
touch CMS001_002_Ra_S1
touch CMS001_003_Ra_S2
touch CMS001_004_Ra_S2
touch CMS001_005_Ra_S3
touch CMS001_006_Ra_S5
touch CMS001_007_Ra_S12
touch CMS001_008_Ra_S3
touch CMS001_009_Ra_S13
touch CMS001_010_Ra_S1
touch CMS001_011_Ra_S4
touch CMS001_012_Ra_S4
touch CMS001_013_Ra_S5
touch CMS001_014_Ra_S5
touch CMS001_015_Ra_S13
touch CMS001_016_Ra_S6
touch CMS001_017_Ra_S6
touch CMS001_018_Ra_S14
touch CMS001_019_Ra_S14
touch CMS001_020_Ra_S15
touch CMS001_021_Ra_S16
touch CMS001_022_Ra_S6
touch CMS001_023_Ra_S17
touch CMS001_024_Ra_S15
touch CMS001_025_Ra_S7
touch CMS001_026_Ra_S18
touch CMS001_027_Ra_S16
touch CMS001_028_Ra_S17
touch CMS001_029_Ra_S18
touch CMS001_030_Ra_S7
touch CMS001_031_Ra_S19
touch CMS001_032_Ra_S7
touch CMS001_033_Ra_S8
touch CMS001_034_Ra_S19
touch CMS001_035_Ra_S20
touch CMS001_036_Ra_S20
touch CMS001_037_Ra_S21
touch CMS001_038_Ra_S22
touch CMS001_039_Ra_S9
touch CMS001_040_Ra_S21
touch CMS001_041_Ra_S10
touch CMS001_042_Ra_S23
touch CMS001_043_Ra_S24
touch CMS001_044_Ra_S25
touch CMS001_045_Ra_S2
touch CMS001_046_Ra_S3
touch CMS001_047_Ra_S4
touch CMS001_048_Ra_S5
touch CMS001_049_Ra_S6
touch CMS001_050_Ra_S23
touch CMS001_051_Ra_S8
touch CMS001_052_Ra_S7
touch CMS001_053_Ra_S8
touch CMS001_054_Ra_S11
touch CMS001_055_Ra_S9
touch CMS001_056_Ra_S10
touch CMS001_057_Ra_S11
touch CMS001_058_Ra_S9
touch CMS001_059_Ra_S10
touch CMS001_060_Ra_S12
touch CMS001_water1_S11
touch CMS001_water2_S24
touch CMS001_water3_Qiagen_S26
touch CMS001_water4_Zymo_S27
touch CMS001_0water5_RNA_A_S12
touch CMS002_001a_Rb_S116_L004
touch CMS002_004a_Rb_S117_L004
touch CMS002_007a_Rb_S118_L004
touch CMS002_010a_Rb_S119_L004
touch CMS002_013a_Rb_S120_L004
touch CMS002_016a_Rb_S121_L004
touch CMS002_017a_Rb_S122_L004
touch CMS002_017b_Rb_S123_L004
touch CMS002_017c_Rb_S124_L004
touch CMS002_017d_Rb_S125_L004
touch CMS002_017e_Rb_S126_L004
touch CMS002_018a_Rb_S128_L004
touch CMS002_018b_Rb_S129_L004
touch CMS002_019a_Rb_S130_L004
touch CMS002_020a_Rb_S131_L004
touch CMS002_020b_Rb_S132_L004
touch CMS002_020c_Rb_S133_L004
touch CMS002_020d_Rb_S134_L004
touch CMS002_020e_Rb_S135_L004
touch CMS002_021a_Rb_S136_L004
touch CMS002_022a_Rb_S137_L004
touch CMS002_023a_Rb_S138_L004
touch CMS002_025a_Rb_S140_L004
touch CMS002_025b_Rb_S141_L004
touch CMS002_025c_Rb_S142_L004
touch CMS002_025d_Rb_S143_L004
touch CMS002_025e_Rb_S144_L004
touch CMS002_025f_Rb_S145_L004
touch CMS002_026a_Rb_S146_L004
touch CMS002_026b_Rb_S147_L004
touch CMS002_026c_Rb_S148_L004
touch CMS002_026d_Rb_S149_L004
touch CMS002_026e_Rb_S150_L004
touch CMS002_027a_Rb_S152_L004
touch CMS002_027b_Rb_S153_L004
touch CMS002_028a_Rb_S154_L004
touch CMS002_028b_Rb_S155_L004
touch CMS002_028c_Rb_S156_L004
touch CMS002_028d_Rb_S157_L004
touch CMS002_028e_Rb_S158_L004
touch CMS002_029a_Rb_S159_L004
touch CMS002_029b_Rb_S160_L004
touch CMS002_029c_Rb_S161_L004
touch CMS002_029d_Rb_S162_L004
touch CMS002_029e_Rb_S164_L004
touch CMS002_031a_Rb_S165_L004
touch CMS002_032a_Rb_S166_L004
touch CMS002_033a_Rb_S167_L004
touch CMS002_034a_Rb_S168_L004
touch CMS002_035a_Rb_S169_L004
touch CMS002_036a_Rb_S170_L004
touch CMS002_037a_Rb_S171_L004
touch CMS002_038a_Rb_S172_L004
touch CMS002_039a_Rb_S173_L004
touch CMS002_040a_Rb_S174_L004
touch CMS002_041a_Rb_S176_L004
touch CMS002_042a_Rb_S177_L004
touch CMS002_044a_Rb_S178_L004
touch CMS002_044b_Rb_S179_L004
touch CMS002_044c_Rb_S180_L004
touch CMS002_044d_Rb_S181_L004
touch CMS002_044e_Rb_S182_L004
touch CMS002_045a_Rb_S183_L004
touch CMS002_045b_Rb_S184_L004
touch CMS002_045c_Rb_S185_L004
touch CMS002_045d_Rb_S186_L004
touch CMS002_045e_Rb_S188_L004
touch CMS002_045f_Rb_S189_L004
touch CMS002_045g_Rb_S190_L004
touch CMS002_046a_Rb_S191_L004
touch CMS002_046b_Rb_S192_L004
touch CMS002_047a_Rb_S193_L004
touch CMS002_047b_Rb_S194_L004
touch CMS002_047c_Rb_S195_L004
touch CMS002_047d_Rb_S196_L004
touch CMS002_047e_Rb_S197_L004
touch CMS002_047f_Rb_S198_L004
touch CMS002_047g_Rb_S200_L004
touch CMS002_047h_Rb_S1_L004
touch CMS002_047i_Rb_S2_L004
touch CMS002_047j_Rb_S3_L004
touch CMS002_049a_Rb_S4_L004
touch CMS002_050a_Rb_S5_L004
touch CMS002_051a_Rb_S6_L004
touch CMS002_053a_Rb_S7_L004
touch CMS002_054a_Rb_S8_L004
touch CMS002_056a_Rb_S9_L004
touch CMS002_057a_Rb_S10_L004
touch CMS002_Water1_Rb_S127_L004
touch CMS002_Water2_Rb_S139_L004
touch CMS002_Water3_Rb_S151_L004
touch CMS002_Water4_Rb_S163_L004
touch CMS002_Water5_Rb_S175_L004
touch CMS002_Water6_Rb_S187_L004
touch CMS002_Water7_Rb_S199_L004
touch CMS002_Water8_Rb_S11_L004
touch CMS002_001a_Rb_S116_L004
touch CMS002_004a_Rb_S117_L004
touch CMS002_007a_Rb_S118_L004
touch CMS002_010a_Rb_S119_L004
touch CMS002_013a_Rb_S120_L004
touch CMS002_016a_Rb_S121_L004
touch CMS002_017a_Rb_S122_L004
touch CMS002_017b_Rb_S123_L004
touch CMS002_017c_Rb_S124_L004
touch CMS002_017d_Rb_S125_L004
touch CMS002_017e_Rb_S126_L004
touch CMS002_018a_Rb_S128_L004
touch CMS002_018b_Rb_S129_L004
touch CMS002_019a_Rb_S130_L004
touch CMS002_020a_Rb_S131_L004
touch CMS002_020b_Rb_S132_L004
touch CMS002_020c_Rb_S133_L004
touch CMS002_020d_Rb_S134_L004
touch CMS002_020e_Rb_S135_L004
touch CMS002_021a_Rb_S136_L004
touch CMS002_022a_Rb_S137_L004
touch CMS002_023a_Rb_S138_L004
touch CMS002_025a_Rb_S140_L004
touch CMS002_025b_Rb_S141_L004
touch CMS002_025c_Rb_S142_L004
touch CMS002_025d_Rb_S143_L004
touch CMS002_025e_Rb_S144_L004
touch CMS002_025f_Rb_S145_L004
touch CMS002_026a_Rb_S146_L004
touch CMS002_026b_Rb_S147_L004
touch CMS002_026c_Rb_S148_L004
touch CMS002_026d_Rb_S149_L004
touch CMS002_026e_Rb_S150_L004
touch CMS002_027a_Rb_S152_L004
touch CMS002_027b_Rb_S153_L004
touch CMS002_028a_Rb_S154_L004
touch CMS002_028b_Rb_S155_L004
touch CMS002_028c_Rb_S156_L004
touch CMS002_028d_Rb_S157_L004
touch CMS002_028e_Rb_S158_L004
touch CMS002_029a_Rb_S159_L004
touch CMS002_029b_Rb_S160_L004
touch CMS002_029c_Rb_S161_L004
touch CMS002_029d_Rb_S162_L004
touch CMS002_029e_Rb_S164_L004
touch CMS002_031a_Rb_S165_L004
touch CMS002_032a_Rb_S166_L004
touch CMS002_033a_Rb_S167_L004
touch CMS002_034a_Rb_S168_L004
touch CMS002_035a_Rb_S169_L004
touch CMS002_036a_Rb_S170_L004
touch CMS002_037a_Rb_S171_L004
touch CMS002_038a_Rb_S172_L004
touch CMS002_039a_Rb_S173_L004
touch CMS002_040a_Rb_S174_L004
touch CMS002_041a_Rb_S176_L004
touch CMS002_042a_Rb_S177_L004
touch CMS002_044a_Rb_S178_L004
touch CMS002_044b_Rb_S179_L004
touch CMS002_044c_Rb_S180_L004
touch CMS002_044d_Rb_S181_L004
touch CMS002_044e_Rb_S182_L004
touch CMS002_045a_Rb_S183_L004
touch CMS002_045b_Rb_S184_L004
touch CMS002_045c_Rb_S185_L004
touch CMS002_045d_Rb_S186_L004
touch CMS002_045e_Rb_S188_L004
touch CMS002_045f_Rb_S189_L004
touch CMS002_045g_Rb_S190_L004
touch CMS002_046a_Rb_S191_L004
touch CMS002_046b_Rb_S192_L004
touch CMS002_047a_Rb_S193_L004
touch CMS002_047b_Rb_S194_L004
touch CMS002_047c_Rb_S195_L004
touch CMS002_047d_Rb_S196_L004
touch CMS002_047e_Rb_S197_L004
touch CMS002_047f_Rb_S198_L004
touch CMS002_047g_Rb_S200_L004
touch CMS002_047h_Rb_S1_L004
touch CMS002_047i_Rb_S2_L004
touch CMS002_047j_Rb_S3_L004
touch CMS002_049a_Rb_S4_L004
touch CMS002_050a_Rb_S5_L004
touch CMS002_051a_Rb_S6_L004
touch CMS002_053a_Rb_S7_L004
touch CMS002_054a_Rb_S8_L004
touch CMS002_056a_Rb_S9_L004
touch CMS002_057a_Rb_S10_L004
touch CMS002_0Water1_Rb_S127_L004
touch CMS002_0Water2_Rb_S139_L004
touch CMS002_0Water3_Rb_S151_L004
touch CMS002_0Water4_Rb_S163_L004
touch CMS002_0Water5_Rb_S175_L004
touch CMS002_0Water6_Rb_S187_L004
touch CMS002_0Water7_Rb_S199_L004
touch CMS002_0Water8_Rb_S11_L004
