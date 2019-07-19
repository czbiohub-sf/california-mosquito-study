# need to create /mnt/data/foo.py file and /mnt/data/names.sh file, files are below
mkdir -p /mnt/data/refseqs /mnt/data/counts /mnt/data/names; cd /mnt/data/names; /mnt/data/names.sh; cd /mnt/data


##FOR REFERENCE ADD-ON
#download reference sequences to refseqs directory
aws s3 cp s3://czbiohub-mosquito/references/mosquito_genomes_20181207/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz /mnt/data/refseqs
aws s3 cp --recursive --exclude "*" --include "*.fasta" s3://czbiohub-mosquito/custom_host_filtering/ /mnt/data/refseqs
cd /mnt/data/refseqs
gunzip GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
cd /mnt/data


#loop through samples, make directories for reads & contigs files, upload reads and contigs from S3 to EC2
contig_align_bowtie () {
    sample=$1
    mkdir /mnt/data/${sample}reads /mnt/data/${sample}contigs /mnt/data/${sample}index /mnt/data/${sample}sam /mnt/data/${sample}bam
    aws s3 cp s3://czbiohub-mosquito/sequences/CMS_cleaned_reads/${sample}_R1_cdh_lzw_trim30_PF.fastq.gz /mnt/data/${sample}reads
    aws s3 cp s3://czbiohub-mosquito/sequences/CMS_cleaned_reads/${sample}_R2_cdh_lzw_trim30_PF.fastq.gz /mnt/data/${sample}reads
    mkdir -p /mnt/data/contigs/${sample}/
    aws s3 cp s3://czbiohub-mosquito/contigs/${sample}/contigs.fasta /mnt/data/contigs/${sample}/contigs.fasta

    ##FOR REFERENCE ADD-ON
    # concatenate contig sequences with refseqs for bowtie2 alignment
    cat /mnt/data/refseqs/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna /mnt/data/refseqs/custom*.fasta /mnt/data/contigs/${sample}/contigs.fasta >> /mnt/data/${sample}contigs/MosRefMosCustom.contigs.fasta

    ##FOR REFERENCE ADD-ON
    #build bowtie index from contigs+reference file
    cd /mnt/data/${sample}index
    echo "Starting index build for ${sample}Mos" >> /mnt/data/log.out
    bowtie2-build /mnt/data/${sample}contigs/MosRefMosCustom.contigs.fasta ${sample}_MosRefMosCustom

    ##FOR REFERENCE ADD-ON
    #run bowtie with insert size 1000 & log alignment start to log.out and stdout from alignment to log_may.out
    echo "Starting ${sample}Mos bowtie ALIGNMENT" >> /mnt/data/log.out
    bowtie2 -p4 -t -x ${sample}_MosRefMosCustom -q -1 /mnt/data/${sample}reads/${sample}_R1_cdh_lzw_trim30_PF.fastq.gz \
    -2 /mnt/data/${sample}reads/${sample}_R2_cdh_lzw_trim30_PF.fastq.gz  \
    --very-sensitive-local -X 1000 -S /mnt/data/${sample}sam/${sample}Mos_1000.sam --no-unal --al-conc /mnt/data/${sample}sam/${sample}Mos_conc_1000.fq \
    2>> /mnt/data/log_may.out

    #generate bam_1000 from sam_1000 output for insert size 1000 bowtie alignment to contigs
    echo "Starting ${sample}Mos.sam to ${sample}Mos.bam" >> /mnt/data/log2.out
    samtools view -S -b /mnt/data/${sample}sam/${sample}Mos_1000.sam > /mnt/data/${sample}bam/${sample}Mos_1000.bam

    #delete files that are extraneous at this stage
    rm -r /mnt/data/${sample}reads /mnt/data/${sample}contigs

    #compile counts from bam file for insert size 1000 with foo.py for concordant paired reads, all reads, and concordant single reads
    echo "Starting counts for ${sample}Mos" >> /mnt/data/log3.out
    samtools view -f 0x2 /mnt/data/${sample}bam/${sample}Mos_1000.bam | grep -v "XS:i:" | /mnt/data/foo.py | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/counts/${sample}Mos_bowtie_csp_1000_counts.txt
    samtools view -F 260 /mnt/data/${sample}bam/${sample}Mos_1000.bam | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}'  > /mnt/data/counts/${sample}Mos_bowtie_all_1000_counts.txt
    samtools view -f 0x2 /mnt/data/${sample}bam/${sample}Mos_1000.bam | grep -v "XS:i:" | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > /mnt/data/counts/${sample}Mos_bowtie_cs_1000_counts.txt

    #upload sam_1000 and bam_1000 files to aws s3 bucket
    aws s3 cp /mnt/data/${sample}bam/${sample}Mos_1000.bam s3://amyk-bucket/bam_1000/${sample}Mos.bam
    aws s3 cp /mnt/data/${sample}sam/${sample}Mos_1000.sam s3://amyk-bucket/sam_1000/${sample}Mos.sam
    #echo "${sample}Mos done" >> /mnt/data/log2.out

    #upload count files to aws s3 bucket
    echo "copying over count data for ${sample}Mos to aws" >> /mnt/data/log3.out
    aws s3 cp /mnt/data/counts/${sample}Mos_bowtie_csp_1000_counts.txt s3://amyk-bucket/${sample}Mos/bowtie_csp_counts_1000.txt
    aws s3 cp /mnt/data/counts/${sample}Mos_bowtie_all_1000_counts.txt s3://amyk-bucket/${sample}Mos/bowtie_all_counts_1000.txt
    aws s3 cp /mnt/data/counts/${sample}Mos_bowtie_cs_1000_counts.txt s3://amyk-bucket/${sample}Mos/bowtie_cs_counts_1000.txt
    echo "${sample}Mos done" >> /mnt/data/log3.out

    #delete files that are extraneous at this stage
    #rm -r /mnt/data/index /mnt/data/bam /mnt/data/sam
    echo "deleted /index /bam and /sam for ${sample}Mos" >> /mnt/data/log3.out
}

export -f contig_align_bowtie

ls /mnt/data/names | parallel -j 1 contig_align_bowtie {}

# #foo.py script
# #!/usr/bin/env python
# import csv
# import sys
#
# f = csv.reader(sys.stdin, dialect="excel-tab")
# of = csv.writer(sys.stdout, dialect="excel-tab")
# last_read = None
# for line in f :
#     #take care of the header
#     if(line[0][0] == "@") :
#         of.writerow(line)
#         continue
#
#     if(last_read == None) :
#         last_read = line
#     else :
#         if(last_read[0] == line[0]) :
#             of.writerow(last_read)
#             of.writerow(line)
#             last_read = None
#         else :
#             last_read = line
#
#
#
# #/mnt/data/names.sh
# touch CMS001_001_Ra_S1
# touch CMS001_002_Ra_S1
# touch CMS001_003_Ra_S2
# touch CMS001_004_Ra_S2
# touch CMS001_005_Ra_S3
# touch CMS001_006_Ra_S5
# touch CMS001_007_Ra_S12
# touch CMS001_008_Ra_S3
# touch CMS001_009_Ra_S13
# touch CMS001_010_Ra_S1
# touch CMS001_011_Ra_S4
# touch CMS001_012_Ra_S4
# touch CMS001_013_Ra_S5
# touch CMS001_014_Ra_S5
# touch CMS001_015_Ra_S13
# touch CMS001_016_Ra_S6
# touch CMS001_017_Ra_S6
# touch CMS001_018_Ra_S14
# touch CMS001_019_Ra_S14
# touch CMS001_020_Ra_S15
# touch CMS001_021_Ra_S16
# touch CMS001_022_Ra_S6
# touch CMS001_023_Ra_S17
# touch CMS001_024_Ra_S15
# touch CMS001_025_Ra_S7
# touch CMS001_026_Ra_S18
# touch CMS001_027_Ra_S16
# touch CMS001_028_Ra_S17
# touch CMS001_029_Ra_S18
# touch CMS001_030_Ra_S7
# touch CMS001_031_Ra_S19
# touch CMS001_032_Ra_S7
# touch CMS001_033_Ra_S8
# touch CMS001_034_Ra_S19
# touch CMS001_035_Ra_S20
# touch CMS001_036_Ra_S20
# touch CMS001_037_Ra_S21
# touch CMS001_038_Ra_S22
# touch CMS001_039_Ra_S9
# touch CMS001_040_Ra_S21
# touch CMS001_041_Ra_S10
# touch CMS001_042_Ra_S23
# touch CMS001_043_Ra_S24
# touch CMS001_044_Ra_S25
# touch CMS001_045_Ra_S2
# touch CMS001_046_Ra_S3
# touch CMS001_047_Ra_S4
# touch CMS001_048_Ra_S5
# touch CMS001_049_Ra_S6
# touch CMS001_050_Ra_S23
# touch CMS001_051_Ra_S8
# touch CMS001_052_Ra_S7
# touch CMS001_053_Ra_S8
# touch CMS001_054_Ra_S11
# touch CMS001_055_Ra_S9
# touch CMS001_056_Ra_S10
# touch CMS001_057_Ra_S11
# touch CMS001_058_Ra_S9
# touch CMS001_059_Ra_S10
# touch CMS001_060_Ra_S12
# touch CMS001_water1_S11
# touch CMS001_water2_S24
# touch CMS001_water3_Qiagen_S26
# touch CMS001_water4_Zymo_S27
# touch CMS001_0water5_RNA_A_S12
# touch CMS002_001a_Rb_S116_L004
# touch CMS002_004a_Rb_S117_L004
# touch CMS002_007a_Rb_S118_L004
# touch CMS002_010a_Rb_S119_L004
# touch CMS002_013a_Rb_S120_L004
# touch CMS002_016a_Rb_S121_L004
# touch CMS002_017a_Rb_S122_L004
# touch CMS002_017b_Rb_S123_L004
# touch CMS002_017c_Rb_S124_L004
# touch CMS002_017d_Rb_S125_L004
# touch CMS002_017e_Rb_S126_L004
# touch CMS002_018a_Rb_S128_L004
# touch CMS002_018b_Rb_S129_L004
# touch CMS002_019a_Rb_S130_L004
# touch CMS002_020a_Rb_S131_L004
# touch CMS002_020b_Rb_S132_L004
# touch CMS002_020c_Rb_S133_L004
# touch CMS002_020d_Rb_S134_L004
# touch CMS002_020e_Rb_S135_L004
# touch CMS002_021a_Rb_S136_L004
# touch CMS002_022a_Rb_S137_L004
# touch CMS002_023a_Rb_S138_L004
# touch CMS002_025a_Rb_S140_L004
# touch CMS002_025b_Rb_S141_L004
# touch CMS002_025c_Rb_S142_L004
# touch CMS002_025d_Rb_S143_L004
# touch CMS002_025e_Rb_S144_L004
# touch CMS002_025f_Rb_S145_L004
# touch CMS002_026a_Rb_S146_L004
# touch CMS002_026b_Rb_S147_L004
# touch CMS002_026c_Rb_S148_L004
# touch CMS002_026d_Rb_S149_L004
# touch CMS002_026e_Rb_S150_L004
# touch CMS002_027a_Rb_S152_L004
# touch CMS002_027b_Rb_S153_L004
# touch CMS002_028a_Rb_S154_L004
# touch CMS002_028b_Rb_S155_L004
# touch CMS002_028c_Rb_S156_L004
# touch CMS002_028d_Rb_S157_L004
# touch CMS002_028e_Rb_S158_L004
# touch CMS002_029a_Rb_S159_L004
# touch CMS002_029b_Rb_S160_L004
# touch CMS002_029c_Rb_S161_L004
# touch CMS002_029d_Rb_S162_L004
# touch CMS002_029e_Rb_S164_L004
# touch CMS002_031a_Rb_S165_L004
# touch CMS002_032a_Rb_S166_L004
# touch CMS002_033a_Rb_S167_L004
# touch CMS002_034a_Rb_S168_L004
# touch CMS002_035a_Rb_S169_L004
# touch CMS002_036a_Rb_S170_L004
# touch CMS002_037a_Rb_S171_L004
# touch CMS002_038a_Rb_S172_L004
# touch CMS002_039a_Rb_S173_L004
# touch CMS002_040a_Rb_S174_L004
# touch CMS002_041a_Rb_S176_L004
# touch CMS002_042a_Rb_S177_L004
# touch CMS002_044a_Rb_S178_L004
# touch CMS002_044b_Rb_S179_L004
# touch CMS002_044c_Rb_S180_L004
# touch CMS002_044d_Rb_S181_L004
# touch CMS002_044e_Rb_S182_L004
# touch CMS002_045a_Rb_S183_L004
# touch CMS002_045b_Rb_S184_L004
# touch CMS002_045c_Rb_S185_L004
# touch CMS002_045d_Rb_S186_L004
# touch CMS002_045e_Rb_S188_L004
# touch CMS002_045f_Rb_S189_L004
# touch CMS002_045g_Rb_S190_L004
# touch CMS002_046a_Rb_S191_L004
# touch CMS002_046b_Rb_S192_L004
# touch CMS002_047a_Rb_S193_L004
# touch CMS002_047b_Rb_S194_L004
# touch CMS002_047c_Rb_S195_L004
# touch CMS002_047d_Rb_S196_L004
# touch CMS002_047e_Rb_S197_L004
# touch CMS002_047f_Rb_S198_L004
# touch CMS002_047g_Rb_S200_L004
# touch CMS002_047h_Rb_S1_L004
# touch CMS002_047i_Rb_S2_L004
# touch CMS002_047j_Rb_S3_L004
# touch CMS002_049a_Rb_S4_L004
# touch CMS002_050a_Rb_S5_L004
# touch CMS002_051a_Rb_S6_L004
# touch CMS002_053a_Rb_S7_L004
# touch CMS002_054a_Rb_S8_L004
# touch CMS002_056a_Rb_S9_L004
# touch CMS002_057a_Rb_S10_L004
# touch CMS002_Water1_Rb_S127_L004
# touch CMS002_Water2_Rb_S139_L004
# touch CMS002_Water3_Rb_S151_L004
# touch CMS002_Water4_Rb_S163_L004
# touch CMS002_Water5_Rb_S175_L004
# touch CMS002_Water6_Rb_S187_L004
# touch CMS002_Water7_Rb_S199_L004
# touch CMS002_Water8_Rb_S11_L004
# touch CMS002_001a_Rb_S116_L004
# touch CMS002_004a_Rb_S117_L004
# touch CMS002_007a_Rb_S118_L004
# touch CMS002_010a_Rb_S119_L004
# touch CMS002_013a_Rb_S120_L004
# touch CMS002_016a_Rb_S121_L004
# touch CMS002_017a_Rb_S122_L004
# touch CMS002_017b_Rb_S123_L004
# touch CMS002_017c_Rb_S124_L004
# touch CMS002_017d_Rb_S125_L004
# touch CMS002_017e_Rb_S126_L004
# touch CMS002_018a_Rb_S128_L004
# touch CMS002_018b_Rb_S129_L004
# touch CMS002_019a_Rb_S130_L004
# touch CMS002_020a_Rb_S131_L004
# touch CMS002_020b_Rb_S132_L004
# touch CMS002_020c_Rb_S133_L004
# touch CMS002_020d_Rb_S134_L004
# touch CMS002_020e_Rb_S135_L004
# touch CMS002_021a_Rb_S136_L004
# touch CMS002_022a_Rb_S137_L004
# touch CMS002_023a_Rb_S138_L004
# touch CMS002_025a_Rb_S140_L004
# touch CMS002_025b_Rb_S141_L004
# touch CMS002_025c_Rb_S142_L004
# touch CMS002_025d_Rb_S143_L004
# touch CMS002_025e_Rb_S144_L004
# touch CMS002_025f_Rb_S145_L004
# touch CMS002_026a_Rb_S146_L004
# touch CMS002_026b_Rb_S147_L004
# touch CMS002_026c_Rb_S148_L004
# touch CMS002_026d_Rb_S149_L004
# touch CMS002_026e_Rb_S150_L004
# touch CMS002_027a_Rb_S152_L004
# touch CMS002_027b_Rb_S153_L004
# touch CMS002_028a_Rb_S154_L004
# touch CMS002_028b_Rb_S155_L004
# touch CMS002_028c_Rb_S156_L004
# touch CMS002_028d_Rb_S157_L004
# touch CMS002_028e_Rb_S158_L004
# touch CMS002_029a_Rb_S159_L004
# touch CMS002_029b_Rb_S160_L004
# touch CMS002_029c_Rb_S161_L004
# touch CMS002_029d_Rb_S162_L004
# touch CMS002_029e_Rb_S164_L004
# touch CMS002_031a_Rb_S165_L004
# touch CMS002_032a_Rb_S166_L004
# touch CMS002_033a_Rb_S167_L004
# touch CMS002_034a_Rb_S168_L004
# touch CMS002_035a_Rb_S169_L004
# touch CMS002_036a_Rb_S170_L004
# touch CMS002_037a_Rb_S171_L004
# touch CMS002_038a_Rb_S172_L004
# touch CMS002_039a_Rb_S173_L004
# touch CMS002_040a_Rb_S174_L004
# touch CMS002_041a_Rb_S176_L004
# touch CMS002_042a_Rb_S177_L004
# touch CMS002_044a_Rb_S178_L004
# touch CMS002_044b_Rb_S179_L004
# touch CMS002_044c_Rb_S180_L004
# touch CMS002_044d_Rb_S181_L004
# touch CMS002_044e_Rb_S182_L004
# touch CMS002_045a_Rb_S183_L004
# touch CMS002_045b_Rb_S184_L004
# touch CMS002_045c_Rb_S185_L004
# touch CMS002_045d_Rb_S186_L004
# touch CMS002_045e_Rb_S188_L004
# touch CMS002_045f_Rb_S189_L004
# touch CMS002_045g_Rb_S190_L004
# touch CMS002_046a_Rb_S191_L004
# touch CMS002_046b_Rb_S192_L004
# touch CMS002_047a_Rb_S193_L004
# touch CMS002_047b_Rb_S194_L004
# touch CMS002_047c_Rb_S195_L004
# touch CMS002_047d_Rb_S196_L004
# touch CMS002_047e_Rb_S197_L004
# touch CMS002_047f_Rb_S198_L004
# touch CMS002_047g_Rb_S200_L004
# touch CMS002_047h_Rb_S1_L004
# touch CMS002_047i_Rb_S2_L004
# touch CMS002_047j_Rb_S3_L004
# touch CMS002_049a_Rb_S4_L004
# touch CMS002_050a_Rb_S5_L004
# touch CMS002_051a_Rb_S6_L004
# touch CMS002_053a_Rb_S7_L004
# touch CMS002_054a_Rb_S8_L004
# touch CMS002_056a_Rb_S9_L004
# touch CMS002_057a_Rb_S10_L004
# touch CMS002_0Water1_Rb_S127_L004
# touch CMS002_0Water2_Rb_S139_L004
# touch CMS002_0Water3_Rb_S151_L004
# touch CMS002_0Water4_Rb_S163_L004
# touch CMS002_0Water5_Rb_S175_L004
# touch CMS002_0Water6_Rb_S187_L004
# touch CMS002_0Water7_Rb_S199_L004
# touch CMS002_0Water8_Rb_S11_L004
