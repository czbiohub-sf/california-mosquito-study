#Download FASTQs from Sequencing Folder
aws s3 cp s3://czbiohub-mosquito/sequences/CMS001_fastq.gz /mnt/data/fastq/
aws s3 cp s3://czbiohub-mosquito/sequences/CMS001_fastq.gz /mnt/data/fastq/

#gunzip all the files
#IN PARALLEL download GNU Parallel through https://www.gnu.org/software/parallel/parallel_tutorial.html
#(wget -O - pi.dk/3 || curl pi.dk/3/ || \ fetch -o - http://pi.dk/3) | bash
cd /mnt/data/fastq
parallel --gnu gunzip  ::: *.gz

#Convert FASTQ to FASTA (parallel to IDseq Pipeline)
#download seqtk: conda install seqtk
mkdir fasta
for sample in `ls /mnt/data/fastq/*R1_001.fastq`
  do
    dir="/mnt/data/fastq"
    base=$(basename $sample "_R1_001.fastq")
    seqtk seq -A  ${dir}/${base}_R1_001.fastq  > /mnt/data/fasta/${base}_R1_001.fasta
    seqtk seq -A  ${dir}/${base}_R2_001.fastq  > /mnt/data/fasta/${base}_R2_001.fasta
  done

#Create file with counts from fastas. will use fastq counts off idseq... takes too much unecessary time otherwise
touch /mnt/data/log_fasta.out
for sample in `ls /mnt/data/fasta/*R1.fasta`
  do
    #dir="/mnt/data/reads"
    dir2="/mnt/data/fasta"
    base=$(basename $sample "_R1.fasta")
    echo "$base R1 FASTA" >> /mnt/data/log_fasta.out
    grep -c "^>" ${dir2}/${base}_R1.fasta >> /mnt/data/log_fasta.out
    #echo "$base R1 FASTQ"
    #awk '{s++}END{print s/4}' ${dir}/${base}_R1_001.fastq
    echo "$base R2 FASTA" >> /mnt/data/log_fasta.out
    grep -c "^>" ${dir2}/${base}_R2.fasta >> /mnt/data/log_fasta.out
    #echo "$base R2 FASTQ"
    #awk '{s++}END{print s/4}' ${dir}/${base}_R2_001.fastq
  done

#cd-hit-dup paired end method
#download cd-hit-dup:
mkdir /mnt/data/cdh
for sample in `ls /mnt/data/fasta/*R1_001.fasta`
  do
    dir="/mnt/data/fasta"
    base=$(basename $sample "_R1_001.fasta")
    cd-hit-dup -i ${dir}/${base}_R1_001.fasta -i2 ${dir}/${base}_R2_001.fasta -o /mnt/data/cdh/${base}_R1_cdh.fasta -o2 /mnt/data/cdh/${base}_R2_cdh.fasta -e 0.05 -u 70
  done

#count cdh files
touch /mnt/data/log_cdh.out
for sample in `ls /mnt/data/cdh2/*R1_cdh.fasta`
  do
    dir="/mnt/data/cdh2"
    base=$(basename $sample "_R1_cdh.fasta")
    echo "$base R1 FASTA CDH" >> /mnt/data/log_cdh.out
    grep -c "^>" ${dir}/${base}_R1_cdh.fasta >> /mnt/data/log_cdh.out
    echo "$base R2 FASTA CDH" >> /mnt/data/log_cdh.out
    grep -c "^>" ${dir}/${base}_R2_cdh.fasta >> /mnt/data/log_cdh.out
  done


#running LZW on _cdh.fasta files..
mkdir lzw
#setting up LZW
git clone https://github.com/chanzuckerberg/idseq-dag.git
cd idseq-dag; pip3 install -e .
nano /mnt/data/idseq-dag/idseq_dag/__main__.py
##ALTER BOTTOM OF MAIN SCRIPT FOR THIS
'''
from idseq_dag.steps.run_lzw import PipelineStepRunLZW

if __name__ == "__main__":
    step_instance = PipelineStepRunLZW("lzw step", [], [],
                              "", "", "",
                              "", "")
    step_instance.run()
'''

nano /mnt/data/idseq-dag/idseq_dag/steps/run_lzw.py
## ALTER THE run(self) entire clause to the following
'''
def run(self):
        input_fas = ["file1.fa", "file2.fa"]
        output_fas = ["output1.fa", "output2.fa"]
        cutoff_fractions = [0.45, 0.42]
        PipelineStepRunLZW.generate_lzw_filtered(input_fas, output_fas, cutoff_fractions,150)
'''
#run LZW
for sample in `ls /mnt/data/cdh/*R1_cdh.fasta`
  do
    dir="/mnt/data/cdh"
    dir2="/mnt/data/lzw"
    base=$(basename $sample "_R1_cdh.fasta")
    mv ${dir}/${base}_R1_cdh.fasta /mnt/data/cdh/file1.fa
    mv ${dir}/${base}_R2_cdh.fasta /mnt/data/cdh/file2.fa
    python3 -m idseq_dag
    mv /mnt/data/cdh/output1.fa ${dir2}/${base}_R1_cdh_lzw.fasta
    mv /mnt/data/cdh/output2.fa ${dir2}/${base}_R2_cdh_lzw.fasta
  done

#count LZW reads
touch /mnt/data/log_lzw.out
for sample in `ls /mnt/data/lzw/*R1_cdh_lzw.fasta`
  do
    dir="/mnt/data/lzw"
    base=$(basename $sample "_R1_cdh_lzw.fasta")
    echo "$base R1 FASTA CDH LZW" >> /mnt/data/log_lzw.out
    grep -c "^>" ${dir}/${base}_R1_cdh_lzw.fasta >> /mnt/data/log_lzw.out
    echo "$base R2 FASTA CDH LZW" >> /mnt/data/log_lzw.out
    grep -c "^>" ${dir}/${base}_R2_cdh_lzw.fasta >> /mnt/data/log_lzw.out
  done

#Find FASTQ from FASTA
conda install seqkit
mkdir /mnt/data/final; touch /mnt/data/log_fastq.out
for sample in `ls /mnt/data/lzw/*R1_cdh_lzw.fasta`
  do
    dir="/mnt/data/lzw"
    dir2="/mnt/data/fastq"
    dir3="/mnt/data/final"
    base=$(basename $sample "_R1_cdh_lzw.fasta")

    grep -e "^>" ${dir}/${base}_R1_cdh_lzw.fasta > ${dir}/R1a.lst
    sed 's|[>,]||g' ${dir}/R1a.lst > ${dir}/R1.lst
    seqkit grep -nf ${dir}/R1.lst ${dir2}/${base}_R1_001.fastq > ${dir3}/${base}_R1_cdh_lzw.fastq
    echo "$base R1 FASTA CDH LZW FASTQ" >> /mnt/data/log_fastq.out
    awk '{s++}END{print s/4}' ${dir3}/${base}_R1_cdh_lzw.fastq >> /mnt/data/log_fastq.out
    rm ${dir2}/${base}_R1_001.fastq
    grep -e "^>" ${dir}/${base}_R2_cdh_lzw.fasta > ${dir}/R2a.lst
    sed 's|[>,]||g' ${dir}/R2a.lst > ${dir}/R2.lst
    seqkit grep -nf ${dir}/R2.lst ${dir2}/${base}_R2_001.fastq > ${dir3}/${base}_R2_cdh_lzw.fastq
    echo "$base R2 FASTA CDH LZW FASTQ" >> /mnt/data/log_fastq.out
    awk '{s++}END{print s/4}' ${dir3}/${base}_R2_cdh_lzw.fastq >> /mnt/data/log_fastq.out
    rm ${dir2}/${base}_R2_001.fastq
  done

#Trim + PriceSeqFilter
#Download trimmomatic: conda install -c bioconda trimmomatic
mkdir /mnt/data/trim_unp /mnt/data/trim /mnt/data/qc; touch /mnt/data/trim_log.out; touch /mnt/data/PF_log.out
for sample in `ls /mnt/data/final/*R1_cdh_lzw.fastq`
  do
    dir="/mnt/data/final"
    base=$(basename $sample "_R1_cdh_lzw.fastq")
    dir2="/mnt/data/trim"
    dir3="/mnt/data/trim_unp"
    dir4="/mnt/data/qc"
    trimmomatic PE -threads 8 -phred33 ${dir}/${base}_R1_cdh_lzw.fastq ${dir}/${base}_R2_cdh_lzw.fastq ${dir2}/${base}_R1_cdh_lzw_trim30.fastq ${dir3}/${base}_R1_cdh_lzw_trim30_unp.fastq ${dir2}/${base}_R2_cdh_lzw_trim30.fastq ${dir3}/${base}_R2_cdh_lzw_trim30_unp.fastq ILLUMINACLIP:/home/ubuntu/anaconda/share/trimmomatic-0.38-1/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30
    echo "$base R1 CDH LZW FASTQ TRIM30" >> /mnt/data/trim_log.out
    awk '{s++}END{print s/4}' ${dir2}/${base}_R1_cdh_lzw_trim30.fastq >> /mnt/data/trim_log.out
    echo "$base R2 CDH LZW FASTQ TRIM30" >> /mnt/data/log_fastq.out >> /mnt/data/trim_log.out
    awk '{s++}END{print s/4}' ${dir2}/${base}_R2_cdh_lzw_trim30.fastq >> /mnt/data/trim_log.out
    PriceSeqFilter -fp ${dir2}/${base}_R1_cdh_lzw_trim30.fastq ${dir2}/${base}_R2_cdh_lzw_trim30.fastq -op ${dir4}/${base}_R1_cdh_lzw_trim30_PF.fastq ${dir4}/${base}_R2_cdh_lzw_trim30_PF.fastq  -a 12 -rnf 90 -log c -rqf 85 0.98
    echo "$base R1 CDH LZW FASTQ TRIM30 PF" >> /mnt/data/PF_log.out
    awk '{s++}END{print s/4}' ${dir4}/${base}_R1_cdh_lzw_trim30_PF.fastq >> /mnt/data/PF_log.out
    echo "$base R2 CDH LZW FASTQ TRIM30 PF" >> /mnt/data/log_fastq.out >> /mnt/data/PF_log.out
    awk '{s++}END{print s/4}' ${dir4}/${base}_R2_cdh_lzw_trim30_PF.fastq >> /mnt/data/PF_log.out
  done

##gzipping files and uploading to s3
#download parallel GNU: https://www.gnu.org/software/parallel/parallel_tutorial.html
parallel gzip ::: *; aws s3 sync . s3://czbiohub-mosquito/sequences/<path>



## CMS001_053a manipulation in kalani-dedup4 because size of files way too large. split into 2
cat CMS_002_53a_Rb_S7_L004_R1_001.fastq | paste - - - - | sort -T /mnt/data/ -k1,1 -t " " | tr "\t" "\n" > CMS_002_53a_R1_sorted.fastq
cat CMS_002_53a_Rb_S7_L004_R2_001.fastq | paste - - - - | sort -T /mnt/data/  -k1,1 -t " " | tr "\t" "\n" > CMS_002_53a_R2_sorted.fastq
mkdir /mnt/data/new/
seqtk seq -A  CMS002_53a_R1_sorted.fastq > CMS002_53a_R1_sorted.fasta
seqtk seq -A  CMS002_53a_R2_sorted.fastq > CMS002_53a_R2_sorted.fasta
head -76159410 CMS002_53a_R1_sorted.fasta > CMS053a_R1a.fasta
head -76159410 CMS002_53a_R2_sorted.fasta > CMS053a_R2a.fasta
tail -76159412 CMS002_53a_R1_sorted.fasta > CMS053a_R1b.fasta
tail -76159412 CMS002_53a_R2_sorted.fasta > CMS053a_R2b.fasta
cd-hit-dup -i CMS053a_R1a.fasta -i2 CMS053a_R2a.fasta -o /mnt/data/new/CMS053a_R1a.fasta -o2 /mnt/data/new/CMS053a_R2a.fasta -e 0.05 -u 70
cd-hit-dup -i CMS053a_R1b.fasta -i2 CMS053a_R2b.fasta -o /mnt/data/new/CMS053a_R1b.fasta -o2 /mnt/data/new/CMS053a_R2b.fasta -e 0.05 -u 70

cat /mnt/data/new/CMS053a_R2b.fasta >> /mnt/data/new/CMS053a_R2a.fasta
cat /mnt/data/new/CMS053a_R1b.fasta >> /mnt/data/new/CMS053a_R1a.fasta
mkdir /mnt/data/cdh/
cd-hit-dup -i CMS053a_R1a.fasta -i2 /mnt/data/new/CMS053a_R2a.fasta -o /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta -o2 /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta -e 0.05 -u 70

mkdir 53a
cp /mnt/data/new/CMS053a_R1a.fasta /mnt/data/53a/CMS053a_R1.fasta
cat /mnt/data/new/CMS053a_R1b.fasta >> /mnt/data/53a/CMS053a_R1.fasta
cp /mnt/data/new/CMS053a_R2a.fasta /mnt/data/53a/CMS053a_R2.fasta
cat /mnt/data/new/CMS053a_R2b.fasta >> /mnt/data/53a/CMS053a_R2.fasta
cd-hit-dup -i /mnt/data/53a/CMS053a_R1.fasta -i2 /mnt/data/53a/CMS053a_R2.fasta -o /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta -o2 /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta -e 0.05 -u 70
cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta
cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta

mkdir /mnt/data/lzw
cd /mnt/data/lzw
cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta
cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta
mv /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta  /mnt/data/lzw/file1.fa
mv /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta /mnt/data/lzw/file2.fa
cd /mnt/data/lzw
python3 -m idseq_dag
mv /mnt/data/lzw/output1.fa /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fasta
mv /mnt/data/lzw/output2.fa /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fasta


#find the fastqs of the fastas
grep -e "^>" /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fasta > /mnt/data/lzw/R1a.lst
sed 's|[>,]||g' /mnt/data/lzw/R1a.lst > /mnt/data/lzw/R1.lst
seqkit grep -nf /mnt/data/lzw/R1.lst /mnt/data/test/CMS_002_53a_Rb_S7_L004_R1_001.fastq > /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fastq

grep -e "^>" /mnt/data/lzw/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fasta > /mnt/data/lzw/R2a.lst
sed 's|[>,]||g' /mnt/data/lzw/R2a.lst > /mnt/data/lzw/R2.lst
seqkit grep -nf /mnt/data/lzw/R2.lst /mnt/data/test/CMS_002_53a_Rb_S7_L004_R2_001.fastq > /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fastq

#trim
trimmomatic PE -threads 8 -phred33 /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fastq /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fastq /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30.fastq /mnt/data/test/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30_unp.fastq /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30.fastq /mnt/data/test/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30_unp.fastq ILLUMINACLIP:/home/ubuntu/anaconda/share/trimmomatic-0.38-0/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30
#PriceSeqFilter
PriceSeqFilter -fp /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30.fastq /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30.fastq -op /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30_PF.fastq /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30_PF.fastq  -a 12 -rnf 90 -log c -rqf 85 0.98

## COUNTING
#cdh
echo "$base R1 FASTA CDH" >> /mnt/data/log_lzw.out >> /mnt/data/053a.out
grep -c "^>" /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta >> /mnt/data/053a.out
echo "$base R2 FASTA CDH" >> /mnt/data/log_lzw.out >> /mnt/data/053a.out
grep -c "^>" /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta  >> /mnt/data/053a.out
#lzw
echo "$base R1 CDH LZW FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fastq   >> /mnt/data/053a.out
echo "$base R2 CDH LZW FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fastq   >> /mnt/data/053a.out
#trim30
echo "$base R1 CDH LZW TRIM30 FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30.fastq   >> /mnt/data/053a.out
echo "$base R2 CDH LZW TRIM30 FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30.fastq   >> /mnt/data/053a.out
#PF
echo "$base R1 CDH LZW TRIM30 FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30_PF.fastq   >> /mnt/data/053a.out
echo "$base R2 CDH LZW TRIM30 FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30_PF.fastq   >> /mnt/data/053a.out

#GZIP and UPLOAD to AWS
gzip /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta.clstr; aws s3 cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta.clstr.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_fasta.gz/
gzip /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta2.clstr; aws s3 cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta2.clstr.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_fasta.gz/
#cdh
gzip /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta; aws s3 cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_fasta.gz/
gzip /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta;  aws s3 cp /mnt/data/cdh/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_fasta.gz/
#lzw
gzip /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fastq; aws s3 cp /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fastq.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_fastq.gz/
gzip /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fastq; aws s3 cp /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fastq.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_fastq.gz/
#trim30
gzip /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30.fastq; aws s3 cp /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30.fastq.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_trim30_fastq.gz/
gzip /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30.fastq; aws s3 cp /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30.fastq.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_trim30_fastq.gz/
#PF
gzip /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30_PF.fastq; aws s3 cp /mnt/data/final/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30_PF.fastq.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_trim30_PF_fastq.gz/
gzip /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30_PF.fastq; aws s3 cp /mnt/data/final/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30_PF.fastq.gz s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_trim30_PF_fastq.gz/

#Pulldown and check off AWS
mkdir test2
aws s3 cp --recursive --exclude "*" --include "CMS_002_53*" s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_fasta.gz/  /mnt/data/test2
aws s3 cp --recursive --exclude "*" --include "CMS_002_53*"  s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_fastq.gz/ /mnt/data/test2
aws s3 cp --recursive --exclude "*" --include "CMS_002_53*"  s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_trim30_fastq.gz/ /mnt/data/test2
aws s3 cp --recursive --exclude "*" --include "CMS_002_53*" s3://czbiohub-mosquito/sequences/CMS002_cdhitdup_lzw_trim30_PF_fastq.gz/  /mnt/data/test2
cd /mnt/data/test2
parallel --gnu gunzip  ::: *.gz
touch /mnt/data/053a.out
#cdh
echo "CMS_002_53a_Rb_S7_L004 R1 FASTA CDH" >> /mnt/data/053a.out
grep -c "^>" /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R1_cdh.fasta >> /mnt/data/053a.out
echo "CMS_002_53a_Rb_S7_L004 R2 FASTA CDH" >> /mnt/data/053a.out
grep -c "^>" /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R2_cdh.fasta  >> /mnt/data/053a.out
#lzw
echo "CMS_002_53a_Rb_S7_L004 R1 CDH LZW FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw.fastq   >> /mnt/data/053a.out
echo "CMS_002_53a_Rb_S7_L004 R2 CDH LZW FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw.fastq   >> /mnt/data/053a.out
#trim30
echo "CMS_002_53a_Rb_S7_L004 R1 CDH LZW TRIM30  FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30.fastq   >> /mnt/data/053a.out
echo "CMS_002_53a_Rb_S7_L004 R2 CDH LZW TRIM30 FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30.fastq   >> /mnt/data/053a.out
#PF
echo "CMS_002_53a_Rb_S7_L004 R1 CDH LZW TRIM30 PF FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R1_cdh_lzw_trim30_PF.fastq   >> /mnt/data/053a.out
echo "CMS_002_53a_Rb_S7_L004 R2 CDH LZW TRIM30 PF FASTQ" >> /mnt/data/053a.out
awk '{s++}END{print s/4}' /mnt/data/test2/CMS_002_53a_Rb_S7_L004_R2_cdh_lzw_trim30_PF.fastq   >> /mnt/data/053a.out
