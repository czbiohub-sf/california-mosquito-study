import os
import urllib
import urllib.request
import time
import re
import logging
import subprocess

logging.basicConfig(level=logging.INFO)

outdir = "scratch/ncbi_download"
os.makedirs(outdir, exist_ok=True)

gc_re = re.compile(r">(GC.)_(\d\d\d)(\d\d\d)(\d\d\d)(.*)_genomic.fna.gz")
nc_re = re.compile(r">(NC_\d+\.\d+)\s")

with open("scratch/mosquito_genomes2_sequences.txt") as f:
    for line in f:
        gc_matched = gc_re.match(line)
        nc_matched = nc_re.match(line)
        if gc_matched:
            groups = [gc_matched.group(i) for i in range(1, 6)]
            sample = "{0}_{1}{2}{3}{4}".format(*groups)
            fname = "{}_genomic.fna.gz".format(sample)
            ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/{sample}/{fname}".format(*groups, sample=sample, fname=fname)
            logging.info(f"Downloading from {ftp_path}")
            urllib.request.urlretrieve(ftp_path, f"{outdir}/{fname}")
        elif nc_matched:
            sample = nc_matched.group(1)
            cmd = f"efetch -db nucleotide -id {sample} -format fasta -mode text | gzip -c > {outdir}/{sample}.fa.gz"
            logging.info(cmd)
            subprocess.call(cmd, shell=True)
        else:
            assert False
        time.sleep(1)
