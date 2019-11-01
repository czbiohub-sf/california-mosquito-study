import boto3
import pandas as pd
from ete3 import NCBITaxa
import subprocess
import itertools
import os
import s3fs
import numpy as np
import matplotlib.pyplot as plt
from sqlalchemy import create_engine
from lca_functions import *
import pdb


import argparse, sys
parser = argparse.ArgumentParser()
parser.add_argument('--fpath', help="Path to file (local or on S3). For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
args = parser.parse_args()


s3 = boto3.resource('s3')
client = boto3.client('s3')
bucket_name = "czbiohub-mosquito"
bucket = s3.Bucket(bucket_name)
contig_folders = [x["Prefix"] for x in client.list_objects(Bucket=bucket_name, Prefix="contigs/", Delimiter="/")["CommonPrefixes"]]
contig_quality_folders = [x["Prefix"] for x in client.list_objects(Bucket=bucket_name, Prefix="contig_quality/", Delimiter="/")["CommonPrefixes"]]

ncbi = NCBITaxa()


read_count_files = [client.list_objects(Bucket=bucket_name, Prefix=x+"bowtie_csp_counts_1000.txt") \
                   for x in contig_quality_folders]
read_count_files = ["s3://"+bucket_name+"/"+x["Prefix"] for x in read_count_files if "Contents" in x.keys()]
read_counts_csp_1000 = pd.concat([pd.read_csv(x, sep="\t", header=None, names=["query", "read_count"]).\
                                  assign(sample=os.path.split(os.path.split(x)[0])[1]) for x in read_count_files])


filtered_contigs_by_read_count = read_counts_csp_1000[read_counts_csp_1000["read_count"]>2]

def create_blast_output (input_file, output_dir, output_fn, read_counts_df, gi2taxid_db, 
                         blast_headings=None):
    default_blast_headings = ["query", "subject", "identity", "align_length", "mismatches", 
        "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "taxid", 
        "sci_name", "common_name", "subject_title"]
    if (blast_headings is not None):
        default_blast_headings = blast_headings
    output = pd.read_csv(input_file, sep="\t", header=None, names=default_blast_headings)
    output["query"] = output["query"].apply(lambda x: x[x.index("~")+1:])
    output = output.assign(gi=output["subject"].str.split('|').apply(lambda x: x[1]))
    output = output.merge(read_counts_df, on="query", how="inner")
    db = create_engine(gi2taxid_db)
    gi2taxid = {}
    gi_list_str = ', '.join(np.sort(output["gi"].unique()))
    sql_str = "SELECT gi,taxid FROM gi2taxid WHERE gi IN ("+gi_list_str+")"
    for row in db.execute(sql_str).fetchall():
        gi2taxid[str(row[0])] = str(row[1])
    output["taxid"] = [gi2taxid[x] if x in gi2taxid else None for x in output["gi"]]
    list_of_samples = list(output["sample"].unique())
    for sample_name in list_of_samples:
        new_results = output[output["sample"] == sample_name].drop("sample", axis=1)
        new_results = new_results.drop(["gi", "read_count"], axis=1)
        filename = output_dir+sample_name+"/"+output_fn
        try:
            previous_output = pd.read_csv(filename, sep="\t", header=None, index_col=False, comment="#")
            previous_output = previous_output.append(new_results)
        except:
            previous_output = new_results
        df_to_s3(previous_output, filename, header=False)
    print (input_file+" ended.")



create_blast_output(args.fpath, output_dir="s3://lucymli/skeeters/blast_nr/", output_fn="blast_nr.m9"+os.path.basename(args.fpath), read_counts_df = filtered_contigs_by_read_count, gi2taxid_db="sqlite:///prot.gi2taxid.db", blast_headings=default_blast_headings)




