# call this from the commandline in the skeeters/analysis directory:
# aws s3 ls s3://czbiohub-mosquito/contigs/ | grep 'PRE' | awk 'NF>1{print $NF}' | parallel python make_summary_stats.py {}


import pandas as pd
import numpy as np
import boto3
import tempfile
import os
import io
import time
import json
import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()




def get_remote_file(s3path):
    try:
        return (pd.read_csv(s3path, sep="\t", header=0))
    except:
        return None

def load_json_from_s3 (s3resource, bucket_name, fn):
    iostring = s3resource.Object(bucket_name, fn).get()['Body'].read().decode('utf-8')
    return (list(map(json.loads, io.StringIO(iostring).readlines()))[0])

def df_to_s3 (s3resource, obj, bucket_name, s3path, header=True):
    out_file = tempfile.NamedTemporaryFile()
    obj.to_csv(out_file.name, sep="\t", index=False, header=header)
    data = open(out_file.name, "rb")
    s3resource.Bucket(bucket_name).put_object(Key=s3path, Body=data)


def in_nt_nr (in_nt, in_nr, hexapoda):
    output_str = ""
    if (hexapoda):
        return (np.nan)
    if in_nt:
        output_str += "nt"
    if in_nr:
        output_str += "nr"
    if (len(output_str)==0):
        output_str = np.nan
    return (output_str)


def choose_nt_nr (contig_name, nt_blast, nr_blast):
    selected_nt = nt_blast[nt_blast["query"]==contig_name]
    if (len(selected_nt)==0):
        return ("nr")
    nt_bitscore = selected_nt["bitscore"].iloc[0]
    selected_nr = nr_blast[nr_blast["query"]==contig_name]
    if (len(selected_nr)==0):
        return ("nt")
    nr_bitscore = selected_nr["bitscore"].iloc[0]
    if (nr_bitscore > nt_bitscore):
        return ("nr")
    else:
        return ("nt")
    

def return_blast_results (row_x, blast_nt, blast_nr, blast_columns):
    if (blast_nt is None) & (blast_nr is None):
        exit()
    elif (blast_nt is None):
        blast = blast_nr
    elif (blast_nr is None):
        blast = blast_nt
    elif (row_x["nt_or_nr"]=="nt"):
        blast = blast_nt
    else:
        blast = blast_nr
    return(blast[blast["query"]==row_x["contig_name"]][blast_columns].head(n=1))


def get_tax_group (taxid):
    taxon_groups = ["Viruses", "Bacteria", "Archaea", "Metazoa", "Eukaryota"]
    taxon_groups_id = [ncbi.get_name_translator([x])[x][0] for x in taxon_groups]
    lineage = ncbi.get_lineage(taxid)
    for i, tax in enumerate(taxon_groups_id):
        if (tax in lineage):
            return (taxon_groups[i])
    return ("Ambiguous")

s3 = boto3.resource('s3')
bucket_name = "czbiohub-mosquito"


contig_folder_name = os.path.join("contigs", sys.argv[1])
contig_quality_folder_name = os.path.join("contig_quality", sys.argv[1])
sample_name = os.path.basename(os.path.dirname(contig_folder_name))

all_contig_fn = "contig_stats_all.tsv"
lca_contig_fn = "contig_stats_lca.tsv"


contig_stats_json = load_json_from_s3(s3, bucket_name, os.path.join(contig_folder_name, "contig_stats.json"))
if ("*" in contig_stats_json):
    contig_stats_json.pop("*")

contig_stats_df = pd.Series(contig_stats_json).reset_index(name="read_count").rename(columns={"index":"contig_name"})
contig_stats_df = contig_stats_df.assign(sample=sample_name)
contig_stats_df = contig_stats_df.assign(contig_length=contig_stats_df["contig_name"].str.split("_").apply(lambda x: int(x[3])))
contig_stats_df = contig_stats_df[["sample", "contig_name", "contig_length", "read_count"]]

if ((contig_stats_df["read_count"]>2).sum()==0):
    exit()

df_to_s3(s3, contig_stats_df, bucket_name, os.path.join(contig_quality_folder_name, all_contig_fn), header=True)



contig_quality_files = ["exclude_contigs_nt.txt", "exclude_contigs_nr.txt", 
                        "blast_lca_nt_filtered.m9", "blast_lca_nr_filtered.m9"]
blast_col_names = ["taxid", "bitscore", "align_length", "identity"]



contig_quality_dfs = dict(zip(contig_quality_files, [get_remote_file(os.path.join("s3://", bucket_name, contig_quality_folder_name, x)) for x in contig_quality_files]))



contig_stats_lca_df = contig_stats_df[contig_stats_df["read_count"]>2]
nt_contigs = contig_quality_dfs["exclude_contigs_nt.txt"][contig_quality_dfs["exclude_contigs_nt.txt"]["taxid_na"] | contig_quality_dfs["exclude_contigs_nt.txt"]["hexapoda"]]["query"].tolist()

if (contig_quality_dfs["blast_lca_nt_filtered.m9"] is not None):
    nt_contigs += contig_quality_dfs["blast_lca_nt_filtered.m9"]["query"].unique().tolist()

    
nr_contigs = contig_quality_dfs["exclude_contigs_nr.txt"][contig_quality_dfs["exclude_contigs_nr.txt"]["taxid_na"] | contig_quality_dfs["exclude_contigs_nr.txt"]["hexapoda"]]["query"].tolist()

    
if (contig_quality_dfs["blast_lca_nr_filtered.m9"] is not None):
    nr_contigs += contig_quality_dfs["blast_lca_nr_filtered.m9"]["query"].unique().tolist()

    
contig_stats_lca_df = contig_stats_lca_df.assign(nt=contig_stats_lca_df["contig_name"].isin(nt_contigs))
contig_stats_lca_df = contig_stats_lca_df.assign(nr=contig_stats_lca_df["contig_name"].isin(nr_contigs))
contig_stats_lca_df = contig_stats_lca_df[contig_stats_lca_df["nt"] | contig_stats_lca_df["nr"]]
contig_stats_lca_df = pd.merge(contig_stats_lca_df, contig_quality_dfs["exclude_contigs_nt.txt"][["query", "hexapoda"]], how="left", left_on="contig_name", right_on="query").drop(columns=["query"], axis=1)
contig_stats_lca_df = pd.merge(contig_stats_lca_df, contig_quality_dfs["exclude_contigs_nr.txt"][["query", "hexapoda"]], how="left", left_on="contig_name", right_on="query", suffixes=["_nt", "_nr"]).drop(columns=["query"], axis=1)
hexapoda_discordance = contig_stats_lca_df.apply(lambda x: ((x["hexapoda_nt"]==True)&(x["hexapoda_nr"]==False))|((x["hexapoda_nr"]==True)&(x["hexapoda_nt"]==False)), axis=1)


contig_stats_lca_df = contig_stats_lca_df.assign(hexapoda=(contig_stats_lca_df["hexapoda_nt"]==True) | (contig_stats_lca_df["hexapoda_nr"]==True))
contig_stats_lca_df = contig_stats_lca_df.drop(columns=["hexapoda_nt", "hexapoda_nr"])



contig_stats_lca_df = contig_stats_lca_df.assign(nt_or_nr=contig_stats_lca_df.apply(lambda x: in_nt_nr(x["nt"], x["nr"], x["hexapoda"]), axis=1))
contig_stats_lca_df = contig_stats_lca_df.assign(taxid=np.nan, bitscore=np.nan, identity=np.nan, align_length=np.nan)
nt_nr_contigs = contig_stats_lca_df[(contig_stats_lca_df["nt_or_nr"]=="ntnr")]["contig_name"]
if (len(nt_nr_contigs)>0):
    contig_stats_lca_df.loc[(contig_stats_lca_df["nt_or_nr"]=="ntnr"), "nt_or_nr"] = nt_nr_contigs.apply(lambda x: choose_nt_nr(x, nt_blast=contig_quality_dfs["blast_lca_nt_filtered.m9"], nr_blast=contig_quality_dfs["blast_lca_nr_filtered.m9"])).tolist()

for i in range(len(contig_stats_lca_df)):
    if not pd.isnull(contig_stats_lca_df.loc[i:i, "nt_or_nr"]).any():
        row_replacement = return_blast_results(contig_stats_lca_df.iloc[i], blast_nt=contig_quality_dfs["blast_lca_nt_filtered.m9"], blast_nr=contig_quality_dfs["blast_lca_nr_filtered.m9"], blast_columns=blast_col_names)
        contig_stats_lca_df.loc[i:i, blast_col_names] = row_replacement.values


contig_stats_lca_df = contig_stats_lca_df.assign(taxon_group=np.nan)
contig_stats_lca_df.loc[~contig_stats_lca_df["taxid"].isnull(), "taxon_group"] = contig_stats_lca_df["taxid"][~contig_stats_lca_df["taxid"].isnull()].apply(get_tax_group)


df_to_s3(s3, contig_stats_lca_df, bucket_name, os.path.join(contig_quality_folder_name, lca_contig_fn), header=True)

