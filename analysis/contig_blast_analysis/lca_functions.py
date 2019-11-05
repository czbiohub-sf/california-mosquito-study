##########
# Determine the Lowest Common Ancestor (LCA) of a set of taxonomic IDs
##########

import pandas as pd
import numpy as np
from ete3 import NCBITaxa
import boto3
import tempfile
import subprocess
import os
import io
import re
import time
import json

import pdb, traceback, sys

# NCBI Entrez functions
from Bio import Entrez
Entrez.email = "lucy.li@czbiohub.org"
api_key = "1a6a75bc7f8a5a3088510eb4f1b35eefa009"

# load ncbi taxonomy database
ncbi = NCBITaxa()
update_tax_database = False
if update_tax_database:
    ncbi.update_taxonomy_database()




default_blast_headings = ["query", "subject", "identity", "align_length", "mismatches", 
        "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "taxid", 
        "sci_name", "common_name", "subject_title"]

##
## For a given accession number, find the corresponding TaxID from NCBI's Taxonomy Database
##
def get_taxid (acc, db):
    result = int(Entrez.read(Entrez.esummary(id=str(acc), db=db, api_key=api_key))[0]["TaxId"])
    time.sleep(0.1)
    return (result)

##
## For a given accession number, return the genbank record
##
def get_gb (acc, db):
    if ("|" in acc):
        acc = acc.split('|')[1]
    result = list(Entrez.efetch(id=str(acc), db=db, rettype="gb", retmode="text", api_key=api_key))
    time.sleep(0.1)
    return (result)

##
## If an accession number is not associated with a TaxId, try to find that information elsewhere, such
##   as in a related NCBI record
##
def find_missing_taxid (acc, db):
    # Strategy 1: check if the record replaced an older record that did contain the TaxId Information
    try:
        match_1 = [get_taxid(re.search("gi:(.*).", x).group(1), db) for x in get_gb(acc, db) if ("replace" in x)]
        if (len(match_1)==1):
            return (int(match_1[0]))
        else:
            return (None)
    except:
        return (None)
    


##
## Find the lowest common ancestor (LCA) for a given set of taxonomic groups
##
def get_lca (taxids, tax_col="taxid", query_col="query"):
    # taxids should be a list of integers, e.g. [1, 12392, 178688] or
    #   a data frame with taxids stored in the `tax_col` column
    # returns an integer if input was a list, or a data frame if input
    #   was a data frame
    if isinstance(taxids, pd.DataFrame):
        lca = taxids[[query_col, tax_col]].iloc[[0]]
        lca.iloc[0, 1] = get_lca(taxids[tax_col].tolist())
    else:
        if (len(set(taxids)) <= 1):
            lca = taxids[0]
        else:
            tree = ncbi.get_topology(taxids)
            lca = tree.get_tree_root().taxid
    return (lca)

##
## Load json from s3 or local
##
def load_json (fpath, colnames):
    if (fpath.startswith("s3://")):
        s3 = boto3.resource('s3')
        s3_bucket_name, s3_path = split_s3_path(fpath)
        data_in_bytes = s3.Object(s3_bucket_name, s3_path).get()["Body"].read().decode('utf-8')
        json_data = list(map(json.loads, io.StringIO(data_in_bytes).readlines()))[0]
    else:
        json_data = json.loads(fpath)
    df = pd.DataFrame(pd.Series(json_data), columns=[colnames[1]]).reset_index(level=0).rename(columns={"index":colnames[0]})
    return (df)


##
## Read in BLAST results stored in tab-delimited files as a pandas data frame
##
def parse_blast_file (fpath, sep="\t", comment=None, blast_type="nt", col_names="auto"):
    # col_names:
    #   'auto' if the column headings should be auto-populated based on the # of columns
    #   '[a, b, c]' if the column headings should be set to a, b, c
    #   'None' if column headings should not be changed after reading in the file
    df = pd.read_csv(fpath, sep=sep, header=None, comment=comment)
    if (col_names=="auto"):
        # blast column headings:
        col_headings = ["query", "subject", "identity", "align_length", "mismatches", 
        "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "taxid", 
        "sci_name", "common_name", "subject_title"]
        df.columns = col_headings[:len(df.columns)]
    elif (col_names is not None):
        df.columns = col_names
    df = df.assign(blast_type=blast_type)
    return (df)


##
## Filter contigs or blast hits for contigs based on matches to a taxonomic id
##
def filter_by_taxid (df, db, taxid):
    check_isin = df["taxid"].apply(lambda x, taxid=taxid: taxid in ncbi.get_lineage(x))
    if (not check_isin.any()):
        return (df)
    if (check_isin.all()):
        return (df[:0])
    qlength = float(df["query"].iloc[0].split("_")[3])
    align_prop = df["align_length"]/qlength
    if (db=="protein"):
        align_prop = align_prop*3
    if (align_prop[check_isin].max() >= 0.8):
        return (df[:0])
    else:
        return (df)

##
## Get HSP for each query-subject pairing
##
def get_single_hsp (df_file, blast_type, col_names, coltypes):
    temp = tempfile.NamedTemporaryFile(delete=False)
    print (temp.name+" tempfile")
    if (df_file.startswith("s3://")):
        download_s3_file(df_file, temp.name)
    process = subprocess.Popen(["python", "parse.py", temp.name], stdout=subprocess.PIPE)
    csv = [line.decode().strip('"\n').split('\t') for line in process.stdout]
    outdf = pd.DataFrame(csv, columns=col_names).assign(blast_type=blast_type).astype(coltypes)
    os.unlink(temp.name)
    return (outdf)


##
## Select taxonomic IDs to perfrom LCA analysis on based on BLAST results
##
def select_taxids_for_lca (df, db="nucleotide", return_taxid_only=True, ident_cutoff=0, align_len_cutoff=0, bitscore_cutoff=0, read_counts=None):
    # df should be a pandas dataframe
    # remove blast hits where identity < ident_cutoff*max(identity) AND
    # align_length < align_len_cutoff*max(align_length) AND
    # bitscore < bitscore_cutoff*max(bitscore)
    if (len(df.index)>1):
        df = df[df["identity"]>=(ident_cutoff*df["identity"].max())]
        df = df[df["align_length"]>=(align_len_cutoff*df["align_length"].max())]
        df = df[df["bitscore"]>=(bitscore_cutoff*df["bitscore"].max())]
#     if (df["taxid"].isnull().any()):
#         df.loc[df["taxid"].isnull(), ["taxid"]] = df[df["taxid"].isnull()]["subject"].apply(find_missing_taxid, db=db)
#         df = df.dropna()
    df["taxid"] = df["taxid"].astype('int64')
    if (return_taxid_only):
        return (list(set(df["taxid"])))
    else:
        return (df)
#     try:
#         original_contigs = read_counts[read_counts["query"].isin(df["query"])]
#         filtered_df = filter_by_taxid(df, db=db, taxid=ncbi.get_name_translator(["Hexapoda"])["Hexapoda"][0])
#         if (len(filtered_df.index)==0):
#             excluded_contigs = original_contigs[~original_contigs["query"].isin(filtered_df["query"])].assign(reason="hexapoda")
#         else:
#             excluded_contigs = pd.DataFrame(columns=list(original_contigs.columns)+["reason"])
#     except:
#         pdb.set_trace()


##
## Split an s3://bucket/path string into the bucket name and the path
##
def split_s3_path (s3_path):
    s3_split = os.path.normpath(s3_path).split(os.sep)
    bucket_name = s3_split[1]
    s3_path = '/'.join(s3_split[2:])
    return bucket_name, s3_path

##
## Upload pandas dataframe to S3
##
def df_to_s3 (obj, s3path, header=True):
    s3 = boto3.resource('s3')
    client = boto3.client('s3')
    out_file = tempfile.NamedTemporaryFile()
    obj.to_csv(out_file.name, sep="\t", index=False, header=header)
    data = open(out_file.name, "rb")
    s3_bucket_name, s3_path = split_s3_path(s3path)
    s3.Bucket(s3_bucket_name).put_object(Key=s3_path, Body=data)

##
## Download file from S3 and return the local path to this file
##
def download_s3_file (s3path, local_path=None):
    s3 = boto3.resource('s3')
    client = boto3.client('s3')
    s3_bucket_name, s3_path = split_s3_path(s3path)
    if (local_path is not None):
        fpath = client.download_file(s3_bucket_name, s3_path, local_path)
    else:
        fpath = client.get_object(Bucket=s3_bucket_name, Key=s3_path)['Body']
        return fpath

##
## Filter blast records belonging to a particular taxid
##
def filter_by_lineage (df, taxid_col, lineage_id):
    ncbi = NCBITaxa()
    if (isinstance(lineage_id, str)):
        lineage_id = ncbi.get_name_translator([lineage_id])[lineage_id][0]
    descendants = ncbi.get_descendant_taxa(lineage_id)
    return (df[df[taxid_col].isin(descendants)])

##
## Print message to STDOUT
##
def print_to_stdout(message, start_time, verbose):
    if verbose:
        elapsed_time = round(time.time() - start_time, 2)
        print (message+"| elapsed time: "+str(elapsed_time)+" seconds")

##
## Combine blast results with LCA analysis
## 
def combine_blast_lca (lca_file_name, blast_file_name, outfile, sample_name, blast_type, output_file_name=None):
    lca_data = pd.read_csv(lca_file_name, sep="\t", header=0)
    blast_data = pd.read_csv(blast_file_name, sep="\t", header=0)
    blast_data_grouped = blast_data.groupby(["query"], as_index=False).\
    apply(lambda x: x[x["bitscore"]==max(x["bitscore"])].head(n=1))
    blast_data_grouped = blast_data_grouped[['query', 'identity', 'align_length', 'mismatches', 'gaps',
                                             'qstart', 'qend', 'sstart', 'send', 'bitscore']]
    blast_data_grouped.columns = blast_data_grouped.columns.get_level_values(0)
    grouped_df = pd.merge(blast_data_grouped, lca_data, how="left", on="query")
    grouped_df.insert(1, "blast_type", value=blast_type)
    grouped_df.insert(2, "sample", value=sample_name)
    df_to_s3(grouped_df, outfile)
    return (outfile)