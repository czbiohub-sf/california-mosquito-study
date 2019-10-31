import s3fs
import sys
import time
import random
import os

import pandas as pd
import numpy as np
import sqlite3
import sqlalchemy


from ete3 import NCBITaxa
ncbi_taxdump_names = ["taxdump_2019-01-01.tar.gz", "taxdump_2019-06-01.tar"]
ncbi_dbs = []
if not os.path.exists(ncbi_taxdump_names[0].replace(".tar.gz", ".sqlite")):
    ncbi_dbs.append(NCBITaxa(ncbi_taxdump_names[0].replace(".tar.gz", ".sqlite"), taxdump_file=ncbi_taxdump_names[0]))
else:
    ncbi_dbs.append(NCBITaxa(ncbi_taxdump_names[0].replace(".tar.gz", ".sqlite")))
if not os.path.exists(ncbi_taxdump_names[1].replace(".tar", ".sqlite")):
    ncbi_dbs.append(NCBITaxa(ncbi_taxdump_names[1].replace(".tar", ".sqlite"), taxdump_file=ncbi_taxdump_names[1]))
else:
    ncbi_dbs.append(NCBITaxa(ncbi_taxdump_names[1].replace(".tar", ".sqlite")))
ncbi_new = NCBITaxa()
#ncbi_new.update_taxonomy_database()
ncbi_dbs.append(ncbi_new)

from lca_functions import *



def get_name (query, ncbi_dbs, sci_name=True):
    for ncbi_db in ncbi_dbs:
        if (sci_name):
            result = ncbi_db.get_taxid_translator([query])
        else:
            result = ncbi_db.get_common_names([query])
        if (len(result)>0):
            break
    if (len(result)==0):
        result = {query:np.nan}
    return (result[query])


def add_taxid(s3_path, sqlite_path):
    print ("starting to process: "+s3_path)
    data_frame = pd.read_csv(s3_path, header=0, sep="\t")
    gi_queries = data_frame["subject"].str.split("|").apply(lambda x: x[1])
    print ("processing "+s3_path+"; reading in sql database.")
    conn = sqlite3.connect(sqlite_path, timeout=300)
    taxids = pd.read_sql_query("select gi, taxid from gi2taxid where gi in ("+','.join(gi_queries)+");", conn).set_index(keys="gi").to_dict()["taxid"]
    conn.close()
    data_frame = data_frame.assign(taxid=[taxids[int(x)] for x in gi_queries])
    data_frame = data_frame.assign(sci_name=data_frame["taxid"].apply(get_name, ncbi_dbs=ncbi_dbs, sci_name=True))
    data_frame = data_frame.assign(common_name=data_frame["taxid"].apply(get_name, ncbi_dbs=ncbi_dbs, sci_name=False))
    df_to_s3(data_frame, s3_path, header=True)
    print ("uploaded dataframe to: "+s3_path)
    return (data_frame)




if (len(sys.argv)>1):
    blast_nr_paths = [sys.argv[1]]
else:
    s3 = boto3.resource('s3')
    client = boto3.client('s3')
    bucket_name = "czbiohub-mosquito"
    bucket = s3.Bucket(bucket_name)
    contig_folders = [x["Prefix"] for x in client.list_objects(Bucket=bucket_name, Prefix="contigs/", Delimiter="/")["CommonPrefixes"]]
    blast_nr_paths = ["s3://"+bucket_name+"/"+x["Prefix"]+"blast_nr.m9" for x in client.list_objects(Bucket=bucket_name, Prefix="contigs/", Delimiter="/")["CommonPrefixes"] if "Mos" not in x["Prefix"]]


[add_taxid(x, sqlite_path="../../data/gi2taxid.db") for x in blast_nr_paths]



