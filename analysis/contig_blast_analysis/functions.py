import pandas as pd
import os
import sys
import json
from ete3 import NCBITaxa
from Bio import Entrez
import time


def get_pmatch_per_contig (df):
    if (len(df.index)==1):
        return (df["pmatch"].iloc[0])
    covered_bases = []
    for row_i in range(len(df.index)):
        start_end = df.iloc[row_i, :][["qend", "qstart"]].sort_values()
        covered_bases += list(range(start_end[0], start_end[1]))
    return(float(len(set(covered_bases))/df.head(1)["qlength"]))


def get_taxid (acc, dir_name):
    fn = dir_name+"/"+acc+".txt"
    if (os.path.isfile(fn)):
        taxid = pd.read_csv(fn, sep=" ", header=None).iloc[0][1]
    else:
        db_name = "nucleotide"
        if ('protein' in dir_name):
            db_name = "protein"
        taxid = int(Entrez.read(Entrez.esummary(id=acc, db=db_name, rettype="gb", retmode="text"))[0]['TaxId'])
        output_string = acc+" "+str(taxid)
        with open (fn, 'w') as f:
            f.write("%s" % output_string)
    return (taxid)

def get_lca (lineages):
    found_lca = False
    for node in reversed(lineages[0]):
        if (all([node in lineage_x for lineage_x in lineages[1:]])):
            found_lca = True
            break
    if (found_lca):
        return (node)
    else:
        return ("no_lca_found")

def get_blast_hits_summary (df):
    blast_type = "protein_matches_unique"
    if (df["blast_type"].iloc[0]=="gsnap"):
        blast_type = "nucleotide_matches_unique"
    if (df["sseqid"].nunique()==1):
        acc = [df["sseqid"].iloc[0]]
        taxid = [get_taxid(acc[0], blast_type)]
        lca = taxid[0]
    else:
        pmatch_table = df.groupby("sseqid").apply(get_pmatch_per_contig)
        selected_acc = pmatch_table[pmatch_table > max(pmatch_table)*0.9].index
        taxid = [get_taxid (acc, blast_type) for acc in selected_acc]
        if (len(set(taxid))==1):
            acc = list(selected_acc)
            lca = taxid[0]
        else:
            lineages = [ncbi.get_lineage(tax) for tax in taxid]
            acc = list(selected_acc)
            lca = get_lca(lineages)
    return [acc], [taxid], lca
