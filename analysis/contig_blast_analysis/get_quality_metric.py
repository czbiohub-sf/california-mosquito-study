import pandas as pd
import os
import sys
import json
from ete3 import NCBITaxa
from Bio import Entrez


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
    return [acc], [taxid], [lca]

x = sys.argv[1]

Entrez.email = "lucy.li@czbiohub.org"

ncbi = NCBITaxa()

contig_coverage = {}
for path, subdirs, files in os.walk("contigs"):
    for name in files:
        fn = os.path.join(path, name)
        sample_name = os.path.basename(path)
        if ((".json" in name) & (sample_name==x)):
            with open (fn) as json_file:
                contig_coverage = json.load(json_file)
                for key in contig_coverage:
                    contig_coverage[key]['len'] = len(contig_coverage[key]['coverage'])


contig_metric_df = {}

metric_column_headings = ["qseqid", # name of contig
                          "sample", # name of sample
                          "pmatch_gsnap", # pcov * pident according to gsnap
                          "pmatch_rapsearch", # ^ rapsearch2
                          "qlength", # size of contig
                          "step_change", # whether there is a dramatic drop in coverage (max coverage > 20, and > 40% of contig has coverage < 3)
                          "other_blast_contigs", # are there other contigs that have blast hits that match any of the blast hits for this contig
                          "blast_hits_gsnap_acc", "blast_hits_gsnap_taxid", "blast_hits_gsnap_lca",
                          "blast_hits_rapsearch_acc", "blast_hits_rapsearch_taxid", "blast_hits_rapsearch_lca",
                          "short", 
                          "low_depth"]
fn = "contig_quality/"+x+".txt"
pd.DataFrame(columns=metric_column_headings).to_csv(fn, sep="\t", header=True, index=False, mode="w")


blast_fn = "/mnt/data/blast_results/"+x+".csv"
if (os.path.isfile(blast_fn)):
    sample_df = pd.read_csv(blast_fn, index_col=False)
else:
    sample_df = pd.DataFrame({'qseqid':[]})

contig_names = list(set(list(sample_df.qseqid.unique()) + list(contig_coverage.keys())))

for y in contig_names:
    qlength = int(y.split("_")[3])
    zipped = zip(metric_column_headings,
                 [[y], [x], [0], [0], [qlength], [False], 
                  [None], [None], [None], [None], [None], [None], 
                  [0], [False], [False]])
    metric_dict = {k:v for k,v in zipped}
    if (sample_df.qseqid.str.contains(y).any()):
        contig_df = sample_df[sample_df.qseqid==y]
        if contig_df.blast_type.str.contains("gsnap").any():
            gsnap_df = contig_df[contig_df["blast_type"]=="gsnap"]
            metric_dict["pmatch_gsnap"] = [get_pmatch_per_contig(gsnap_df)]
            metric_dict["blast_hits_gsnap_acc"], metric_dict["blast_hits_gsnap_taxid"], metric_dict["blast_hits_gsnap_lca"] = get_blast_hits_summary(gsnap_df)
        if contig_df.blast_type.str.contains("rapsearch2").any():
            rapsearch_df = contig_df[contig_df["blast_type"]=="rapsearch2"]
            metric_dict["pmatch_rapsearch"] = [get_pmatch_per_contig(rapsearch_df)]
            metric_dict["blast_hits_rapsearch_acc"], metric_dict["blast_hits_rapsearch_taxid"], metric_dict["blast_hits_rapsearch_lca"] = get_blast_hits_summary(rapsearch_df)
        other_blast_contigs = sample_df[sample_df["qseqid"]!=y]["sseqid"].isin(contig_df["sseqid"])
        if (other_blast_contigs.any()):
            metric_dict["other_blast_contigs"] = [sum(other_blast_contigs)]
    if (y in contig_coverage.keys()):
        coverage_map = contig_coverage[y]['coverage']
        if (len(coverage_map) <= 300):
            metric_dict["short"] = [True]
        if (max(coverage_map) < 3):
            metric_dict["low_depth"] = [True]
        else:
            if (max(coverage_map) > 20):
                if (sum([pos < 3 for pos in coverage_map])/len(coverage_map) > 0.4):
                    metric_dict["step_change"] = [True]
    metric = pd.DataFrame(metric_dict)[metric_column_headings]
    with open(fn, 'a') as f:
        metric.to_csv(f, sep="\t", header=False, index=False)

 