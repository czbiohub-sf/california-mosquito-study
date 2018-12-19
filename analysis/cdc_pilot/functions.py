import pandas as pd
from ete3 import NCBITaxa
import os
import seaborn
import matplotlib.pyplot as plt
import math
import re
import numpy as np
from Bio import Entrez


def read_concat_csv (dirname, sample_table_fn, study_name):
    reports_fn = list(filter(lambda y: y.endswith(".csv"), os.listdir(dirname)))
    reports = pd.concat([pd.read_csv(dirname+"/"+x).assign(sample_name=x.replace(".csv", "").replace("_report", "")) for x in reports_fn], sort=False)
    sample_table = pd.read_csv(sample_table_fn)
    reports = pd.merge(reports, sample_table, how='left',
                       left_on=reports.sample_name.str.lower(),
                       right_on=sample_table.sample_name.str.lower())
    reports = reports.assign(study=study_name)
    reports.key_0 = reports.sample_name_y
    reports = reports.drop(["sample_name_x", "sample_name_y"], axis=1)
    reports = reports.rename(columns={'key_0':'sample_name'})
    selection = [(reports.category_name=="Uncategorized"), (reports.tax_level==1), reports.name.str.contains('virus'),
                 reports.name.str.contains('bacter'), reports.name.str.contains('rchae'), reports.name.str.contains("synthetic")]
    reports.loc[selection[0]&selection[1], "category_name"] = "Eukaryota"
    reports.loc[selection[0]&selection[1]&selection[2], "category_name"] = "Viruses"
    reports.loc[selection[0]&selection[1]&selection[3], "category_name"] = "Bacteria"
    reports.loc[selection[0]&selection[1]&selection[4], "category_name"] = "Archaea"
    reports.loc[selection[0]&selection[1]&selection[5], "category_name"] = "Uncategorized"
    return reports

def add_cdc_metadata (df, csv):
    cdc_pilot_metadata = pd.read_csv(csv)
    cdc_pilot_metadata = cdc_pilot_metadata.rename(
    columns={"Sample_type ":"sample_type", "Species":"mosquito_species", "Mosquito_number":"n_mosquito", "Blood Fed?":"bloodfed_label"})
    cdc_pilot_metadata = cdc_pilot_metadata[cdc_pilot_metadata.iloc[:, 0].notnull()]
    left_idx = df.sample_name.str.replace("a", "").str.split('-').apply(lambda x: int(x[1]))
    right_idx = cdc_pilot_metadata.iloc[:, 0].str.split("-").apply(lambda x: int(x[1]))
    df = pd.merge(df,
                     cdc_pilot_metadata[["sample_type", "mosquito_species", "n_mosquito", "bloodfed_label"]],
                     how="left", left_on=left_idx, right_on=right_idx)
    df.bloodfed_label = df.bloodfed_label.str.lower().replace({"y": True, "n": False})
    df.sample_type = df.sample_type.fillna("water")
    df = df.drop(["key_0"], axis=1)
    df = df.replace({'sample_type':{'aliquot of whole insect homogenate':'homogenate', 'homegenate':'homogenate'}})
    df.loc[df.sample_name.str.contains("-hom-"), "sample_type"] = "homogenate"
    return(df)

def add_extra_info (df, species_csv):
    df = df.assign(sample_type="whole")
    cms_species_calls_df = pd.read_csv(species_csv)
    cms_species_calls = cms_species_calls_df[["corrected genus", "corrected species"]].apply(lambda x: ' '.join(x) if pd.notnull(x[0]) else x[0], axis=1)
    cms_species_calls_df = cms_species_calls_df.assign(mosquito_species=cms_species_calls)[[cms_species_calls_df.columns[0], "mosquito_species"]]
    cms_species_calls_df = cms_species_calls_df.rename(columns={cms_species_calls_df.columns[0]:"sample_name"})
    cms_species_calls_df["sample_name"] = cms_species_calls_df["sample_name"].str.lower()
    df = pd.merge(df, cms_species_calls_df, how="left", on="sample_name")
    df = df.assign(n_mosquito=1).assign(bloodfed_label=True)
    df.loc[df["sample_name"].str.startswith("CMS_002_"), "bloodfed_label"] = False
    df.loc[:, "location"] = "wild-caught"
    return (df)

def find_common_name(input_row):
    Entrez.email = "lucy.li@czbiohub.org"
    ncbi = NCBITaxa()
    tax_id = input_row[0]
    name = input_row[1]
    common_name = input_row[2]
    if (name=='All taxa with neither family nor genus classification'):
        input_row[2] = name
    elif (pd.notnull(common_name)):
        input_row[2] = common_name
    else:
        if (tax_id > 0):
            search_id = tax_id
        else:
            search_term = name.split(' ')[-1]
            search_id = Entrez.read(Entrez.esearch(db="Taxonomy", term=search_term, retmode="xml"))['IdList'][0]
        search_result = list(Entrez.parse(Entrez.esummary(db="taxonomy", id=search_id, retmode="xml")))[0]
        if (7742 in ncbi.get_lineage(str(search_id))):
            input_row[3] = True
        if (len(search_result['CommonName'])>0):
            input_row[2] = search_result['CommonName']
        else:
            input_row[2] = search_result['Division']
    return (input_row)

def infer_bloodmeal (report):
    ncbi = NCBITaxa()
    conditions = [(report.tax_level==2), (report.category_name=="Eukaryota"), (report.NT_r>=50)]
    subreport = report[conditions[0]&conditions[1]&conditions[2]]
    if (len(subreport.index)==0):
        result = report[report.tax_level==2].head(1)
        result = result.assign(bloodfed_idseq=False)
    else: 
        subreport = subreport.sort_values(by="NT_r", ascending=False)
        most = list(subreport.head(1).NT_r)[0]
        subreport = subreport[subreport.NT_r >= (most * 0.2)]
        subreport = subreport.sort_values(by="NT_percentidentity")
        subreport.tax_id = subreport.tax_id.replace(
            to_replace=[133898, 39952, 122377], 
            value=[133894, 91493, 133894]
        )
        selection = subreport.tax_id.apply(lambda x: (7742 in ncbi.get_lineage(str(x))) if x > 0 else False)
        if (selection.any()):
            subreport = subreport[selection]
            result = subreport.tail(1)
            result = result.assign(bloodfed_idseq=True)
        else:
            result = subreport.tail(1)
            result = result.assign(bloodfed_idseq=False)
    result = result[["tax_id", "name", "common_name", "bloodfed_idseq", "bloodfed_label", "location"]]
    result = result.apply(find_common_name, axis=1)
    result["bloodfed_label"] = result["bloodfed_label"].astype('bool')
    return (result)

def calc_prop_reads (idseq_data):
    prop_reads_df = pd.pivot_table(idseq_data[idseq_data["tax_level"]==1] \
                               .groupby(["sample_name", "category_name"]) \
                               .agg({'NT_r': 'sum'}).groupby(level=0) \
                               .apply(lambda x: x/x.sum()), 
                               index="sample_name", columns="category_name", values="NT_r", fill_value = 0.0)
    prop_reads_df = pd.merge(prop_reads_df, idseq_data[['sample_name', 'location', 'sample_type', 'study']], how='left',
                             left_index=True, right_on='sample_name').groupby(['sample_name']).first().reset_index()
    return (prop_reads_df)


def calc_viral_diversity (idseq_data, log=False):
    result = idseq_data[idseq_data["category_name"].str.match("Viruses") & (idseq_data["tax_level"]==1)].pivot("sample_name", columns="name", values="NT_r").fillna(0)
    if (log):
        result = result.transform(lambda x: np.log10(x+0.00001))
    return (result)


def plot_heatmap (df, title='', vmin=None, vmax=None, xticklabels=True, yticklabels=True):
    seaborn.set(style="whitegrid", font_scale=2)
    ax = seaborn.heatmap(data=df, vmin=vmin or df.min().min(), vmax=vmax or df.max().max(),
                         cmap="RdPu", xticklabels=xticklabels, yticklabels=yticklabels)
    ax.set_title(title)
    return (ax)

