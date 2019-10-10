from glob import glob
import pandas as pd
from Bio import Entrez
from ete3 import NCBITaxa
from collections import defaultdict
import numpy as np

Entrez.email = "joshua.batson@czbiohub.org"

# load ncbi taxonomy database
ncbi = NCBITaxa()

def load_lca():
    nr_df = load_blast_output(sorted(glob('/Users/josh/src/skeeters/data/s3/'
                                          'contig_quality/*/blast_lca_nr_filtered.m9')))
    nt_df = load_blast_output(sorted(glob('/Users/josh/src/skeeters/data/s3/'
                                          'contig_quality/*/blast_lca_nt_filtered.m9')))

    # Get all taxids present in or ancestors of our dataset
    taxids = set(nr_df.taxid).union(nt_df.taxid)
    taxids = {t for taxid in taxids for t in ncbi.get_lineage(taxid)}

    taxid2name = ncbi.get_taxid_translator(taxids)
    taxid2name[-1] = 'Ambiguous'
    taxid2rank = ncbi.get_rank(taxids)
    def taxid2kingdom(taxid):
        lineage = ncbi.get_rank(ncbi.get_lineage(taxid))
        for k, v in lineage.items():
            if v == 'superkingdom':
                kingdom = taxid2name.get(k)
                return kingdom

    for df in [nt_df, nr_df]:
        df['name'] = df['taxid'].map(taxid2name)
        df['rank'] = df['taxid'].map(taxid2rank)
        df['kingdom'] = df['taxid'].map(taxid2kingdom)

    lca_df = merge_nr_nt(nr_df, nt_df)

    return lca_df

def merge_nr_nt(nr_df, nt_df):
    lca_df = nr_df.merge(nt_df, how='outer', on='contig_key', suffixes= ('_nr', '_nt'))

    lca_df.loc[lca_df['bitscore_nr'].isna(), 'bitscore_nr'] = 0
    lca_df.loc[lca_df['bitscore_nt'].isna(), 'bitscore_nt'] = 0

    selector = lca_df['bitscore_nr'] < lca_df['bitscore_nt']

    lca_df['db'] = np.where(selector, 'nt', 'nr')

    for field in ['contig', 'sample', 'name', 'rank', 'taxid', 'kingdom', 'align_length', 'identity', 'bitscore']:
        lca_df[field] = np.where(selector,
                                   lca_df[field + '_' + 'nt'],
                                   lca_df[field + '_' + 'nr'])

    lca_df['align_length'] = np.where(selector, lca_df['align_length'], lca_df['align_length'] * 3)

    return lca_df

def load_blast_output(files):
    def get_sample_from_path(path):
        return path.split('contig_quality/')[1].split('/')[0]

    dfs = []

    for file in files:
        df = pd.read_csv(file, sep='\t')
        df['sample'] = get_sample_from_path(file)
        dfs.append(df)

    df = pd.concat(dfs, axis=0)
    df = df.rename(columns={'query': 'contig'})
    df['contig_key'] = df['sample'] + '~' + df['contig']
    return df

def get_kingdom(taxid, taxid2name):
    if taxid is None:
        return None
    lineage = ncbi.get_rank(ncbi.get_lineage(taxid))
    for k, v in lineage.items():
        if v == 'superkingdom':
            kingdom = taxid2name.get(k)
            return kingdom

rank_list = ['superkingdom',
            'kingdom',
            'subkingdom',
            'superphylum', 'phylum', 'subphylum'
            'superclass', 'class', 'subclass', 'infraclass',
            'cohort', 'subcohort',
            'superorder', 'order', 'parvorder', 'suborder', 'infraorder',
            'section',
            'superfamily', 'family', 'subfamily',
            'tribe', 'subtribe',
            'genus', 'subgenus',
            'section',
            'species group', 'species subgroup', 'species', 'subspecies',
            'varietas', 'forma']

def rank_to_int(rank):
    """Gives an integer to each rank, such that i < j means j is more specific."""
    if rank in rank_list:
        return rank_list.index(rank)
    else:
        return -1

def member_key(member):
    return member['sample'] + '~' + member['contig']

def parse_cdhit_row(row):
    if '*' in row:
        index, length, name, percent_id = row.split()
        percent_id_sign, percent_id = '0', 100
        is_ref = True
    else:
        index, length, name, _, percent_id = row.split()
        is_ref = False
    length = int(length.strip(',nt'))
    name = name.strip('>').strip('.')
    sample, contig = name.split('~')
    coverage = float(contig.split('_')[-1])
    if not is_ref:
        percent_id_sign, percent_id = percent_id.strip('%').split('/')
        percent_id = float(percent_id)
    return {'contig':contig,
            'sample':sample,
            'length':length,
            'percent_id':percent_id,
            'percent_id_sign':percent_id_sign,
            'coverage':coverage,
            'is_ref':is_ref}

def load_cdhit_clusters(filename):
    clusters = defaultdict(list)
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>Cluster'):
                cluster_id = line.split()[-1]
            else:
                member = parse_cdhit_row(line)
                clusters[cluster_id].append(member)
    return clusters

def get_cluster_rep(cluster):
    for member in cluster:
        if member['is_ref']:
            return member

def get_cluster_kingdom(cluster):
    taxid_list = [member.taxid for member in cluster]
    kingdoms = [get_kingdom(taxid) for taxid in taxid_list]
    kingdoms = [k for k in kingdoms if k is not None]
    if len(kingdoms) > 0:
        return mode(kingdoms)
    else:
        return None

def merge_clusters_lca(clusters, lca_df):
    cluster_ids = [int(id) for id in clusters]
    cluster_reps = [get_cluster_rep(clusters[id]) for id in clusters]
    cluster_lengths = [rep['length'] for rep in cluster_reps]
    cluster_sizes = [len(clusters[id]) for id in clusters]
    cluster_contig_keys = [member_key(rep) for rep in cluster_reps]

    df = pd.DataFrame({'contig_key': cluster_contig_keys,
                             'cluster': cluster_ids,
                             'cluster_size': cluster_sizes,
                             'contig_length': cluster_lengths})

    df = df.merge(lca_df, on='contig_key', how='inner')

    df['align_percent'] = np.round(100*df['align_length']/df['contig_length'], 1)
    df['percent_identity'] = np.round(df['identity'], 1)

    df = df[['cluster', 'cluster_size', 'contig_length',
                         'name', 'rank', 'kingdom', 'taxid', 'db', 'align_length',
                         'align_percent', 'percent_identity', 'bitscore',
                        'contig', 'sample', 'contig_key']]
    df = df.sort_values('cluster')

    return df

def add_lca_to_clusters(clusters, lca_df):
    contig_to_lca = dict(zip(lca_df['contig_key'], lca_df['taxid']))
    for id, cluster in clusters.items():
        for member in cluster:
            taxid = contig_to_lca.get(member_key(member))
            member['taxid'] = taxid
            member['name'] = taxid2name.get()
            member['name'] = taxid2name.get()
