import pandas as pd
import json as json

from collections import defaultdict, namedtuple

Member = namedtuple('Member', ['contig', 'length', 'percent_id', 'percent_id_sign', 'sample', 'coverage'])

def parse_cdhit_row(row):
    if '*' in row:
        index, length, name, percent_id = row.split()
        percent_id_sign, percent_id = '0', 100
    else:
        index, length, name, _, percent_id = row.split()
    length = int(length.strip(',nt'))
    name = name.strip('>').strip('.')
    sample, contig = name.split('~')
    coverage = float(contig.split('_')[-1])

    if percent_id != 100:
        percent_id_sign, percent_id = percent_id.strip('%').split('/')
        percent_id = float(percent_id)
    return Member(contig=contig, sample=sample, length=length,
                  percent_id=percent_id, percent_id_sign=percent_id_sign, coverage=coverage)

clusters = defaultdict(list)
with open('/Users/josh/src/skeeters/data/500_contigs_cluster.clstr', 'r') as file:
    for line in file:
        if line.startswith('>Cluster'):
            cluster_id = line.split()[-1]
        else:
            member = parse_cdhit_row(line)
            clusters[cluster_id].append(member)

def cluster_to_contigs(clusters, out_file):

    cluster_contig_df = [{'sample': member.sample, 'contig': member.contig,
                          'cluster': cluster_id}
                         for cluster_id in clusters for member in clusters[cluster_id]]
    cluster_contig_df = pd.DataFrame(cluster_contig_df)
    cluster_contig_df.to_csv(out_file, index=False)

def cluster_to_virus(json_file, out_file):

    with open(json_file) as fp:
        virus_json = json.load(fp)

    cluster_to_virus = []
    for poly_group in virus_json:
        virus = virus_json[poly_group]
        taxid = virus['submission_taxid']
        clusters = []
        for segment in virus['segments']:
            for cluster in virus['segments'][segment]['clusters']:
                cluster_to_virus.append({'cluster': cluster,
                                     'submission_taxid': virus['submission_taxid'],
                                     'poly_group': poly_group,
                                     'name': virus['name'],
                                     'provisional_name': (
        virus['name'] if 'provisional_name' not in virus else virus['provisional_name']),
                                     'segment': segment})
    pd.DataFrame(cluster_to_virus).to_csv(out_file, index=None)

cluster_to_contigs(clusters, '/Users/josh/src/skeeters/data/darkmatter/cluster_contig.csv')
cluster_to_virus('/Users/josh/src/skeeters/data/darkmatter/virus.json',
    '/Users/josh/src/skeeters/data/darkmatter/clusters_to_virus.csv')
