from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Extract contig from sample')
parser.add_argument('sample', help='sample name from idseq')
parser.add_argument('node', help='name of contig (node)')

args = parser.parse_args()

def add_sample(record, sample):
    record.id = sample + '~' + record.id
    record.name = record.id
    record.description = ''
    return record

def print_contig(sample, node):
    filename = '../data/s3/contigs/' + sample + '/contigs.fasta'
    records = [add_sample(record, sample) for record in SeqIO.parse(filename, 'fasta') if record.id == node]
    if len(records) ==1:
        print(records[0].format("fasta"))
    else:
        print("Record not found")
        
print_contig(args.sample, args.node)
