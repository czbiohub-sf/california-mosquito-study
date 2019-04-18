##########
# Perform lca anlaysis on nt and nr hits
##########

import argparse, sys, os, tempfile
parser = argparse.ArgumentParser()
parser.add_argument('--fpath', help="Path to file (local or on S3). For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--outpath', help="Path to output file (local or on S3).  For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--blast_type', help="Either 'nt' or 'nr'")
parser.add_argument('--filtered_blast_path', default=None, help="Path to filtered blast results file (local or on S3). For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--ident_cutoff', default=0, help="Remove blast hits if their ident < ident_cutoff*max(ident) for a given contig.")
parser.add_argument('--align_len_cutoff', default=0, help="Remove blast hits if their alignment length < align_len_cutoff*max(align) for a given contig.")
parser.add_argument('--bitscore_cutoff', default=0, help="Remove blast hits if their bitscore < bitscore_cutoff*max(bitscore) for a given contig.")
args = parser.parse_args()

# import functions and variables from lca_functions.py
exec(open("lca_functions.py").read())

def split_s3_path (s3_path):
    s3_split = os.path.normpath(s3_path).split(os.sep)
    bucket_name = s3_split[1]
    s3_path = '/'.join(s3_split[2:])
    return bucket_name, s3_path

# Set up aws s3 access
if (args.fpath.startswith("s3://")):
    import boto3
    s3 = boto3.resource('s3')
    client = boto3.client('s3')
    s3_bucket_name, s3_path = split_s3_path(args.fpath)
    blast_file = client.get_object(Bucket=s3_bucket_name, Key=s3_path)['Body']
else:
    blast_file = args.fpath

# read in objects
blast_results = parse_blast_file(blast_file, sep="\t", comment="#", blast_type=args.blast_type, col_names="auto")

# lca analysis
filtered_blast_results = blast_results.groupby(["query"]).apply(select_taxids_for_lca, return_taxid_only=False, ident_cutoff=args.ident_cutoff, align_len_cutoff=args.align_len_cutoff, bitscore_cutoff=args.bitscore_cutoff)

if (args.filtered_blast_path is not None):
    s3_bucket_name, s3_path = split_s3_path(args.filtered_blast_path)
    s3.Bucket(s3_bucket_name).put_object(Key=s3_path, Body=data)

lca_results = filtered_blast_results.groupby(["query"]).apply(get_lca)

# write results to file
if (args.outpath.startswith("s3://")):
    out_file = tempfile.NamedTemporaryFile()
    lca_results.to_csv(out_file.name, sep="\t", index=False)
    s3_bucket_name, s3_path = split_s3_path(args.outpath)
    data = open(out_file.name, "rb")
    s3.Bucket(s3_bucket_name).put_object(Key=s3_path, Body=data)
else:
    lca_results.to_csv(args.outpath, sep="\t", index=False)
