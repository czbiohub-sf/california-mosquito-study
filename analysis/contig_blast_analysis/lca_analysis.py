##########
# Perform lca anlaysis on nt and nr hits
##########

import argparse, sys
parser = argparse.ArgumentParser()
parser.add_argument('--fpath', help="Path to file (local or on S3). For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--outpath', help="Path to output file (local or on S3).  For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--blast_type', help="Either 'nt' or 'nr'")
parser.add_argument('--filtered_blast_path', default=None, help="Path to filtered blast results file (local or on S3). For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--excluded_contigs_path', default=None, help="Path to the file containing excluded contigs and the reason for exclusion (local or on S3). For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--read_count_path', default=None, help="Path to read counts file (local or on S3). For S3 files, the path should be in the form s3://bucket_name/path/to/file.")
parser.add_argument('--ident_cutoff', default=0, help="Remove blast hits if their ident < ident_cutoff*max(ident) for a given contig.")
parser.add_argument('--align_len_cutoff', default=0, help="Remove blast hits if their alignment length < align_len_cutoff*max(align) for a given contig.")
parser.add_argument('--bitscore_cutoff', default=0, help="Remove blast hits if their bitscore < bitscore_cutoff*max(bitscore) for a given contig.")
args = parser.parse_args()

# import functions and variables from lca_functions.py
exec(open("lca_functions.py").read())

read_counts = pd.read_csv(args.read_count_path, sep="\t", header=None, names=["query", "read_count"])
filtered_contigs_by_read_count = read_counts[read_counts["read_count"]>2]

if (len(filtered_contigs_by_read_count.index)==0):
    print(args.fpath+": no contig had more than 2 reads.")
    exit()

if args.blast_type =="nt":
    db = "nucleotide"
elif args.blast_type=="nr":
    db = "protein"

# Set up aws s3 access
blast_file = args.fpath
if (args.fpath.startswith("s3://")):
    blast_file = download_s3_file(args.fpath)    

# read in objects
blast_results = parse_blast_file(blast_file, sep="\t", comment="#", blast_type=args.blast_type, col_names="auto")

# lca analysis
filtered_blast_results, excluded_contigs = blast_results.groupby(["query"], as_index=False).apply(
    select_taxids_for_lca, db=db,
    return_taxid_only=False,
    ident_cutoff=float(args.ident_cutoff),
    align_len_cutoff=float(args.align_len_cutoff),
    bitscore_cutoff=float(args.bitscore_cutoff),
    read_counts=filtered_contigs_by_read_count
)

excluded_contigs = pd.concat([excluded_contigs, read_counts[read_counts["read_count"]<=2].assign(reason="low_read_count")])

if (args.filtered_blast_path is not None):
    df_to_s3(filtered_blast_results, args.filtered_blast_path)


lca_results = filtered_blast_results.groupby(["query"]).apply(get_lca)

# write results to file
if (args.outpath.startswith("s3://")):
    df_to_s3(lca_results, args.outpath)
    df_to_s3(excluded_contigs, args.excluded_contigs_path)
else:
    lca_results.to_csv(args.outpath, sep="\t", index=False)
    excluded_contigs.to_csv(args.excluded_contigs_path, sep="\t", index=False)
