##########
# Perform lca anlaysis on nt and nr hits
##########

import argparse, sys, time
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
parser.add_argument('--verbose', default="False", help="If true, print messages after each step of the script.")
args = parser.parse_args()

verbose = args.verbose=="True"

# import functions and variables from lca_functions.py
exec(open("lca_functions.py").read())

start_time = time.time()

# if contigs have 2 or fewer reads, then exclude them from LCA analysis
read_counts = load_json(args.read_count_path, colnames=["query", "read_count"])
filtered_contigs_by_read_count = read_counts[read_counts["read_count"]>2]

print_to_stdout("Read counts have been loaded: "+args.read_count_path, start_time, verbose)

if (len(filtered_contigs_by_read_count.index)==0):
    print_to_stdout(args.fpath+": no contig had more than 2 reads.", start_time, verbose)
    exit()

if args.blast_type =="nt":
    db = "nucleotide"
    col_names = ["query", "subject", "identity", "align_length", "mismatches", 
        "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "taxid", 
        "sci_name", "common_name", "subject_title", "qcov", "hsp_count"]
elif args.blast_type=="nr":
    db = "protein"
    col_names = ["query", "subject", "identity", "align_length", "mismatches", 
        "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "taxid", 
        "sci_name", "common_name", "subject_title", "qcov", "hsp_count"]

# # Set up aws s3 access
# blast_file = args.fpath
        
# # read in objects
# blast_results = parse_blast_file(blast_file, sep="\t", comment="#", blast_type=args.blast_type, col_names="auto")
    
# obtain a single HSP for each query-subject pairing and read in blast results as a pandas data frame
blast_results = get_single_hsp(args.fpath, col_names) 
print_to_stdout("Loaded blast file: "+args.fpath, start_time, verbose)


# data frame: whether or not each contig was included or excluded from blast analysis, and reason for exclusion
excluded_contigs = pd.DataFrame({"query":blast_results["query"].unique()})
excluded_contigs = excluded_contigs.assign(contig_length=excluded_contigs["query"].str.split("_").apply(lambda x: int(x[3])))
excluded_contigs = pd.merge(excluded_contigs, read_counts, how="left", on="query").fillna(0)
excluded_contigs = excluded_contigs.assign(low_read_count=~excluded_contigs["query"].isin(filtered_contigs_by_read_count["query"]))

# find missing taxids
blast_results["taxid"] = blast_results["taxid"].replace(to_replace=0, value=np.nan) # some synthetic constructs have taxid 0
if (blast_results["taxid"].isnull().any()):
    subjects_to_search = list(blast_results[blast_results["taxid"].isnull()]["subject"].unique())
    print_to_stdout(str(blast_results["taxid"].isnull().sum())+" blast hits corresponding to "+str(len(subjects_to_search))+" accession numbers have taxid 'NA'. Trying to find the taxid for these hits on NCBI.", start_time, verbose)
    subjects_taxids = [find_missing_taxid(x, db=db) for x in subjects_to_search]
    subjects_taxid_dict = dict(zip(subjects_to_search, subjects_taxids))
    blast_results.loc[blast_results["taxid"].isnull(), ["taxid"]] = blast_results[blast_results["taxid"].isnull()]["subject"].apply(lambda x: subjects_taxid_dict[x])
    blast_results = blast_results[~blast_results["taxid"].isnull()]

excluded_contigs = excluded_contigs.assign(taxid_na=~excluded_contigs["query"].isin(blast_results["query"]))
print_to_stdout(str(excluded_contigs["taxid_na"].sum())+" contigs were excluded because none of the subject taxids could be found.", start_time, verbose)        

# exclude contigs with hits to mosquito
all_hits_queries = blast_results["query"].unique()
blast_results = blast_results.groupby(["query"], as_index=False).apply(
        filter_by_taxid, db=db, taxid=ncbi.get_name_translator(["Hexapoda"])["Hexapoda"][0]
)
hexa_contigs = blast_results["query"][~blast_results["query"].isin(all_hits_queries)]
excluded_contigs = excluded_contigs.assign(hexapoda=excluded_contigs["query"].isin(hexa_contigs))

if (len(hexa_contigs.index)>0):
    print_to_stdout(str(len(hexa_contigs.index))+" contigs were likely hexapoda.", start_time, verbose)

# remove contigs with too few reads from blast_results
blast_results = blast_results[~blast_results["query"].isin(excluded_contigs["query"][excluded_contigs["low_read_count"]])]

print_to_stdout (str(len(blast_results["query"].unique())) + " out of " + str(len(excluded_contigs["query"])) + " contigs passed the filters and will be processed for LCA analysis.", start_time, verbose)

# lca analysis
filtered_blast_results = blast_results.groupby(["query"], as_index=False).apply(
    select_taxids_for_lca, db=db,
    return_taxid_only=False,
    ident_cutoff=float(args.ident_cutoff),
    align_len_cutoff=float(args.align_len_cutoff),
    bitscore_cutoff=float(args.bitscore_cutoff),
    read_counts=filtered_contigs_by_read_count
)

print_to_stdout("BLAST results have been filtered: "+args.filtered_blast_path, start_time, verbose)

if (args.filtered_blast_path is not None):
    df_to_s3(filtered_blast_results, args.filtered_blast_path)
    print_to_stdout("BLAST results have been uploaded to : "+args.filtered_blast_path, start_time, verbose)


lca_results = filtered_blast_results.groupby(["query"]).apply(get_lca)

print_to_stdout("LCA analysis complete: "+args.outpath, start_time, verbose)

# write results to file
if (args.outpath.startswith("s3://")):
    df_to_s3(lca_results, args.outpath)
    df_to_s3(excluded_contigs, args.excluded_contigs_path)
else:
    lca_results.to_csv(args.outpath, sep="\t", index=False)
    excluded_contigs.to_csv(args.excluded_contigs_path, sep="\t", index=False)

    
    
print_to_stdout("LCA analysis saved to file: "+args.outpath, start_time, verbose)
 

