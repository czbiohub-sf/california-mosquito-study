# California Mosquito Study

Description TK.

## Data Manifest

Data for this project are stored in:

### The Short Read Archive
The Short Read Archive (SRA), hosted by NCBI, contains the raw sequencing reads (`fastq` files).
Accession PRJNA605178.

### Figshare
Large derived data files and data underlying manuscript figures. See [dx.doi.org/10.6084/m9.figshare.11832999](dx.doi.org/10.6084/m9.figshare.11832999).

* Assembled contig files (`contigs.zip`)
* CD-HIT-EST clusters of contigs > 500 base pairs long. (`500_contigs_cluster.clstr`)
* Read counts per contig (`contig_stats_all.tsv`)
* Least Common Ancestor (LCA) assignments for each contig (`contig_stats_lca.tsv`)
* Sample Demographics (`sample_demographics.xlsx`)
* Sequence yields per mosquito (`read_stats_per_sample.xlsx`)
* SNP distance matrices between samples (`snp_distances.tsv`)
* Underlying data: Treemap (`treemap.tsv`)
* Viral taxa detected (`virus_details.tsv`)
* Underlying data: peribunya figure (`source_data_for_perbunya_figure.zip`)
* Branch length contributions: alignments and phylogenetic trees (`branch_length.zip`)
* Branch length contributions: figures of each tree (`branch_length_trees.zip`)
* Per-mosquito nonhost reads for select taxa (`select_taxa_nonhost_reads.zip`)
* Viral segments found, grouped by genome, including sample and contig information. (`virus.json`)
* Blood meal species calls and blood meal parasites (`bloodmeal.zip`)
* Output of co-occurrence analysis. (`co-occurrence.xlsx`)
* Narnavirus analysis (`narnavirus.zip`)
* Wuhan Mosquito Virus 6 phylogenetic trees (`wmv6_trees.zip`)


(Internally, Amazon s3 was used for data storage. If you would like access to 
some intermediate files not released above, please raise and issue and we will 
add them to the public release.)

## Work in Progress

As improvements are being made to the manuscript, filenames and locations are subject to change. If you have trouble finding or using anything, please raise an issue and we will be happy to help.
