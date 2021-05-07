# California Mosquito Study

Mosquitoes are a disease vector with a complex ecology involving interactions between transmissible pathogens, endogenous microbiota, and human and animal blood meal sources. Unbiased metatranscriptomic sequencing of individual mosquitoes offers a straightforward and rapid way to characterize these dynamics. Here, we profile 148 diverse wild-caught mosquitoes collected in California, detecting sequences from eukaryotes, prokaryotes, and over 70 known and novel viral species. Because we sequenced singletons, it was possible to compute the prevalence of each microbe and recognize a high frequency of viral co-infection. By analyzing the pattern of co-occurrence of sequences across samples, we associated "dark matter" sequences with recognizable viral polymerases, and animal pathogens with specific blood meal sources. We were also able to detect frequent genetic reassortment events in a highly prevalent quaranjavirus undergoing a recent intercontinental sweep. In the context of an emerging disease, where knowledge about vectors, pathogens, and reservoirs is lacking, the approaches described here can provide actionable information for public health surveillance and intervention decisions.

# Citation 
Joshua Batson, Gytis Dudas, Eric Haas-Stapleton, Amy L Kistler, Lucy M Li, Phoenix Logan, Kalani Ratnasiri, and Hanna Retallack. Single mosquito metatranscriptomics identifies vectors, emerging pathogens and reservoirs in one assay. eLife 2021;10:e68353 DOI: 10.7554/eLife.68353


## Code

All code used to generate the figures and analyses in this manuscript are available in this repo, across the `analysis` and `scripts` repos.

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
