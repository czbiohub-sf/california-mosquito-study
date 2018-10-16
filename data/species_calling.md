Based on genetic evidence, we have reclassified the species of a few of the samples.

There were two sources of evidence:

1. kmer similarity of RNA seq (with abundances)
2. kmer similarity of RNA seq (without abundances)
3. quality of read alignment to the COI gene of each species

## CMS_045_RNA_A_S2
Culex	tarsalis to	Culex	erythrothorax based on k51 and COI

## CMS_026_RNA_A_S18
Culex	tarsalis to	Culex	indicens based on k51 and k31
#20181015 AK noticed typo in genus, updated genus from Culex to Culiseta

## CMS_058_RNA_A_S9
NaN to Culex erythrothorax based on k51 and COI

## CMS_025_RNA_A_S7
Culex tarsalis to Culex particeps based on k21,31,51.
#20181015 AK noticed typo in genus, updated genus from Culex to Culiseta


A few samples appear to be ambiguous:

## CMS_044_RNA_A_S22
The sample `CMS_044_RNA_A_S22` is an extreme outlier for no-abundance (unique kmers). Its nearest neighbor is an order of magnitude less close than any other sample’s nearest neighbor.
#20181015 AK update - this is a failed sequencing run for this sample, which may explain outlier. Replaced with successfully run CMS_044_RNA_A_S25

## CMS_040_RNA_A_S21
This sample sometimes clusters with erythrothorax (including for all no-abun). Perhaps inbreeding?

## Culex pipiens/quinquefasciatus
These species are quite closely related genetically, and may actually
have some interbreeding. They are distinguished completely by unique kmers. They get mixed up when using abundances, where different values of k give different
subgroupings, and certain samples appear to be outliers
in terms of similarity. These include

* CMS_002_27b_Rb_S153_L004
* CMS_029_RNA_A_S18
* CMS_002_10a_Rb_S119_L004
