from ete3 import NCBITaxa
from Bio import Entrez
import sys

ncbi = NCBITaxa()
Entrez.email = "lucy.li@czbiohub.org"

acc = sys.argv[1]
taxid = str(list(Entrez.parse(Entrez.esummary(db="nucleotide", id=acc)))[0]['TaxId'])
lineage = ncbi.get_lineage(taxid)

handle = open("contig_lineages/"+acc, "w")
for x in lineage:
    handle.write(str(x)+"\n")
handle.close()