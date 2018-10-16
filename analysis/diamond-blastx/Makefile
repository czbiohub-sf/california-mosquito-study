.PHONY: all

all: blast/contigs.combined.bloodfed.curated.diamond_blastx_sensitive.m8.gz blast/contigs_combined_CMS02.withName.diamond_blastx_sensitive.m8.gz

blast/%.diamond_blastx_sensitive.m8.gz: nr.dmnd | contig_miner_output
	diamond blastx --sensitive -d nr -q contig_miner_output/$*.fasta -o $@ --compress 1 -c 1 -b 10

nr.dmnd: nr.gz
	diamond makedb --in $< -d nr

contig_miner_output:
	aws s3 cp --recursive s3://czbiohub-mosquito/contig_miner_output $@

nr.gz:
	aws s3 cp s3://czbiohub-ncbi-store/blast/db/FASTA/nr.gz $@
