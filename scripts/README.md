# Scripts

Scripts used for processing data, and running external tools.



After filtering to contigs with >2 reads, we ran aligned against the `nt` (nucleotide) and `nr` (protein) databases from NCBI using versions of BLAST. In particular, we used discontinuous megablast for `nt` and PLAST for `nr`. The precise commands were:

```BLASTDB=/mnt/data/blast blastn -task dc-megablast -db nt -evalue 1e-2 -num_threads 48 \
    -query /mnt/data/contigs/admissable_contigs.fasta \
    -outfmt "7 std staxid ssciname scomname stitle"
    -out /mnt/data/admissable_contigs_dc_megablast.m9```

```!/home/ubuntu/plastbinary_linux_20160121/plast -p plastx \
    -i /mnt/data/contigs/all.fasta \
    -d /mnt/data/blast/nr.pal \
    -o /mnt/data/plast_output.tab \
    -e 1e-2 -a 48 -max-hit-per-query 30 -outfmt 1 \
    -bargraph -verbose \
    -max-database-size 200000000```
