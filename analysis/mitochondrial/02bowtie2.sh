for x in `ls COIref*.fasta`; do bowtie2-build $x ${x%.*}; done

parallel < bowtie2_commands.txt 

for x in `ls */*consensus.fasta`; do var='>'; y=$(basename $x); y=${y%.*}; var=$var$y; sed -i "1s/.*/$var/" $x; done
cat */*consensus.fasta > consensus.fasta

muscle -maxiters 1 -diags -in consensus.fasta -out consensus_aligned.fasta

modeltest-ng -c 4 -d nt -p 2 -i consensus_aligned.fasta

raxml-ng --msa consensus_outgroup_aligned.fasta --model TPM1uf+I+G4 --outgroup MG384714_Aedes_aegypti