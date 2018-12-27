load("aedes_coi.RData")
samples <- aedes_samples

fastq_loc <- lapply(1:2, function (i) {
  prefix <- paste0("/mnt/data/s3/sequences/CMS00", i, "_fastq.gz/")
  lapply(samples[[i]], function (x) data.frame(x, filename=paste0(prefix, x$id), stringsAsFactors=FALSE))
}) %>% unlist(recursive=FALSE)

lapply(fastq_loc, function (x) {
  refname <- paste0("/mnt/data/analysis/coi/", x$genus[1], "_", x$species[1],"/COIref_", x$genus[1], "_", x$species[1])
  samname <- paste0("/mnt/data/analysis/coi/", x$genus[1], "_", x$species[1], "/", basename(x$filename), ".sam")
  check_command <- paste0("if [ ! -f ", x$filename, "_R1_001.fastq.gz ]; then continue; fi")
  bowtie_command <- paste0("bowtie2 --threads 1 --qupto 1000000 -x ", refname, 
         " -1 ", x$filename, "_R1_001.fastq.gz -2 ", x$filename, "_R2_001.fastq.gz -S ", samname)
  bam_command <- paste0("samtools view -Sb ", samname, " | samtools sort > ", gsub(".sam", ".bam", samname))
  vcf_command1 <- paste0("bcftools mpileup -Ou -f ", refname, ".fasta ", gsub(".sam", ".bam", samname), " | bcftools call -mv -Oz -o ", gsub(".sam", ".vcf.gz", samname))
  vcf_command2 <- paste("tabix", gsub(".sam", ".vcf.gz", samname))
  consensus_command <- paste0("cat ", refname, ".fasta | bcftools consensus ", gsub(".sam", ".vcf.gz", samname), " > ", gsub(".sam", "_consensus.fasta", samname))
  apply(data.frame(check_command, bowtie_command, bam_command, vcf_command1, vcf_command2, consensus_command), 1, paste, collapse="; ")
}) %>%
  unlist() %>%
  cat(sep="\n", file="bowtie2_commands_coi.txt")

gsub(".sam", "_.sam", temp, fixed=TRUE) %>% gsub(".bam", "_.bam", ., fixed=TRUE) %>% gsub(".vcf", "_.vcf", ., fixed=TRUE) %>% gsub("sus.f", "sus_.f", ., fixed=TRUE) %>% gsub("COIref_Culex_erythrothorax", "COIref_Culex_brami", .) %>% cat(sep="\n", file="brami.txt")