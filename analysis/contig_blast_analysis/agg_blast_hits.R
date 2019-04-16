library(dplyr)
library(magrittr)
library(readr)

agg_blast_hits <- function (query_x) {
  output <- query_x[1, ] %>% select(qseqid, sseqid, pident, length, blast_type, sample, pmatch, qlength, taxid, lineage)
  if (nrow(query_x)==1) return (output)
  pmatch <- apply(query_x, 1, function (row_x) {
    output <- rep(0, row_x["qlength"])
    output[row_x["qstart"]:row_x["qend"]] <- as.numeric(row_x["pident"])/100
    return (output)
  }) %>% apply(1, max)
  output$pident <- round(mean(pmatch[pmatch>0])*100, 1)
  output$length <- sum(pmatch>0)
  output$pmatch <- mean(pmatch)
  output
}

args = commandArgs(trailingOnly=TRUE)

if (file.exists(args[1])) {
  raw_file <- read_tsv(args[1])
  filtered_contigs <- raw_file %>%
    group_by(sample, qseqid, sseqid, blast_type) %>%
    do(agg_blast_hits(.)) %>%
    group_by(sample, qseqid, taxid, blast_type) %>%
    filter(pmatch==max(pmatch))
  write.table(filtered_contigs, gsub("blast_results", "blast_agg_results", args[1]), sep="\t", row.names=FALSE, quote=FALSE)
}