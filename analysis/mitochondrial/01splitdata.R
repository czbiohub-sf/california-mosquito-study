library(dplyr)
library(reutils)
library(XML)
library(magrittr)
library(parallel)
sp_table <- read.csv("../../data/180806_project-mosquito_sample-table_pruned.csv", header=FALSE)
names(sp_table) <- c("id", "genus", "species")
cms1 <- filter(sp_table, !grepl("CMS_002", id), !grepl("[wW]ater", id)) %>%
  split(., .$genus)
cms2 <- filter(sp_table, grepl("CMS_002", id), !grepl("[wW]ater", id)) %>%
  split(., .$genus)

sapply(names(cms1), function (x) {
  if (nrow(cms1[[x]])==0) return (NULL)
  writeLines(as.character(cms1[[x]]$id), paste0("../../data/mito_cms1_", x, ".txt"))
})
sapply(names(cms2), function (x) {
  writeLines(as.character(cms2[[x]]$id), paste0("../../data/mito_cms2_", x, ".txt"))
})


# one CMS001 sample did not have genus and species info:
filter(sp_table, !grepl("[wW]ater", id), is.na(genus))
blast_results <- read.csv("../../data/CMS_058_RNA_A_S9_blast_results.csv", header=FALSE, stringsAsFactors = FALSE)
colnames(blast_results) <-  c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                              "qstart", "qend", "sstart", "send", "evalue", "bitscore")
freq_table <- table(blast_results$sseqid)
system.time(sp_names <- mclapply(names(freq_table), function (x) {
  reutils::esummary(x, db="nuccore") %>%
    content() %>%
    getNodeSet(., "//Organism") %>%
    xmlToDataFrame() %>%
    `[[`("text") %>%
    as.character()
}, mc.cores=4))

possible <- split(blast_results$sseqid, blast_results$qseqid) %>%
  lapply(table) %>%
  lapply(sort) %>%
  lapply(tail, 20) %>%
  lapply(names) %>%
  unlist() %>%
  unname() %>%
  `[`(freq_table, .) %>%
  `[`(sp_names, .) %>%
  unlist() %>%
  table() %>%
  sort()
