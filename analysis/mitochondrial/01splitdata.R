library(dplyr)
library(reutils)
library(XML)
library(magrittr)
library(parallel)
sp_table <- read.csv("../../data/sample_genus_and_species.csv", header=FALSE)
names(sp_table) <- c("id", "genus", "species")
sp_table %<>% mutate(genus=gsub(" ", "", genus), species=gsub(" ", "", species))
cms1 <- filter(sp_table, !grepl("CMS_002_[0-9]", id), !grepl("[wW]ater", id)) %>%
  split(., paste(.$genus, .$species))
cms2 <- filter(sp_table, grepl("CMS_002_[0-9]", id), !grepl("[wW]ater", id)) %>%
  split(., paste(.$genus, .$species))

sapply(names(cms1), function (x) {
  if (nrow(cms1[[x]])==0) return (NULL)
  writeLines(as.character(cms1[[x]]$id), paste0("../../data/mito_cms1_", x, ".txt"))
})
sapply(names(cms2), function (x) {
  writeLines(as.character(cms2[[x]]$id), paste0("../../data/mito_cms2_", x, ".txt"))
})



# Culex spp ---------------------------------------------------------------
culex_species <- unique(gsub("  ", " ", grep("Culex", c(names(cms1), names(cms2)), value=TRUE)))

search_results <- paste0(culex_species, "[Organism] AND COI [Title]") %>%
  lapply(reutils::esearch, db="nuccore") %>%
  lapply(reutils::esummary) %>%
  lapply(content, "parsed")

selected_culex_coi_ref <- lapply(search_results, slice, which.max(Slen)) %>%
  lapply(., `[[`, "Caption") %>%
  lapply(., efetch, db="nuccore", rettype="fasta", retmode="text") %>%
  lapply(content)
mapply(cat, selected_culex_coi_ref, file=paste0("../../data/COIref_", gsub(" ", "_", culex_species), ".fasta"))

culex_samples <- list(cms1=cms1 %>% `[`(., grep("Culex", names(.))),
                      cms2=cms2 %>% `[`(., grep("Culex", names(.))))

cat(paste("mkdir", gsub(" ", "_", culex_species)), sep="\n")

save(culex_species, culex_samples, file='culex_coi.RData')


# Aedes spp ---------------------------------------------------------------
aedes_species <- unique(gsub("  ", " ", grep("Aedes", c(names(cms1), names(cms2)), value=TRUE)))

selected_aedes_coi_ref <- paste0(aedes_species, "[Organism] AND COI [Title]") %>%
  lapply(reutils::esearch, db="nuccore") %>%
  lapply(reutils::esummary) %>%
  lapply(content, "parsed") %>%
  lapply(., slice, which.max(Slen)) %>%
  lapply(., `[[`, "Caption") %>%
  lapply(., efetch, db="nuccore", rettype="fasta", retmode="text") %>%
  lapply(content)
mapply(cat, selected_aedes_coi_ref, file=paste0("../../data/COIref_", gsub(" ", "_", aedes_species), ".fasta"))

aedes_samples <- list(cms1=cms1 %>% `[`(., grep("Aedes", names(.))),
                      cms2=cms2 %>% `[`(., grep("Aedes", names(.))))

cat(paste("mkdir", gsub(" ", "_", aedes_species)), sep="\n")

save(aedes_species, aedes_samples, file='aedes_coi.RData')


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
