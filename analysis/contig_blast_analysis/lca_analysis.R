library(dplyr)
library(magrittr)
library(readr)
library(rjson)

get_lca <- function (x) {
  lca <- x[[1]]
  for (i in seq(x)[-1]) {
    lca <- intersect(lca, x[[i]])
  }
  as.character(tail(lca, 1))
}

summarize_hits <- function (query_x) {
    if (max(query_x$pmatch)==1) {
        threshold <- 1
    } else if (max(query_x$pmatch)>=0.9) {
        threshold <- max(query_x$pmatch)-0.05
    } else if (max(query_x$pmatch)>=0.8) {
        threshold <- max(query_x$pmatch)-0.1
    } else {
        threshold <- max(query_x$pmatch)-0.2
    }
    filtered_table <- filter(query_x, pmatch >= threshold)
    output <- filtered_table[1, ]
    for (colname in c("sseqid", "pident", "length", "taxid", "lineage")) {
        output[[colname]] <- rjson::toJSON(as.list(filtered_table[[colname]])) %>% gsub('"', "", .)
    }
    lca_list <- tryCatch({filtered_table$lineage %>% lapply(rjson::fromJSON)}, error=function (e) browser())
    if (nrow(filtered_table)==1) {
        output$lca <- as.character(output$taxid)
        output$lca_lineage <- lca_list[[1]] %>% unlist() %>% rjson::toJSON() %>% gsub('"', "", .)
        return (output)
    }
    output$lca <- get_lca(lca_list)
    output$lca_lineage <- lca_list[[1]] %>% `[`(., 1:which(output$lca==.)) %>% unlist() %>% rjson::toJSON() %>% gsub('"', "", .)
    return (output)
}

args = commandArgs(trailingOnly=TRUE)

if (file.exists(args[1])) {
  raw_file <- read_tsv(args[1])
  filtered_contigs <- raw_file %>%
    group_by(sample, qseqid, blast_type)
  filtered_contigs <- filtered_contigs %>%
    do(summarize_hits(.))
  write.table(filtered_contigs, gsub("blast_agg_results", "blast_lca_results", args[1]), sep="\t", row.names=FALSE, quote=FALSE)
}