



calc_eucl_dist <- function (v1, v2) {
  sqrt(sum((v2-v1)^2))
}

capitalize_first <- function (x) {
  A <- toupper(substr(x, 1, 1))
  B <- tolower(substr(x, 2, nchar(x)))
  paste0(A, B)
}



get_ancestor <- function (input_df, taxonomic_rank=NULL) {
  # Return the taxonomic ID 
  if (is.null(taxonomic_rank)) return (input_df)
  output_df <- input_df
  updated_df <- lapply(1:nrow(input_df), function (i) {
    output <- select(input_df, taxid, taxid_rank, sci_name, lineage, lineage_rank, lineage_sci_name)[i, ]
    ranks <- fromJSON(output$lineage_rank)
    pos <- NA
    if (taxonomic_rank %in% ranks) {
      # check to see if the taxonomic_rank exists in the lineage leading to a taxid
      pos <- which(ranks==taxonomic_rank)
    } else {
      # else check for the closest rank to the taxonomic_rank. e.g. if taxonomic_rank=="family", then 
      # the lines below check to see if "genus" exists in the lineage, and if it does exist then the
      # pointer (pos) moves up one taxonomic level from "genus" 
      remaining_ranks <- taxize::rank_ref$ranks[-1:-which(taxonomic_rank==taxize::rank_ref$ranks)]
      if (any(ranks %in% remaining_ranks)) {
        pos <- which(ranks %in% remaining_ranks)[1]
        if (pos > 1) pos <- pos - 1
      } else {
        pos <- length(fromJSON(output$lineage))
      }
    }
    while (grepl("unclassified", fromJSON(output$lineage_sci_name)[pos])) {
      pos <- pos - 1
    }
    output[, c("lineage", "lineage_rank", "lineage_sci_name")] %<>% 
      lapply(fromJSON) %>% lapply(`[`, 1:pos) %>% lapply(toJSON) %>%
      data.frame(stringsAsFactors=FALSE)
    output[, c("taxid", "taxid_rank", "sci_name")] <- 
      output[, c("lineage", "lineage_rank", "lineage_sci_name")] %>%
      lapply(fromJSON) %>% lapply(`[`, pos) %>%
      data.frame(stringsAsFactors=FALSE)
    output$taxid %<>% as.numeric()
    return (output)
  }) %>%
    do.call(what=rbind)
  output_df[, names(updated_df)] <- updated_df
  return (output_df)
}




