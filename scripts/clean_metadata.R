library(dplyr)
library(stringr)
library(magrittr)

project_dir <- ""
for (i in 1:4) {
  project_dir <- paste(rep("..", i), collapse="/")
  if (any(list.files(project_dir) == "skeeters")) {
    project_dir <- normalizePath(paste0(project_dir, "/skeeters"))
    break
  }
}

meta_df <- read.csv(paste0(project_dir, "/data/project-mosquito_sample-table.csv"),
                    stringsAsFactors=FALSE)

meta_df <- strsplit(meta_df$location, ", ") %>%
  lapply(function (x) {
    if (is.na(x[1])) return (data.frame(lat=NA, long=NA))
    if (length(x)==0) return (data.frame(lat=NA, long=NA))
    if (length(x)==2) {
      x <- as.numeric(x)
      return (data.frame(lat=x[1], long=x[2]))
    }
    if (length(x)==1) {
      x <- try(as.numeric(strsplit(x, "_")[[1]]))
      if (!is.numeric(x)) return (data.frame(lat=NA, long=NA))
      return (data.frame(lat=x[1], long=x[2]))
    }
    return (data.frame(lat=NA, long=NA))
  }) %>%
  do.call(what=rbind) %>%
  data.frame(meta_df, .)
meta_df$long[which(meta_df$long < -180)] %<>%
  sapply(function (x) {
    as.numeric(paste0(substr(x, 1, 4), ".", substr(x, 5, nchar(x))))
  })

write.csv(meta_df, paste0(project_dir, "/data/project-mosquito_sample-table_cleaned.csv"))
