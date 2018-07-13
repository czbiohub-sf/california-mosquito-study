library(dplyr)

meta_df <- read.csv("~/Downloads/project-mosquito_sample-table.csv")

library(stringr)

meta_df$location %>%
  str_match("(\\d+\\.\\d+)_-(\\d+\\.\\d+)") ->
  locations_str_matched

meta_df$lat <- as.numeric(locations_str_matched[,2])
meta_df$long <- -as.numeric(locations_str_matched[,3])
