library(ggmap)
library(dplyr)
library(magrittr)

project_dir <- ""
for (i in 1:4) {
  project_dir <- paste(rep("..", i), collapse="/")
  if (any(list.files(project_dir) == "skeeters")) {
    project_dir <- normalizePath(paste0(project_dir, "/skeeters"))
    break
  }
}

metadata <- read.csv(paste0(project_dir, "/data/project-mosquito_sample-table_cleaned.csv"))

metadata_clean <- filter(metadata, !is.na(lat))

create_map <- function (input_data, zoom=5) {
  freq_table <- group_by(input_data, location) %>%
    summarise(n(), long=mean(long), lat=mean(lat)) %>%
    rename(Number=`n()`)
  long_range <- range(input_data$long)
  lat_range <- range(input_data$lat)
  cali_map <- get_map(location=bbox(SpatialPoints(input_data[, c("long", "lat")])), zoom=zoom)
  map <- ggmap(cali_map, maptype="roadmap", color="color")
  map.points <- geom_point(data=freq_table,
                           aes(x=long, y=lat, size=Number),
                           color="red", alpha=.7)
  initial_map <- map + map.points
}

cms001_map <- create_map(filter(metadata_clean, grepl("^CMS_[0-9][0-9][0-9]_RNA", metadata_clean$sample_name)),
                         zoom=10) + 
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) +
  ggtitle("CMS001")
cms002_map <- create_map(filter(metadata_clean, grepl("^CMS_002_[0-9]", metadata_clean$sample_name)),
                         zoom=6) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) +
  ggtitle("CMS002")

print(cms001_map)
print(cms002_map)


ggsave(paste0(project_dir, "/figures/initial_cms001_map.pdf"), cms001_map, width=5, height=5)
ggsave(paste0(project_dir, "/figures/initial_cms002_map.pdf"), cms002_map, width=5, height=5)