library(ggmap)
library(dplyr)
library(magrittr)
metadata <- read.csv("../data/project-mosquito_sample-table_cleaned.csv")

metadata_clean <- filter(metadata, !is.na(lat))
freq_table <- group_by(metadata_clean, location) %>%
  summarise(n(), long=mean(long), lat=mean(lat)) %>%
  rename(Number=`n()`)

cali_map <- get_map(location=c(lon=mean(range(metadata_clean$long, na.rm=TRUE)),
                               lat=mean(range(metadata_clean$lat, na.rm=TRUE))))

map <- ggmap(cali_map, zoom=6, maptype="roadmap", color="color")
map.points <- geom_point(data=freq_table,
                             aes(x=long, y=lat, size=Number),
                             color="red")
initial_map <- map + map.points

ggsave("../../figures/initial_map.pdf", initial_map, width=5, height=5)
