
```{r get-shapefiles, cache=TRUE}
list_of_counties <- c("Alameda", "Placer", "West Valley", "Coachella Valley", "San Diego")
# https://data.ca.gov/dataset/ca-geographic-boundaries
shapefile_names <- c(paste0("CA_", c("State", "Counties", "Places")), "WestValley")
ca_shp <- file.path(data_dir, "shapefiles", shapefile_names) %>%
  lapply(., readOGR, stringsAsFactors=FALSE) %>%
  lapply(., spTransform, CRS("+proj=longlat +datum=WGS84")) %>%
  setNames(., gsub("CA_", "", shapefile_names))
coachellavalley_names <- ca_shp$Places@data %>% 
  select(INTPTLON, INTPTLAT) %>%
  apply(1, as.numeric) %>% 
  apply(2, calc_eucl_dist, c(-116.0825, 33.51680)) %>%
  order(.) %>% 
  slice(ca_shp$Places@data, .) %>%
  slice(., 1:(which(NAME=="Coachella")-1)) %>%
  `[[`("NAME")
ca_shp$Areas <- ca_shp$Counties %>%
  `[`(., .$NAME %in% c("Alameda", "Placer", "San Diego"), )
ca_shp$Areas@plotOrder <- 1:5
ca_shp$Areas@data %<>% rbind(., .[1:2, ])
places_data <- slice(ca_shp$Places@data, match(c("Ontario", "Coachella"), NAME)) %>%
  mutate(NAME=c("West Valley", "Coachella Valley"))
common_columns <- names(ca_shp$Areas@data) %>% `[`(., . %in% names(places_data))
ca_shp$Areas@data[4:5, common_columns] <- places_data[, common_columns]
ca_shp$Areas@data[4:5, c("COUNTYFP", "COUNTYNS")] <- places_data[, c("PLACEFP", "PLACENS")]
ca_shp$Areas@data %<>% 
  mutate(colours=c("#6b2500", "#de2d6a", "#565fd5", "#ff6f53", "#22b14e"),
         order=match(NAME, list_of_counties))
gpclibPermit() # https://stackoverflow.com/questions/30790036/error-istruegpclibpermitstatus-is-not-true
ca_shp$Areas@polygons %<>% 
  c(list(unionSpatialPolygons(ca_shp$WestValley, IDs=rep(1, length(ca_shp$WestValley@polygons))),
         unionSpatialPolygons(ca_shp$Places[ca_shp$Places$NAME %in% coachellavalley_names, ], rep(1, length(coachellavalley_names)))) %>%
      lapply(., spTransform, CRS("+proj=longlat +datum=WGS84")) %>%
      lapply(., slot, "polygons") %>% 
      lapply(., `[[`, 1)
  )
for (i in seq(ca_shp$Areas@polygons)) {
  slot(slot(ca_shp$Areas, "polygons")[[i]], "ID") <- as.character(i)
}
ca_shp$Areas@data %<>% mutate(., id=sapply(ca_shp$Areas@polygons, function (x) x@ID))
zoom_options <- get_zoom_options (long=ca_shp$State@bbox[1, ], lat=ca_shp$State@bbox[2, ])
map <- openmap(c(zoom_options$lat_range[2], zoom_options$lon_range[1]),
               c(zoom_options$lat_range[1], zoom_options$lon_range[2]), 
               zoom = zoom_options$zoom, type = "osm", mergeTiles = TRUE) %>%
  openproj(CRS("+proj=longlat +datum=WGS84"))
base_cali_map <- autoplot(map) +
  geom_polygon(data=ca_shp$State, aes(x=long, y=lat, group=group), colour = "black", fill = NA) +
  geom_polygon(data=ca_shp$Areas, aes(x=long, y=lat, group=group, fill=factor(id, levels=with(ca_shp$Areas@data, id[order(order)]))), colour = "black", alpha=.75) +
  scale_fill_manual("Collection site", labels=with(ca_shp$Areas@data, NAME[order(order)]), values=with(ca_shp$Areas@data, colours[order(order)])) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="top", legend.direction="vertical",
        plot.margin = unit(c(1,1,1,1), "mm"))
```

```{r infer_place_names, cache=TRUE}
# cms_all_data <- filter(metadata, grepl("^CMS", sample_name)) %>%
#     mutate(OldIDseqName=substr(sample_name, 1, 11) %>% sapply(function (x) {
#     merged_metadata$OldIDseqName[c(grep(x, merged_metadata$czbiohub.mosquito_sequences_id),
#       grep(x, merged_metadata$OldIDseqName))[1]]
# })%>%factor()) %>%
#   left_join(., select(merged_metadata, contains("IDseqName")), by="OldIDseqName") %>%
#   rename(sample=NewIDseqName)
# 
# cms_all_data[is.na(cms_all_data$genus), c("genus", "species", "corrected.genus", "corrected.species")] <- 
#   data.frame(genus="Culex", species=c("erythrothorax", "tarsalis"),
#              corrected.genus="Culex", corrected.species=c("erythrothorax", "tarsalis"),
#              stringsAsFactors = FALSE)
# cms_all_data$location_name
merged_metadata <- read.csv(file.path(data_dir, "CMS001_CMS002_MergedAnnotations_181218.csv"))
merged_metadata$location_name <- apply(select(merged_metadata, paste0("collection_", c("long", "lat"))), 1, function (x) {
  place_name <- lapply(ca_shp$Areas@data[, c("INTPTLON", "INTPTLAT")], as.numeric) %>% 
    data.frame() %>%
    apply(., 1, calc_eucl_dist, x) %>%
    which.min() %>%
    `[`(ca_shp$Areas@data$NAME, .)
  if (length(place_name)==0) return (NA)
  return (place_name)
})
merged_metadata %<>% mutate(location=paste(collection_lat, collection_long, sep=",")) %>%
  rename(long=collection_long, lat=collection_lat) %>%
  mutate(long=sapply(long, function (x) {
    if (is.na(x)) return (x)
    ifelse(x< -1000, 100*(prettyNum(x)%>%stringr::str_split(., "e", simplify=TRUE)%>%`[`(1)%>%as.numeric()), x)
  }))
```

```{r plot-shapefiles, cache=TRUE}
polygon_plots <- lapply(with(ca_shp$Areas@data, id[order(order)]), function (x) {
  plot_data <- filter(tidy(ca_shp$Areas), id==x)
  width_height <- plot_data[, c("long", "lat")]%>%lapply(range)%>%sapply(diff)
  place_name <- filter(ca_shp$Areas@data, id==x)$NAME
  sites <- filter(merged_metadata, location_name==place_name) %>% 
    mutate(lat_lon=paste0(lat, "_", long)) %>% `[[`("lat_lon") %>% table()
  n <- sum(sites)
  create_map(filter(merged_metadata, location_name==place_name), polygon=plot_data, polygon_col=filter(ca_shp$Areas@data, id==x)$colours) +
    ggtitle(paste0(filter(ca_shp$Areas@data, id==x)$NAME, "\nNumber of sites: ", length(sites), "     n: ", n)) +
    theme(plot.title=element_text(hjust=0.5))
})
grouped_polygons <- group_polygons(polygon_plots)

full_demo_map <- arrangeGrob(grobs=list(base_cali_map, grouped_polygons), ncol=2, widths=c(3, 5))
ggsave("01_demographics/demographic_map.pdf", full_demo_map, width=18, height=8)

polygon_plots_species <-  lapply(with(ca_shp$Areas@data, id[order(order)]), function (x) {
  plot_data <- filter(tidy(ca_shp$Areas), id==x)
  width_height <- plot_data[, c("long", "lat")]%>%lapply(range)%>%sapply(diff)
  place_name <- filter(ca_shp$Areas@data, id==x)$NAME
  sites <- filter(merged_metadata, location_name==place_name) %>% mutate(lat_lon=paste0(lat, "_", long)) %>% `[[`("lat_lon") %>% table()
  n <- sum(sites)
  title_text <- paste0(filter(ca_shp$Areas@data, id==x)$NAME, "\nNumber of sites: ", length(sites), "     n: ", n)
  filter(merged_metadata, location_name==place_name) %>%
    split(., paste(.$corrected.genus, .$corrected.species)) %>%
    lapply(function (input_data) {
      create_map(filter(merged_metadata, location_name==place_name),
                 polygon=plot_data, 
                 polygon_col=filter(ca_shp$Areas@data, id==x)$colours,
                 use_species=TRUE) + 
        ggtitle(paste(input_data$corrected.genus, input_data$corrected.species)[1])
    }) %>%
    arrangeGrob(grobs=., nrow=1, top=title_text)
})

grouped_polygons_species <- group_polygons(polygon_plots_species)
full_demo_map_species <- arrangeGrob(grobs=list(base_cali_map, grouped_polygons_species), ncol=2, widths=c(3, 5))
ggsave("01_demographics/demographic_map_species.pdf", full_demo_map_species, width=18, height=8)

```


```{r demographics-tables, cache=TRUE}
species_by_location <- merged_metadata %>%
  group_by(., location_name, compute_genus, compute_species) %>%
  summarize(., n=n())
species_by_location_plot <- 
  ggplot(species_by_location, 
         aes(x=paste0(compute_genus, "\n", compute_species), 
             y=factor(location_name, levels=with(ca_shp$Areas@data, NAME[order(order, decreasing = TRUE)])), 
             fill=n)) +
  theme_light() +
  geom_tile() +
  facet_grid(.~compute_genus, scales="free_x") +
  xlab("Species") +
  ylab("Location") +
  scale_fill_continuous("Number of\nSamples")
ggsave("01_demographics/species_by_location.pdf", species_by_location_plot, width=12, height=4)

species_by_location %>%
  mutate(species_name=paste(compute_genus, compute_species)) %>%
  spread(., key=species_name, value=n, fill=0) %>%
  `[`(., , -2:-3) %>%
  write_csv(., "01_demographics/demographics.csv")
```

# 2. mNGS analysis

```{r mNGS}
# mNGS analysis pipeline
# Placeholder diagram for now
DiagrammeR::mermaid("
                    graph TD
                    id1[raw reads]--filtering--> reads 
                    reads --assembly --> id2[contigs] 
                    ", height=200)
```


# 3. Microbiota overview

```{bash}
/Users/lucy.li/anaconda3/bin/aws s3 cp s3://czbiohub-mosquito/contig_quality/contig_quality_df.tsv ../../../data/
  ```


```{r, message=FALSE, warning=FALSE}
contig_hits <- read_tsv(file.path(data_dir, "contig_quality_df.tsv"))
```

```{r}
table(contig_hits$kingdom)
```


```{r contig_hits_summary}
presence_df <- lapply(na.omit(unique(contig_hits$kingdom)), function (x) {
  tab <- filter(contig_hits, kingdom==x, !((qlength<500) & (pmatch_rapsearch<0.1) & (pmatch_gsnap<0.1)))
  tab$CommonName <- coalesce(!!!select(tab, CommonName, scientific_name)) %>%
    gsub("[", "", ., fixed=TRUE) %>% gsub("]", "", ., fixed=TRUE)
  all_species <- levels(factor(tab$CommonName))
  all_samples <- levels(factor(contig_hits$sample))
  lapply(all_species, function (species_x) {
    subdata <- filter(tab, CommonName==species_x)
    division <- subdata %>% select("Division") %>% `[`(1, 1) %>% unlist() %>% unname()
    data.frame(sample=all_samples, 
               present=all_samples%in%subdata$sample,
               taxa=species_x, 
               division=division,
               rank=subdata$rank[1])
  }) %>%
    do.call(what=rbind, .) %>%
    data.frame(., kingdom=x)
}) %>%
  do.call(what=rbind) %>%
  left_join(., merged_metadata %>%
              mutate(sample_mosquito_species=paste(compute_genus, compute_species) %>% gsub("NA NA", NA, .)) %>%
              rename(sample=NewIDseqName) %>%
              select(sample, location_name, sample_mosquito_species, blood_fed), by="sample")
write_tsv(presence_df, file.path(data_dir, "taxa_presence_table.tsv"), col_names=TRUE, quote_escape=FALSE)
```





```{r}
taxa_rank <- c("superkingdom", "kingdom", "subkingdom", "infrakingdom", "superphylum", "phylum", "phylum", "subphylum", "infraphylum", "microphylum", "superclass", "no rank", "class", "subclass", "infraclass", "parvclass", "cohort", "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "infratribe", "genus", "subgenus", "section", "series", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")
presence_by_kingdom <- lapply(levels(presence_df$kingdom), function (kingdom_x) {
  if (kingdom_x=="Metazoa") {
    kingdom_x <- c(kingdom_x, "Chordata", "Arthropoda")
  } else if (kingdom_x=="Eukaryota") {
    kingdom_x <- c(kingdom_x, "Fungi", "Viridiplantae", "Chordata", "Metazoa", "Arthropoda")
  }
  filter(presence_df, kingdom%in%kingdom_x) %>% 
    group_by(sample, division) %>% 
    summarize(ntaxa=sum(present), location_name=location_name[1], rank=rank[1], sample_mosquito_species=sample_mosquito_species[1])
}) %>%
  lapply(function (x) {
    division_order <- x %>%
      group_by(division) %>% 
      summarize(num_samples=sum(ntaxa>0), rank=rank[1]) %>% 
      arrange(match(rank, taxa_rank), num_samples) %>%
      slice(rev(seq(nrow(.))))
    x %<>% mutate(division=factor(division, levels=division_order$division))
    x
  }) %>%
  mapply(data.frame, ., kingdom=as.list(levels(presence_df$kingdom)), SIMPLIFY = FALSE) %>%
  do.call(what=rbind, .) %>%
  mutate(., sample=factor(sample, unique(sample[order(sample_mosquito_species)])))
write.table(presence_by_kingdom, file=file.path(data_dir, "presence_of_contigs_by_kingdoms.tsv"),
            sep="\t", row.names=TRUE, quote=FALSE)

```


```{r kingdom_numbers_plot, fig.width=10, fig.height=15}
presence_by_kingdom_plot <- 
  ggplot(presence_by_kingdom, aes(x=division, y=sample, alpha=ntaxa>0, fill=sample_mosquito_species)) +
  theme_bw() +
  geom_tile() +
  facet_grid(location_name~kingdom, scales="free", space="free") +
  scale_alpha_manual("", values=c(0, 1)) +
  guides(alpha="none") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text.y=element_text(angle=0, hjust=0))



ggsave("03_microbiota_overview/taxa_by_location.pdf", presence_by_kingdom_plot, width=15, height=15)
```

```{r prokaryotes}
filtered_bacteria_contigs <- filter(contig_hits, kingdom=="Bacteria", !duplicated(taxid), (pmatch_rapsearch>0.8) | (pmatch_gsnap>0.8) ) %>%
  filter(match(rank, taxa_rank)>match("genus", taxa_rank))
filtered_bacteria_contigs_ncbi = vector(mode="list", length=0)
for (x in 1:ceiling(length(filtered_bacteria_contigs$taxid)/500)) {
  start <- (x-1)*500+1
  end <- min(length(filtered_bacteria_contigs$taxid), x*500)
  filtered_bacteria_contigs_ncbi[[x]] <- filtered_bacteria_contigs$taxid[start:end] %>% reutils::efetch(., db="taxonomy") %>% content() %>% XML::xmlToList()
}
filtered_bacteria_contigs_ncbi_parsed <- lapply(filtered_bacteria_contigs_ncbi, function (x) {
  lapply(x, `[[`, "LineageEx") %>%
    lapply(lapply, data.frame) %>%
    lapply(function (y) {
      do.call(what=rbind, y)
    })
}) %>%
  do.call(what=c, .)
filtered_bacteria_results <- mapply(function (i, results) {
  taxa_ranks <- c("genus", "family", "order", "class", "phylum")
  c(filtered_bacteria_contigs$scientific_name[i],  as.character(results[match(taxa_ranks, results$Rank), "ScientificName"])) %>%
    as.list() %>%
    setNames(., c("species", taxa_ranks)) %>%
    data.frame() %>%
    data.frame(taxid=filtered_bacteria_contigs$taxid[i], .)
}, as.list(1:nrow(filtered_bacteria_contigs)), filtered_bacteria_contigs_ncbi_parsed, SIMPLIFY=FALSE) %>%
  do.call(what=rbind, .)
filtered_bacteria_unique <- contig_hits %>% filter(kingdom=="Bacteria", match(rank, taxa_rank)>match("genus", taxa_rank)) %>%
  split(., .$sample) %>% 
  lapply(function (x) filter(x, !duplicated(taxid)))
find_top_hits <- function (results, grouping, n) {
  results %>%
    split(., .[[grouping]]) %>% 
    sapply(function (x) sum(mapply(function(a, b) any(a%in%b), list(x$taxid), filtered_bacteria_unique %>% lapply(`[[`, "taxid")))) %>% 
    sort(decreasing=TRUE) %>% 
    data.frame() %>% 
    setNames("num_samples") %>% 
    tibble::rownames_to_column(grouping) %>% 
    mutate(prop_samples=num_samples/n)
}
bacteria_n <- filter(contig_hits, kingdom=="Bacteria")$sample %>% unique %>% length()
popular_phyla <- find_top_hits(filtered_bacteria_results, "phylum", bacteria_n) %>% filter(prop_samples>0.1)
popular_genera <- popular_phyla %>%
  apply(1, function (x) {
    find_top_hits(filter(filtered_bacteria_results, phylum==x[1]), "genus", bacteria_n) %>% filter(prop_samples>0.1)
  }) %>%
  setNames(., popular_phyla$phylum[seq(.)])
```




```{r prokaryotic_plot}
popular_phyla_plot <- ggplot(popular_phyla, aes(x=factor(phylum, levels=phylum), y=prop_samples)) +
  theme_bw() +
  geom_bar(fill="royalblue3", stat="identity") +
  scale_y_continuous(labels=scales::percent) +
  xlab("Phylum") + ylab("% of sample") +
  theme(text=element_text(size=12))
paste0("03_microbiota_overview/bacteria_phyla.", c("pdf", "png")) %>%
  sapply(ggsave, popular_phyla_plot, width=6.5, height=3)

popular_bacteria_genera_plot <- ggplot(popular_genera[1:2] %>% mapply(data.frame, phylum=as.list(names(.)), ., SIMPLIFY=FALSE) %>% do.call(what=rbind)) +
  theme_bw() +
  geom_bar(aes(x=factor(genus, levels=genus), y=prop_samples), fill="salmon1", stat="identity") +
  facet_grid(.~phylum, scales="free_x", space="free_x") +
  scale_y_continuous(labels=scales::percent) +
  xlab("Genus") + ylab("% of sample") +
  theme(text=element_text(size=12))
paste0("03_microbiota_overview/bacteria_genera.", c("pdf", "png")) %>%
  sapply(ggsave, popular_bacteria_genera_plot, width=12.5, height=3)

```


```{r viruses}
virus_plot <-
  ggplot(filter(presence_df, kingdom=="Viruses") %>%
           mutate(., sample=factor(sample, unique(sample[order(sample_mosquito_species)]))), aes(x=taxa, y=sample, alpha=present, fill=sample_mosquito_species)) +
  theme_bw() +
  geom_tile() +
  facet_grid(location_name~., scales="free", space="free") +
  scale_alpha_manual("", values=c(0, 1)) +
  guides(alpha="none") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=90, hjust=0))

ggsave("03_microbiota_overview/viruses.pdf", width=25, height=15)
```

```{r single_mosquito}
ggplot(filter(contig_hits, kingdom=="Viruses", pmatch_rapsearch>0.8 | pmatch_gsnap>0.8)) +
  
  ```


```{r}
ggplot(filter(presence_df, kingdom=="Fungi") %>%
         group_by(sample) %>%
         summarize(present=sum(present), location_name=location_name[1]),
       aes(x=location_name, y=present)) +
  theme_bw() +
  geom_boxplot() +
  geom_jitter()
```


```{r fungi}
fungi_plot <-
  ggplot(filter(presence_df, kingdom=="Fungi") %>%
           mutate(., sample=factor(sample, unique(sample[order(sample_mosquito_species)]))), aes(x=taxa, y=sample, alpha=present, fill=sample_mosquito_species)) +
  theme_bw() +
  geom_tile() +
  facet_grid(location_name~., scales="free", space="free") +
  scale_alpha_manual("", values=c(0, 1), guide=FALSE) +
  guides(alpha="none", fill=guide_legend("Mosquito species:", title.position="top")) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=90, hjust=0)) +
  xlab("Taxonomic groups") + ylab("Samples")


ggsave("03_microbiota_overview/fungi.pdf", width=25, height=15)
```

```{r}
dense_sampling <- merged_metadata %>%
  mutate(species_name=paste(compute_genus, compute_species)) %>%
  mutate(lat_lon=paste(lat, long)) %>%
  group_by(species_name, lat_lon) %>%
  summarize(n=n()) %>% 
  `[`(., which.max(.$n), )
```

```{r}
culex_erythrothorax_samples_df <- merged_metadata %>%
  mutate(species_name=paste(compute_genus, compute_species)) %>%
  mutate(lat_lon=paste(lat, long)) %>%
  filter(species_name %in% dense_sampling$species_name,
         lat_lon %in% dense_sampling$lat_lon)
culex_erythrothorax_samples <- culex_erythrothorax_samples_df %>%
  select(., NewIDseqName) %>% 
  unlist() %>%
  as.character()
culex_erythrothorax_samples_plot_df <- filter(presence_df, sample %in% culex_erythrothorax_samples) %>%
  mutate(kingdom=sapply(kingdom, function (x) ifelse(x%in%c("Bacteria", "Viruses", "Archaea"), as.character(x), "Eukaryotes"))) %>%
  mutate(sample=strsplit(as.character(sample), "_")%>%lapply(head, 2)%>%sapply(paste, collapse="_"))
culex_erythrothorax_samples_plot <- ggplot(culex_erythrothorax_samples_plot_df, aes(x=sample, y=taxa, alpha=present)) +
  theme_bw() +
  geom_tile() +
  facet_grid(kingdom~., space="free", scale="free") +
  theme(axis.text.y=element_blank(), strip.text.y=element_text(angle=0),
        axis.text.x=element_text(angle=90),
        text=element_text(size=20)) + 
  scale_alpha_manual("", values=c(0, 1), guide=FALSE) +
  guides(alpha="none") +
  xlab("Samples") + 
  ylab("Taxonomic groups")

culex_erythrothorax_samples_df %>% `[[`("collection_date") %>% as.character() %>% substr(., 1, 10) %>% as.Date() %>% range  

paste0("03_microbiota_overview/culex_erythrothorax_samples_plot.", c("pdf", "png")) %>%
  sapply(ggsave, culex_erythrothorax_samples_plot, width=8, height=8)

```

# 4. Bloodmeal

```{r}
bloodmeal_plot <-
  ggplot(filter(presence_df, kingdom=="Chordata") %>%
           mutate(., sample=factor(sample, unique(sample[order(sample_mosquito_species)]))), aes(x=taxa, y=sample, alpha=present, fill=sample_mosquito_species)) +
  theme_bw() +
  geom_tile() +
  facet_grid(location_name~division, scales="free", space="free") +
  scale_alpha_manual("", values=c(0, 1)) +
  guides(alpha="none") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=90, hjust=0))

ggsave("04_bloodmeal/bloodmeal_species.pdf", width=35, height=30)
```

```{r focused_bloodmeal_plot}
focused_bloodmeal_groups <- list(divisions=c("birds", "even-toed ungulates", "rabbits & hares"),
                                 species=c("domestic cat", "human", "house mouse"))
focused_bloodmeal_df <- filter(presence_df, kingdom=="Chordata") %>%
  mutate(., sample=factor(sample, unique(sample[order(sample_mosquito_species)])),
         taxa=droplevels(taxa)) %>%
  mutate(., group=apply(., 1, function (x) {
    if (x["taxa"] %in% focused_bloodmeal_groups$species) {
      return (unname(x["taxa"]))
    } else if (x["division"] %in% focused_bloodmeal_groups$divisions) {
      return (unname(x["division"]))
    } else {
      return ("other chordates")
    }
  })) %>%
  filter(., !(grepl("^CMS", sample) & is.na(location_name) & (!grepl("ater", sample)))) %>%
  mutate(location_name = gsub(" ", "\n", location_name)) %>%
  replace_na(list(location_name="Water", blood_fed="no")) %>%
  mutate(sample=strsplit(as.character(sample), "_")%>%lapply(head, 2)%>%sapply(paste, collapse="_")) %>%
  mutate(sample=factor(sample, unique(sample[order(sample_mosquito_species)]))) %>%
  group_by(sample, group) %>% 
  summarize(present=any(present), sample_mosquito_species=sample_mosquito_species[1], location_name=location_name[1], blood_fed=blood_fed[1]) %>% 
  mutate(location_name=factor(location_name, c(list_of_counties, "Water")%>%gsub(" ", "\n", .))) %>%
  mutate(group=factor(group, rev(unname(c(unlist(focused_bloodmeal_groups), "other chordates")))))
# bloodfed_face <- arrange(focused_bloodmeal_df, location_name, sample_mosquito_species) %>% filter(!duplicated(sample)) %>% replace_na(list(blood_fed="no"))
# bloodfed_face$new_name <- sapply(1:nrow(bloodfed_face), function (i) ifelse(bloodfed_face$blood_fed[i]=="yes", paste0("BLOODFED ", bloodfed_face$sample[i]), as.character(bloodfed_face$sample[i])))
# focused_bloodmeal_df$new_sample_name <- slice(bloodfed_face, match(focused_bloodmeal_df$sample, sample)) %>% mutate(new_name=factor(new_name, levels=slice(bloodfed_face, match(levels(focused_bloodmeal_df$sample), sample)))) %>% `[[`("new_name")



focused_bloodmeal_plot <- ggplot(focused_bloodmeal_df, aes(y=group, x=sample, alpha=present, fill=sample_mosquito_species)) +
  theme_bw() +
  geom_tile() +
  facet_grid(.~location_name, scales="free", space="free") + 
  scale_alpha_manual("", values=c(0, 1), guide=FALSE) +
  guides(alpha="none", fill=guide_legend("Mosquito species:", title.position="top")) +
  theme(strip.text.y=element_text(angle=0, hjust=0),
        text=element_text(size=40),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="top") +
  ylab("Taxonomic groups")
focused_bloodmeal_bloodfed_plot <- ggplot(filter(focused_bloodmeal_df, group=="human") %>% mutate(group="blood_fed"), aes(y=group, x=sample, alpha=blood_fed), fill="black") +
  theme_bw() +
  geom_tile() +
  facet_grid(.~location_name, scales="free", space="free") +
  scale_alpha_manual("", values=c(0, 1), guide=FALSE) +
  guides(alpha="none") +
  theme(strip.text.y=element_text(angle=0, hjust=0),
        text=element_text(size=30),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=90, size=10, vjust=0.5),
        strip.background.x=element_blank(),
        strip.text=element_blank()) +
  xlab("Samples")  

gA <- ggplotGrob(focused_bloodmeal_plot)
gB <- ggplotGrob(focused_bloodmeal_bloodfed_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
focused_bloodmeal_combined_plot <- arrangeGrob(gA, gB, ncol=1, heights=c(5,2))
#grid.arrange(focused_bloodmeal_plot, focused_bloodmeal_bloodfed_plot, ncol=1)

c("png", "pdf") %>%
  paste0("04_bloodmeal/bloodmeal_summary.", .) %>%
  sapply(., ggsave, focused_bloodmeal_combined_plot, width=35, height=15)
```
