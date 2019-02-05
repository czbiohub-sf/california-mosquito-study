calc_eucl_dist <- function (v1, v2) {
  sqrt(sum((v2-v1)^2))
}

capitalize_first <- function (x) {
  A <- toupper(substr(x, 1, 1))
  B <- tolower(substr(x, 2, nchar(x)))
  paste0(A, B)
}

get_zoom_options <- function (zoom=NULL, long, lat, multiplier=0.1) {
  require(ggmap)
  output <- list(zoom=zoom)
  if (is.null(zoom)) {
    output$lon_range <- extendrange(long, f=multiplier)
    output$lat_range <- extendrange(lat, f=multiplier)
    output$zoom <- calc_zoom(output$lon_range, output$lat_range)
  } else {
    output$lon_range <- range(long)
    output$lat_range <- range(lat)
  }
  return (output)
}


get_summary_data <- function (x) {
  freq_table <- group_by(x, location) %>%
    summarise(n(), long=mean(long), lat=mean(lat)) %>%
    rename(Number=`n()`) %>%
    left_join(., x, "location") %>% 
    group_by(location) %>% 
    summarise_all(first) %>%
    select(-long.y, -lat.y) %>%
    rename(long=long.x, lat=lat.x)
  return (freq_table)
}

split_by_location <- function (input_data) {
  list_names <- c(unique(input_data$location_name), "Northern California", "Southern California")
  lapply(list_names, function (x) {
    if (x=="Northern California") {
      dataframe <- filter(input_data, location_name %in% c("Alameda", "Sacramento"))
    } else if (x == "Southern California") {
      dataframe <- filter(input_data, !(location_name %in% c("Alameda", "Sacramento")))
    } else {
      dataframe <- filter(input_data, location_name==x)
    }
    dataframe$location <- x
    dataframe %>% 
      mutate(location=paste(round(long, 2), round(lat, 2)))
  }) %>%
    setNames(., list_names)
}


create_map <- function (input_data, zoom=NULL, polygon=NULL, polygon_col=NULL, use_species=FALSE) {
  freq_table <- get_summary_data(input_data)
  if (use_species) {
    freq_table_species <- split(input_data, with(input_data, paste(corrected.genus, corrected.species))) %>%
      lapply(get_summary_data) %>%
      do.call(what=rbind, .)
  }
  zoom_options <- get_zoom_options (zoom,
                                    c(input_data$long, polygon$long), 
                                    c(input_data$lat, polygon$lat))
  map <- openmap(c(zoom_options$lat_range[2], zoom_options$lon_range[1]),
                 c(zoom_options$lat_range[1], zoom_options$lon_range[2]), 
                 zoom = zoom_options$zoom, type = "osm", mergeTiles = TRUE) %>%
    openproj(CRS("+proj=longlat +datum=WGS84"))
  if (use_species) {
    map.points <- geom_point(data=freq_table_species, 
                             aes(x=long, y=lat, size=Number, 
                                 color=paste(corrected.genus, corrected.species)))
  } else {
    map.points <- geom_point(data=freq_table,
                             aes(x=long, y=lat, size=Number),
                             alpha=.8)
  }
  map_plot <- autoplot(map) + 
    theme(axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          axis.title=element_blank()) +
    scale_size_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
  if (!is.null(polygon)) {
    map_plot %<>% `+`(geom_polygon(data=polygon, aes(x=long, y=lat, group=group), fill=polygon_col, alpha=.5))
  }
  map_plot + map.points
}

create_summary_map <- function (input_data, zoom=NULL) {
  freq_table <- get_summary_data(input_data) %>%
    mutate(label=paste0(location_name, "\n", Number, " samples")) %>%
    mutate(long_label=long-0.55*ifelse(long>mean(long), 1, -1),
           lat_label=lat-0.55*ifelse(lat>mean(lat), 1, -1))
  zoom_options <- get_zoom_options (zoom, input_data$long, input_data$lat)
  map <- openmap(c(zoom_options$lat_range[2], zoom_options$lon_range[1]),
                 c(zoom_options$lat_range[1], zoom_options$lon_range[2]), 
                 zoom = zoom_options$zoom,
                 type = "osm",
                 mergeTiles = TRUE)
  map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  map.points <- geom_point(data=freq_table,
                           aes(x=long, y=lat),
                           alpha=.8)
  map.segments <- geom_segment(data=freq_table, 
                               aes(x=long, xend=long_label, y=lat, yend=lat_label))
  map.labels <- geom_label(data=freq_table,
                           aes(x=long_label, y=lat_label, label=label),
                           position=position_jitter())
  autoplot(map.latlon) + 
    map.points +
    map.segments +
    map.labels +
    theme(axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          axis.title=element_blank())
}

group_polygons <- function (input_plots) {
  split(input_plots, c(1, 1, 2, 2, 2)) %>%
    mapply(function (x, y) {
      arrangeGrob(grobs=x, nrow=1)
    }, ., split(seq(polygons), c(1, 1, 2, 2, 2)), SIMPLIFY=FALSE) %>%
    arrangeGrob(grobs=., ncol=1)
}
