

# Libraries and functions -------------------------------------------------

library(tidyverse)
library(raster)
library(sf)
library(stars)
library(fasterize)
library(rgeos)
library(gdalUtils)

# Create gridded and tiles version of disturbance map ---------------------

countries <- list.files("../TimeSync/data/countries/", pattern = ".shp$", full.names = TRUE) %>%
  map(read_sf) %>%
  purrr::reduce(., sf:::rbind.sf)

countries <- lwgeom::st_make_valid(countries)

# grid of 30 m reference grid (pixel-level)
grid_30m <- raster(ext = extent(countries), resolution = 30, crs = CRS(st_crs(countries)[[2]]))

# grid of 9000 m reference grid (attribution-level)
grid_9000m <- aggregate(grid_30m, fact = 300)
values(grid_9000m) <- 1:ncell(grid_9000m)
names(grid_9000m) <- "grid_id"
grid_9000m_shp <- st_as_stars(grid_9000m) %>% st_as_sf() # Raster to sf

for (i in sample(1:nrow(grid_9000m_shp), 100)) {
  print(i)
  t <- grid_30m %>% crop(., grid_9000m_shp[i, ])
  values(t) <- 1:ncell(t)
  names(t) <- "pixel_id"
  writeRaster(t, paste0("data/reference_tiles/tile_9000m2_", i, ".tif"), datatype = "INT4U")
}

# Resample disturbance maps to tiles --------------------------------------

countries_master <- read_csv("../TimeSync/data/countries_master.csv")

disturbances <- list.files("../Mapping/results/maps/", pattern = glob2rx("*metrics*tif"), recursive = TRUE, full.names = TRUE) %>%
  grep(glob2rx("*unmasked*"), ., invert = TRUE, value = TRUE)

forests <- list.files("../Mapping/results/maps/", pattern = glob2rx("*forest*tif"), recursive = TRUE, full.names = TRUE) %>%
  grep(glob2rx("*unmasked*"), ., invert = TRUE, value = TRUE)

cntr <- strsplit(disturbances, "_") %>% map(., ~ strsplit(.[3], "\\.")[[1]][1]) %>% unlist(.)

tile_ids <- list.files("data/reference_tiles", ".tif$") %>%
  str_split(., "_") %>%
  map(., ~ str_split(.[3], "\\.")[[1]][1])

for (tile_id in tile_ids[c(1, 2, 3)]) {
  
  print(paste0("Processing ", tile_id))
  
  tile <- raster(paste0("data/reference_tiles/tile_9000m2_", tile_id, ".tif"))
  c_tile <- st_crop(countries, tile)
  cntr_sel <- countries_master[countries_master$iso_code %in% c_tile$ISO_CC, "country_name_short"][[1]]
  
  disturbances_sel <- disturbances[which(cntr %in% cntr_sel)]
  disturbances_sel_ras <- disturbances_sel %>% map(stack)
  disturbances_sel_ras_aligned <- disturbances_sel_ras %>% 
    map(., ~ resample(., tile, method = "ngb"))
  
  forests_sel <- forests[which(cntr %in% cntr_sel)]
  forests_sel_ras <- forests_sel %>% map(raster)
  forests_sel_ras_aligned <- forests_sel_ras %>% 
    map(., ~ resample(., tile, method = "ngb"))
  
  if (length(disturbances_sel_ras_aligned) > 1) {
    disturbances_sel_ras_aligned <- do.call("merge", disturbances_sel_ras_aligned)
    forests_sel_ras_aligned <- do.call("merge", forests_sel_ras_aligned)
  } else {
    disturbances_sel_ras_aligned <- disturbances_sel_ras_aligned[[1]]
    forests_sel_ras_aligned <- forests_sel_ras_aligned[[1]]
  }
  
  out <- stack(disturbances_sel_ras_aligned, forests_sel_ras_aligned)
  
  writeRaster(out, paste0("data/metrics_tiles/metrics_9000m2_", tile_id, ".tif"))
  
}




#

### AFTER attribution!

grid_attributes <- read_sf("data/referencegrid/referencegrid_10km.shp")

grid_attributes <- grid_attributes %>% filter(., !is.na(agent))

table(grid_attributes$agent)

disturbances <- list.files("../Mapping/results/maps/", pattern = glob2rx("*metrics*tif"), recursive = TRUE, full.names = TRUE) %>%
  grep(glob2rx("*unmasked*"), ., invert = TRUE, value = TRUE)

forests <- list.files("../Mapping/results/maps/", pattern = glob2rx("*forest*tif"), recursive = TRUE, full.names = TRUE) %>%
  grep(glob2rx("*unmasked*"), ., invert = TRUE, value = TRUE)

for (i in 1:length(disturbances)) {
  
  cntr <- strsplit(strsplit(disturbances[[i]], "_")[[1]][3], "\\.")[[1]][1]
  
  print(paste0("Extracting data for ", cntr))
  
  dist <- stack(disturbances[[i]])
  fore <- stack(forests[[i]])
  grid_attributes_sel <- grid_attributes %>% st_crop(extent(dist))
  
  if (nrow(grid_attributes_sel) > 0) {
    
    for (j in 1:nrow(grid_attributes_sel)) {
      
      print(paste0("...polygon ", j, " of ", nrow(grid_attributes_sel)))
      
      dist_sel <- dist %>% crop(grid_attributes_sel[j, ], snap = "near")
      
      if (mean(is.na(values(subset(dist_sel, 1)))) < 1) {
        
        fore_sel <- fore %>% crop(grid_attributes_sel[j, ])
        
        out <- stack(dist_sel, fore_sel)
        
        writeRaster(out, paste0("data/training_nnet/training_", grid_attributes_sel[j, "gridid"][[1]], "_", grid_attributes_sel[j, "agent"][[1]], "_", cntr, ".tif"))
        
      }
       
    }
    
  }
  
}

# Combine grid-cells intersecting with two countries

files <- list.files("data/training_nnet", pattern = glob2rx("training*tif"), full.names = TRUE)

duplicates <- files %>%
  str_split("_") %>%
  map(., ~ .[3]) %>%
  unlist(.) %>%
  .[duplicated(.)]

for (d in duplicates) {
  
  agent <- grep(d, files, value = TRUE) %>%
    str_split("_") %>%
    map(., ~ .[4]) %>%
    unlist(.) %>%
    unique(.)
  
  cntr <- grep(d, files, value = TRUE) %>%
    str_split("_") %>%
    map(., ~ str_split(.[5], "\\.")[[1]][1]) %>%
    unlist(.) %>%
    paste(., collapse = "-")
  
  files_d <- grep(d, files, value = TRUE) %>%
    map(stack)
  
  # Resample if clipping led to slightly smaller/larger extent (i.e., 334 or 332 pixels). This also fixes origin issues during mosaicing.
  dims <- files_d %>% map(., ~ sum(dim(.)[1:2] == 333))
  ref <- files_d[unlist(dims) == 2][[1]]
  files_d_rsp <- map(files_d[unlist(dims) != 2], ~ resample(., y = ref, method = "ngb"))
  
  # Mosaic
  files_d <- c(files_d[unlist(dims) == 2], files_d_rsp, fun = "mean", tolerance = 0.1)
  files_d <- do.call(mosaic, files_d)
  
  writeRaster(files_d, paste0("data/training_nnet/training_", d, "_", agent, "_", cntr, ".tif"))
  
}

rm()
