
# Libraries ---------------------------------------------------------------

library(raster)
library(sf)
library(tidyverse)

# Get data ----------------------------------------------------------------

ref_pts <- list.files("data/referencepoint/by_country", pattern = ".shp", full.names = TRUE)

years <- 1986:2016

countries <- strsplit(basename(ref_pts), "_") %>%  
  map(., ~ substring(.[3], 1, (nchar(.[3]) - 4))) %>%
  unlist()

# Run through extraction by country ---------------------------------------

process <- c("ukraine", "unitedkingdom")

for (j in 1:length(countries)) {
  
  cntr <- countries[j]
  
  if (cntr %in% process) {
    
    outfile <- paste0("data/model_input/attribution_model_input_", cntr, ".csv")
    
    if (!file.exists(outfile)) {
      
      print(paste0("Processing ", cntr))
      
      ref <- read_sf(ref_pts[j])
      
      extr <- vector("list", length(years))
      
      for (i in 1:length(years)) {
        
        print(years[i])
        
        patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))  
        extr_tmp <- raster::extract(patches, ref, df = TRUE)
        extr_tmp$year <- years[i]
        extr_tmp$ID <- NULL
        names(extr_tmp) <- c("patch", "year")
        extr_tmp <- extr_tmp %>%
          mutate(agent = ref$agent,
                 id = ref$id)
        extr[[i]] <- extr_tmp
        
      }
      
      extr_sum <- extr %>%
        bind_rows() %>%
        na.omit() %>%
        group_by(patch, year) %>%
        summarise(agent = unique(agent),
                  id = id[1]) %>%
        ungroup()
      
      predictors <- read_csv(paste0("data/predictors/attribution_predictors_", cntr, ".csv"))
      
      dat_mod <- predictors %>%
        left_join(extr_sum, by = c("year", "patch")) %>%
        mutate(agent = ifelse(is.na(agent), "background", agent)) %>%
        mutate(id = ifelse(is.na(id), 0, id)) %>%
        mutate(country = cntr)
      
      write_csv(dat_mod, outfile)
      
    }
    
  }
  
}

### Add Moldova, for which no reference data is available

inp <- read_csv("data/predictors/attribution_predictors_moldova.csv")
inp$id <- 0
inp$agent <- "background"
inp$country <- "moldova"
write_csv(inp, "data/model_input/attribution_model_input_moldova.csv")
