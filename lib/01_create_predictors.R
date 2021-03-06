
# Libraries ---------------------------------------------------------------

library(raster)
library(tidyverse)
library(landscapemetrics)
library(foreach)
library(doParallel)
library(data.table)
library(rstudioapi)

rasterOptions(tmptime = 2, tmpdir = "temp/")
write("TMP = temp", file = file.path('~/.Renviron'))
source("lib/raster_as_data_table.R")
Rcpp::sourceCpp('lib/ngb_rcpp.cpp')

distmap_path <- "/home/csenf/Projects/mapping/results/version1.1"

# Settings ----------------------------------------------------------------

#cntr <- "czechia"
countries <- list.files(distmap_path)
years <- 1986:2016
radius <- 5000
grouping <- TRUE

### Check if country was already processed, otherwise switch to next

cntr <- 0
k <- 0
while (cntr == 0) {
  k <- k + 1
  if (!file.exists(paste0("data/patches/", countries[k]))) cntr <- countries[k]
}

print(cntr)

# Group neighboring patches -----------------------------------------------

if (grouping) {
  
  print("Grouping patches")
  
  tictoc::tic()
  
  dist_map <- raster(paste0(distmap_path, "/", cntr, "/disturbance_year_", cntr, ".tif"))
  
  years_reclass <- raster::unique(dist_map)
  
  # Loop through years and group neighboring patches of consecutive years
  
  tmp_folder <- paste0("temp/temp_grouping_", cntr)
  dir.create(tmp_folder, recursive = TRUE, showWarnings = FALSE)
  
  for (y in years) {
    
    print(y)
    
    if (y > min(years)) dist_map <- raster(paste0(tmp_folder, "/disturbance_", y - 1, "_grouped_", cntr, ".tif"))
    
    reclass_matrix <- matrix(c(years, years %in% c(y, y + 1)), ncol = 2)
    selection <- reclassify(dist_map, reclass_matrix)
    
    patches <- clump(selection)
    
    years_dt <- as.data.table.raster(dist_map) # Convert years into data.table
    names(years_dt) <- "year"
    patches_dt <- as.data.table.raster(patches) # Convert binary disturbance patches into data.table
    names(patches_dt) <- "patch"
    dt <- cbind(years_dt, patches_dt) # Combine both in new data.table
    dt[ , n := length(unique(na.omit(year))), by = patch] # Count number of unique years in each patch (can be 1 or 2, with two identifying neighboring patches of consecutive years)
    dt[ , year_new := as.integer(modal(year)), by = patch] # Calculate 'new' year based on modal value (most frequent value)
    dt[n == 2, year := year_new] # Overwirte 'old' year with 'new' year if patch consists of two years (n=2)
    
    values(dist_map) <- dt$year # Overwrite raster values with 'new' year
    
    writeRaster(dist_map, paste0(tmp_folder, "/disturbance_", y, "_grouped_", cntr, ".tif"), overwrite = TRUE)
    
    if (y < max(years)) rm(dist_map, selection, patches)
    
    gc()
    
  }
  
  writeRaster(dist_map, paste0("data/disturbances/disturbance_year_grouped_", cntr, ".tif"), overwrite = TRUE, type = "INT2U")
  #dist_map <- raster(paste0("data/disturbances/disturbance_year_grouped_", cntr, ".tif"))
  
  f <- list.files(tmp_folder, include.dirs = TRUE, full.names = TRUE, recursive = TRUE)
  file.remove(f)
  
  rm(dist_map, selection, patches)
  
  gc()
  
  print("Finished grouping patches")
  
  tictoc::toc()
  
}

# Identify patches --------------------------------------------------------

print("Identifying final patches")

tictoc::tic()

dist_map <- raster(paste0("data/disturbances/disturbance_year_grouped_", cntr, ".tif"))

years_reclass <- raster::unique(dist_map)

inp_size <- file.size(paste0("data/disturbances/disturbance_year_grouped_", cntr, ".tif")) * 1e-6
cores <- ifelse(inp_size < 20, 20, ifelse(inp_size < 40, 10, 1))

if (cores == 1) {
  
  patches_centroids <- foreach(i = 1:length(years)) %do% {
    
    patch_dir <- paste0("data/patches/", cntr)
    dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)
    
    dist_map_tmp <- reclassify(dist_map, matrix(c(years_reclass, years_reclass == years[i]), ncol = 2))
    
    patches <- clump(dist_map_tmp, gaps = FALSE, 
                     filename = paste0(patch_dir, "/patches_", years[i], "_", cntr, ".tif"), 
                     datatype = "INT4U", overwrite = TRUE)
    
    index_values <- expand.grid(col = 1:ncol(patches), row = 1:nrow(patches))
    index_values <- cbind(index_values, patch = getValues(patches))
    
    patches_centroids <- index_values %>%
      as.data.frame(.) %>%
      na.omit(.) %>%
      group_by(patch) %>%
      summarize(row = round(mean(row), 0),
                col = round(mean(col), 0)) %>%
      mutate(year = years[i])
    
    patches_centroids <- cbind(patches_centroids, xyFromCell(dist_map, cellFromRowCol(dist_map, patches_centroids$row, patches_centroids$col)))
    
    rm(patches)
    
    patches_centroids
    
  }
  
} else {

  cl <- makeCluster(cores)
  
  registerDoParallel(cl)
  
  patches_centroids <- foreach(i = 1:length(years), .packages = c("raster", "landscapemetrics", "tidyverse")) %dopar% {
    
    patch_dir <- paste0("data/patches/", cntr)
    dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)
    
    dist_map_tmp <- reclassify(dist_map, matrix(c(years_reclass, years_reclass == years[i]), ncol = 2))
    
    patches <- clump(dist_map_tmp, gaps = FALSE, 
                     filename = paste0(patch_dir, "/patches_", years[i], "_", cntr, ".tif"), 
                     datatype = "INT4U", overwrite = TRUE)
    
    index_values <- expand.grid(col = 1:ncol(patches), row = 1:nrow(patches))
    index_values <- cbind(index_values, patch = getValues(patches))
    
    patches_centroids <- index_values %>%
      as.data.frame(.) %>%
      na.omit(.) %>%
      group_by(patch) %>%
      summarize(row = round(mean(row), 0),
                col = round(mean(col), 0)) %>%
      mutate(year = years[i])
    
    patches_centroids <- cbind(patches_centroids, xyFromCell(dist_map, cellFromRowCol(dist_map, patches_centroids$row, patches_centroids$col)))
    
    rm(patches)
    
    patches_centroids
    
  }
  
  stopCluster(cl)
  
}

patches_centroids <- patches_centroids %>%
  bind_rows() %>%
  as_tibble()

write_csv(patches_centroids, paste0("data/patch_centroids/", cntr, "_patch_centroids.csv"))
#patches_centroids <- read_csv(paste0("data/patch_centroids/", cntr, "_patch_centroids.csv"))

gc()

print("Finished identifying final patches")

tictoc::toc()

# Calculate ngb metrics ---------------------------------------------------

print("Calculating neighborhood metrics")

tictoc::tic()

### Create kernel

maxd <- round(radius / 30)

dmat <- matrix(0, nrow= 2*maxd+1, ncol=2*maxd+1)

for (i in 1:(2 * maxd + 1))
  for (j in 1:(2 * maxd + 1))
    dmat[j, i] <- sqrt((maxd - i) * (maxd - i) + (maxd - j) * (maxd - j) ) * 30

dtab <- data.frame()

for (i in 1:(2 * maxd + 1))
  for (j in 1:(2 * maxd + 1))
    if (dmat[j, i] < radius)
      dtab <- rbind(dtab, data.frame(ix = j - maxd, iy = i - maxd))

### Loop through years and extract metrics

dist_map_matrix <- as.matrix(dist_map)

ngb_out <- vector("list", length(years))

for (i in 1:length(years)) {
  
  print(years[i])
  
  patches_centroids_tmp <- patches_centroids %>% filter(year == years[i])
  
  patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
  
  ngb <- countPxCentroidMT(dist_map_matrix, 
                           patches_centroids_tmp$col, 
                           patches_centroids_tmp$row,
                           patches_centroids_tmp$year,
                           dtab$ix, 
                           dtab$iy,
                           getValues(patches))
  
  colnames(ngb) <- c("n_t0", "n_total", "n_tminus1", "n_tplus1", "n_patches")
  ngb <- cbind(patches_centroids_tmp, ngb)
  ngb_out[[i]] <- as.data.frame(ngb)
  
}

### Write into final dataframe and save

ngb_metrics <- ngb_out %>%
  bind_rows() %>%
  mutate(perc_t0 = n_t0 / n_total,
         perc_tminus1 = n_tminus1 / n_total,
         perc_tplus1 = n_tplus1 / n_total) %>%
  as_tibble()

write_csv(ngb_metrics, paste0("data/predictors/ngb/nbg_", cntr, ".csv"))
#ngb_metrics <- read_csv(paste0("data/predictors/ngb/nbg_", cntr, ".csv"))

rm(patches)
rm(ngb)
rm(ngb_out)
rm(dist_map_matrix)
rm(dtab)
rm(dmat)

gc()

print("Finished calculating neighborhood metrics")

tictoc::toc()

# Patch metrics -----------------------------------------------------------

print("Calculating landscape metrics")

tictoc::tic()

### Calculate landscape metrics and save to disc

inp_size <- file.size(paste0("data/disturbances/disturbance_year_grouped_", cntr, ".tif")) * 1e-6
cores <- ifelse(inp_size < 20, 20, ifelse(inp_size < 40, 10, 1))

if (cores == 1) {
  
  
  ls_metrics <- foreach(i = 1:length(years)) %do% {
    
    patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
    edges <- boundaries(patches)
    
    lsm <- data.frame(patch = getValues(patches),
                      edges = getValues(edges)) %>%
      filter(!is.na(patch)) %>%
      group_by(patch) %>%
      summarize(area = n() * 900,
                perimeter = sum(edges) * 30) %>%
      mutate(frac = 2 * log(0.25 * perimeter) / log(area)) %>%
      mutate(year = years[i])
    
    gc()
    
    lsm
    
  }
  
} else {
  
  cl <- makeCluster(cores)
  
  registerDoParallel(cl)
  
  ls_metrics <- foreach(i = 1:length(years), .packages = c("raster", "tidyverse")) %dopar% {
    
    patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
    edges <- boundaries(patches)
    
    lsm <- data.frame(patch = getValues(patches),
                      edges = getValues(edges)) %>%
      filter(!is.na(patch)) %>%
      group_by(patch) %>%
      summarize(area = n() * 900,
                perimeter = sum(edges) * 30) %>%
      mutate(frac = 2 * log(0.25 * perimeter) / log(area)) %>%
      mutate(year = years[i])
    
    gc()
    
    lsm
    
  }
  
  stopCluster(cl)
  
}

ls_metrics <- ls_metrics %>%
  bind_rows()

write_csv(ls_metrics, paste0("data/predictors/lsm/lsm_", cntr, ".csv"))
#ls_metrics <- read_csv(paste0("data/predictors/lsm/lsm_", cntr, ".csv"))

gc()

print("Finished calculating landscape metrics")

tictoc::toc()

# Spectral data -----------------------------------------------------------

print("Calculating spectral metrics")

tictoc::tic()

landtrendr_metrics <- stack(paste0("/home/csenf/Projects/mapping/results/maps/", cntr, "/disturbance_metrics_", cntr, ".tif"))
names(landtrendr_metrics) <- as.vector(outer(c("year", "magnitude", "duration", "pre", "rate", "dsnr"), c("B5", "B7", "NBR", "TCW"), paste, sep = "."))
landtrendr_metrics <- subset(landtrendr_metrics, c(14, 16, 17))

landtrendr_metrics_df <- as.data.table.raster(landtrendr_metrics)

inp_size <- file.size(paste0("data/disturbances/disturbance_year_grouped_", cntr, ".tif")) * 1e-6
cores <- ifelse(inp_size < 20, 20, ifelse(inp_size < 40, 10, 1))

if (cores == 1) {
  
  landtrendr_metrics <- foreach(i = 1:length(years)) %do% {
    
    patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
    
    patches_df <- cbind(landtrendr_metrics_df, patch = getValues(patches)) %>%
      mutate(year = years[i]) %>%
      filter(!is.na(patch)) %>%
      as_tibble() %>%
      group_by(patch, year) %>%
      summarize_at(.vars = vars(magnitude.NBR, pre.NBR, rate.NBR), mean)
    
    patches_df
    
  }
  
} else {
  
  cl <- makeCluster(cores)
  
  registerDoParallel(cl)
  
  landtrendr_metrics <- foreach(i = 1:length(years), .packages = c("raster", "tidyverse")) %dopar% {
    
    patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
    
    patches_df <- cbind(landtrendr_metrics_df, patch = getValues(patches)) %>%
      mutate(year = years[i]) %>%
      filter(!is.na(patch)) %>%
      as_tibble() %>%
      group_by(patch, year) %>%
      summarize_at(.vars = vars(magnitude.NBR, pre.NBR, rate.NBR), mean)
    
    patches_df
    
  }
  
  stopCluster(cl)
  
}

landtrendr_metrics <- landtrendr_metrics %>%
  bind_rows()

write_csv(landtrendr_metrics, paste0("data/predictors/spectral/spectral_", cntr, ".csv"))
#landtrendr_patches <- read_csv(paste0("data/predictors/spectral/spectral_", cntr, ".csv"))

rm(landtrendr_metrics_df)

gc()

print("Finished calculating spectral metrics")

tictoc::tic()

# Combine everything and save to disc -------------------------------------

dat <- ngb_metrics %>%
  left_join(ls_metrics) %>%
  left_join(landtrendr_metrics)

write_csv(dat, paste0("data/predictors/attribution_predictors_", cntr, ".csv"))

print(paste0("Finished ", cntr, "!"))

# Remove everything and restart R session ---------------------------------

rm(list = ls())

gc(reset = TRUE)

removeTmpFiles(h=0)

restartSession(command = 'source("lib/01_create_predictors.R")')
