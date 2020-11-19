
library(tidyverse)
library(raster)
library(landscapemetrics)
library(data.table)

### Functions

# -> This functions converts a raster into a data.table, which is a faster and more RAM-safe version of a data.frame

as.data.table.raster <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = canProcessInMemory(x, 2), ...) {
  stopifnot(require("data.table"))
  if(inmem) {
    v <- as.data.table(as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
  } else {
    tr <- blockSize(x, n=2)
    l <- lapply(1:tr$n, function(i) 
      as.data.table(as.data.frame(getValues(x, 
                                            row=tr$row[i], 
                                            nrows=tr$nrows[i]), 
                                  row.names=row.names, optional=optional, xy=xy, ...)))
    v <- rbindlist(l)
  }
  coln <- names(x)
  if(xy) coln <- c("x", "y", coln)
  setnames(v, coln)
  v
}

### Loop through years and group neighboring patches of consecutive years

rasterOptions(maxmemory = 1e+09) # I increase the number of pixels stored in RAM to increase speed

years <- raster("../Mapping/results/maps/slovenia/disturbance_year_filtered_slovenia.tif") # Get raster of disturbance year

year_range <- 1986:2018 # Define year range

# Loop through years

for (y in year_range) {
  
  print(y)
  
  # First step: Reclassify consecutive years into binary raster
  reclass_matrix <- matrix(c(year_range, year_range %in% c(y, y + 1)), ncol = 2) # Defines the reclassification matrix
  selection <- reclassify(years, reclass_matrix)
  
  # Second step: Identify neighboring patches
  patches <- get_patches(selection, class = 1, directions = 4)
  
  # Third step: Assign modal year to combined patches
  years_dt <- as.data.table.raster(years) # Convert years into data.table
  names(years_dt) <- "year"
  patches_dt <- as.data.table.raster(patches$`1`) # Convert binary disturbance patches into data.table
  names(patches_dt) <- "patch"
  dt <- cbind(years_dt, patches_dt) # Combine both in new data.table
  dt[ , n := length(unique(na.omit(year))), by = patch] # Count number of unique years in each patch (can be 1 or 2, with two identifying neighboring patches of consecutive years)
  dt[ , size := length(n), by = list(patch, year)] # Calculate size of neighboring patches
  dt[ , size := min(size), by = patch] # Identify minimum size of neighboring patches (for applying the MMU)
  dt[ , year_new := as.integer(modal(year)), by = patch] # Calculate 'new' year based on modal value (most frequent value)
  dt[n == 2, year := year_new] # Overwirte 'old' year with 'new' year if patch consists of two years (n=2)
  
  values(years) <- dt$year # Overwrite raster values with 'new' year

}

writeRaster(years, "temp/patch_aggregation_test_slovenia.tif", overwrite = TRUE)


