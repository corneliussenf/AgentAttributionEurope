
library(randomForest)
library(raster)

# Load files and predict agent --------------------------------------------

inputs <- list.files("data/model_input", pattern = ".csv$", full.names = TRUE)

load("temp/fit_rf.RData")

years <- 1986:2016

for (inp in inputs) {

  cntr <- strsplit((strsplit(basename(inp), "_")[[1]][4]), "\\.")[[1]][1]
  
  if (!file.exists(paste0("results/maps/agent_classes_", cntr, ".tif"))) {
    
    print(cntr)
    
    dat_pred <- read_csv(paste0("data/predictors/attribution_predictors_", cntr, ".csv"))
    pred <- predict(fit.rf, newdata = dat_pred, type = "class")
    dat_pred$prediction <- as.integer(pred)
    dat_pred$prediction_class <- pred
    
    write_csv(dat_pred, paste0("results/predictions/prediction_", cntr, ".csv"))
    
    rcls <- vector("list", length(years))
    
    for (i in 1:length(years)) {
      
      print(years[i])
      
      patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))  
      
      rcls[[i]] <- reclassify(patches, rcl = as.matrix(dat_pred[dat_pred$year == years[i], c("patch", "prediction")]))
      
    }
    
    rcls_stack <- stack(rcls)
    
    rcls_stack_sum <- sum(rcls_stack, na.rm = TRUE)
    
    writeRaster(rcls_stack_sum, paste0("results/maps/agent_classes_", cntr, ".tif"), overwrite = TRUE)
    
    
  }
  
}

# Remove everything, restart R and rerun

rm(list = ls())

gc(reset = TRUE)

removeTmpFiles(h=0)

restartSession(command = 'source("lib/05_predictions.R")')
