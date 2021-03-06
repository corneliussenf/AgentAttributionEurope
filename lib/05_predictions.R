
library(randomForest)
library(raster)
library(tidyverse)

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
    pred_prob <- predict(fit.rf, newdata = dat_pred, type = "prob")
    dat_pred$prediction <- as.integer(pred)
    dat_pred$prediction_class <- pred
    dat_pred$prediction_breakage_prob <- pred_prob[, "breakage"]
    dat_pred$prediction_fire_prob <- pred_prob[, "fire"]
    dat_pred$prediction_background_prob <- pred_prob[, "background"]
    
    nrow_before <- nrow(dat_pred)
    dat_pred <- dat_pred[!is.na(dat_pred$prediction), ]
    nrow_after <- nrow(dat_pred)
    
    nrow_before - nrow_after
    
    dat_pred$certainty_background <- NA
    dat_pred$certainty_confusion <- NA
    
    dat_pred[dat_pred$prediction_class == "breakage", "certainty_background"] <- abs((dat_pred[dat_pred$prediction_class == "breakage", "prediction_background_prob"] - dat_pred[dat_pred$prediction_class == "breakage", "prediction_breakage_prob"])[[1]])
    dat_pred[dat_pred$prediction_class == "breakage", "certainty_confusion"] <- abs((dat_pred[dat_pred$prediction_class == "breakage", "prediction_fire_prob"] - dat_pred[dat_pred$prediction_class == "breakage", "prediction_breakage_prob"])[[1]])
    
    dat_pred[dat_pred$prediction_class == "fire", "certainty_background"] <- abs((dat_pred[dat_pred$prediction_class == "fire", "prediction_background_prob"] - dat_pred[dat_pred$prediction_class == "fire", "prediction_fire_prob"])[[1]])
    dat_pred[dat_pred$prediction_class == "fire", "certainty_confusion"] <- abs((dat_pred[dat_pred$prediction_class == "fire", "prediction_breakage_prob"] - dat_pred[dat_pred$prediction_class == "fire", "prediction_fire_prob"])[[1]])
    
    write_csv(dat_pred, paste0("results/predictions/prediction_", cntr, ".csv"))
    
    rcls <- vector("list", length(years))
    rcls_certainty_background <- vector("list", length(years))
    rcls_certainty_confusion <- vector("list", length(years))
    
    for (i in 1:length(years)) {
      
      print(years[i])
      
      patches <- raster(paste0("data/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))  
      
      rcls[[i]] <- reclassify(patches, rcl = as.matrix(dat_pred[dat_pred$year == years[i], c("patch", "prediction")]))
      rcls_certainty_background[[i]] <- reclassify(patches, rcl = as.matrix(dat_pred[dat_pred$year == years[i], c("patch", "certainty_background")]))
      rcls_certainty_confusion[[i]] <- reclassify(patches, rcl = as.matrix(dat_pred[dat_pred$year == years[i], c("patch", "certainty_confusion")]))
      
    }
    
    rcls_stack <- stack(rcls)
    rcls_stack_certainty_background <- stack(rcls_certainty_background)
    rcls_stack_certainty_confusion <- stack(rcls_certainty_confusion)
    
    rcls_stack_sum <- sum(rcls_stack, na.rm = TRUE)
    rcls_stack_certainty_background_sum <- sum(rcls_stack_certainty_background, na.rm = TRUE)
    rcls_stack_certainty_confusion_sum <- sum(rcls_stack_certainty_confusion, na.rm = TRUE)
    
    writeRaster(rcls_stack_sum, paste0("results/maps/agent_classes_", cntr, ".tif"), overwrite = TRUE)
    writeRaster(rcls_stack_certainty_background_sum, paste0("results/maps/agent_certainty_background_", cntr, ".tif"), overwrite = TRUE)
    writeRaster(rcls_stack_certainty_confusion_sum, paste0("results/maps/agent_certainty_confusion_", cntr, ".tif"), overwrite = TRUE)
    
  }
  
}

# Remove everything, restart R and rerun

rm(list = ls())

gc(reset = TRUE)

removeTmpFiles(h=0)

restartSession(command = 'source("lib/05_predictions.R")')
