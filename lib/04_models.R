
# Libraries ---------------------------------------------------------------

library(randomForest)
library(tidyverse)
library(blockCV)
library(sf)
library(pROC)

rasterOptions(tmptime = 2, tmpdir = "temp/")

years <- 1986:2016

# Get reference data ------------------------------------------------------

dat_mod <- list.files("data/model_input", pattern = ".csv$", full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

# Sample pseudo absences and split into cal/val ---------------------------

dat_mod <- dat_mod %>%
  filter(agent != "biotic")

### Equal number of pseudo-absences

dat_mod_luc <- dat_mod %>%
  filter(agent == "landusechange") %>%
  mutate(agent = "background")

dat_mod_no_luc <- dat_mod %>%
  filter(agent != "landusechange")

sampsizes <- table(dat_mod_no_luc$agent)

sampsizes["background"] <- sum(sampsizes[c("breakage", "fire")], na.rm = TRUE)

dat_mod_sample <- dat_mod_no_luc %>%
  split(.$agent) %>%
  map2(.y = sampsizes,
       ., ~ sample_n(., .y)) %>%
  bind_rows()

dat_mod_sample <- list(dat_mod_sample, dat_mod_luc) %>%
  bind_rows()

save(dat_mod_sample, file = "temp/dat_mod_sample.RData")

### Plots

dat_mod_sample %>% 
  gather(key = predictor, value = value, -(patch:y), -(agent:country)) %>%
  ggplot(., aes(x = agent, y = value, fill = agent)) +
  geom_boxplot() +
  facet_wrap(~predictor, scales = "free") +
  scale_y_log10() +
  scale_fill_manual(values = c("grey", "#004488", "#BB5566"))

dat_mod_sample %>%
  gather(key = predictor, value = value, -(patch:y), -(agent:country)) %>%
  filter(predictor %in% c("area", "frac", "perimeter")) %>%
  mutate(predictor = case_when(
    predictor == "area" ~ "Size (ha)",
    predictor == "frac" ~ "Fractional dimension",
    predictor == "perimeter" ~ "Perimenter (m)"
  )) %>%
  mutate(agent = case_when(
    agent == "background" ~ "Background",
    agent == "breakage" ~ "Breakage",
    agent == "fire" ~ "Fire"
  )) %>%
  ggplot(., aes(x = value, fill = agent)) +
  geom_density(alpha = 0.3, adjust = 2) +
  facet_wrap(~predictor, scales = "free", ncol = 3) +
  scale_x_log10(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("grey", "#004488", "#BB5566")) +
  theme_classic() +
  labs(x = NULL, 
       y = "Density",
       fill = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        strip.background = element_blank())

dat_mod_sample %>%
  gather(key = predictor, value = value, -(patch:y), -(agent:country)) %>%
  filter(predictor %in% c("magnitude.NBR", "rate.NBR", "pre.NBR")) %>%
  mutate(predictor = case_when(
    predictor == "magnitude.NBR" ~ "Spectral magnitude",
    predictor == "rate.NBR" ~ "Spectral recovery",
    predictor == "pre.NBR" ~ "Pre-dist. spectral value"
  )) %>%
  mutate(agent = case_when(
    agent == "background" ~ "Background",
    agent == "breakage" ~ "Breakage",
    agent == "fire" ~ "Fire"
  )) %>%
  ggplot(., aes(x = value, fill = agent)) +
  geom_density(alpha = 0.3, adjust = 2) +
  facet_wrap(~predictor, scales = "free", ncol = 3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("grey", "#004488", "#BB5566")) +
  theme_classic() +
  labs(x = NULL, 
       y = "Density",
       fill = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        strip.background = element_blank())

dat_mod_sample %>%
  gather(key = predictor, value = value, -(patch:y), -(agent:country)) %>%
  filter(predictor %in% c("n_patches", "perc_t0")) %>%
  mutate(predictor = case_when(
    predictor == "n_patches" ~ "Number of patches",
    predictor == "perc_t0" ~ "Percent of disturbance area"
  )) %>%
  mutate(agent = case_when(
    agent == "background" ~ "Background",
    agent == "breakage" ~ "Breakage",
    agent == "fire" ~ "Fire"
  )) %>%
  ggplot(., aes(x = value, fill = agent)) +
  geom_density(alpha = 0.3, adjust = 2) +
  facet_wrap(~predictor, scales = "free", ncol = 2) +
  scale_x_log10(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("grey", "#004488", "#BB5566")) +
  theme_classic() +
  labs(x = NULL, 
       y = "Density",
       fill = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        strip.background = element_blank())


# Write out for control
write_csv(dat_mod_sample, "temp/dat_mod_sample.csv")

# Fit models --------------------------------------------------------------

sampsizes_train <- min(sampsizes)

dat_mod_sample <- dat_mod_sample %>%
  mutate(agent = factor(agent))

fit.rf <- randomForest(agent ~ n_patches + perc_t0 + perc_tminus1 + perc_tplus1 + 
                         area + frac + 
                         magnitude.NBR + rate.NBR + pre.NBR + 
                         x + y, 
                       data = dat_mod_sample)

plot(fit.rf)

fit.rf

save(fit.rf, file = "temp/fit_rf.RData")

# Model performance --------------------------------------------------------

load(file = "temp/dat_mod_sample.RData")

### Create fold for cross-validation

presence_sf <- st_as_sf(x = dat_mod_sample, 
         coords = c("x", "y"),
         crs = st_crs(3035))

sb <- spatialBlock(speciesData = presence_sf,
                   species = "agent",
                   theRange = 100000, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 100)

folds <- sb$folds

### Run cross-validation

testdat <- dat_mod_sample
testdat$agent_pred <- NA
testdat$breakage_prob_pred <- NA
testdat$fire_prob_pred <- NA

for(k in seq_len(length(folds))){
  
  cal <- unlist(folds[[k]][1]) 
  val <- unlist(folds[[k]][2]) 
  
  dat_mod_sample$agent <- factor(dat_mod_sample$agent)
  
  rf <- randomForest(agent ~ n_patches + perc_t0 + perc_tminus1 + perc_tplus1 + 
                       area + frac + 
                       magnitude.NBR + rate.NBR + pre.NBR + 
                       x + y, 
                     data = dat_mod_sample[cal, ])
  
  
  testdat[val, "agent_pred"] <- predict(rf, dat_mod_sample[val, ], type = "class")
  testdat[val, "breakage_prob_pred"] <- predict(rf, dat_mod_sample[val, ], type = "prob")[, 2]
  testdat[val, "fire_prob_pred"] <- predict(rf, dat_mod_sample[val, ], type = "prob")[, 3]
  
}

### AUC

testdat <- testdat %>%
  mutate(observed_breakage = ifelse(agent == "breakage", 1, 0),
         observed_fire = ifelse(agent == "fire", 1, 0))

roc_breakage <- roc(testdat$observed_breakage, testdat$breakage_prob_pred)
auc(roc_breakage)

roc_fire <- roc(testdat$observed_fire, testdat$fire_prob_pred)
auc(roc_fire)

### Accuracies

testdat <- testdat %>%
  mutate(agent_pred = factor(agent_pred, labels = c("background", "breakage", "fire"))) %>%
  mutate(agent_pred_threshold = case_when(
    agent_pred == "breakage" & breakage_prob_pred > 0.5 ~ "breakage",
    agent_pred == "fire" & fire_prob_pred > 0.5 ~ "fire",
    TRUE ~ "background"
  ))

conf <- table(testdat$agent_pred_threshold, testdat$agent)

1 - (sum(diag(conf)) / sum(conf)) # overall error
1 - (diag(conf) / rowSums(conf)) # commission error
1 - (diag(conf) / colSums(conf)) # omission error


