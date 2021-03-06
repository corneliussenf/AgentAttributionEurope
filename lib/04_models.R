
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

p <- dat_mod_sample %>% 
  gather(key = predictor, value = value, -(patch:y), -(agent:country)) %>%
  filter(predictor %in% c("area", "frac", "magnitude.NBR", "n_patches", "perc_t0", "perc_tminus1", "perc_tplus1", "perimenter", "pre.NBR", "rate.NBR")) %>%
  mutate(agnt = ifelse(agent == "breakage", "Storm", ifelse(agent == "background", "Other", "Fire"))) %>%
  ggplot(., aes(x = agent, y = value, fill = agent)) +
  geom_boxplot() +
  facet_wrap(~predictor, scales = "free") +
  scale_y_log10() +
  scale_fill_manual(values = c("grey", "#004488", "#BB5566")) +
  labs(x = NULL, y = "Value", fill = "Agent") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        plot.title = element_text(size = 11),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))

ggsave("figures/predictors.pdf", p, width = 7.5, height = 7)

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

# Variable importance -----------------------------------------------------

load(file = "temp/fit_rf.RData")

p <- importance(fit.rf) %>%
  as.data.frame() %>%
  rownames_to_column(var = "predictor") %>%
  ggplot(., aes(x = MeanDecreaseGini, y = reorder(predictor, MeanDecreaseGini))) +
  geom_segment(aes(x = 0, 
                   xend = MeanDecreaseGini, 
                   y = reorder(predictor, MeanDecreaseGini),
                   yend = reorder(predictor, MeanDecreaseGini))) +
  geom_point(col = "darkgrey", size = 3) +
  labs(x = "Variable importance", y = NULL) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11),
        strip.background = element_blank())

ggsave("figures/variable_importance.pdf", p, width = 5.5, height = 5.5)
