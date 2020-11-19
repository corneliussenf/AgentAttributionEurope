
library(tidyverse)
library(patchwork)

### Global functions

country_renamer <- function(cntr) {
  
  cntr <- str_to_title(cntr)
  
  cntr <- gsub("Bosniaandherzegovina", "Bosnia and Herzegovina", cntr)
  cntr <- gsub("Macedonia", "North Macedonia", cntr)
  cntr <- gsub("Unitedkingdom", "United Kingdom", cntr)
  
  return(cntr)
  
}

# Get data ----------------------------------------------------------------

load(file = "temp/out_temporal.RData")

wind <- read_csv("data/FORWIND/forwind_country_totals.csv")

dat <- out_temporal %>% 
  filter(agent == "breakage") %>%
  mutate(area = area / 10000) %>%
  dplyr::select(-agent) %>%
  group_by(country) %>%
  dplyr::rename(landsat_breakage_ha = area) %>%
  summarize(landsat_breakage_ha = sum(landsat_breakage_ha, na.rm = TRUE)) %>%
  ungroup() %>%
  right_join(wind %>% filter(country %in% unique(out_temporal$country)), by = c("country"))

