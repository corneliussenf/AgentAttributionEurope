
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

# Create summary files ----------------------------------------------------

cntrs <- list.files("results/predictions") %>%
  strsplit(., "_") %>%
  map(., ~ .[[2]]) %>%
  gsub(".csv", "", .)

out <- vector("list", length(cntrs))

k <- 0

for (cntr in cntrs) {
  
  print(cntr)
  
  k <- k + 1
  
  patches <- read_csv(paste0("results/predictions/prediction_", cntr, ".csv"))
  
  patches <- patches %>%
    mutate(prediction_minimal = case_when(
      prediction == 2 & prediction_breakage_prob > 0.5 ~ 2,
      prediction == 3 & prediction_fire_prob > 0.5 ~ 3,
      TRUE ~ 1
    ))
  
  labels <- c("background", "breakage", "fire")[sort(unique(patches$prediction_minimal))]  
  
  patches <- patches %>%
    mutate(prediction_class_minimal = factor(prediction_minimal, labels = labels))
  
  out[[k]] <- patches %>%
    group_by(year, agent = prediction_class_minimal) %>%
    summarise(area = sum(area)) %>%
    ungroup() %>%
    mutate(country = cntr)
  
}

out_temporal <- out  %>%
  bind_rows()

save(out_temporal, file = "temp/out_temporal.RData")
load(file = "temp/out_temporal.RData")

# Prevalence by country ---------------------------------------------------

p <- out_temporal %>%
  filter(!is.na(agent)) %>%
  group_by(agent, country) %>%
  summarize(area = sum(area)) %>%
  group_by(country) %>%
  mutate(share = area / sum(area)) %>%
  ungroup() %>%
  filter(agent != "background") %>%
  mutate(country = str_to_title(country_renamer(country))) %>%
  mutate(agent = ifelse(agent == "breakage", "Storm", "Fire")) %>%
  ggplot(., aes(x = reorder(country, share, sum), y = share * 100, fill = agent)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme_linedraw() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, 
       y = "Prevalence for 1986-2016 (%)",
       fill = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11)) +
  scale_fill_manual(values = c("#BB5566", "#004488"))

ggsave("figures/prevalence_country_1986-2016.pdf", p, width = 7.5, height = 3.5)

# Temporal analysis ----------------------------------------------------------

p1 <- out_temporal %>%
  filter(!is.na(agent)) %>%
  group_by(agent, year) %>%
  summarize(area = sum(area)) %>%
  group_by(year) %>%
  mutate(share = area / sum(area)) %>%
  ungroup() %>%
  filter(agent != "background") %>%
  mutate(agent = ifelse(agent == "breakage", "Storm", "Fire")) %>%
  ggplot(., aes(x = year, y = share * 100, col = agent)) +
  geom_line(size = 1.25) +
  theme_classic() +
  theme_linedraw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Year", 
       y = "Annual prevalence (%)",
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11)) +
  scale_color_manual(values = c("#BB5566", "#004488"))

out_temporal %>%
  filter(!is.na(agent)) %>%
  group_by(agent, year) %>%
  summarize(area = sum(area)) %>%
  group_by(year) %>%
  mutate(share = area / sum(area)) %>%
  ungroup() %>%
  filter(agent != "background") %>%
  mutate(agent = ifelse(agent == "breakage", "Storm", "Fire")) %>%
  split(.$agent) %>%
  map(., ~ trend::sens.slope(.$share, conf.level = 0.95))

trends <- out_temporal %>%
  filter(!is.na(agent)) %>%
  group_by(country, agent, year) %>%
  summarize(area = sum(area)) %>%
  group_by(year, country) %>%
  mutate(share = area / sum(area)) %>%
  ungroup() %>%
  filter(agent != "background") %>%
  right_join(expand_grid(year = 1986:2016,
                         agent = c("breakage", "fire"),
                         country = cntrs)) %>%
  mutate(share = ifelse(is.na(share), 0, share)) %>%
  mutate(agent = ifelse(agent == "breakage", "Storm", "Fire")) %>%
  split(list(.$agent, .$country)) %>%
  map(., ~ trend::sens.slope(.$share, conf.level = 0.95)) %>%
  map(., ~ data.frame(estimate = .$estimates, z = .$statistic, p = .$p.value)) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("agent", "country"), "\\.") %>%
  mutate(country = country_renamer(country)) %>%
  mutate(trend_label = ifelse(p < 0.01, round(estimate * 100, 2), NA))

write_csv(trends, "results/prevalence_trends_country.csv")

p2 <- out_temporal %>%
  filter(!is.na(agent)) %>%
  group_by(agent, year, country) %>%
  summarize(area = sum(area)) %>%
  group_by(year, country) %>%
  mutate(share = area / sum(area)) %>%
  ungroup() %>%
  filter(agent != "background") %>%
  right_join(expand_grid(year = 1986:2016, country = unique(out_temporal$country), agent = c("breakage", "fire"))) %>%
  mutate(share = ifelse(is.na(share), 0, share)) %>%
  mutate(country = country_renamer(country)) %>%
  mutate(agent = ifelse(agent == "breakage", "Storm", "Fire")) %>%
  ggplot(.) +
  geom_line(aes(x = year, y = share * 100, col = agent, group = interaction(agent, country)), size = 1.25) +
  theme_linedraw() +
  scale_x_continuous(expand = c(0, 0), breaks = c(1990, 2000, 2010)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Year", 
       y = "Annual prevalence (%)",
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11)) +
  scale_color_manual(values = c("#BB5566", "#004488")) +
  facet_wrap(~country, ncol = 5) +
  geom_text(aes(x = 1986.5, y = 80, label = country), show.legend = FALSE, check_overlap = TRUE, hjust = 0, size = 2.5)
  # geom_text(data = trends,
  #           aes(x = 1986.5, y = 70, label = trend_label, col = agent), 
  #           show.legend = FALSE, check_overlap = TRUE, hjust = 0, size = 2.5)

ggsave("figures/annual_prevalence_europe.pdf", p1, width = 3.5, height = 3.5)
ggsave("figures/annual_prevalence_country.pdf", p2, width = 7.5, height = 7)

annual_prevalences <- out_temporal %>%
  filter(!is.na(agent)) %>%
  group_by(country, agent, year) %>%
  summarize(area = sum(area)) %>%
  group_by(year, country) %>%
  mutate(share = area / sum(area)) %>%
  ungroup() %>%
  filter(agent != "background") %>%
  right_join(expand_grid(year = 1986:2016,
                         agent = c("breakage", "fire"),
                         country = cntrs)) %>%
  mutate(share = ifelse(is.na(share), 0, share)) %>%
  mutate(agent = ifelse(agent == "breakage", "Storm", "Fire")) %>%
  mutate(country = country_renamer(country))

annual_prevalences_summary <- annual_prevalences %>% 
  mutate(period = ifelse(year <= 2001, "early (1986-2001)", "late (2002-2016)")) %>%
  group_by(period, country, agent) %>%
  summarize(median = median(share)) %>%
  spread(key = period, value = median) %>%
  right_join(annual_prevalences %>% 
              group_by(country, agent) %>%
              summarise(`full (1986-2016)` = median(share)))

write_csv(annual_prevalences_summary, "results/prevalence_country.csv")

