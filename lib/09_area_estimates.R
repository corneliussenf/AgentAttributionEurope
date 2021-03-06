
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(patchwork)

# Get data ----------------------------------------------------------------

load(file = "temp/out_temporal.RData")

forest_disturbance_rates_eur <- read_csv("data/senfetal2021/annual_rates_europe_summary.csv")

forest_area_europe <-  227000000 # Accroding to State of Europe's forest 2020

# Calculate rate and absolute disturbance area ----------------------------

rates_and_areas <- out_temporal %>%
  filter(!is.na(agent)) %>%
  group_by(agent, year) %>%
  summarize(area = sum(area)) %>%
  group_by(year) %>%
  mutate(share = area / sum(area)) %>%
  ungroup() %>%
  filter(agent != "background") %>%
  mutate(agent = ifelse(agent == "breakage", "Storm", "Fire")) %>%
  left_join(forest_disturbance_rates_eur) %>%
  mutate(rate_est = share * mean,
         rate_lower = share * (mean - sd),
         rate_upper = share * (mean + sd)) %>%
  dplyr::select(agent, year, rate_est, rate_lower, rate_upper) %>%
  mutate(area_est = (rate_est / 100) * forest_area_europe,
         area_lower = (rate_lower / 100) * forest_area_europe,
         area_upper = (rate_upper / 100) * forest_area_europe)

rates_and_areas %>%
  group_by(agent) %>%
  summarize(rate_total_est = sum(rate_est),
            rate_total_lower = sum(rate_lower),
            rate_total_upper = sum(rate_upper),
            area_total_est = sum(area_est),
            area_total_lower = sum(area_lower),
            area_total_upper = sum(area_upper),
            rate_mean_est = mean(rate_est),
            rate_mean_lower = mean(rate_lower),
            rate_mean_upper = mean(rate_upper),
            area_mean_est = mean(area_est),
            area_mean_lower = mean(area_lower),
            area_mean_upper = mean(area_upper)) %>%
  View(.)

p <- rates_and_areas %>%
  mutate(period = ifelse(year < 2002, "1986-2001", "2002-2016")) %>%
  group_by(agent, period) %>%
  summarize(rate_total_est = sum(rate_est),
            rate_total_lower = sum(rate_lower),
            rate_total_upper = sum(rate_upper),
            area_total_est = sum(area_est),
            area_total_lower = sum(area_lower),
            area_total_upper = sum(area_upper),
            rate_mean_est = mean(rate_est),
            rate_mean_lower = mean(rate_lower),
            rate_mean_upper = mean(rate_upper),
            area_mean_est = mean(area_est),
            area_mean_lower = mean(area_lower),
            area_mean_upper = mean(area_upper)) %>%
  ggplot(.) +
  geom_point(aes(x = period, 
                 y = area_total_est / 10^6, 
                 col = agent),
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(x = period, 
                    ymin = area_total_lower / 10^6, 
                    ymax = area_total_upper / 10^6, 
                    col = agent),
                width = 0.1, position = position_dodge(width = 0.2)) +
  geom_line(aes(x = period, 
                y = area_total_est / 10^6, 
                col = agent, 
                group = agent),
            position = position_dodge(width = 0.2)) +
  theme_linedraw() +
  labs(x = "Observation period", 
       y = bquote("Disturbance area ("*10^6~"ha)"),
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11)) +
  scale_color_manual(values = c("#BB5566", "#004488"))

ggsave("figures/are_estimate.pdf", p, width = 3.5, height = 3.5)
