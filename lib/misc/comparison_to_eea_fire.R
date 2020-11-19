
library(tidyverse)
library(patchwork)

load(file = "temp/out_temporal.RData")

fire <- read_csv("data/eea_fire/eea_burned_area.csv")
fire <- fire %>%
  dplyr::rename(total_eumed5 = `Total EUMED5`,
                other_countries = `Other countries`,
                year = Year) %>%
  gather(key = country, value = fire_area, -year, -total_eumed5, -other_countries) %>%
  mutate(country = tolower(country))

dat <- fire %>%
  left_join(out_temporal %>% filter(agent == "fire"), by = c("year", "country"))

ggplot(dat, aes(x = fire_area, y = area / 10000, col = country)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1)

ggplot(dat) +
  geom_point(aes(x = fire_area, y = area / 10000, col = country)) +
  # geom_point(aes(x = fire_area, 
  #                y = zoo::rollmean(area * 0.0001, k = 3, na.pad = TRUE), col = country), 
  #            shape = 2) +
  geom_smooth(aes(x = fire_area, y = area / 10000, col = country), 
              method = "lm", se = FALSE) +
  # geom_smooth(aes(x = fire_area, 
  #                 y = zoo::rollmean(area * 0.0001, k = 3, na.pad = TRUE), col = country), 
  #             method = "lm", se = FALSE, linetype = 2) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~country, scales = "free", ncol = 5)

dat %>%
  mutate(country = str_to_title(country)) %>%
  filter(year %in% 1986:2016) %>%
  ggplot() +
  geom_line(aes(x = year, y = fire_area), col = "#BB5566") +
  geom_line(aes(x = year, y = area / 10000), col = "#BB5566", linetype = "dashed") +
  facet_wrap(~country, scales = "free", ncol = 5) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Year", 
       y = "Burned area (ha)",
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11),
        strip.background = element_blank())

dat %>%
  group_by(country) %>%
  summarize(effis_buntarea_ha = sum(fire_area, na.rm = TRUE),
            landsat_burntarea_ha = sum(area / 10000, na.rm = TRUE)) %>%
  ggplot(., aes(x = effis_buntarea_ha, y = landsat_burntarea_ha)) +
  geom_text(aes(label = country)) +
  geom_abline(intercept = 0, slope = 1)

dat %>%
  group_by(country) %>%
  summarize(effis_buntarea_ha = sum(fire_area, na.rm = TRUE),
            landsat_burntarea_ha = sum(area / 10000, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(pearson = cor(effis_buntarea_ha, landsat_burntarea_ha, method = "pearson"),
            spearman = cor(effis_buntarea_ha, landsat_burntarea_ha, method = "spearman"))

p <- dat %>%
  mutate(country = str_to_title(country)) %>%
  ggplot() +
  geom_line(aes(x = year, y = scale(fire_area)), col = "#BB5566") +
  geom_line(aes(x = year, y = scale(area)), col = "#BB5566", linetype = "dashed") +
  facet_wrap(~country, scales = "free", ncol = 5) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Year", 
       y = "Burned area (ha)",
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11),
        strip.background = element_blank())

ggsave("figures/comparison_fire_europe.png", p, width = 9, height = 2.5)

p <- dat %>%
  mutate(country = str_to_title(country)) %>%
  filter(year >= 1986 & year < 2017) %>%
  group_by(country) %>%
  mutate(fire_area_cum = cumsum(fire_area),
         area_cum = cumsum(area)) %>%
  mutate(fire_area_cum_p = fire_area_cum / sum(fire_area),
         area_cum_p = area_cum / sum(area)) %>%
  ggplot() +
  geom_step(aes(x = year, y = fire_area_cum_p), col = "#BB5566") +
  geom_step(aes(x = year, y = area_cum_p), col = "#BB5566", linetype = "dashed") +
  geom_bar(aes(x = year, y = fire_area_cum_p - area_cum_p), stat = "identity", fill = "grey", col = NA) +
  facet_wrap(~country, scales = "free", ncol = 5) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Year", 
       y = "Proportion of total area burned",
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11),
        strip.background = element_blank())

ggsave("figures/comparison_fire_europe_cumsum.png", p, width = 11, height = 2.5)

dat %>%
  group_by(year) %>%
  summarise(area = sum(area / 10000, na.rm = TRUE),
            fire_area = sum(fire_area, na.rm = TRUE)) %>%
  filter(year >= 1986 & year < 2017) %>%
  ggplot(., aes(x = area, y = fire_area)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0, 800000) + ylim(0, 800000)


dat %>%
  group_by(country) %>%
  mutate(area = area / 10000,
         area_mean = zoo::rollmean(area, k = 3, na.pad = TRUE),
         area_max = zoo::rollmax(area, k = 3, na.pad = TRUE),
         fire_area_mean = zoo::rollmean(fire_area, k = 3, na.pad = TRUE),
         fire_area_max = zoo::rollmax(fire_area, k = 3, na.pad = TRUE)) %>%
  summarise(r = cor(area, fire_area, use = "complete.obs", method = "spearman"),
            r_mean = cor(area_mean, fire_area, use = "complete.obs", method = "spearman"),
            r_max = cor(area_max, fire_area, use = "complete.obs", method = "spearman"),
            r_mean_mean = cor(area_mean, fire_area_mean, use = "complete.obs", method = "spearman"),
            r_max_max = cor(area_max, fire_area_max, use = "complete.obs", method = "spearman"))

