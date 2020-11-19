
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

fire <- read_csv("data/EFFIS/eudb-ba.csv") %>%
  gather(key = country, value = effis_buntarea_ha, -year)

dat <- out_temporal %>% 
  filter(agent == "fire") %>%
  mutate(area = area / 10000) %>%
  dplyr::select(-agent) %>%
  dplyr::rename(landsat_burntarea_ha = area) %>%
  right_join(fire %>% filter(country %in% unique(out_temporal$country)), by = c("year", "country")) %>%
  mutate(landsat_burntarea_ha = ifelse(is.na(landsat_burntarea_ha), 0, landsat_burntarea_ha)) %>%
  filter(year %in% 1986:2016)

# Comparison --------------------------------------------------------------

dat %>%
  group_by(country) %>%
  summarize(n = sum(!is.na(effis_buntarea_ha)),
            effis_buntarea_ha = sum(effis_buntarea_ha, na.rm = TRUE),
            landsat_burntarea_ha = sum(landsat_burntarea_ha, na.rm = TRUE)) %>%
  mutate(effis_buntarea_ha = effis_buntarea_ha / n,
         landsat_burntarea_ha = landsat_burntarea_ha / n) %>%
  filter(landsat_burntarea_ha > 0) %>%
  lm(landsat_burntarea_ha~effis_buntarea_ha, data = .) %>%
  summary(.)

dat %>%
  group_by(country) %>%
  summarize(n = sum(!is.na(effis_buntarea_ha)),
            effis_buntarea_ha = sum(effis_buntarea_ha, na.rm = TRUE),
            landsat_burntarea_ha = sum(landsat_burntarea_ha, na.rm = TRUE)) %>%
  mutate(effis_buntarea_ha = effis_buntarea_ha / n,
         landsat_burntarea_ha = landsat_burntarea_ha / n) %>%
  filter(landsat_burntarea_ha > 0) %>%
  summarize(mean(landsat_burntarea_ha - effis_buntarea_ha))

p <- dat %>%
  group_by(country) %>%
  summarize(n = sum(!is.na(effis_buntarea_ha)),
            effis_buntarea_ha = sum(effis_buntarea_ha, na.rm = TRUE),
            landsat_burntarea_ha = sum(landsat_burntarea_ha, na.rm = TRUE)) %>%
  mutate(effis_buntarea_ha = effis_buntarea_ha / n,
            landsat_burntarea_ha = landsat_burntarea_ha / n) %>%
  filter(landsat_burntarea_ha > 0) %>%
  ggplot(., aes(x = effis_buntarea_ha, y = landsat_burntarea_ha)) +
  geom_smooth(method = "lm", se = FALSE, col = "grey",
              size = 0.75) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = country_renamer(country)), size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_log10(limits = c(1, 10000000)) +
  scale_y_log10(limits = c(1, 10000000)) +
  labs(x = "EFFIS annual area burnt (ha)",
       y = "Landsat annual area burnt (ha)") +
  theme_linedraw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11)) +
  annotate("text", x = 1, y = 1e07, size = 2.5,
           label = "y = -66.16 + 0.17 * x",
           hjust = 0) +
  annotate("text", x = 1, y = 4e06, size = 2.5,
           label = bquote("R^2==0.72"),
           hjust = 0, parse = TRUE) +
  annotate("text", x = 1, y = 1.5e06, size = 2.5,
           label = "Bias = -19,186 ha",
           hjust = 0)

ggsave("figures/effis_vs_landsat_burnt_area.pdf", p, height = 3.5, width = 3.5)

dat %>%
  group_by(country) %>%
  summarize(effis_buntarea_ha = sum(effis_buntarea_ha, na.rm = TRUE),
            landsat_burntarea_ha = sum(landsat_burntarea_ha, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(pearson = cor(effis_buntarea_ha, landsat_burntarea_ha, method = "pearson"),
            spearman = cor(effis_buntarea_ha, landsat_burntarea_ha, method = "spearman"))

dat %>%
  mutate(country = str_to_title(country)) %>%
  ggplot() +
  geom_line(aes(x = year, y = effis_buntarea_ha), col = "#BB5566") +
  geom_line(aes(x = year, y = landsat_burntarea_ha), col = "#BB5566", linetype = "dashed") +
  facet_wrap(~country, ncol = 4, scales = "free") +
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
