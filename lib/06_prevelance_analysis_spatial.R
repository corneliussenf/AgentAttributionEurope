
library(tidyverse)
library(data.table)
library(patchwork)
library(sf)

# Create summary files ----------------------------------------------------

grid <- read_sf("data/grids/hexagon_50km.shp")

grid <- grid %>%
  dplyr::rename(gridid = id)

cntrs <- list.files("results/predictions/") %>%
  strsplit(., "_") %>%
  map(., ~ .[[2]]) %>%
  map(., ~strsplit(., "\\.")[[1]]) %>%
  map(., ~ .[[1]]) %>%
  unlist()

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
  
  patch_centers_sf <- st_as_sf(sp::SpatialPoints(coords = as.matrix(patches[, c("x", "y")])), dim = "XY")
  patch_centers_sf$patch <- patches$patch
  patch_centers_sf$agent <- patches$prediction_class_minimal
  patch_centers_sf$area <- patches$area
  patch_centers_sf$year <- patches$year
  
  st_crs(patch_centers_sf) <- st_crs(grid)
  
  grid_at_patches <- st_intersection(grid, patch_centers_sf)
  
  out[[k]] <- grid_at_patches %>%
    st_drop_geometry() %>%
    group_by(gridid, agent, year) %>%
    summarise(area = sum(area)) %>%
    ungroup()
  
}

out_grid <- out  %>%
  bind_rows()

save(out_grid, file = "temp/out_grid.RData")
load(file = "temp/out_grid.RData")

# Map prevalance ----------------------------------------------------------

selector <- read_csv("data/grid_selector.csv") # Select hexagons based on land area (similar to Senf et al. 2020)

prevalence_annual <- expand.grid(gridid = unique(out_grid$gridid),
                          year = 1986:2016,
                          agent = unique(out_grid$agent)) %>%
  as_tibble() %>%
  left_join(out_grid) %>%
  mutate(area = ifelse(is.na(area), 0, area)) %>%
  group_by(gridid, agent, year) %>%
  summarize(area = sum(area, na.rm = TRUE)) %>%
  group_by(gridid, year) %>%
  mutate(prevalence = area / sum(area, na.rm = TRUE)) %>%
  dplyr::select(-area) %>%
  mutate(prevalence = ifelse(is.nan(prevalence), 0, prevalence)) %>%
  spread(key = agent, value = prevalence)

prevalence <- expand.grid(gridid = unique(out_grid$gridid),
                          year = 1986:2016,
                          agent = unique(out_grid$agent)) %>%
  as_tibble() %>%
  left_join(out_grid) %>%
  mutate(area = ifelse(is.na(area), 0, area)) %>%
  group_by(gridid, agent) %>%
  summarize(area = sum(area, na.rm = TRUE)) %>%
  group_by(gridid) %>%
  mutate(prevalence = area / sum(area, na.rm = TRUE)) %>%
  dplyr::select(-area) %>%
  spread(key = agent, value = prevalence)

library(cowplot)

countries_sf <- read_sf("../Recovery/data/admin/countries_europe_simplyfied.shp")

### Overall

prevelance_map <- grid %>%
  right_join(prevalence) %>%
  filter(gridid %in% selector$grid_id)

data <- prevelance_map %>%
  mutate(fire_cut = cut(fire, c(0, 0.25, 0.5, 1), labels = 1:3, include.lowest = TRUE),
         breakage_cut = cut(breakage, c(0, 0.25, 0.5, 1), labels = 1:3, include.lowest = TRUE)) %>%
  mutate(bi_class = paste(fire_cut, breakage_cut, sep = "-"))

map <- ggplot() +
  geom_point(aes(x = mean(st_coordinates(data)[, 1]), 
                 y = mean(st_coordinates(data)[, 2]), 
                 col = "Combination does not exist")) +
  geom_sf(data = data, 
          mapping = aes(fill = bi_class, geometry = geometry), 
          show.legend = FALSE, 
          color = NA) +
  geom_sf(data = countries_sf, aes(geometry = geometry), 
          show.legend = FALSE, 
          fill = NA) +
  scale_color_manual(values = "grey80") +
  scale_fill_manual(values = c("#BBBBBB", 
                               "#00448880", 
                               "#004488", 
                               "#BB556680", 
                               "#BB5566"),
                    na.value = "#BBBBBB") +
  theme_void() +
  theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        legend.position = c(0.085, 0.7325),
        legend.title = element_blank(),
        legend.spacing.x = unit(0, "cm"),
        legend.spacing.y = unit(0, "cm"),
        legend.text = element_text(size = 7),
        legend.justification = c(0, 1)) +
  coord_sf(expand = FALSE)

legend <- ggplot(data = expand.grid(x = c("<25%", "25-50%", ">50%"), 
                          y = c(">50%", "25-50%", "<25%")),
       aes(x = x, y = y, 
           fill = factor(c(5, 6, 6, 4, 6, 6, 1, 2, 3)))) +
  geom_tile(col = "grey20") +
  theme_void() +
  scale_fill_manual(values = c("#BBBBBB", 
                               "#00448880", 
                               "#004488", 
                               "#BB556680", 
                               "#BB5566", 
                               "grey90")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = -1),
        axis.text = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9, angle = 90, vjust = 0.5)) +
  labs(x = "Storm\n", y = "Fire\n", fill = "Prevalence") +
  scale_x_discrete(position = "top") +
  coord_equal()

plot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.025, 0.7, 0.3, 0.3)

ggsave("figures/prevalence_map.pdf", plot, width = 7.5, height = 7.5)

### Annual

prevelance_annual_map <- grid %>%
  right_join(prevalence_annual) %>%
  filter(gridid %in% selector$grid_id)

map_annual_breakage <- ggplot() +
  geom_sf(data = prevelance_annual_map, 
          mapping = aes(fill = breakage, geometry = geometry), 
          show.legend = TRUE, color = NA) +
  geom_sf(data = countries_sf, aes(geometry = geometry), fill = NA) +
  scale_fill_gradientn(colors = c("#BBBBBB", "#00448880", "#004488")) +
  theme_void() +
  theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        strip.text = element_text(hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.margin = ggplot2::margin(0.75, 0.75, 0.75, 0.75),
        legend.box.margin = ggplot2::margin(-10, 0.15, 0.15, 0.15)) +
  coord_sf(expand = FALSE) +
  facet_wrap(~year) +
  labs(title = "Prevalence of storm-related disturbances")

ggsave("figures/prevalence_annual_breakage_map.pdf", map_annual_breakage, width = 7.5, height = 8.5)

map_annual_fire <- ggplot() +
  geom_sf(data = prevelance_annual_map, 
          mapping = aes(fill = fire, geometry = geometry), 
          show.legend = TRUE, color = NA) +
  geom_sf(data = countries_sf, aes(geometry = geometry), fill = NA) +
  scale_fill_gradientn(colors = c("#BBBBBB", "#BB556680", "#BB5566")) +
  theme_void() +
  theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        strip.text = element_text(hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.margin = ggplot2::margin(0.75, 0.75, 0.75, 0.75),
        legend.box.margin = ggplot2::margin(-10, 0.15, 0.15, 0.15)) +
  coord_sf(expand = FALSE) +
  facet_wrap(~year) +
  labs(title = "Prevalence of fire-related disturbances")

ggsave("figures/prevalence_annual_fire_map.pdf", map_annual_fire, width = 7.5, height = 8.5)


