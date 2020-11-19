
# Libraries ---------------------------------------------------------------

library(raster)
library(sf)
library(tidyverse)

# Load reference data -----------------------------------------------------

ref_pts <- list.files("data/referencepoint", pattern = ".shp$", full.names = TRUE) %>%
  map(read_sf) %>%
  map2(.x = ., .y = c("biotic", "breakage", "fire", "landusechange"), 
       ~ mutate(., 
                id = 1:length(id),
                agent = .y)) %>%
  reduce(., sf:::rbind.sf)

# Map reference data ------------------------------------------------------

countries <- read_sf("data/admin/countries_europe_simplyfied.shp")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, projection(countries))
st_crs(world) <- st_crs(countries)
world <- st_crop(world, st_bbox(countries) + c(-0.05, -0.01, 0.01, 0.01) * as.double(st_bbox(countries)))

# p <- ggplot() +
#   geom_sf(data = world, color = NA, fill = "lightgray") +
#   geom_sf(data = ref_pts %>% filter(!agent %in% c("biotic", "landusechange")), 
#           aes(col = str_to_title(agent)),
#           size = 0.75) +
#   geom_sf(data = world, color = "black", fill = NA) +
#   theme_linedraw() +
#   labs(x = NULL, y = NULL, col = "Disturbance agent") +
#   theme(panel.spacing = unit(0, "cm"),
#         panel.background = element_rect(fill = "#d1e5f0"),
#         legend.position = c(0, 1),
#         legend.justification = c(-0.05, 1.1),
#         legend.box.background = element_rect(colour = "black"),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10, face = "bold"),
#         legend.key.height = unit(0.3, "cm"),
#         legend.key.width = unit(0.2, "cm"),
#         strip.background = element_blank(),
#         strip.text = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         plot.margin = unit(rep(0.5, 4), "cm")) +
#   coord_sf(expand = FALSE) +
#   scale_color_manual(values = c("#004488", "#BB5566"))

p <- ggplot() +
  geom_sf(data = countries, color = "black", fill = "grey") +
  geom_sf(data = ref_pts %>% 
            filter(!agent %in% c("biotic", "landusechange")) %>%
            mutate(agent = ifelse(agent == "breakage", "storm", agent)), 
          aes(col = str_to_title(agent)),
          size = 0.55) +
  labs(x = NULL, y = NULL, col = NULL) +
  theme_void() +
  theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        plot.title = element_text(size = 9),
        legend.key.size = unit(0.1, "cm"),
        legend.position = c(0, 1),
        legend.justification = c(0, 1)) +
  coord_sf(expand = FALSE) +
  scale_color_manual(values = c("#BB5566", "#004488"))

ggsave("figures/map_reference_points.pdf", p, width = 3.5, height = 3.5)

# Export by country ------------------------------------------------------

countries_master <- read_csv("data/countries_master.csv")

for (cntr_iso in unique(countries$ISO_CC)) {
  
  cntr <- countries_master[countries_master$iso_code ==  cntr_iso, "country_name"][[1]]
  
  print(cntr)
  
  cntr_tmp <- countries %>% filter(ISO_CC == cntr_iso)
  ref_pts_tmp <- st_intersection(ref_pts, cntr_tmp)
  ref_pts_tmp <- dplyr::select(ref_pts_tmp, -AREA_km2, -COUNTRY, -ISO_CC)
  
  if (nrow(ref_pts_tmp) > 0) {
    ref_pts_tmp$country <- cntr
    write_sf(ref_pts_tmp, paste0("data/referencepoint/by_country/ref_points_", cntr, ".shp"), delete_layer = TRUE)
  }
  
}
