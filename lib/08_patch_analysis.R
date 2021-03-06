
library(tidyverse)
library(patchwork)

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
  
  out[[k]] <- patches %>%
    sample_frac(., 0.01)
  
}

out_patches <- out  %>%
  bind_rows()

save(out_patches, file = "temp/out_patches.RData")
load(file = "temp/out_patches.RData")

# Patch analysis ----------------------------------------------------------

patch_size_distribution <- out_patches %>%
  mutate(area_cut = cut(area * 0.0001, seq(0, 3535, 1), labels = seq(0.5, 3534.5, 1))) %>%
  mutate(area_cut = as.integer(area_cut)) %>%
  group_by(patchsize = area_cut, agent = prediction_class) %>%
  summarise(n = n()) %>%
  group_by(agent) %>%
  mutate(p = n / sum(n),
         pcum = cumsum(p)) %>%
  ungroup()

fits <- patch_size_distribution %>%
  split(.$agent) %>%
  map(., ~ lm(log(p) ~ log(patchsize), data = .))

power_fits <- data.frame(patchsize = seq(1, max(patch_size_distribution$patchsize), length.out = 100))
power_fits$background <- predict(fits$background, newdata = power_fits)
power_fits$breakage <- predict(fits$breakage, newdata = power_fits)
power_fits$fire <- predict(fits$fire, newdata = power_fits)
power_fits <- power_fits %>%
  gather(key = agent, value = value, -patchsize)

p <- ggplot(patch_size_distribution, 
            aes(x = log(patchsize), y = log(p), col = str_to_title(agent))) +
  # geom_line(data = power_fits, 
  #           aes(x = log(patchsize), y = value, col = str_to_title(agent))) +
  see::geom_point2(alpha = 0.75) +
  theme_linedraw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Patch size (log-scaled)", 
       y = "Relative frequency (log-scaled)",
       title = "Patch size distributions",
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11)) +
  scale_color_manual(values = c("grey", "#004488", "#BB5566"))

ggsave("figures/figure05.pdf", p, width = 3.5, height = 3.5)

ggplot(patch_size_distribution, aes(x = (patchsize), y = pcum, col = str_to_title(agent))) +
  see::geom_point2(alpha = 0.5) +
  theme_linedraw() +
  scale_x_log10(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Patch size", 
       y = "Empirical CDF",
       title = "Patch size distributions",
       col = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 11)) +
  scale_color_manual(values = c("grey", "#004488", "#BB5566"))
