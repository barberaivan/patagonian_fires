library(tidyverse); theme_set(theme_bw())
library(viridis)
# Topography --------------------------------------------------------------

tmarg <- read.csv("topo_marginal.csv")
tburn <- read.csv("topo_burned.csv")

for(i in 2:4) {
  plot(density(tmarg[, i]), main = names(tmarg)[i])
  lines(density(tburn[, i]), col = 2)
}


# Vegetation and NDVI -----------------------------------------------------

dndvi <- read.csv("vegetation_ndvi_data.csv")
head(dndvi)
enes <- aggregate(ndvi_max ~ burned + vegetation, dndvi, length) # OK
dndvi <- dndvi[dndvi$vegetation != 1, ]
dndvi$fire <- factor(dndvi$burned, levels = c(0, 1), labels = c("unburned", "burned"))
dndvi$class <- factor(dndvi$vegetation, labels = c(
  "Subalpine forest", "Wet forest", "Dry forest", 
  "Shrubland", "Grassland", 
  "Anthropogenic prairie and shrubland", "Plantation"
), levels = 2:8)

ggplot(dndvi, aes(x = ndvi_max, colour = fire, fill = fire)) +
  geom_density(alpha = 0.5) + 
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  facet_wrap(vars(class), ncol = 2) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.118)) 
ggsave("figure_ndvi by vegetation and burned.jpeg", width = 14, height = 22,
       units = "cm")

boxplot(ndvi_max ~ class, data = dndvi)

# Falta hacer el marginal, mezclando todas las vegetaciones.
# Pensar más qué conviene mostrar, si las distris marginales o las burned y unburned. 