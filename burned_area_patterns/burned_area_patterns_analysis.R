# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(circular) # density.circular, for aspect

library(rgdal)       # readOGR
library(rgeos)       # gIntersection, gArea...
library(maptools)    # UnionSpatialPolygons


# Data --------------------------------------------------------------------

burned_veg <- read.csv("data_burned area by vegetation type_valdivian ecorregion.csv")
spatial_vars <- read.csv("data_spatial_variables.csv")
ndvi_raw <- read.csv("data_ndvi_by_vegetation.csv")
colnames(ndvi_raw)[3] <- "vegetation_code"
ndvi <- left_join(ndvi_raw, 
                  burned_veg[, c("vegetation_code", "vegetation_class")], 
                  by = "vegetation_code")
head(ndvi)

veg_levels <- c("Wet forest",
                "Subalpine forest",
                "Plantation",
                "Dry forest",
                "Shrubland",
                "Anthropogenic prairie and shrubland",
                "Steppe and grassland",
                "Non burnable")


# Burned area by vegetation type ------------------------------------------

burned_veg_long <- pivot_longer(burned_veg[, c("vegetation_class", 
                                               "percent_available_area",
                                               "percent_burned_area")],
                                cols = 2:3, names_to = "class",
                                values_to = "percent_area")
burned_veg_long$class <- plyr::revalue(
  burned_veg_long$class,
  c("percent_available_area" = "Available",
    "percent_burned_area" = "Burned")
)
burned_veg_long <- burned_veg_long[complete.cases(burned_veg_long), ]
burned_veg_long$vegetation_class <- factor(burned_veg_long$vegetation_class,
                                           levels = veg_levels)

ggplot(burned_veg_long, aes(x = vegetation_class, y = percent_area, 
                            colour = class, fill = class)) + 
  geom_bar(stat = "identity", position = position_dodge2(width = 0.5,
                                                         padding = 0.05)) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.5))
ggsave("figure_burned_area_by_veg.png")




# Distribution of spatial variables in burned and available area ----------



spatial_vars_long <- pivot_longer(
  spatial_vars, cols = which(names(spatial_vars) != "class"),
  names_to = "variable", values_to = "value"
)

spatial_vars_long$variable <- factor(spatial_vars_long$variable,
                                     levels = c("ndvi_max", "distance_humans_m",
                                                "distance_roads_m", 
                                                "elevation_m", "slope",
                                                "aspect"),
                                     labels = c("NDVI max", "Distance from settlements (m)",
                                                "Distance from roads (m)",
                                                "Elevation (masl)", "Slope (°)",
                                                "Aspect (°)"))


ggplot(spatial_vars_long, aes(x = value, colour = class, fill = class)) + 
  geom_density(alpha = 0.2) + 
  facet_wrap(vars(variable), scales = "free") + 
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  xlab("variable")
ggsave("figure_spatial_variables.png")



# Aspect, circular density -----------------------------------------------

filter_b <- spatial_vars$class == "burned"
filter_ub <- spatial_vars$class == "available"

aspect_burned <- circular(spatial_vars$aspect[filter_b], 
                          type = "angles", units = "degrees",
                          template = "geographic")
aspect_unburned <- circular(spatial_vars$aspect[filter_ub], 
                            type = "angles", units = "degrees",
                            template = "geographic")

bw_circ <- 30 # smaller = smoother

d_burned <- density.circular(aspect_burned, bw = bw_circ)
d_unburned <- density.circular(aspect_unburned, bw = bw_circ)

x_mine <- seq(0, 360, length.out = length(d_unburned$x))

aspect_plot <- rbind(
  data.frame(aspect = x_mine, density = d_unburned$y),
  data.frame(aspect = x_mine, density = d_burned$y)
)
aspect_plot$class <- rep(c("available", "burned"), each = length(d_unburned$x))

# order (not necessary)
# aspect_plot <- aspect_plot[order(aspect_plot$type, 
#                                  aspect_plot$aspect, 
#                                  aspect_plot$density), ]

ggplot(aspect_plot, aes(x = aspect, y = density, colour = class,
                        fill = class)) + 
  geom_line() + 
  coord_polar() + 
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_x_continuous(breaks = c(0, 90, 180, 270), 
                     labels = c("N", "E", "S", "W")) + 
  theme(legend.title = element_blank())
ggsave("figure_aspect.png")




# NDVI distribution by vegetation class -----------------------------------

ndvi$vegetation_class <- factor(ndvi$vegetation_class, levels = veg_levels)

ggplot(ndvi[ndvi$vegetation_code > 1, ], 
       aes(x = ndvi_max, colour = class, fill = class)) +
  geom_density(alpha = 0.5) + 
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  facet_wrap(vars(vegetation_class), ncol = 2) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.118)) 
ggsave("figure_ndvi by vegetation and burned.jpeg", width = 14, height = 22,
       units = "cm")

boxplot(ndvi_max ~ vegetation_class, data = ndvi)




# Climate analysis --------------------------------------------------------

# Clip fires to study area to compute burned area only in the region of interest.

fires <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/patagonian_fires/patagonian_fires.shp")
study_area <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/study_area/study_area.shp")

plot(fires[fires@data$year == 1999, ])

# Hay problemas en la base de datos
# San Ramon figura como siendo de 2015
