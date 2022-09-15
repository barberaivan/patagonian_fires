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
# ggsave("figure_burned_area_by_veg.png")




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
# ggsave("figure_spatial_variables.png")



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
# ggsave("figure_aspect.png")




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
# ggsave("figure_ndvi by vegetation and burned.jpeg", width = 14, height = 22,
#        units = "cm")

# boxplot(ndvi_max ~ vegetation_class, data = ndvi)




# Fire-based analysis -----------------------------------------------------

# Clip fires to study area to compute burned area only in the region of interest.

fires_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/patagonian_fires/patagonian_fires.shp")
study_area_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/study_area/study_area.shp")

proj_posgar <- "+proj=tmerc +lat_0=-90 +lon_0=-72 +k=1 +x_0=1500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
fires_posgar <- spTransform(fires_wgs, proj_posgar)
sa_posgar <- spTransform(study_area_wgs, proj_posgar)

fires_clipped <- do.call("rbind", lapply(1:nrow(fires), function(i) {
  print(i)
  fire_clip <- gIntersection(polygons(fires_posgar[i, ]), polygons(sa_posgar))
  row.names(fire_clip) <- row.names(fires_posgar[i, ])
  fire_df <- SpatialPolygonsDataFrame(fire_clip, 
                                      data = fires_posgar@data[i, , drop = F])
  # recompute area
  fire_df$area_ha <- gArea(polygons(fire_clip)) * 0.0001 # m2 to ha
  
  return(fire_df)
}))

# plot(sa_posgar)
# plot(fires_clipped, add = TRUE)


# Burned area and number of fires by year ----------------------------------

sum_length <- function(x) c(sum = sum(x), lenght = length(x))
burned_annual <- do.call("data.frame",
                         aggregate(area_ha ~ year, data = fires_clipped@data,
                                   FUN = sum_length))
colnames(burned_annual) <- c("year", "area_ha", "fires")


# Burned area by year and vegetation class ---------------------------------

veg_map <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/vegetation_valdivian_reclassified/vegetation_valdivian_reclassified.shp")
vegmap_posgar <- spTransform(veg_map, proj_posgar)
vegmap_posgar <- vegmap_posgar[!is.na(vegmap_posgar$clss_nm), ]

burned_annual_veg <- expand.grid(
  year = 1999:2022,
  vegetation_code = 2:8, # without non-burnable
  area_ha = 0
)

# for(i in 1:nrow(burned_annual_veg)) {
#   print(i)
#   # i = 3
#   fires_year <- fires_clipped[fires_clipped$year == burned_annual_veg$year[i], ]
#   veg_local <- vegmap_posgar[vegmap_posgar$clss_nm == burned_annual_veg$vegetation_code[i], ]
#   
#   burned_veg <- gIntersection(polygons(fires_year), polygons(veg_local))
#   if(is.null(burned_veg)) {burned_annual_veg$area_ha[i] <- 0} else {
#     burned_annual_veg$area_ha[i] <- gArea(burned_veg) * 0.0001 # m2 to ha
#   }
# } # takes long time

# bring vegetation data
burned_annual_veg <- left_join(
  burned_annual_veg[, c("year", "vegetation_code", "area_ha")], 
  burned_veg[, c("vegetation_code", "vegetation_class", "area_ha_available",
                 "area_ha_burned")],
  by = "vegetation_code"
)

# add burned area by year
burned_annual_sum <- aggregate(area_ha ~ year, burned_annual_veg, sum)
names(burned_annual_sum)[2] <- "area_ha_annual"

burned_annual_veg <- left_join(
  burned_annual_veg, 
  burned_annual_sum[, c("year", "area_ha_annual")],
  by = "year"
)

# burned area / available area by year and vegetation
burned_annual_veg$proportion_burned_annual <- burned_annual_veg$area_ha / 
                                                burned_annual_veg$area_ha_available

# burned area / total burned area in a given year, by year and vegetation
burned_annual_veg$burn_distribution <- burned_annual_veg$area_ha / 
                                         burned_annual_veg$area_ha_annual
# sum to 1?
# aggregate(burn_distribution ~ year, burned_annual_veg, sum)


# write.csv(burned_annual_veg, "data_burned_area_year_veg.csv")
burned_annual_veg <- read.csv("data_burned_area_year_veg.csv")[, -1]

ggplot(burned_annual_veg, aes(year, y = proportion_burned_annual)) +
  geom_point() + geom_line() +
  facet_wrap(vars(vegetation_class))

ggplot(burned_annual_veg, aes(year, y = burn_distribution)) +
  geom_point() + geom_line() +
  facet_wrap(vars(vegetation_class))



# Burn distribution as a function of climate ------------------------------

climate_long <- read.csv("data_climate_interannual.csv")
climate <- pivot_wider(climate_long[, c("variable", "value", "year")], 
                       names_from = "variable", 
                       values_from = "value")

burned_annual_veg_clim <- left_join(
  burned_annual_veg, climate, by = "year"  
)


ggplot(burned_annual_veg_clim, aes(x = vpd, y = burn_distribution, 
                              colour = vegetation_class)) +
  geom_point() + geom_line()



ggplot(burned_annual_veg_clim, aes(x = pp, y = proportion_burned_annual)) +
  geom_point() + geom_smooth() + 
  facet_wrap(vars(vegetation_class), scales = "free")



# Burned area as a function of climate ------------------------------------

burned_annual_clim <- left_join(
  burned_annual, climate, by = "year"  
)
burned_annual_clim_long <- pivot_longer(
  burned_annual_clim, which(names(burned_annual_clim) %in% c("pp", "temp", 
                                                             "vpd", "wind")),
  names_to = "clim_var", values_to = "clim_value"
)

ggplot(burned_annual_clim_long, 
       aes(x = clim_value, y = area_ha_annual)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              method.args = list(family = Gamma(link = "log"))) +
  geom_point() + 
  facet_wrap(vars(clim_var), scales = "free")

  