# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(grid); library(gtable)
# library(gridExtra) # conflictúa con gtable?

library(circular) # density.circular, for aspect

library(rgdal)       # readOGR
library(rgeos)       # gIntersection, gArea...
library(maptools)    # UnionSpatialPolygons

# function to share legend --------------------------------------------------
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

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

veg_levels2 <- c("Wet forest",
                 "Subalpine\nforest",
                 "Plantation",
                 "Dry forest",
                 "Shrubland",
                 "Anthropogenic prairie\nand shrubland",
                 "Steppe and\ngrassland",
                 "Non burnable")

# Distribution of spatial variables in burned and available area ----------


# Vegetation type


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
                                           levels = veg_levels,
                                           labels = veg_levels2)
burned_veg_long$variable <- "Vegetation type"

veg_dist <- 
ggplot(burned_veg_long, 
       aes(x = vegetation_class, y = percent_area / 100, 
           colour = class, fill = class)) + 
  geom_bar(stat = "identity", position = position_dodge2(width = 0.5,
                                                         padding = 0.05),
           alpha = 0.2) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.67),
        axis.title.x = element_blank()) +
  facet_wrap(vars(variable)) +
  ylab("Probability mass")
veg_dist
# ggsave("figure_burned_area_by_veg.png")




# Compute density separately for the next plots, to make a facet_wrap later.
l <- 512 * 2 * 5
dlist <- vector(mode = "list", length = 5)
names(dlist) <- c("elevation_m", "slope", "aspect", 
                  "distance_humans_m", "distance_roads_m")

# rescale distances to km
spatial_vars[, c("distance_humans_m", "distance_roads_m")] <-
  spatial_vars[, c("distance_humans_m", "distance_roads_m")] / 1000  


for(i in 1:length(dlist)) {
  
  # i = 1
  v <- names(dlist)[i]
  d_burned <- density(spatial_vars[spatial_vars$class == "burned", v], from = 0)
  d_av <- density(spatial_vars[spatial_vars$class == "available", v], from = 0)
  
  dlist[[i]] <- data.frame(variable = v,
                           value = c(d_burned$x, d_av$x),
                           density = c(d_burned$y, d_av$y),
                           class = factor(rep(c("Burned", "Available"), 
                                              each = length(d_av$x)),
                                          levels = c("Burned", "Available"))
  )
}


# Compute density for aspect (circular distribution)

filter_b <- spatial_vars$class == "burned"
filter_ub <- spatial_vars$class == "available"

aspect_burned <- circular(spatial_vars$aspect[filter_b], 
                          type = "angles", units = "degrees",
                          template = "geographic")
aspect_unburned <- circular(spatial_vars$aspect[filter_ub], 
                            type = "angles", units = "degrees",
                            template = "geographic")

bw_circ <- 20 # smaller = smoother

d_burned <- density.circular(aspect_burned, bw = bw_circ)
d_unburned <- density.circular(aspect_unburned, bw = bw_circ)

x_mine <- seq(0, 360, length.out = length(d_unburned$x))

aspect_data <- rbind(
  data.frame(variable = "aspect", value = x_mine, density = d_unburned$y, class = "Available"),
  data.frame(variable = "aspect", value = x_mine, density = d_burned$y, class = "Burned")
)
aspect_data$class <- factor(aspect_data$class, levels = c("Burned", "Available"))

dlist$aspect <- aspect_data



# merge variables in df
d_data <- do.call("rbind", dlist)
d_data$variable <- plyr::revalue(d_data$variable,
                         c("elevation_m" = "Elevation (m a.s.l.)", 
                           "slope" = "Slope (°)",
                           "aspect" = "Aspect",
                           "distance_humans_m" = "Distance from human\nsettlements (km)", 
                           "distance_roads_m" = "Distance from\nroads (km)"))
d_data$variable <- factor(d_data$variable, levels = 
                          c("Elevation (m a.s.l.)", 
                            "Aspect",
                            "Slope (°)",
                            "Distance from human\nsettlements (km)", 
                            "Distance from\nroads (km)"))


var5plot <- ggplot(d_data, 
                   aes(x = value, y = density, ymin = 0, ymax = density,
                       colour = class, fill = class)) + 
  # geom_line(show.legend = F) + 
  # geom_ribbon(alpha = 0.2, colour = NA) +
  geom_ribbon(alpha = 0.2, size = 0.4) +  
  geom_hline(yintercept = 0, colour = "white", size = 0.45, alpha = 1) +
  facet_wrap(vars(variable), ncol = 3, scales = "free") +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.25)) +
  xlab("Variable") + 
  ylab("Probability density")
var5plot




# Merge this and vegetation type
veg_dist2 <- veg_dist + theme(axis.title.y = element_text(vjust = 7),
                              plot.margin = unit(c(2,2,2,7), "mm"))
# veg_dist2
ggpubr::ggarrange(veg_dist2, 
                  var5plot, 
                  nrow = 2, heights = c(1.3, 2))
ggsave("Figure 03 - spatial variables distributions.jpeg",
       width = 17, height = 21, units = "cm", dpi = 500)






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

  