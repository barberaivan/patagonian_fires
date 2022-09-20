# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(ggpubr)
library(grid)
library(scales)   # log scale 
library(circular) # density.circular, for aspect

library(rgdal)       # readOGR
library(rgeos)       # gIntersection, gArea...
library(maptools)    # UnionSpatialPolygons

# functions ---------------------------------------------------------------

# Normalize
normalize <- function(x) x / sum(x)

# To shaer legend (not used I think)
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# custom ggplot theme -----------------------------------------------------

# from https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_mine <- function() { 
  font <- "Arial"   #assign font family up front
  marg <- 2 # figure margin in mm
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 16,                #set font size
        #face = 'bold',            #bold typeface
        hjust = -0.1,                #left align
        vjust = 1),               
      
      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 12),               
      
      # para separar el eje y de los nros
      axis.title.y = element_text(             
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),
      
      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size
      
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9, family = font),
      
      strip.text = element_text(size = 12, family = font, color = "white"),
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),
      
      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

theme_set(theme_mine())



# Data --------------------------------------------------------------------

burned_veg <- read.csv("data_burned_and_available_area_by_vegetation.csv")
spatial_vars <- read.csv("data_spatial_variables.csv")
ndvi_raw <- read.csv("data_ndvi_by_vegetation.csv")
colnames(ndvi_raw)[colnames(ndvi_raw) == "vegetation_valdivian"] <- "vegetation_code"
ndvi_raw <- ndvi_raw[, c("class", "ndvi_max", "vegetation_code")]
ndvi <- left_join(ndvi_raw, 
                  burned_veg[, c("vegetation_code", "vegetation_class")], 
                  by = "vegetation_code")

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

burned_veg$prob_av <- c(NA, normalize(burned_veg$area_available_ha[-1]))
burned_veg$prob_b <- c(NA, normalize(burned_veg$area_burned_ha[-1]))
veg_dist_data <- data.frame(
  vegetation_class = factor(
    rep(burned_veg$vegetation_class[-1], 2),
    levels = veg_levels, labels = veg_levels2
  ),
  prob = c(burned_veg$prob_av[-1], burned_veg$prob_b[-1]),
  class = factor(
    rep(c("Available", "Burned"), each = nrow(burned_veg) - 1),
    levels = c("Available", "Burned")
    ),
  variable = "Vegetation type"
)

veg_dist <- 
ggplot(veg_dist_data, 
       aes(x = vegetation_class, y = prob, 
           colour = class, fill = class)) + 
  geom_bar(stat = "identity", position = position_dodge2(width = 0.5,
                                                         padding = 0.05),
           alpha = 0.2) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  facet_wrap(vars(variable)) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.67),
        axis.title.x = element_blank()) +
  ylab("Probability mass")
veg_dist
# ggsave("figure_burned_area_by_veg.png")




# Compute density separately for the next plots, to make a facet_wrap later.
l <- 512 * 2 * 5
dlist <- vector(mode = "list", length = 5)
names(dlist) <- c("elevation", "slope", "aspect", 
                  "distance_humans_km", "distance_roads_km")

# rescale distances to km
spatial_vars$distance_humans_km <- spatial_vars$distance_humans / 1000  
spatial_vars$distance_roads_km <- spatial_vars$distance_roads / 1000  

for(i in 1:length(dlist)) {
  
  # i = 1
  v <- names(dlist)[i]
  d_burned <- density(spatial_vars[spatial_vars$class == "Burned", v], from = 0)
  d_av <- density(spatial_vars[spatial_vars$class == "Available", v], from = 0)
  
  dlist[[i]] <- data.frame(variable = v,
                           value = c(d_burned$x, d_av$x),
                           density = c(d_burned$y, d_av$y),
                           class = factor(rep(c("Burned", "Available"), 
                                              each = length(d_av$x)),
                                          levels = c("Burned", "Available"))
  )
}


# Compute density for aspect (circular distribution)

filter_b <- spatial_vars$class == "Burned"
filter_ub <- spatial_vars$class == "Available"

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
                         c("elevation" = "Elevation (m a.s.l.)", 
                           "slope" = "Slope (°)",
                           "aspect" = "Aspect",
                           "distance_humans_km" = "Distance from human\nsettlements (km)", 
                           "distance_roads_km" = "Distance from\nroads (km)"))
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

ndvi_veg <- read.csv("data_ndvi_by_vegetation.csv")
names(ndvi_veg)[names(ndvi_veg) == "vegetation_valdivian"] <- "vegetation_code"
ndvi_veg <- left_join(ndvi_veg,
                      burned_veg[, c("vegetation_code", "vegetation_class")], 
                      by = "vegetation_code")

# ndvi marginal to vegetation type
ndvi_marginal <- spatial_vars[, c("class", "ndvi_max")]
ndvi_marginal$vegetation_class <- "All vegetation types"

ndvi_data <- rbind(ndvi_marginal, 
                   ndvi_veg[ndvi_veg$vegetation_code > 1, 
                            c("class", "ndvi_max", "vegetation_class")])

ndvi_data$vegetation_class <- factor(
  ndvi_data$vegetation_class,
  levels = c(veg_levels, "All vegetation types"),
  labels = c(veg_levels2, "All vegetation\ntypes")
)

ggplot(ndvi_data, 
       aes(x = ndvi_max, colour = class, fill = class)) +
  geom_density(alpha = 0.2) + 
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  facet_wrap(vars(vegetation_class), ncol = 3) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.15),
        legend.title = element_blank()) +
  ylab("Probability density") +
  xlab("NDVI_mean_max")
ggsave("figure_ndvi by vegetation and burned.jpeg", width = 17, height = 15,
       units = "cm")

# boxplot(ndvi_max ~ vegetation_class, data = ndvi)




# Fire-based analysis -----------------------------------------------------

# Clip fires to study area to compute burned area only in the region of interest.

fires_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/patagonian_fires/patagonian_fires.shp")
study_area_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/study_area/study_area.shp")

proj_posgar <- "+proj=tmerc +lat_0=-90 +lon_0=-72 +k=1 +x_0=1500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
fires_posgar <- spTransform(fires_wgs, proj_posgar)
sa_posgar <- spTransform(study_area_wgs, proj_posgar)

fires_clipped <- do.call("rbind", lapply(1:nrow(fires_posgar), function(i) {
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
# write.csv(burned_annual, "data_burned_area_by_year.csv")


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


# bring fwi data
fwi_data <- read.csv("data_climate_interannual_fwi.csv")
fwi_data$date <- as.Date(fwi_data$date, format = "%Y-%m-%d")
fwi_data$year <- format(fwi_data$date, format = "%Y") %>% as.numeric
fwi_data$month <- format(fwi_data$date, format = "%m") %>% as.numeric
fwi_data$fseason <- NA
fwi_data$fseason <- fwi_data$year
fwi_data$fseason[fwi_data$month == 12] <- fwi_data$year[fwi_data$month == 12] + 1

fwi_agg <- aggregate(fwi ~ fseason, fwi_data[fwi_data$month %in% c(12, 1:3), ],
                     FUN = mean)
names(fwi_agg) <- c("year", "fwi")

# remaining climatic variables
climate_long <- read.csv("data_climate_interannual.csv")
climate <- pivot_wider(climate_long[, c("variable", "value", "year")], 
                       names_from = "variable", 
                       values_from = "value")

# merge with fwi
climate <- cbind(climate[climate$year > 1998, ], fwi = fwi_agg$fwi)

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
                                                             "vpd", "wind", "fwi")),
  names_to = "clim_var", values_to = "clim_value"
)

clim_area <- 
ggplot(burned_annual_clim_long, 
       aes(x = clim_value, y = area_ha)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              method.args = list(family = Gamma(link = "log")),
              color = viridis(1, begin = 0.2, option = "B"),
              fill = viridis(1, begin = 0.2, option = "B"),
              alpha = 0.2, size = 0.8) +
  geom_point(shape = 19, alpha = 0.7, size = 2.5) + 
  facet_wrap(vars(clim_var), scales = "free_x", ncol = 1) +
  theme(panel.grid.minor = element_blank(),
        # remove strip
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  xlab("Climatic variable") + 
  ylab("Annual burned area (ha)") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                labels = trans_format("log10", math_format(10 ^ .x)))
clim_area

clim_fires <- 
  ggplot(burned_annual_clim_long, 
         aes(x = clim_value, y = fires)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              method.args = list(family = mgcv::nb(link = "log")),
              color = viridis(1, begin = 0.2, option = "B"),
              fill = viridis(1, begin = 0.2, option = "B"),
              alpha = 0.2, size = 0.8) +
  geom_point(shape = 19, alpha = 0.7, size = 2.5) + 
  facet_wrap(vars(clim_var), scales = "free_x", ncol = 1,
             strip.position = "right") +  
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 270)) +
  xlab("Climatic variable") + 
  ylab("Number of fires")
clim_fires

# con gg arrange (no anda bien el eje)
clim_fig <-
ggarrange(clim_area + theme(axis.title.x = element_blank()), 
          clim_fires + theme(axis.title.x = element_blank()),
          ncol = 2)
annotate_figure(clim_fig, 
                bottom = textGrob("Climatic variable", 
                                  gp = gpar(cex = 1),
                                  hjust = 0.3,
                                  vjust = -0.25))
ggsave("figure_climate and annual burned area.png", width = 15, height = 20,
       units = "cm")
# with number of fires (Problema con facet_wrap)
# even_longer <- pivot_longer(
#   burned_annual_clim_long, 
#   which(names(burned_annual_clim_long) %in% c("area_ha", "fires")),
#   names_to = "fire_var", values_to = "fire_value"
# )
# 
# ggplot(even_longer, 
#        aes(x = clim_value, y = fire_value)) +
#   geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
#               method.args = list(family = Gamma(link = "log"))) +
#   geom_point() + 
#   facet_grid(rows = vars(clim_var), cols = vars(fire_var), 
#              scales = "free") +
#   theme(panel.grid.minor = element_blank())






# Intraannual fire activity -----------------------------------------------

remove <- which(fires_clipped@data$obs == "uncertain_year")
dmonth <- fires_clipped@data[-remove, ]
dmonth$month_num <- format(as.Date(dmonth$date, format = "%Y-%m-%d"), 
                           format = "%m") %>% as.numeric
dmonth$month <- factor(dmonth$month_num, levels = c(7:12, 1:6), 
                       labels = month.name[c(7:12, 1:6)])


### Plot with annual values (failure)


# aggregate by month and year (sum area, count fires)
burned_my <- do.call("data.frame",
                     aggregate(area_ha ~ month_num + month + year,
                               data = dmonth,
                               FUN = sum_length))
colnames(burned_my) <- c("month_num", "month", "year", "area", "fires")

# Average area and fire number by month
burned_monthly <- do.call("data.frame",
                          aggregate(cbind(area, fires) ~ month_num + month, 
                                    data = burned_my,
                                    FUN = mean))
colnames(burned_monthly) <- c("month_num", "month", "area", "fires")
# write.csv(burned_monthly, "data_burned_area_by_month.csv")


# categorical numerical month
burned_my$month_num_cat <- factor(burned_my$month_num, 
                                  levels = c(7:12, 1:6))
burned_monthly$month_num_cat <- factor(burned_monthly$month_num, 
                                       levels = c(7:12, 1:6))

# Longanize monthly and month-yearly
burned_monthly2 <- burned_monthly
burned_my2 <- burned_my
coeff <- max(burned_my$area) / max(burned_my$fires)

burned_monthly2$fires <- burned_monthly$fires * coeff
burned_my2$fires <- burned_my$fires * coeff


burned_monthly_long <- pivot_longer(burned_monthly2, cols = 3:4, 
                                    names_to = "var", values_to = "y")
burned_my_long <- pivot_longer(burned_my2, cols = 3:4, 
                               names_to = "var", values_to = "y")


ggplot(mapping = aes(x = month_num_cat, y = y, colour = var, shape = var, 
                     group = var)) +
  geom_point(data = burned_my_long, size = 3, alpha = 0.7) + 
  geom_line(data = burned_monthly_long) +
  geom_point(data = burned_monthly_long, size = 3, alpha = 0.7) + 
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  scale_y_continuous(
    # Features of the first axis
    name = "Mean monthly burned area (ha)",
    trans = "log10",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . / coeff, name = "Mean number of fires",
                        trans = NULL)
  ) +
  theme(legend.position = "right")


### Plot without annual values (failure)

# Average area and fire number by month
burned_monthly <- do.call("data.frame",
                          aggregate(area_ha ~ month_num + month, 
                                    data = dmonth,
                                    FUN = sum_length))
colnames(burned_monthly) <- c("month_num", "month", "area", "fires")
# categorical numerical month
burned_monthly$month_num_cat <- factor(burned_monthly$month_num, 
                                       levels = c(7:12, 1:6))

# Rescale for common axis and longanize  
burned_monthly2 <- burned_monthly
coeff <- max(burned_monthly$area) / max(burned_monthly$fires)
burned_monthly2$fires <- burned_monthly$fires * coeff
burned_monthly_long <- pivot_longer(burned_monthly2, cols = 3:4, 
                                    names_to = "var", values_to = "y")
intra_fire <-
ggplot(data = burned_monthly_long,
       mapping = aes(x = month_num_cat, y = y / 24, colour = var, shape = var, 
                     group = var)) +
  geom_line() +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  theme(legend.position = "right") + 
  scale_y_continuous(
    name = "Mean burned area (ha)",
    sec.axis = sec_axis(~ . / coeff, name = "Mean number of fires")
  ) 
intra_fire

# climate variables

clim_intra <- read.csv("data_climate_monthly.csv")
names(clim_intra)[1] <- "month_num"
clim_intra <- left_join(
  clim_intra,
  burned_monthly[, c("month", "month_num", "month_num_cat")],
  by = "month_num"
)

# scale
coef_clim <- max(clim_intra$value[clim_intra$variable == "prec"]) /
               max(clim_intra$value[clim_intra$variable == "tavg"])

clim_intra$y_scaled <- clim_intra$value
clim_intra$y_scaled[clim_intra$variable == "tavg"] <- 
  clim_intra$value[clim_intra$variable == "tavg"] * coef_clim

clim_intra$variable <- plyr::revalue(
  clim_intra$variable,
  c("prec" = "pp", "tavg" = "temp")
)

intra_clim <-
ggplot(data = clim_intra,
       mapping = aes(x = month_num_cat, y = y_scaled, 
                     colour = variable, shape = variable, 
                     group = variable)) +
  geom_line() +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_viridis(discrete = TRUE, option = "D", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.5, direction = -1) +
  theme(legend.position = "right",
        axis.title.y.left = element_text(vjust = 5),
        plot.margin = unit(c(2, 2, 2, 4), "mm")) + 
  xlab("Month") +
  scale_y_continuous(
    name = "Precipitation (mm)",
    sec.axis = sec_axis(~ . / coef_clim, name = "Temperature (°C)")
  )

intra_clim


# merge

ggarrange(intra_fire + theme(axis.title.x = element_blank()), 
          intra_clim, ncol = 1, heights = c(0.93, 1))
ggsave("figure_intraannual fire activity.png", width = 10, height = 13, 
       units = "cm")






# Fire size distribution --------------------------------------------------







