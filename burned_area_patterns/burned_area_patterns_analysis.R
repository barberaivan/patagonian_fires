# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

# library(ggpubr)
# library(gridExtra)
# library(gtable)   # merge density plots
library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr

library(scales)   # log scale 
library(circular) # density.circular, for aspect
library(brms)     # Dirichlet model

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

# compute overlap
ol_num <- 1 - sum(na.omit(abs(burned_veg$prob_av - burned_veg$prob_b))) / 2
(ol_text <- paste(round(ol_num * 100, 2), "%"))

data_ol <- data.frame(prob = 0.3, 
                      vegetation_class = 6.5,
                      lab = ol_text)

# plot
veg_dist <- 
ggplot(veg_dist_data, 
       aes(x = vegetation_class, y = prob, 
           colour = class, fill = class)) + 
  geom_bar(stat = "identity", 
           # position = "identity",
           position = position_dodge2(width = 2, padding = 0.05),
           alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  #facet_wrap(vars(variable)) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", #c(0.85, 0.8),
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.67),
        axis.title.x = element_text(vjust = 7),
        axis.title.y = element_blank(),
        strip.text.x = element_text(margin = margin(3,0,3,0, "mm"))) +
  ylab("Probability mass") + 
  xlab("Vegetation type") + 
  geom_text(data = data_ol, inherit.aes = FALSE,
            mapping = aes(x = vegetation_class, y = prob, label = lab), 
            colour = "black")
veg_dist
# ggsave("figure_burned_area_by_veg.png")

# Densities and overlap

n_var <- 6
l <- 512 * 2 * n_var

# list to save densities
dlist <- vector(mode = "list", length = n_var)
names(dlist) <- c("elevation", "slope", "aspect", "solar",
                  "distance_humans_km", "distance_roads_km")
# list to save overlaps
ovlist <- dlist

# rescale distances to km
spatial_vars$distance_humans_km <- spatial_vars$distance_humans / 1000  
spatial_vars$distance_roads_km <- spatial_vars$distance_roads / 1000  

# Compute densities and overlap
for(i in 1:length(dlist)) {
  
  # i = 3
  v <- names(dlist)[i]
  
  # densities
  d_burned <- density(spatial_vars[spatial_vars$class == "Burned", v], from = 0)
  d_av <- density(spatial_vars[spatial_vars$class == "Available", v], from = 0)
  
  if(v == "aspect") {
    filter_b <- spatial_vars$class == "Burned"
    filter_av <- spatial_vars$class == "Available"
    
    aspect_burned <- circular(spatial_vars$aspect[filter_b], 
                              type = "angles", units = "degrees",
                              template = "geographic")
    aspect_av <- circular(spatial_vars$aspect[filter_av], 
                                type = "angles", units = "degrees",
                                template = "geographic")
    
    bw_circ <- 20 # smaller = smoother
    
    d_burned <- density.circular(aspect_burned, bw = bw_circ)
    d_av <- density.circular(aspect_av, bw = bw_circ)
    
  }
  
  # densities df
  dlist[[i]] <- data.frame(variable = v,
                           value = c(d_burned$x, d_av$x),
                           density = c(d_burned$y, d_av$y),
                           class = factor(rep(c("Burned", "Available"), 
                                              each = length(d_av$x)),
                                          levels = c("Burned", "Available")))
  
  # Compute delta and overlap
  
  # approximate densities on a common x sequence
  xseq <- seq(min(dlist[[i]]$value), max(dlist[[i]]$value), length.out = 200)
  if(v == "aspect") {
    xseq <- seq(max(d_av$x), min(d_av$x), length.out = 200)
  }
  diff_size <- abs(unique(diff(xseq))[1])
  
  d_burned_pred <- approx(d_burned$x, d_burned$y, xseq, method = "linear",
                          yleft = 0, yright = 0)
  d_av_pred <- approx(d_av$x, d_av$y, xseq, method = "linear",
                      yleft = 0, yright = 0)

  # circular densities are not normalized, so we get the normalizing factor
  if(v == "aspect") {
    auc_b <- sum(d_burned_pred$y * diff_size)
    auc_av <- sum(d_av_pred$y * diff_size)
    norm_factor <- 1/ mean(auc_b, auc_av)
  }
  
  den_diff <- (d_burned_pred$y - d_av_pred$y) 
  delta <-  ((abs(den_diff) * diff_size) %>% sum) / 2
  if(v == "aspect") {delta <- delta * norm_factor}
  overlap <- 1 - delta
  
  # delta and overlap df
  ovlist[[i]] <- data.frame(variable = v, 
                            value = quantile(xseq, 0.8),
                            density = max(c(d_burned$y, d_av$y)) * 0.8,
                            overlap = paste(round(overlap * 100, 2), "%"),
                            delta = paste(round(delta * 100, 2), "%"))
  
}

# good variable names:
var_names <- 
c("Elevation (m a.s.l.)", 
  "Slope (°)",
  "Aspect (°)",
  "Mean daily solar\nirradiation (kWh / m2)",
  "Distance from human\nsettlements (km)", 
  "Distance from\nroads (km)")

plot_list <- list()
for(i in 1:n_var) {
  
  # i <- 1
  plot_list[[i]] <- 
  ggplot(dlist[[i]], 
         aes(x = value, y = density, ymin = 0, ymax = density,
             colour = class, fill = class)) + 
    geom_line(show.legend = F) + 
    geom_ribbon(alpha = 0.2, colour = NA) +
    # geom_ribbon(alpha = 0.2, size = 0.4) +  
    # geom_hline(yintercept = 0, colour = "white", size = 0.45, alpha = 1) +
    # facet_wrap(vars(variable), ncol = 3, scales = "free") +
    scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
    scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
    geom_text(data = ovlist[[i]], mapping = aes(x = value, y = density, 
                                                label = overlap),
              size = 4, inherit.aes = FALSE) +
    theme(panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.position = "none", 
          axis.title.y = element_blank()) +
    xlab(var_names[i]) + 
    ylab("Probability density")
  
  # print(plot_list[[i]])
}

# Relabel aspect plot axis
plot_list[[3]] <- plot_list[[3]] + 
                  scale_x_continuous(breaks = c(-270, -180, -90, 0, 90),
                                     labels = c("E", "S", "W", "N", "E")) + 
                  xlab("Aspect")

# Merge plots using egg approach
# https://cran.r-project.org/web/packages/egg/vignettes/Overview.html

joined_dens <- 
ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],
          plot_list[[4]], plot_list[[5]], plot_list[[6]],
          nrow = 2)

join_all <- 
grid.arrange(
  # veg dist
  arrangeGrob(
    textGrob("Probability mass", 
             gp = gpar(fontsize = 12), rot = 90,
             hjust = 0.2), 
    ggplotGrob(veg_dist + 
               theme(plot.margin = margin(l = 7.5, r = 3, t = 2, unit = "mm"))), 
    nrow = 1, widths = c(1, 20)
  ),
  
  # 6 paneles
  arrangeGrob(
    textGrob("Probability density", gp = gpar(fontsize = 12), 
             rot = 90,
             hjust = 0.4), 
    joined_dens, 
    nrow = 1, widths = c(1, 20)
  ),
  nrow = 2, heights = c(1, 1.5)
)
ggsave("figures/trial6.png", join_all,
       height = 17, width = 17, units = "cm")









# # merge variables in df
# d_data <- do.call("rbind", dlist)
# d_data$variable <- plyr::revalue(d_data$variable,
#                          c("elevation" = "Elevation (m a.s.l.)", 
#                            "slope" = "Slope (°)",
#                            "aspect" = "Aspect (°)",
#                            "solar" = "Mean daily solar\nirradiation (kWh / m2)",
#                            "distance_humans_km" = "Distance from human\nsettlements (km)", 
#                            "distance_roads_km" = "Distance from\nroads (km)"))
# d_data$variable <- factor(d_data$variable, levels = 
#                           c("Elevation (m a.s.l.)", 
#                             "Aspect (°)",
#                             "Slope (°)",
#                             "Mean daily solar\nirradiation (kWh / m2)",
#                             "Distance from human\nsettlements (km)", 
#                             "Distance from\nroads (km)"))




# Plot


var6plot <- ggplot(d_data, 
                   aes(x = value, y = density, ymin = 0, ymax = density,
                       colour = class, fill = class)) + 
  geom_line(show.legend = F) + 
  geom_ribbon(alpha = 0.2, colour = NA) +
  # geom_ribbon(alpha = 0.2, size = 0.4) +  
  # geom_hline(yintercept = 0, colour = "white", size = 0.45, alpha = 1) +
  facet_wrap(vars(variable), ncol = 3, scales = "free") +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  xlab("Variable") + 
  ylab("Probability density")
var6plot




# Merge this and vegetation type
veg_dist2 <- veg_dist + theme(axis.title.y = element_text(vjust = -10),
                              plot.margin = unit(c(2,2,2,5.8), "mm"))
# veg_dist2
ggpubr::ggarrange(veg_dist2, 
                  var6plot, 
                  nrow = 2, heights = c(1.3, 2))
# ggsave("figures/Figure - spatial variables distributions.jpeg",
#        width = 17, height = 21, units = "cm", dpi = 500)






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

# ggsave("figures/Figure - ndvi distribution by vegetation type.jpeg", 
#        width = 17, height = 15, dpi = 500, units = "cm")

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

# sum_length <- function(x) c(sum = sum(x), lenght = length(x))
# burned_annual <- do.call("data.frame",
#                          aggregate(area_ha ~ year, data = fires_clipped@data,
#                                    FUN = sum_length))
# colnames(burned_annual) <- c("year", "area_ha", "fires")
# write.csv(burned_annual, "data_burned_area_by_year.csv")
burned_annual <- read.csv("data_burned_area_by_year.csv")


# Burned area as a function of climate ------------------------------------


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



# merge with fire data
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

# con gg arrange (no anda bien el eje), con gridExtra no anda ggsave.

x_axis <- textGrob("Climatic variable", gp = gpar(cex = 1), hjust = 0.3,
                   vjust = -0.25)
clim_fig <-
grid.arrange(arrangeGrob(clim_area + theme(axis.title.x = element_blank()), 
                         clim_fires + theme(axis.title.x = element_blank()),
                         ncol = 2, widths = c(1, 1)),
             x_axis, nrow = 2, heights = c(10, 1))

# jpeg(file = "figures/Figure - fire activity and climate interannual.jpeg",
#      width = 12, height = 20, units = "cm", quality = 100,
#      res = 900)
# plot(clim_fig)
# dev.off()



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



# Burned area by year and vegetation class -----------------------------


# veg_map <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/vegetation_valdivian_reclassified/vegetation_valdivian_reclassified.shp")
# vegmap_posgar <- spTransform(veg_map, proj_posgar)
# vegmap_posgar <- vegmap_posgar[!is.na(vegmap_posgar$clss_nm), ]
# 
# burned_annual_veg <- expand.grid(
#   year = 1999:2022,
#   vegetation_code = 2:8, # without non-burnable
#   area_ha = 0
# )
# 
# # for(i in 1:nrow(burned_annual_veg)) {
# #   print(i)
# #   # i = 3
# #   fires_year <- fires_clipped[fires_clipped$year == burned_annual_veg$year[i], ]
# #   veg_local <- vegmap_posgar[vegmap_posgar$clss_nm == burned_annual_veg$vegetation_code[i], ]
# #   
# #   burned_veg <- gIntersection(polygons(fires_year), polygons(veg_local))
# #   if(is.null(burned_veg)) {burned_annual_veg$area_ha[i] <- 0} else {
# #     burned_annual_veg$area_ha[i] <- gArea(burned_veg) * 0.0001 # m2 to ha
# #   }
# # } # takes long time
# 
# # bring vegetation data
# burned_annual_veg <- left_join(
#   burned_annual_veg[, c("year", "vegetation_code", "area_ha")], 
#   burned_veg[, c("vegetation_code", "vegetation_class", "area_ha_available",
#                  "area_ha_burned")],
#   by = "vegetation_code"
# )
# 
# # add burned area by year
# burned_annual_sum <- aggregate(area_ha ~ year, burned_annual_veg, sum)
# names(burned_annual_sum)[2] <- "area_ha_annual"
# 
# burned_annual_veg <- left_join(
#   burned_annual_veg, 
#   burned_annual_sum[, c("year", "area_ha_annual")],
#   by = "year"
# )
# 
# # burned area / available area by year and vegetation
# burned_annual_veg$proportion_burned_annual <- burned_annual_veg$area_ha / 
#   burned_annual_veg$area_ha_available
# 
# # burned area / total burned area in a given year, by year and vegetation
# burned_annual_veg$burn_distribution <- burned_annual_veg$area_ha / 
#   burned_annual_veg$area_ha_annual
# # sum to 1?
# # aggregate(burn_distribution ~ year, burned_annual_veg, sum)


# write.csv(burned_annual_veg, "data_burned_area_year_veg.csv")
burned_annual_veg <- read.csv("data_burned_area_year_veg.csv")[, -1]


# Burn distribution as a function of climate ------------------------------


burned_annual_veg_clim <- left_join(
  burned_annual_veg, climate, by = "year"  
)

# relabel vegetation
burned_annual_veg_clim$vegetation_class <- factor(
  burned_annual_veg_clim$vegetation_class, 
  levels = veg_levels, labels = veg_levels2
)

# convert zeroes
burned_annual_veg_clim$burn_distribution2 <-
  burned_annual_veg_clim$burn_distribution
burned_annual_veg_clim$burn_distribution2[burned_annual_veg_clim$burn_distribution == 0] <- 1e-5
burned_annual_veg_clim$burn_distribution2[burned_annual_veg_clim$burn_distribution == 1] <- 1 - 1e-5


# Burn distribution (veg) ~ c(fwi, pp, temp) ------------------------------


brms_data <- pivot_wider(burned_annual_veg_clim[, c("vegetation_class", 
                                                    "burn_distribution2",
                                                    "year", "fwi", "pp", "temp")], 
                         names_from = "vegetation_class", 
                         values_from = "burn_distribution2")
y_matrix <- brms_data[, veg_levels2[-8]] %>% as.matrix
y_matrix <- apply(y_matrix, 1, normalize) %>% t
brms_data$y_matrix <- y_matrix

# add area
brms_data <- left_join(brms_data, burned_annual[, c("year", "area_ha")],
                       by = "year")


# modelling dispersion to account for higher dispersion in years with low
# burned area (the fwi effect decreases in shrublands compared to the previous
# model because this one downweights extreme values in years with low burned 
# area; the opposite occurs for subalpine forests).

var_names <- c("fwi", "pp", "temp")
var_labels <- c("Fire weather index (FWI)", "Precipitation (mm)", "Temperature (°C)")

# for(i in 1:3) {
i <- 2 # variar en 1 a 3

var_name <- var_names[i]
var_label <- var_labels[i]
var_breaks <- list(fwi = seq(8, 20, length.out = 4),
                   pp = seq(20, 260, length.out = 4),
                   temp = 18:22)
var <- brms_data[, var_name, drop = TRUE]
brms_data$var <- var

pmodel_brms_d <- brm(bf(y_matrix ~ var, phi ~ log10(area_ha)),
                     data = brms_data, 
                     family = dirichlet(link = "logit", link_phi = "log"),
                     chains = 4, cores = 4)

print(paste("min ESS:",
            min(c(sss$fixed$Bulk_ESS, sss$fixed$Tail_ESS))))
print(paste("max Rhat:", max(sss$fixed$Rhat)))

nseq <- 150
pdata_wide <- data.frame(row = 1:nseq,
                         area_ha = rep(log10(median(burned_annual$area_ha)), # mean(log10(burned_annual$area_ha)), 
                                       nseq),
                         var = seq(min(var), max(var), length.out = nseq))

preds2 <- fitted(pmodel_brms_d, pdata_wide)

dimnames(preds2) <- list(row = 1:nseq,
                         summary = c("p_mean", "Est.Error", "p_lower", "p_upper"),
                         vegetation_class = veg_levels2[-8])

preds3 <- pivot_wider(as.data.frame.table(preds2), 
                      names_from = "summary",
                      values_from = "Freq")

preds3$row <- as.numeric(preds3$row)
preds4 <- left_join(preds3, pdata_wide, by = "row")

# plot dirichlet model

dtemp <- burned_annual_clim[, c("year", "area_ha")]
names(dtemp) <- c("year", "area_ha_total")

d <- left_join(burned_annual_veg_clim,
               dtemp[, c("year", "area_ha_total")],
               by = "year")
d$var <- burned_annual_veg_clim[, var_name]
  
vdist2 <- 
  ggplot() +
  geom_ribbon(data = preds4, 
              mapping = aes(x = var, y = p_mean, ymin = p_lower, ymax = p_upper,
                            fill = vegetation_class, color = vegetation_class), 
              alpha = 0.3, size = 0.1) +
  geom_line(data = preds4, mapping = aes(x = var, y = p_mean),
            alpha = 1) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.8) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.8) +
  geom_point(data = d, 
             mapping = aes(x = var, y = burn_distribution2, 
                           color = vegetation_class, fill = vegetation_class,
                           alpha = area_ha_total),
             size = 3) +
  scale_alpha_continuous(trans = "log10", name = "Annual burned\narea (ha)",
                         breaks = c(100, 1000, 10000, 30000)) +
  scale_x_continuous(breaks = var_breaks[[i]]) +
  facet_wrap(vars(vegetation_class), nrow = 2) +
  theme(legend.position = c(0.88, 0.25), 
        axis.text = element_text(size = 9),
        strip.text = element_text(size = 9),
        legend.title = element_text(size = 9, vjust = 2),
        legend.text = element_text(size = 8),
        legend.spacing.y = unit(0.01, "mm")) +
  guides(color = "none", fill = "none", alpha = guide_legend(byrow = TRUE)) + 
  ylab("Proportion of annual burned area") + 
  xlab(var_label)
vdist2

# Stacked proportions

# ord_veg <- veg_levels2[c(1, 2, 5, 4, 3, 7, 6)]

means_wide <- preds2[, "p_mean", ]
means_wide <- means_wide[, ncol(means_wide):1]
means_cum <- apply(means_wide, 1, cumsum) %>% t
means_cum <- cbind(zero = rep(0, nrow(means_cum)), means_cum)

uppers <- as.data.frame(means_cum[, -1])
uppers$var <- pdata_wide$var
uppers$row <- pdata_wide$row

lowers <- as.data.frame(means_cum[, -8])
colnames(lowers) <- colnames(means_cum[, -1])
lowers$row <- pdata_wide$row

lowers_long <- pivot_longer(lowers, 1:7, values_to = "p_lower", 
                            names_to = "vegetation_class")
uppers_long <- pivot_longer(uppers, 1:7, values_to = "p_upper", 
                            names_to = "vegetation_class")
# plot(lowers_long$row ~ uppers_long$row)

# p_ribbons <- left_join(uppers_long, 
#                        lowers_long[, c("row", "p_lower")], 
#                        by = "row")
p_ribbons <- uppers_long
p_ribbons$p_lower <- lowers_long$p_lower
p_ribbons$vegetation_class <- factor(p_ribbons$vegetation_class,
                                     levels = veg_levels2)

veg_stack <- 
ggplot(p_ribbons, aes(x = var, ymin = p_lower, ymax = p_upper, 
                      colour = vegetation_class, fill = vegetation_class,
                      group = vegetation_class)) +
  geom_ribbon(color = NA, alpha = 0.85) + 
  geom_line(mapping = aes(x = var, y = p_lower), size = 0.1,
            color = "black") +#, linetype = "dotted") +
  scale_x_continuous(breaks = var_breaks[[i]]) + # , limits = range(var_breaks[[i]]
  theme(legend.position = "right") +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.8) +
  ylab("Cumulative probability") + 
  xlab(var_label) +
  theme(legend.spacing.y = unit(1, "mm"),
        legend.text = element_text(size = 8)) + 
  guides(fill = guide_legend(byrow = TRUE))
  # facet_wrap(vars(vegetation_class))
veg_stack

merged <- 
ggarrange(veg_stack, vdist2, nrow = 2, heights = c(0.8, 1),
          labels = c("A", "B"))
merged

filename <- paste("figures/Figure - Proportion burned area by veg and ", 
                  var_name, ".jpeg", sep = "")
ggsave(filename, merged,
       units = "cm", dpi = 300, width = 17, height = 16)

# }
  


# Burned proportion (of available) ~ fwi ----------------------------------

# convert zeroes
burned_annual_veg_clim$proportion_burned_annual2 <-
  burned_annual_veg_clim$proportion_burned_annual
burned_annual_veg_clim$proportion_burned_annual2[burned_annual_veg_clim$proportion_burned_annual == 0] <- 1e-5
burned_annual_veg_clim$proportion_burned_annual2[burned_annual_veg_clim$proportion_burned_annual == 1] <- 1 - 1e-5

# Se podría ajustar un modelo Dirichlet, pero primero vayamos con un logitnormal.
# (stat_smooth se caga.)

burned_annual_veg_clim$logit_p <- qlogis(burned_annual_veg_clim$proportion_burned_annual2)
pmodel <- lm(logit_p ~ fwi * vegetation_class, data = burned_annual_veg_clim)
pdata <- expand.grid(
  fwi = seq(min(burned_annual_veg_clim$fwi), max(burned_annual_veg_clim$fwi), 
            length.out = 150),
  vegetation_class = veg_levels2[-8]
)
preds <- predict(pmodel, pdata, se.fit = TRUE)

# Compute predicted probs
pdata$p_lower <- NA
pdata$p_upper <- NA
pdata$p_mle <- NA

for(i in 1:nrow(pdata)) {
  pdata$p_mle[i] <- momentsLogitnorm(preds$fit[i], preds$residual.scale)["mean"]
  pdata$p_lower[i] <- momentsLogitnorm(preds$fit[i] - qnorm(0.975) * preds$se.fit[i], 
                                       preds$residual.scale)["mean"]
  pdata$p_upper[i] <- momentsLogitnorm(preds$fit[i] + qnorm(0.975) * preds$se.fit[i], 
                                       preds$residual.scale)["mean"]
  
}

dtemp <- burned_annual_clim[, c("year", "area_ha")]
names(dtemp) <- c("year", "area_ha_total")
d2 <- left_join(burned_annual_veg_clim,
               dtemp[, c("year", "area_ha_total")],
               by = "year")
ggplot(d2, aes(x = fwi, y = proportion_burned_annual2, color = area_ha_total)) +
  # geom_ribbon(data = pdata, mapping = aes(x = fwi, y = p_mle, ymin = p_lower,
  #                                          ymax = p_upper),
  #             color = NA, alpha = 0.3, inherit.aes = F) +
  geom_line(data = pdata, mapping = aes(x = fwi, y = p_mle),
            alpha = 1, inherit.aes = F) +
  scale_color_viridis(option = "B", begin = 0.2, end = 0.8, direction = -1,
                      name = "Annual burned area (ha)",
                      trans = "log10") +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(vars(vegetation_class)) +
  # guides(color = "Annual burned area (ha)") +
  theme(legend.position = c(0.7, 0.15), 
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 9)) +
  guides(colour = guide_colorbar(title.position = "right")) +
  ylab("Burned / available area") + 
  xlab("Fire Weather Index")
  # xlab("Mean maximum temperature (C°)")

ggsave("figures/Figure - Proportion burned over available fwi.jpeg", 
       units = "cm", width = 17, height = 15, dpi = 500)


# Burned proportion (of available) ~ temp ----------------------------------

pmodel <- lm(logit_p ~ temp * vegetation_class, data = burned_annual_veg_clim)
pdata <- expand.grid(
  temp = seq(min(burned_annual_veg_clim$temp), max(burned_annual_veg_clim$temp), 
            length.out = 150),
  vegetation_class = veg_levels2[-8]
)
preds <- predict(pmodel, pdata, se.fit = TRUE)

# Compute predicted probs
pdata$p_lower <- NA
pdata$p_upper <- NA
pdata$p_mle <- NA

for(i in 1:nrow(pdata)) {
  pdata$p_mle[i] <- momentsLogitnorm(preds$fit[i], preds$residual.scale)["mean"]
  pdata$p_lower[i] <- momentsLogitnorm(preds$fit[i] - qnorm(0.975) * preds$se.fit[i], 
                                       preds$residual.scale)["mean"]
  pdata$p_upper[i] <- momentsLogitnorm(preds$fit[i] + qnorm(0.975) * preds$se.fit[i], 
                                       preds$residual.scale)["mean"]
  
}

ggplot(d2, aes(x = temp, y = proportion_burned_annual2, color = area_ha_total)) +
  # geom_ribbon(data = pdata, mapping = aes(x = temp, y = p_mle, ymin = p_lower,
  #                                         ymax = p_upper),
  #             color = NA, alpha = 0.3, inherit.aes = F) +
  geom_line(data = pdata, mapping = aes(x = temp, y = p_mle),
            alpha = 1, inherit.aes = F) +
  scale_color_viridis(option = "B", begin = 0.2, end = 0.8, direction = -1,
                      name = "Annual burned area (ha)",
                      trans = "log10") +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(vars(vegetation_class)) +
  # guides(color = "Annual burned area (ha)") +
  theme(legend.position = c(0.7, 0.15), 
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 9)) +
  guides(colour = guide_colorbar(title.position = "right")) +
  ylab("Burned / available area") + 
  xlab("Mean maximum temperature (C°)")

ggsave("figures/Figure - Proportion burned over available temp.jpeg", 
       units = "cm", width = 17, height = 15, dpi = 500)














ggplot(burned_annual_veg_clim, aes(x = fwi, y = burn_distribution2)) +
  geom_ribbon(data = pdata, mapping = aes(x = fwi, y = p_mle, ymin = p_lower,
                                          ymax = p_upper),
              color = NA, alpha = 0.3,
              fill = viridis(1, begin = 0.2, option = "B")) +
  geom_line(data = pdata, mapping = aes(x = fwi, y = p_mle),
            alpha = 1,
            color = viridis(1, begin = 0.2, option = "B")) +
  geom_point() +
  facet_wrap(vars(vegetation_class))#, scales = "free"

# how did the sum-to-one go?
hist(aggregate(p_mle ~ fwi, pdata, sum)[, 2])
# bad. Fit dirichlet model


ggplot(burned_annual_veg_clim, aes(x = fwi, y = proportion_burned_annual)) +
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              method.args = list(family = Gamma(link = "log")),
              color = viridis(1, begin = 0.2, option = "B"),
              fill = viridis(1, begin = 0.2, option = "B"),
              alpha = 0.2, size = 0.8) + 
  facet_wrap(vars(vegetation_class))#, scales = "free"










# Density difference trials -----------------------------------------------

# Compute density separately for the next plots, to make a facet_wrap later.
n_var <- 6
l <- 512 * 2 * n_var

# list to save densities
dlist <- vector(mode = "list", length = n_var)
names(dlist) <- c("elevation", "slope", "aspect", "solar",
                  "distance_humans_km", "distance_roads_km")
# list to save density differences
diff_list <- dlist

# rescale distances to km
spatial_vars$distance_humans_km <- spatial_vars$distance_humans / 1000  
spatial_vars$distance_roads_km <- spatial_vars$distance_roads / 1000  

# rescale to have max = 1
spat_subset <- spatial_vars[, c("class", names(dlist))]
for(i in 2:7) spat_subset[, i] <- spat_subset[, i] / max(spat_subset[, i])

# Now densities live in similar scales.

# Compute densities
for(i in 1:length(dlist)) {
  
  # i = 1
  v <- names(dlist)[i]
  d_burned <- density(spat_subset[spat_subset$class == "Burned", v], from = 0)
  d_av <- density(spat_subset[spat_subset$class == "Available", v], from = 0)
  
  dlist[[i]] <- data.frame(variable = v,
                           value = c(d_burned$x, d_av$x),
                           value_original = c(d_burned$x, d_av$x) * max(spat_subset[, v]),
                           density_scaled = c(d_burned$y, d_av$y),
                           class = factor(rep(c("Burned", "Available"), 
                                              each = length(d_av$x)),
                                          levels = c("Burned", "Available"))
  )
  
  # approximate densities on a common x sequence
  xseq <- seq(min(dlist[[i]]$value), max(dlist[[i]]$value), length.out = 200)
  d_burned_pred <- approx(d_burned$x, d_burned$y, xseq, method = "linear",
                         yleft = 0, yright = 0)
  d_av_pred <- approx(d_av$x, d_av$y, xseq, method = "linear",
                      yleft = 0, yright = 0)
  
  # normalize on the sequence
  diff_size <- unique(diff(xseq))[1]
  d_burned_pred$y <- d_burned_pred$y / sum(d_burned_pred$y * diff_size)
  # (d_burned_pred$y * diff_size) %>% sum  # it's normalized
  d_av_pred$y <- d_av_pred$y / sum(d_av_pred$y * diff_size)
   
  # compute density difference
  den_diff <- (d_burned_pred$y - d_av_pred$y) 
  
  # compute delta and overlap by hand
  delta <-  ((abs(den_diff) * diff_size) %>% sum) / 2
  overlap <- 1 - delta
  
  # get positive and negative values for ribbon
  positive <- den_diff
  positive[positive < 0] <- 0
  negative <- den_diff
  negative[negative >= 0] <- 0
  
  diff_list[[i]] <- data.frame(
    variable = v,
    value = xseq,
    value_original = xseq * max(spatial_vars[, v]), 
    density_diff = den_diff,
    sign = factor(sign(den_diff), levels = c(-1, 1)),
    positive = positive,
    negative = negative,
    delta = delta, 
    overlap = overlap
  )
  
  

}


# Compute density for aspect (circular distribution)

# filter_b <- spatial_vars$class == "Burned"
# filter_ub <- spatial_vars$class == "Available"
# 
# aspect_burned <- circular(spatial_vars$aspect[filter_b], 
#                           type = "angles", units = "degrees",
#                           template = "geographic")
# aspect_unburned <- circular(spatial_vars$aspect[filter_ub], 
#                             type = "angles", units = "degrees",
#                             template = "geographic")
# 
# bw_circ <- 20 # smaller = smoother
# 
# d_burned <- density.circular(aspect_burned, bw = bw_circ)
# d_unburned <- density.circular(aspect_unburned, bw = bw_circ)
# 
# x_mine <- seq(0, 360, length.out = length(d_unburned$x))
# 
# aspect_data <- rbind(
#   data.frame(variable = "aspect", value = x_mine, density = d_unburned$y, class = "Available"),
#   data.frame(variable = "aspect", value = x_mine, density = d_burned$y, class = "Burned")
# )
# aspect_data$class <- factor(aspect_data$class, levels = c("Burned", "Available"))
# 
# dlist$aspect <- aspect_data



# merge variables in df
diff_data <- do.call("rbind", diff_list)
diff_data$variable <- plyr::revalue(diff_data$variable,
                                 c("elevation" = "Elevation (m a.s.l.)", 
                                   "slope" = "Slope (°)",
                                   "aspect" = "Aspect (°)",
                                   "solar" = "Mean daily solar\nirradiation (kWh / m2)",
                                   "distance_humans_km" = "Distance from human\nsettlements (km)", 
                                   "distance_roads_km" = "Distance from\nroads (km)"))
diff_data$variable <- factor(diff_data$variable, levels = 
                            c("Elevation (m a.s.l.)", 
                              "Aspect (°)",
                              "Slope (°)",
                              "Mean daily solar\nirradiation (kWh / m2)",
                              "Distance from human\nsettlements (km)", 
                              "Distance from\nroads (km)"))


# var6plot <- ggplot(diff_data, 
#                    aes(x = value, y = density, ymin = 0, ymax = density,
#                        colour = class, fill = class)) + 
#   geom_line(show.legend = F) + 
#   geom_ribbon(alpha = 0.2, colour = NA) +
#   # geom_ribbon(alpha = 0.2, size = 0.4) +  
#   # geom_hline(yintercept = 0, colour = "white", size = 0.45, alpha = 1) +
#   facet_wrap(vars(variable), ncol = 3, scales = "free") +
#   scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
#   scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
#   theme(panel.grid.minor = element_blank(),
#         legend.title = element_blank(),
#         legend.position = "none") +
#   xlab("Variable") + 
#   ylab("Probability density")
# var6plot



# plot diff
data_fill <- data.frame(x = c(0, 0, 1, 1), y = c(0, -2, 0, -2))
index_data <- aggregate(cbind(delta, overlap) ~ variable, diff_data, mean)
index_data$x_original <- aggregate(value_original ~ variable, diff_data, 
                                   quantile, probs = c(0.80))[, 2]
index_data$delta_text <- paste(round(index_data$delta * 100, 1), "%", sep = "")


# with predictors scaled at [0, 1]
ggplot(diff_data, aes(x = value, y = density_diff)) +
  # reference area
  geom_ribbon(data = data_fill, mapping = aes(x = x, ymin = -2, ymax = 0),
              inherit.aes = F, alpha = 0.2) +
  geom_line() + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-2, 2), alpha = 0.5, size = 0.3) + 
  ylim(-2, 2) +
  geom_ribbon(mapping = aes(ymin = 0, ymax = positive), colour = NA,
              fill = viridis(option = "B", end = 0.5, direction = -1, n = 1), 
              alpha = 0.3) +
  geom_ribbon(mapping = aes(ymin = negative, ymax = 0), colour = NA,
              fill = viridis(option = "B", end = 0.1, direction = 1, n = 1), 
              alpha = 0.3) +
  geom_text(data = index_data,
            mapping = aes(x = 0.75, y = 1, label = round(delta, 3)),
            size = 4) +
  facet_wrap(vars(variable), ncol = 3) +
  ylab("Probability density difference") +
  xlab("Variable")
  


vircols <- viridis(2, begin = 0.2, end = 0.6, option = "B")
ggplot(diff_data, aes(x = value_original, y = density_diff)) +
  # reference area
  geom_ribbon(data = diff_data, mapping = aes(x = value_original, 
                                              ymin = -2, ymax = 0),
              inherit.aes = F, alpha = 0.2) +
  geom_line() + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-2, 2), alpha = 0.5, size = 0.3) + 
  ylim(-2, 2) +
  geom_ribbon(mapping = aes(ymin = 0, ymax = positive), colour = NA,
              #fill = viridis(option = "B", end = 0.5, direction = -1, n = 1), 
              fill = vircols[2],
              alpha = 0.3) +
  geom_ribbon(mapping = aes(ymin = negative, ymax = 0), colour = NA,
              # fill = viridis(option = "B", end = 0.1, direction = 1, n = 1), 
              fill = vircols[1],
              alpha = 0.3) +
  geom_text(data = index_data,
            mapping = aes(x = x_original, y = 1.2, label = delta_text),
            size = 3.5,
            inherit.aes = F) +
  facet_wrap(vars(variable), ncol = 3, scales = "free_x") +
  ylab("Fire selectivity index") +
  #ylab("Probability density difference") +
  xlab("Variable")

ggsave("figures/Figure - spatial variables selectivity index.jpeg",
       width = 17, height = 12, units = "cm", dpi = 500)




# Intraannual fire activity -----------------------------------------------

remove <- which(fires_clipped@data$obs == "uncertain_year")
dmonth <- fires_clipped@data[-remove, ]
dmonth$month_num <- format(as.Date(dmonth$date, format = "%Y-%m-%d"), 
                           format = "%m") %>% as.numeric
dmonth$month <- factor(dmonth$month_num, levels = c(7:12, 1:6), 
                       labels = month.name[c(7:12, 1:6)])


# ### Plot with annual values (failure)
# 
# 
# # aggregate by month and year (sum area, count fires)
# burned_my <- do.call("data.frame",
#                      aggregate(area_ha ~ month_num + month + year,
#                                data = dmonth,
#                                FUN = sum_length))
# colnames(burned_my) <- c("month_num", "month", "year", "area", "fires")
# 
# # Average area and fire number by month
# burned_monthly <- do.call("data.frame",
#                           aggregate(cbind(area, fires) ~ month_num + month, 
#                                     data = burned_my,
#                                     FUN = mean))
# colnames(burned_monthly) <- c("month_num", "month", "area", "fires")
# # write.csv(burned_monthly, "data_burned_area_by_month.csv")
# 
# 
# # categorical numerical month
# burned_my$month_num_cat <- factor(burned_my$month_num, 
#                                   levels = c(7:12, 1:6))
# burned_monthly$month_num_cat <- factor(burned_monthly$month_num, 
#                                        levels = c(7:12, 1:6))
# 
# # Longanize monthly and month-yearly
# burned_monthly2 <- burned_monthly
# burned_my2 <- burned_my
# coeff <- max(burned_my$area) / max(burned_my$fires)
# 
# burned_monthly2$fires <- burned_monthly$fires * coeff
# burned_my2$fires <- burned_my$fires * coeff
# 
# 
# burned_monthly_long <- pivot_longer(burned_monthly2, cols = 3:4, 
#                                     names_to = "var", values_to = "y")
# burned_my_long <- pivot_longer(burned_my2, cols = 3:4, 
#                                names_to = "var", values_to = "y")
# 
# 
# ggplot(mapping = aes(x = month_num_cat, y = y, colour = var, shape = var, 
#                      group = var)) +
#   geom_point(data = burned_my_long, size = 3, alpha = 0.7) + 
#   geom_line(data = burned_monthly_long) +
#   geom_point(data = burned_monthly_long, size = 3, alpha = 0.7) + 
#   scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
#   scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
#   scale_y_continuous(
#     # Features of the first axis
#     name = "Mean monthly burned area (ha)",
#     trans = "log10",
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~ . / coeff, name = "Mean number of fires",
#                         trans = NULL)
#   ) +
#   theme(legend.position = "right")


### Plot without annual values 

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
# ggsave("figures/Figure - intraannual fire activity.jpeg",
#        width = 12, height = 13, dpi = 500, units = "cm")




# Fire size distribution --------------------------------------------------


size_data <- data.frame("area_abs" = fires_posgar$area_ha)

size_data <- size_data[order(size_data$area_abs, decreasing = TRUE), , drop = F]
size_data$area_prop <- cumsum(size_data$area_abs) / sum(size_data$area_abs)
size_data$number_abs <- 1:nrow(size_data)
size_data$number_prop <- cumsum(size_data$number_abs) / sum(size_data$number_abs)

size_props <- 
  ggplot(size_data, aes(y = area_prop, x = number_prop)) + 
  geom_point(
    color = "black",#viridis(1, begin = 0.2, option = "B"),
    size = 2, alpha = 0.6
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
  ylab("Cumulative proportion of burned area") + 
  xlab("Cumulative proportion of fires")
size_props

# bin to compute log-log relationship
# k <- 15
# size <- ceiling(nrow(size_data) / k)
# size_data$group <- rep(1:k, each = size)[1:nrow(size_data)]

size_data$area_abs_log10 <- log10(size_data$area_abs)
size_data$group <- cut(size_data$area_abs_log10, breaks = 15)

mean_length <- function(x) c(mean = mean(x), length = length(x))

size_bins <- do.call("data.frame", 
                     aggregate(area_abs ~ group, size_data, mean_length))

# acá falta agregar un bin que tenga freq = 13-5.
# y además, no debería usar las medias si no las intermedias de cada clase. 
size_bins$group

size_freq <- 
  ggplot(size_bins, aes(x = area_abs.mean, y = area_abs.length / nrow(size_data))) + 
  # geom_smooth(method = "glm", formula = y ~ log10(x), 
  #             method.args	= list(family = Gamma(link = "log"))) +
  geom_smooth(method = "lm", se = TRUE,
              color = viridis(1, begin = 0.2, option = "B"),
              fill = viridis(1, begin = 0.2, option = "B"),
              alpha = 0.2, size = 0.8) +
  geom_point(size = 2.5, alpha = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  ylab("Relative frequency") + 
  xlab("Area class (ha)")
size_freq

ggarrange(size_props, size_freq, nrow = 1)
ggsave("figures/Figure - fire size distribution.jpeg", 
       width = 17, height = 8.5, units = "cm", dpi = 300)
