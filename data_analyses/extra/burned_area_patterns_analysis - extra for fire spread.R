# Analyses of the relationship between static climatic and topographic variables
# and NDVI_max to decide which variables to include in the fire spread model.
# The analysis is done conditioning on vegetation type, so results are more
# similar to what the spread model is going to estimate.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr

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

d <- read.csv("data_ndvi-pp-temp-elev-aspect-solar_by_vegetation.csv")
d$temp <- d$temp * 0.1 # scale to °C
varnames <- c("ndvi", "pp", "temperature", "elevation", "aspect", "radiation",
              "vegetation_code", "class")
names(d) <- varnames

# bring name of each vegetation code from another dataset
burned_veg <- read.csv("data_burned_and_available_area_by_vegetation.csv")

d <- left_join(d, burned_veg[, c("vegetation_code", "vegetation_class")],
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

d$vegetation_class <- factor(d$vegetation_class, levels = veg_levels,
                             labels = veg_levels2)

# View(d)

# Densities for all variables by veg type for burned and unburned -----------

npred <- 6

plot_list <- vector(mode = "list", length = npred)

for(i in 1:npred) {
  # i <- 1
  plot_list[[i]] <- ggplot(d[d$vegetation_code > 1, ],
                           aes_string(x = varnames[i],
                                      fill = "class")) +
    geom_density(alpha = 0.4, size = 0.2) +
    facet_wrap(vars(vegetation_class)) +
    scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
    scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
    theme(legend.position = c(0.66, 0.15))

  # plot_list[[i]]

  ggsave(filename = paste("figures/fire_spread_patterns_densities_", varnames[i], ".png", sep = ""),
         plot = plot_list[[i]],
         width = 17, height = 15, dpi = 300, units = "cm")
}



# Dispersion: NDVI as a function of ...  ----------------------------------

dsub <- d[d$vegetation_code > 1, ]
disp_list <- vector(mode = "list", length = npred)

# precipitation
p <- ggplot(dsub,#[sample(1:nrow(dsub), size = 100000), ],
            aes(x = pp, y = ndvi)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  facet_wrap(vars(vegetation_class)) +
  theme(legend.position = c(0.66, 0.15))
p
ggsave("figures/fire_spread_patterns_ndvi-response_pp.png", plot = p,
       width = 17, height = 15, dpi = 300, units = "cm")


# aspect
p <- ggplot(dsub,
            aes(x = aspect, y = ndvi)) +
  geom_point(alpha = 0.03) +
  geom_smooth() +
  facet_wrap(vars(vegetation_class)) +
  theme(legend.position = c(0.66, 0.15))
p
ggsave("figures/fire_spread_patterns_ndvi-response_aspect.png", plot = p,
       width = 17, height = 15, dpi = 300, units = "cm")

# elevation
p <- ggplot(dsub,
            aes(x = elevation, y = ndvi)) +
  geom_point(alpha = 0.03) +
  geom_smooth() +
  facet_wrap(vars(vegetation_class)) +
  theme(legend.position = c(0.66, 0.15))
p
ggsave("figures/fire_spread_patterns_ndvi-response_elevation.png", plot = p,
       width = 17, height = 15, dpi = 300, units = "cm")

# radiation
p <- ggplot(dsub,
            aes(x = radiation, y = ndvi)) +
  geom_point(alpha = 0.03) +
  geom_smooth() +
  facet_wrap(vars(vegetation_class)) +
  theme(legend.position = c(0.66, 0.15))
ggsave("figures/fire_spread_patterns_ndvi-response_radiation.png", plot = p,
       width = 17, height = 15, dpi = 300, units = "cm")

# temperature
p <- ggplot(dsub,
            aes(x = temperature, y = ndvi)) +
  geom_point(alpha = 0.03) +
  geom_smooth() +
  facet_wrap(vars(vegetation_class)) +
  theme(legend.position = c(0.66, 0.15))
ggsave("figures/fire_spread_patterns_ndvi-response_temperature.png", plot = p,
       width = 17, height = 15, dpi = 300, units = "cm")




# Dispersion: elevation and temperature -----------------------------------

dsub <- d[d$vegetation_code > 1, ]

a <-
ggplot(dsub, aes(x = elevation, y = temperature)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  facet_wrap(vars(vegetation_class)) +
  theme(legend.position = c(0.66, 0.15))

ggsave("figures/fire_spread_patterns_elevation-temperature.png",
       plot = a,
       width = 17, height = 15, dpi = 300, units = "cm")


# Burned area as a function of climate (dynamic) by veg type --------------

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

# Veg burned by year
burned_annual_veg <- read.csv("data_burned_area_year_veg.csv")[, -1]
burned_annual_veg$vegetation_class <- factor(burned_annual_veg$vegetation_class,
                                             levels = veg_levels,
                                             labels = veg_levels2)

# View(burned_annual_veg)


# merge
bclim <- left_join(burned_annual_veg,
                   climate[, c("year", "pp", "temp", "fwi", "wind", "vpd")],
                   by = "year")
clim_vars <- c("pp", "temp", "fwi", "vpd")

# duplicate dataset for geom_smooth with and without 2015
bclim_with <- bclim
bclim_with$dataset <- "with 2015"
bclim_without <- bclim_with[bclim_with$year != 2015, ]
bclim_without$dataset <- "without 2015"

bclim2 <- rbind(bclim_with, bclim_without)
bclim2$area_ha[bclim2$area_ha == 0] <- 0.001

# Plots in original scale

for(i in 1:length(clim_vars)) {
  # i <- 1

  # theme_set(theme_bw()) # problema del fondo no es por mi theme.
  p <-
  ggplot(bclim2, aes_string(x = clim_vars[i], y = "area_ha",
                            color = "dataset", fill = "dataset")) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
                method.args = list(family = Gamma(link = "log")),
                size = 0.3, alpha = 0.3) +
    geom_point(data = bclim, mapping = aes_string(x = clim_vars[i],
                                                  y = "area_ha"),
               inherit.aes = FALSE) +
    geom_point(data = bclim[bclim$year == 2015, ],
               mapping = aes_string(x = clim_vars[i], y = "area_ha"),
               color = "red", size = 1.5,
               inherit.aes = FALSE) +
    facet_wrap(vars(vegetation_class), scales = "free_y") +
    scale_color_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
    scale_fill_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
    theme(legend.position = c(0.66, 0.15)) #+
    # scale_y_continuous(trans = "log10")

  # print(p)
  ggsave(paste("figures/fire_spread_patterns_clim_dynamic_ori_",
               clim_vars[i], ".png", sep = ""),
         plot = p,
         width = 17, height = 15, dpi = 300, units = "cm")
}


# without ci

for(i in 1:length(clim_vars)) {
  # i <- 1
  p <-
  ggplot(bclim2, aes_string(x = clim_vars[i], y = "area_ha",
                            color = "dataset", fill = "dataset")) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
                method.args = list(family = Gamma(link = "log")),
                size = 0.3, alpha = 0.3, se = FALSE) +
    geom_point(data = bclim, mapping = aes_string(x = clim_vars[i],
                                                  y = "area_ha"),
               inherit.aes = FALSE) +
    geom_point(data = bclim[bclim$year == 2015, ],
               mapping = aes_string(x = clim_vars[i], y = "area_ha"),
               color = "red", size = 1.5,
               inherit.aes = FALSE) +
    facet_wrap(vars(vegetation_class), scales = "free_y") +
    scale_color_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
    scale_fill_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
    theme(legend.position = c(0.66, 0.15)) #+
  # scale_y_continuous(trans = "log10")

  # print(p)
  ggsave(paste("figures/fire_spread_patterns_clim_dynamic_ori_",
               clim_vars[i], "_noci.png", sep = ""),
         plot = p,
         width = 17, height = 15, dpi = 300, units = "cm")
}

#
#
# # Plots at log scale
#
# for(i in 1:length(clim_vars)) {
#   # i <- 1
#   p <-
#     ggplot(bclim2, aes_string(x = clim_vars[i], y = "area_ha_annual",
#                               color = "dataset", fill = "dataset")) +
#     geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
#                 method.args = list(family = Gamma(link = "log")),
#                 size = 0.3, alpha = 0.3) +
#     geom_point(data = bclim, mapping = aes_string(x = clim_vars[i],
#                                                   y = "area_ha_annual"),
#                inherit.aes = FALSE) +
#     geom_point(data = bclim[bclim$year == 2015, ],
#                mapping = aes_string(x = clim_vars[i], y = "area_ha_annual"),
#                color = "red", size = 1.5,
#                inherit.aes = FALSE) +
#     facet_wrap(vars(vegetation_class)) +
#     scale_color_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
#     scale_fill_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
#     theme(legend.position = c(0.66, 0.15)) +
#     scale_y_continuous(trans = "log10")
#
#   print(p)
#   ggsave(paste("figures/fire_spread_patterns_clim_dynamic_log_",
#                clim_vars[i], ".png", sep = ""),
#          plot = p,
#          width = 17, height = 15, dpi = 300, units = "cm")
# }

# Plot altogether

# scale predictors to [0, 1] for facet_grid to show a good plot
bclim3 <- bclim2
for(v in c("pp", "fwi", "temp")) {
  var <- bclim2[, v]
  var_0 <- var - min(var)
  bclim3[, v] <- var_0 / max(var_0)
}

bclim2_long <- pivot_longer(
  bclim3,
  cols = all_of(which(names(bclim2) %in% c("pp", "fwi", "temp"))),
  names_to = "clim_name", values_to = "clim_value"
  )


ggplot(bclim2_long, aes(x = clim_value, y = area_ha,
                        color = dataset, fill = dataset)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              method.args = list(family = Gamma(link = "log")),
              size = 0.3, alpha = 0.3) +
  geom_point(data = bclim2_long, mapping = aes(x = clim_value,
                                               y = area_ha),
             inherit.aes = FALSE) +
  geom_point(data = bclim2_long[bclim2_long$year == 2015, ],
             mapping = aes(x = clim_value, y = area_ha),
             color = "red", size = 1.5,
             inherit.aes = FALSE) +
  facet_grid(rows = vars(vegetation_class), cols = vars(clim_name),
             scales = "free_y") +
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 270))

ggsave("figures/fire_spread_patterns_clim_dynamic_all_variables.png",
       width = 17, height = 25, dpi = 300, units = "cm")


# without ci
ggplot(bclim2_long, aes(x = clim_value, y = area_ha,
                        color = dataset, fill = dataset)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              method.args = list(family = Gamma(link = "log")),
              size = 0.3, alpha = 0.3, se = FALSE) +
  geom_point(data = bclim2_long, mapping = aes(x = clim_value,
                                               y = area_ha),
             inherit.aes = FALSE) +
  geom_point(data = bclim2_long[bclim2_long$year == 2015, ],
             mapping = aes(x = clim_value, y = area_ha),
             color = "red", size = 1.5,
             inherit.aes = FALSE) +
  facet_grid(rows = vars(vegetation_class), cols = vars(clim_name),
             scales = "free_y") +
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.5, direction = -1) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 270))

ggsave("figures/fire_spread_patterns_clim_dynamic_all_variables_noci.png",
       width = 17, height = 25, dpi = 300, units = "cm")




# Correlation between climatic (dynamic) variables ------------------------

filt_pairs <- which(names(bclim2) %in% c("pp", "temp", "fwi", "wind", "vpd"))

ppp <- GGally::ggpairs(bclim2[bclim2$vegetation_code == 2, ],
                filt_pairs) +
  theme(strip.text.y = element_text(angle = 270))

ggsave("figures/fire_spread_patterns_clim_dynamic_pairs.png",
       plot = ppp,
       width = 17, height = 17, dpi = 300, units = "cm")
