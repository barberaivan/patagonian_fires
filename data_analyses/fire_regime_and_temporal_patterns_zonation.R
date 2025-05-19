# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(dplyr)

library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr
library(ggh4x)    # varying strip theme for veg_types and all together

library(scales)   # log scale
library(lubridate)

library(terra)
# library(rgdal)    # readOGR
# library(rgeos)    # gIntersection, gArea...
# library(maptools) # UnionSpatialPolygons     ## all replaced by terra

library(mgcv)
library(DHARMa)

# custom ggplot theme -----------------------------------------------------

# from https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_mine <- function() {
  # font <- "Liberation Sans"   #assign font family up front
  marg <- 2 # figure margin in mm
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
        # family = font,            #set font family
        size = 16,                #set font size
        #face = 'bold',            #bold typeface
        hjust = -0.1,                #left align
        vjust = 1),
      
      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        # family = font,            #font family
        size = 12),
      
      # para separar el eje y de los nros
      axis.title.y = element_text(
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),
      
      axis.text = element_text(              #axis text
        # family = font,            #axis family
        size = 9),                #font size
      
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9), #, family = font
      
      strip.text = element_text(size = 12, color = "white"), # family = font, 
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),
      
      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

theme_set(theme_mine())

# Data --------------------------------------------------------------------

# fires
fires <- vect(file.path("..", "patagonian_fires/patagonian_fires.shp"))

# Regions
regions <- vect(file.path(
  "data", 
  "study_area_partitions_data",
  "study_area_avialable_partitioned.shp"
))
regions$region <- c("Sur", "Centro-Sur", "Centro-Norte", "Norte")

# Regions data
regions_data <- read.csv(
  file.path(
    "data", 
    "study_area_partitions_data",
    "fire_freq_mean_area_by_quadrant.csv"
  )
)
regions_data$region <- factor(
  c("Sur", "Centro-Sur", "Centro-Norte", "Norte"),
  levels = c("Norte", "Centro-Norte", "Centro-Sur", "Sur")
)
regions_data <- regions_data[order(regions_data$region), ]

# Annual burned area (ha) by region --------------------------------------

# Ensure CRS match
fires <- project(fires, crs(regions))

# Define target years
years <- 1999:2022

nreg <- nrow(regions)

# Initialize list to hold results
burned_by_region_list <- vector("list", nreg)

# Loop over each region
for (i in 1:nreg) {
  
  reg <- regions[i, ]
  reg_name <- reg$region
  
  # Intersect fires with this region
  fires_in_reg <- intersect(fires, reg)
  
  if (nrow(fires_in_reg) > 0) {
    # Compute area in ha for each polygon inside region
    fires_in_reg$area_ha <- expanse(fires_in_reg, unit = "ha")
    
    # Aggregate by year
    tab <- aggregate(area_ha ~ year, as.data.frame(fires_in_reg), sum)
    
  } else {
    # No fires in this region: create empty table
    tab <- data.frame(year = integer(0), area_ha = numeric(0))
  }
  
  # Ensure all years are present (fill 0s for missing)
  full_tab <- data.frame(year = years) %>%
    left_join(tab, by = "year") %>%
    mutate(
      region = reg_name,
      area_ha = ifelse(is.na(area_ha), 0, area_ha)
    )
  
  burned_by_region_list[[i]] <- full_tab
}

# Combine all into one data.frame
burned_area_table <- do.call(rbind, burned_by_region_list)

# sort
burned_area_table$region <- factor(
  burned_area_table$region,
  levels = c("Norte", "Centro-Norte", "Centro-Sur", "Sur")
)
oo <- order(burned_area_table$region, burned_area_table$year)
burned_area_table <- burned_area_table[oo, ]


# Burned proportion -------------------------------------------------------

bprop <- left_join(
  burned_area_table, 
  regions_data[, c("region", "area_ha")],
  by = "region"
)

names(bprop) <- c("year", "burned_ha", "region", "available_ha")
bprop$prop <- bprop$burned_ha / bprop$available_ha * 100


# Time series -------------------------------------------------------------

# replace zeroes
bprop$prop2 <- bprop$prop
bprop$prop2[bprop$prop == 0] <- 1e-4

ggplot(bprop, aes(year, prop2)) + 
  geom_line() +
  geom_point() +
  facet_wrap(vars(region), nrow = 4, scales = "free") +
  # scale_y_continuous(trans = "log10")
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")))



# Fire rotation period ----------------------------------------------------

regions_data$frp <- 1 / (regions_data$mean / 24)
barplot(frp ~ region, data = regions_data,
        ylab = "Período de rotación de fuego (años)",
        xlab = "Región")

regions_data$annual_prob <- regions_data$mean / 24 * 100
barplot(annual_prob ~ region, data = regions_data,
        ylab = "Probabilidad de quema anual (%)",
        xlab = "Región")

write.csv(
  regions_data,
  file.path(
    "data", 
    "study_area_partitions_data",
    "fire_freq_mean_area_by_quadrant_computed.csv"
  ),
  row.names = F
)

knitr::kable(
  regions_data[, c("region", "frp", "annual_prob")],
  format = "latex",
  booktabs = F,
  digits = 2,
  caption = "Incidencia del fuego por región."
)

# \begin{table}
# \caption{Incidencia del fuego por región.}
# \centering
# \begin{tabular}[t]{l|r|r}
# \hline
# Región & Tiempo de rotación (años) & Probabilidad de quema anual (\%) \\
# \hline
# Norte & 1259 & 0.08\\
# Centro-Norte & 722 & 0.14\\
# Centro-Sur & 129 & 0.78\\
# Sur & 557 & 0.18\\
# \hline
# \end{tabular}
# \end{table}