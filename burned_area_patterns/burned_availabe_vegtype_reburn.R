library(rgdal)       # readOGR
library(rgeos)       # gIntersection, gArea...
library(maptools)    # UnionSpatialPolygons
library(tidyverse); theme_set(theme_bw())
library(viridis)

options(scipen = 999)


# Reburn analysis ---------------------------------------------------------

reburn_wgs <- readOGR("reburn_frequency_polygons.shp")
proj_posgar <- "+proj=tmerc +lat_0=-90 +lon_0=-72 +k=1 +x_0=1500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
reburn <- spTransform(reburn_wgs, proj_posgar)

reburn_freq <- data.frame(fire_freq = 1:5, area_ha = NA)
for(i in 1:nrow(reburn_freq)) {
  reburn_freq$area_ha[i] <- gArea(reb_vec[reb_vec$fire_freq == i, ]) * 1e-4
}
reburn_freq$prob <- reburn_freq$area_ha / sum(reburn_freq$area_ha)
# Same result as GEE.
write.csv(reburn_freq, "data_reburn_frequency.csv")



# Vegetation in available and burn areas ----------------------------------

fires_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/patagonian_fires/patagonian_fires.shp")
veg_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/vegetation_valdivian_reclassified/vegetation_valdivian_reclassified.shp")


fires <- spTransform(fires_wgs, proj_posgar)
veg <- spTransform(veg_wgs, proj_posgar)
veg <- veg[!is.na(veg$clss_nm), ]

# gArea(veg[65, ]) * 1e-4
# gArea(veg[1719, ]) * 1e-4 #    los na son re peques.
# gArea(veg[1731, ]) * 1e-4



# sa <- spTransform(veg_wgs, proj_posgar)
veg$class
veg_data <- data.frame(vegetation_code = 1:8,
                       vegetation_class = NA, 
                       area_available_ha = NA,
                       area_burned_ha = NA)

for(i in 1:nrow(veg_data)) {
  print(i)
  # subset vegetation
  veg_temp <- veg[veg$clss_nm == i, ]
  
  # available area
  veg_data$area_available_ha[i] <- gArea(veg_temp) * 1e-4
  
  # burned area
  v_int <- gIntersection(polygons(veg_temp), polygons(fires))
  veg_data$area_burned_ha[i] <- gArea(v_int) * 1e-4
}
# write.csv(veg_data, "data_burned_and_available_area_by_vegetation.csv")


sum(veg_data$area_burned_ha); sum(reburn_freq$area_ha) 
# veg tiene como 10000 ha menos. El problema no era GEE entonces?

# Intersección entre fuegos y area de estudio?
sa_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/study_area/study_area.shp")
sa <- spTransform(sa_wgs, proj_posgar)

gArea(gIntersection(sa, fires)) * 1e-4 # 126593.8, exactamente igual que con veg...

# Será que el reburn mira algo mal?
gArea(reburn) * 1e-4 # 136581.6
# En GEE, fires_clipped mide 129383.858 ha, más cercano que veg acá.
# Debe ser que al rasterizar a escala de 30 m muchos "pixeles" con casi cero
# fuego se convierten en burned, y ahí aumenta el área.