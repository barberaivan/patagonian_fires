library(terra)
library(ggplot2)
library(sf)

vegim <- rast("comarca_vegetation.tif")
poly <- vect("comarca_polygon.shp")
fires <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")


poly <- project(poly, crs(vegim))
fires <- project(fires, crs(vegim))

# plot(vegim)
# plot(poly, add = T)
# plot(fires, add = T, col = "red")

## write polygon in kml
# writeVector(poly, "comarca_polygon.kml")


# Burnable area in region of interest -------------------------------------

vegim2 <- mask(vegim, poly)

mask_cond <- vegim2 < 8
mask_cond[is.na(mask_cond)] <- FALSE  # Asegurarse de que los NA se conviertan en FALSE

# Especificamos unit = "km" para obtener el área en km^2.
areas <- cellSize(vegim2, unit = "km")

# 3. Multiplicar la máscara (convertida a numérico: TRUE = 1, FALSE = 0) por el raster de áreas.
areas_sel <- areas * as.numeric(mask_cond)

# 4. Sumar las áreas de las celdas seleccionadas.
# global() devuelve una tabla; tomamos el primer valor.
burnable <- global(areas_sel, sum, na.rm = TRUE)[1,1]
burnable # 2044.749 km2 quemable


mask_cond <- vegim2 < 20 # all, not only burnable.
mask_cond[is.na(mask_cond)] <- FALSE  # Asegurarse de que los NA se conviertan en FALSE
areas_sel <- areas * as.numeric(mask_cond)
total_area <- global(areas_sel, sum, na.rm = TRUE)[1,1]
total_area # total 2535.798 km2

burnable / total_area # 80.63533 % quemable

# Proportion burned in the whole period -----------------------------------

fires_subset <- fires[poly, ]
fires_subset$dummy <- 1
fires_dissolved <- aggregate(fires_subset, by = "dummy")

# Crop to the polygon extent (this limits the features to a similar area)
fires_cropped <- crop(fires_dissolved, poly)

# Now compute the geometric intersection (this should clip the geometries)
fires_clipped <- intersect(fires_cropped, poly)

# Compute the area, for example:
clipped_area <- expanse(fires_clipped, unit = "km")
clipped_area # 360.3203 km2
clipped_area / burnable * 100 # 17.62173 %

# Proportion burned by year -----------------------------------------------

burned_time <- data.frame(
  year = 1999:2022,
  burned_ha = 0,
  perc = 0,
  nfires = 0
)
nyears <- nrow(burned_time)

for(i in 1:nyears) {
  fires_local <- fires_subset[fires_subset$year == burned_time$year[i], ]
  if(nrow(fires_local) > 0) {
    fires_diss <- aggregate(fires_local, by = "dummy")
    fires_crop <- crop(fires_diss, poly)
    fires_clip <- intersect(fires_crop, poly)

    burned_time$burned_ha[i] <- expanse(fires_clip, unit = "ha")  # By default, units will be those of the CRS (often m²)
    burned_time$perc[i] <- burned_time$burned_ha[i] / (burnable * 100) * 100

    burned_time$nfires[i] <- nrow(fires_local)
  }
}
sum(burned_time$perc) # 18.62977 %

ggplot(burned_time, aes(year, perc)) +
  geom_point()

write.csv(burned_time, "buned_area_ts.csv", row.names = F)

summary(burned_time$perc)

plot(fires_subset)
as.data.frame(fires_subset[order(fires_subset$area_ha), ])[, c("year", "area_ha", "fire_id")]

aggregate

nrow(fires_subset)

mean(burned_time$nfires)
