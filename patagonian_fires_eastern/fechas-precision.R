library(terra)
v <- vect("patagonian_fires_eastern.shp")
dlow <- as.Date(v$date_l, format = "%Y-%m-%d")
dup <- as.Date(v$date_u, format = "%Y-%m-%d")

