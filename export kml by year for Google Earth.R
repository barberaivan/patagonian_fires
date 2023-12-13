library(tidyverse)
library(rgdal)       # readOGR

fires <- readOGR("patagonian_fires/patagonian_fires.shp")

# save by year
years <- na.omit(unique(fires@data$year))

# For a correct display of names in Google Earth:

# https://stackoverflow.com/questions/57415464/label-kml-features-using-st-write
# https://michaelminn.net/tutorials/r-kml/index.html

# (from the latter)
# # Write
# writeOGR(utm_zones, "utm_zones.kml", "Universal Transverse Mercator Zones", driver="KML",
#          layer_options=c(NameField = "NAME", DescriptionField = "NAME"))

fires2 <- fires
names(fires2@data)

# rename
# if name is NA, make it fire_id
for(f in 1:nrow(fires)) {
  fires2$name[f] <- ifelse(is.na(fires2$name[f]), fires2$fire_id[f], fires2$name[f])
}

# View(fires2@data) # OK

# order by date
fires2 <- fires2[order(fires2$year, fires2$date), ]

# write csv
# write.csv(fires2@data, "patagonian_fires/google_earth_files/fires_data.csv")

for(y in years) {
  #y = 1999
  writeOGR(fires2[fires2@data$year == y, ],
           dsn = paste("patagonian_fires/google_earth_files/",
                       y, ".kml", sep = ""),
           layer = as.character(y),
           layer_options = c(NameField = "name2"),
           driver = "KML")
}
# En Google Earth fire_id se llama "name", para que los nombres se vean en
# la lista deplegable.



# Do the same for patagonian_fires_eastern --------------------------------

fires <- readOGR("patagonian_fires_eastern/patagonian_fires_eastern.shp")

fires2 <- fires

# rename
# if name is NA, make it fire_id
for(f in 1:nrow(fires)) {
  fires2$name[f] <- ifelse(is.na(fires2$name[f]), fires2$fire_id[f], fires2$name[f])
}

# write.csv(fires2@data, "patagonian_fires_eastern/data to edit.csv")

# save by year, only the ones to review date
fires3 <- fires2[fires2$obs_date == "revisar_con_modis", ]
years <- na.omit(unique(fires3@data$year))

for(y in years) {
  #y = 1999
  writeOGR(fires3[fires3@data$year == y, ],
           dsn = paste("patagonian_fires_eastern/google_earth_files/",
                       y, ".kml", sep = ""),
           layer = as.character(y),
           layer_options = c(NameField = "name2"),
           driver = "KML")
}




