library(terra)
library(tidyverse)
library(lubridate)
library(readxl)

# load excel with edited dates
dnew <- read_excel("google_earth_files/fires_data.xls")

# remove one fire, which was duplicated
dnew <- dnew[dnew$obs != "delete", ]

# merge with old database
v0 <- vect("patagonian_fires.shp")

# subset
v1 <- v0[v0$fire_id %in% dnew$fire_id, ]
length(v1) == nrow(dnew) # OK

# order both by fire_id
v1 <- v1[order(v1$fire_id), ]
dnew <- dnew[order(dnew$fire_id), ]
all.equal(v1$fire_id, dnew$fire_id) # OK

# Editar estos fire_id:
# "1999_1546963766" a "2000_..."
# "2014_-1075171770" a "2016_..."
# lo demás ya fue editado en el excel
dnew$fire_id[dnew$fire_id == "1999_1546963766"] <- "2000_1546963766"
dnew$fire_id[dnew$fire_id == "2014_-1075171770"] <- "2016_-1075171770"

# Replace all columns that were changed
v1$fire_id <- dnew$fire_id
v1$name <- dnew$name
v1$year <- dnew$year

v1$date_l <- dnew$date_l
v1$date_u <- dnew$date_u
v1$date <- dnew$date
v1$datesrc <- dnew$datesrc

v1$obs <- dnew$obs

unique(v1$datesrc)

# Check date_l < date_u
v1$date_l <- as.Date(v1$date_l)
v1$date_u <- as.Date(v1$date_u)

date_diff <- v1$date_u - v1$date_l
range(date_diff)
v1$fire_id[date_diff < 0]

v1$fire_id[date_diff > 30]
v1$fire_id[date_diff > 60]
v1$fire_id[date_diff > 200]

v1$fire_id[dnew$obs == "uncertain_year"]

summary(date_diff %>% as.numeric)
table(v1$datesrc)
plot(ecdf(date_diff))
abline(v = 35)

# maybe backward 130 d, forward = 40 d for FWI analyses



# edit patagonian eastern: ------------------------------------------------

# remove
# "1999_1546963766"
# "2014_-1075171770"
# which fire_id were changed.

e0 <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires_eastern/patagonian_fires_eastern.shp")

# keep those not present at v1:
e0 <- e0[!(e0$fire_id %in% v1$fire_id), ]
length(e0)

# remove the ones with names changed
e0 <- e0[!(e0$fire_id %in% c("1999_1546963766", "2014_-1075171770")), ]

names(e0) == names(v1)

v2 <- v1
v2$obs_date <- "revisada_el_2023-09-04"
e1 <- e0
e1$obs_date <- "revisar_con_modis"
e2 <- rbind(v2, e1)
plot(e2)
length(e2)


# Escribir ----------------------------------------------------------------

plot(v1)
plot(e2)
length(v1)
View(as.data.frame(v1))
View(as.data.frame(e2))

# writeVector(v1, "patagonian_fires.shp", overwrite = TRUE)
# writeVector(e2, "/home/ivan/Insync/patagonian_fires/patagonian_fires_eastern/patagonian_fires_eastern.shp",
#             overwrite = TRUE)

# Done, overwritten on 2023-09-04.


# Resolving problems ------------------------------------------------------

# Hubo problemas. Las fechas parece que se escribieron como millis(), y se borraron
# en la data de eastern

v0 <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")

dlnum <- v0$date_l %>% as.numeric
dunum <- v0$date_u %>% as.numeric
# están en segundos desde 1970
sec_day <- 60*60*24

dl <- as.Date(dlnum / sec_day, origin = "1970-01-01", format = "%Y-%m-%d")
du <- as.Date(dunum / sec_day, origin = "1970-01-01", format = "%Y-%m-%d")
dm <- as.Date(v0$date, format = "%Y-%m-%d")

# reasignar fechas como character
v0$date_l <- as.character(dl)
v0$date_u <- as.character(du)

# check differences
ddd <- du - dl
summary(as.numeric(ddd))

# detect problems in dates
wrong <- which(dm < dl | dm > du)
v0$date[wrong] <- as.character(dl[wrong] + (du[wrong] - dl[wrong]) / 2)

# dates are character now
summary(as.data.frame(v0))

# writeVector(v0, "patagonian_fires/tmp/patagonian_fires.shp")


# Merge again with eastern fires ------------------------------------------

e0 <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires_eastern/patagonian_fires_eastern.shp")

# keep those not present at v1:
e0 <- e0[!(e0$fire_id %in% v0$fire_id), ]
length(e0)

# remove the ones with names changed
e0 <- e0[!(e0$fire_id %in% c("1999_1546963766", "2014_-1075171770")), ]

names(e0) == names(v0)

v2 <- v0
v2$obs_date <- "revisada_el_2023-09-04"
e1 <- e0
e1$obs_date <- "revisar_con_modis"

# check dates are character
summary(as.data.frame(e1))
summary(as.data.frame(v2))

e2 <- rbind(v2, e1)
plot(e2)
length(e2)
view(as.data.frame(e2))

# write.csv(as.data.frame(e2), "patagonian_fires_eastern/data_to_edit.csv",
#           row.names = F)
# writeVector(e2, "patagonian_fires_eastern/patagonian_fires_eastern.shp")


# Add thomas edits --------------------------------------------------------

# a few fires were removed or edited in qgis.
# load excel with edited dates
dt <- read_excel("google_earth_files/fires_data_thomas.xls")
v0 <- vect("patagonian_fires.shp")

# Does Thomas info conflict modis hotspots?
fires_th_1 <- dt$fire_id[which(dt$notas_ivan == "use_thomas")]
fires_modis <- v0$fire_id[grep("modis", v0$datesrc)]

fires_conflict <- fires_modis[fires_modis %in% fires_th_1]


# Export dataset to see conflicts
fires_conf <- as.data.frame(v0[v0$fire_id %in% fires_conflict, ])
fires_conf <- fires_conf[, c("year", "fire_id", "date", "date_l", "date_u",
                             "obs", "datesrc")]
fires_conf <- fires_conf[order(fires_conf$fire_id), ]


fires_conf2 <- dt[dt$fire_id %in% fires_conflict, ]
fires_conf2 <- fires_conf2[, c("fire_id", "date_l", "date_u", "obs")]
fires_conf2 <- fires_conf2[order(fires_conf2$fire_id), ]

names(fires_conf2) <- c("fire_id", "date_l2", "date_u2", "obs2")

# fires_both <- cbind(fires_conf, fires_conf2)
# fires_both$diff_l <- as.Date(fires_both$date_l) - as.Date(fires_both$date_l2)
# fires_both$diff_u <- as.Date(fires_both$date_u) - as.Date(fires_both$date_u2)

# write.csv(fires_both, "data_compare_thomas_modis.csv", row.names = F)

# added column named use_thomas
fires_both <- read.csv("data_compare_thomas_modis.csv")

# in these fires, where use_thomas == 1, his date_l will be "date",
# the datesrc will be records/modis, and date_info will have the link provided
# by thomas.
# date_lr and date_ur will be the upper and lower provided by thomas.

# Also add info in fires from thomas that dont conflict modis.
# loop over all fires

v0$date_lr <- NA
v0$date_ur <- NA
v0$date_info <- NA

for(f in 1:nrow(v0)) {
  # f = 1
  print(f)
  fire <- v0$fire_id[f]

  # was this fire detected by thomas?
  if(fire %in% dt$fire_id) {

    # copy info
    v0$date_lr[f] <- dt$date_l[dt$fire_id == fire]
    v0$date_ur[f] <- dt$date_u[dt$fire_id == fire]
    v0$date_info[f] <- dt$date_l[dt$fire_id == fire]  ## this is a bug???
    v0$name[f] <- dt$name[dt$fire_id == fire]

    # if there is no modis info, replace datescr and date
    # (most likely date is the first)
    if(!(fire %in% fires_both$fire_id)) {
      v0$date[f] <- dt$date_l[dt$fire_id == fire]
      v0$datesrc[f] <- paste("records", v0$datesrc[f], sep = "/")
    }

    # if there is conflict, but records date is preferred, change date
    # and datesrc
    if(fire %in% fires_both$fire_id) {
      if(fires_both$use_thomas[fires_both$fire_id == fire] == 1) {
        v0$date[f] <- v0$date_lr[f]
        v0$datesrc[f] <- paste("records", v0$datesrc[f], sep = "/")
      }
    }
  }
}

View(as.data.frame(v0))


# guardado el 2023-09-08
writeVector(v0, "patagonian_fires.shp")






for(f in 1:nrow(v0)) {
  # f = 1
  print(f)
  fire <- v0$fire_id[f]

  # was this fire detected by thomas?
  if(fire %in% dt$fire_id) {

    # copy info
    v0$date_lr[f] <- dt$date_l[dt$fire_id == fire]
    v0$date_ur[f] <- dt$date_u[dt$fire_id == fire]
    v0$date_info[f] <- dt$date_l[dt$fire_id == fire]    ### bug here???
    v0$name[f] <- dt$name[dt$fire_id == fire]

    # if there is no modis info, replace datescr and date
    # (most likely date is the first)
    if(!(fire %in% fires_both$fire_id)) {
      v0$date[f] <- dt$date_l[dt$fire_id == fire]
      v0$datesrc[f] <- paste("records", v0$datesrc[f], sep = "/")
    }

    # if there is conflict, but records date is preferred, change date
    # and datesrc
    if(fire %in% fires_both$fire_id) {
      if(fires_both$use_thomas[fires_both$fire_id == fire] == 1) {
        v0$date[f] <- v0$date_lr[f]
        v0$datesrc[f] <- paste("records", v0$datesrc[f], sep = "/")
      }
    }
  }
}


dt <- read_excel("google_earth_files/fires_data_thomas.xls")
v0 <- vect("patagonian_fires.shp")

v0$date_info <- NA
for(f in 1:nrow(v0)) {
  fire <- v0$fire_id[f]
  if(fire %in% dt$fire_id) {
    v0$date_info[f] <- dt$date_info[dt$fire_id == fire]
  }
}


writeVector(v0, "patagonian_fires.shp", overwrite = TRUE)




# En patagonian_fires faltan algunos incendios ----------------------------

# como el de san ramón

v1 <- vect("patagonian_fires.shp")
v2 <- vect("backup files/patagonian_fires.shp")

ids_not <- which(!(v2$fire_id %in% v1$fire_id))

plot(v2[ids_not, ])
v2[ids_not, ]$year
v2[ids_not, ]$fire_id

names(v2) %in% names(v1) # all names in v2 are in v1

# get names from v1 absent from v2
add_names <- names(v1)[!(names(v1) %in% names(v2))]
v2$date_lr <- NA
v2$date_ur <- NA
v2$date_info <- NA

# order names in v2
v2 <- v2[, names(v1)]
all(names(v2) == names(v1))

# add missing fires to v1
v3 <- rbind(v1, v2[ids_not, ])
plot(v3)
writeVector(v3, "patagonian_fires.shp", overwrite = T)
