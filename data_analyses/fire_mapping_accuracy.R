library(terra)
library(sf)
library(aod)

# Extract points from point-line kmz file
get_points <- function(sf_data) {
  if(any(st_geometry_type(sf_data) == "POINT")) {
    sf_points <- sf_data[st_geometry_type(sf_data) == "POINT", ]
  } else {
    sf_points <- NULL
  }

  if(any(st_geometry_type(sf_data) == "LINESTRING")) {
    sf_lines <- sf_data[st_geometry_type(sf_data) == "LINESTRING", ]
    sf_points_from_lines <- st_cast(sf_lines, "POINT") # Extraer vértices como puntos
    sf_points <- rbind(sf_points, sf_points_from_lines)
  }

  return(vect(sf_points))
}

# directory for test points
sdir <- "accuracy_files"

# Load fire polygons
fires <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")

# Get fires names

bb <- list.files(sdir, pattern = "-burned")
uu <- list.files(sdir, pattern = "-unburned")

bnames <- sapply(bb, function(x) {
  strsplit(x, split = "-burned.kmz")[[1]][1]
}) |> unname()
unames <- sapply(uu, function(x) {
  strsplit(x, split = "-unburned.kmz")[[1]][1]
}) |> unname()

all.equal(bnames, unames) # OK
all(bnames %in% fires$fire_id) # OK

# Make table with data
etab <- data.frame(
  fire_id = bnames,
  n_burned_ge = NA,
  n_burned_in = NA,    # inside target polygon (OK)
  n_unburned_ge = NA,
  n_unburned_out = NA  # outside target polygon (OK)
)

nf <- nrow(etab)

for(i in 1:nf) {
  print(i)

  focal_fire <- fires[fires$fire_id == etab$fire_id[i], ]

  # burned points
  burned <- st_read(file.path(sdir, bb[i]))
  b_points <- get_points(burned)
  etab$n_burned_ge[i] <- nrow(b_points)
  etab$n_burned_in[i] <- relate(focal_fire, b_points, "intersects") |> sum()

  # unburned points
  unburned <- st_read(file.path(sdir, uu[i]))
  u_points <- get_points(unburned)
  etab$n_unburned_ge[i] <- nrow(u_points)
  etab$n_unburned_out[i] <- (!relate(focal_fire, u_points, "intersects")) |> sum()
}

etab$acc <- with(
  etab,
  (n_burned_in + n_unburned_out) / (n_burned_ge + n_unburned_ge)
)

etab$commission <- with(etab, (n_unburned_ge - n_unburned_out) / n_unburned_ge)
etab$omission <- with(etab, (n_burned_ge - n_burned_in) / n_burned_ge)

etab$n_error_com <- with(etab, n_unburned_ge - n_unburned_out)
etab$n_error_om <- with(etab, n_burned_ge - n_burned_in)

# Hierarchical model for global error -------------------------------------

m_com <- betabin(cbind(n_error_com, n_unburned_ge - n_error_com) ~ 1,
                 random = ~ 1, data = etab, link = "logit")
m_om <- betabin(cbind(n_error_om, n_burned_ge - n_error_om) ~ 1,
                random = ~ 1, data = etab, link = "logit")

pcom <- predict(m_com, se.fit = T, type = "link")

error_com <- c(
  "mle" = plogis(pcom$fit[1]),
  "lower" = plogis(pcom$fit[1] - qnorm(0.975) * pcom$se.fit[1]),
  "upper" = plogis(pcom$fit[1] + qnorm(0.975) * pcom$se.fit[1])
)
round(error_com * 100, 2)
# mle    lower   upper
# 0.64    0.32    1.30 (%)

pom <- predict(m_om, se.fit = T, type = "link")

error_om <- c(
  "mle" = plogis(pom$fit[1]),
  "lower" = plogis(pom$fit[1] - qnorm(0.975) * pom$se.fit[1]),
  "upper" = plogis(pom$fit[1] + qnorm(0.975) * pom$se.fit[1])
)
round(error_om * 100, 2)
# mle    lower   upper
# 3.36    2.02    5.52 (%)


# accuracy
etab$n_good <- etab$n_burned_in + etab$n_unburned_out
etab$n_total <- etab$n_burned_ge + etab$n_unburned_ge
etab$n_bad <- etab$n_total - etab$n_good

m_acc <- betabin(cbind(n_good, n_bad) ~ 1,
                 random = ~ 1, data = etab, link = "logit")

pacc <- predict(m_acc, se.fit = T, type = "link")

acc <- c(
  "mle" = plogis(pacc$fit[1]),
  "lower" = plogis(pacc$fit[1] - qnorm(0.975) * pacc$se.fit[1]),
  "upper" = plogis(pacc$fit[1] + qnorm(0.975) * pacc$se.fit[1])
)
round(acc * 100, 2)

# mle     lower   upper
# 98.06   96.92   98.78 (%)

# Export table ------------------------------------------------------------

write.csv(etab, "exports/mapping_accuracy_table.csv", row.names = F)

