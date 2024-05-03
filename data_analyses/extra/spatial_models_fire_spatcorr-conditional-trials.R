# in 2, predictions for topographic variables co-vary with other topographic
# variables.
# Here, I would do the same, but first, models for separate vegetation types
# will be fitted as a mixed model.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(terra)
library(mgcv)
# library(INLA)
# library(fmesher)
library(HDInterval) # HDI

library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr
library(ggh4x)    # varying strip theme for veg_types and all together
library(GGally)
library(tidyterra)

library(circular)

# Functions ---------------------------------------------------------------

softmax <- function(x) exp(x) / sum(exp(x)) # inverse multinomial-logit link

# to compute mean and 95 % CI from samples.
mean_ci <- function(x, name = "p_") {
  qs <- quantile(x, probs = c(0.025, 0.975), method = 8)
  tmp <- c(mean(x), qs) %>% unname
  names(tmp) <- paste(name, c("mean", "lower", "upper"), sep = "")
  return(tmp)
}

# just the ci
ci <- function(x, name = "p_") {
  qs <- quantile(x, probs = c(0.025, 0.975), method = 8)
  names(qs) <- paste(name, c("lower", "upper"), sep = "")
  return(qs)
}

# cicular mean for aspect
# from
# https://stackoverflow.com/questions/32404222/circular-mean-in-r
mean_circular_raw <- function (x) { # takes angle in radians
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  return(circmean)
}

mean_circular_deg <- function (x) { # takes angle in degrees
  # x <- c(180, 160, 170)
  conv <- 2 * pi / 360 # degrees to radians factor
  mm_raw <- mean_circular_raw(conv * x) / conv
  mm <- (mm_raw + 360) %% 360
  return(mm)
}

# function to make new data varying only one predictor.
make_newdata <- function(varying = "elevation", data, l = 150, HDI = TRUE,
                         prob = 0.95) {

  var_names <- colnames(data)

  values_list_mean <- lapply(var_names, function(v) {
    if(v != varying) {
      res <- mean(data[, v])
      if(v == "aspect") {
        res <- mean_circular_deg(data[, v])
      }
    } else {
      if(HDI) {
        hh <- hdi(data[, v], prob)
        res <- seq(hh["lower"], hh["upper"], length.out = l)
      } else {
        res <- seq(min(data[, v]), max(data[, v]), length.out = l)
      }
    }
    return(res)
  })

  names(values_list_mean) <- var_names

  new_data <- expand.grid(values_list_mean)

  # add columns indicating which is the varying predictor and which are its
  # values, useful to plot later
  new_data$varying_var <- varying
  new_data$varying_val <- new_data[, varying]

  return(new_data)
}

r2bern <- function(p) {
  var_mu <- var(p)
  var_y <- mean(p * (1 - p))
  return(var_mu / (var_mu + var_y))
}

normalize <- function(x) x / sum(x)

# function to compute something similar to a continuous weighted derivative,
# but for a categorical predictor. It computes the average of the pairwise
# differences in the response (y) across all pairs of given levels, but weighted
# by the each pairs' probability. The weights of pairs are computed from the
# marginal -unnormalized- weights (w), as a product. If weights are not given,
# they are assumed constant.
categorical_effect <- function(y, w = NULL) {
  # TEST
  # y <- rnorm(6, sd = 3); w <- normalize(runif(6))

  if(is.null(w)) w <- rep(1 / length(y), length(y))
  w <- normalize(w)

  # pairwise diff
  ydiff <- outer(y, y, "-") %>% abs
  wpairs <- outer(w, w, "*")

  dvec <- ydiff[lower.tri(ydiff)]
  wvec <- wpairs[lower.tri(wpairs)] %>% normalize

  return(sum(dvec * wvec))
}

# simpler categorical effect: absolute difference between predicted values
# and the mean predicted value, averaged over categories. The mean prediction
# and the average over categories are weighted by the abundance of each
# category.
categorical_effect2 <- function(y, w = NULL) {
  # TEST
  # y <- rnorm(6, sd = 3); w <- normalize(runif(6))

  if(is.null(w)) w <- rep(1 / length(y), length(y))
  w <- normalize(w)

  # mean prediction
  m <- sum(y * w)

  # mean difference
  d <- sum(abs(y - m) * w)

  return(d)
}


# Data --------------------------------------------------------------------

v0 <- vect("data_and_files/data_spatial_variables_mindist_800m.shp")
pp <- rast("data_and_files/pp_atlas-climatico/geonode_precipitacion_anual.tif")
tt <- rast("data_and_files/temp_atlas-climatico/13_tempera_anual.tif")
wb <- rast("data_and_files/pp_atlas-climatico/geonode_p_ep_balanceh.tif")

# rename worldClim temperature, so temp is used for the atlas
v0$temp_wc <- v0$temp

v2 <- extract(pp, v0)
v0$pp <- v2$geonode_precipitacion_anual

t2 <- extract(tt, v0)
v0$temp <- t2$`13_tempera_anual`

v3 <- extract(wb, v0)
v0$wb <- v3$geonode_p_ep_balanceh

v <- project(v0, "EPSG:5343")
out_ids <- is.na(v$pp) #| is.na(v$temp)
v <- v[!out_ids, ]
# veg_type distribution
burned_veg0 <- read.csv("data_and_files/data_burned_and_available_area_by_vegetation_dryforest2.csv")


# Tidy data ---------------------------------------------------------------

# inputs to assess spatial corr

# apply a 1500 m buffer, find neighbours, count them and count the number
# of burned ones. Then, use the count or the proportion of burned neighbours
# as predictor.

vneighs <- nearby(v, distance = 2000)

# compute number of burned neighbors and the number of them
bb <- matrix(NA, nrow(v), 2)
colnames(bb) <- c("burned", "total")

for(p in 1:nrow(v)) {
  if(p %% 100 == 0) print(p)
  neighs <- c(vneighs[vneighs[, 1] == p, ], vneighs[vneighs[, 2] == p, ])
  neighs <- unique(neighs)
  neighs <- neighs[neighs != p]
  bb[p, "total"] <- length(neighs)
  bb[p, "burned"] <- v$burned[neighs] %>% sum
}

barplot(table(bb[, "total"]))
barplot(table(bb[, "burned"]))

v$neighs_n <- bb[, "total"]
v$neighs_burned <- bb[, "burned"]
v$neighs_bprop <- bb[, "burned"] / bb[, "total"]

d <- as.data.frame(values(v))
coords <- crds(v)
d$lat <- coords[, "y"]
d$long <- coords[, "x"]
d <- rename(d, ndvi_mean = ndvi_mean_, ndvi_max = ndvi_max_m)
names(d)

# remove unburnable
d <- d[d$vegetation > 1, ]
nrow(d) # 8278

# turn distances into km
d$dist_human <- d$dist_human / 1000
d$dist_roads <- d$dist_roads / 1000

# recode vegetation and remove unuseful categories
burned_veg <- read.csv("data_and_files/data_burned_and_available_area_by_vegetation_dryforest2.csv")
d$vegetation_code <- d$vegetation

# code dry forest A as subalpine forest
d$vegetation_code[d$vegetation_code == 4] <- 2

## Do not transform plantation and prairie into other stuff.
## let them be NA

data <- left_join(d, burned_veg[, c("vegetation_code", "vegetation_class")],
                  by = "vegetation_code")

data$vegetation_class[data$vegetation_class == "Dry forest B"] <- "Dry forest"
data$vegetation_class[data$vegetation_class == "Steppe and grassland"] <- "Grassland"
data$vegetation_class[data$vegetation_class == "Anthropogenic prairie and shrubland"] <- "Anthropogenic prairie"
unique(data$vegetation_class)

veg_levels <- c(
  "Wet forest",
  "Subalpine forest",
  "Dry forest",
  "Shrubland",
  "Grassland",
  "Anthropogenic prairie",
  "Plantation"
)

veg_labels <- c(
  "Wet forest",
  "Subalpine\nforest",
  "Dry forest",
  "Shrubland",
  "Grassland",
  "Anthropogenic\nprairie",
  "Plantation"
)

# levels without low-abundance classes
veg_levels_sub <- c(
  "Wet forest",
  "Subalpine forest",
  "Dry forest",
  "Shrubland",
  "Grassland"
)

veg_labels_sub <- c(
  "Wet forest",
  "Subalpine\nforest",
  "Dry forest",
  "Shrubland",
  "Grassland"
)

# classes for which fire models are to be estimated
veg_model <- c(
  "All vegetation\ntypes",
  "Wet forest",
  "Subalpine\nforest",
  "Dry forest",
  "Shrubland",
  "Grassland"
)

veg_model2 <- c(
  "All vegetation types",
  "Wet forest",
  "Subalpine forest",
  "Dry forest",
  "Shrubland",
  "Grassland"
)

# with numbers for each plot
veg_model_num <- c(
  "(1) All vegetation\ntypes",
  "(2) Wet forest",
  "(3) Subalpine\nforest",
  "(4) Dry forest",
  "(5) Shrubland",
  "(6) Grassland"
)


data$vegetation_class <- factor(data$vegetation_class, levels = veg_levels,
                                labels = veg_labels)

flat_coords <- crds(v)


predictors <- c("elevation", "slope", "aspect", "TPI2k",
                "ndvi_mean", "pp", "dist_human", "dist_roads")

names_frame <- data.frame(variable = predictors,
                          var_name = c("(A) Elevation (m a.s.l.)",
                                       "(B) Slope (°)",
                                       "(C) Aspect",
                                       "(D) Topographic\nposition",
                                       "(B) NDVI max",
                                       "(A) Precipitation\n(mm / year)",
                                       "(C) Distance from\nhuman settlements (km)",
                                       "(D) Distance from\nroads (km)"))

names_frame$name2 <- c("Elevation (m a.s.l.)",
                       "Slope (°)",
                       "Aspect",
                       "Topographic\nposition",
                       "NDVI max",
                       "Precipitation\n(mm / year)",
                       "Distance from\nhuman settlements (km)",
                       "Distance from\nroads (km)")

names_frame$name3 <- c("Elevation",
                       "Slope",
                       "Aspect",
                       "Topographic\nposition",
                       "NDVI max",
                       "Precipitation",
                       "Distance from\nhuman settlements",
                       "Distance from\nroads")

names_frame$name4 <- c("Elevation",
                       "Slope",
                       "Aspect",
                       "Top. position",
                       "NDVI max",
                       "Precipitation",
                       "Dist. settlements",
                       "Dist. roads")


# data for vegetation model (removes low-abundance classes)
data_veg <- data[data$vegetation_class %in% veg_labels_sub, ]
data_veg$vegetation_class <- factor(as.character(data_veg$vegetation_class),
                                    levels = veg_labels_sub)
data_veg$veg_num <- as.numeric(data_veg$vegetation_class) - 1
# unique(data_veg$veg_num)
# numeric categories start at zero because mgcv likes it this way.


P <- length(predictors)
V <- length(veg_model)


# spatial data

crds_data <- as.data.frame(crds(v))
data_coords <- cbind(data, crds_data)


# Load randomized fires

# ff <- list.files("data_and_files")
# ff <- ff[grep("data_randomized_fires_mindist_800m_", ff)]
# dsim0 <- do.call("rbind", lapply(ff, function(x) {
#   read.csv(paste("data_and_files", x, sep = "/"))
# }))
# nrow(dsim0) / length(unique(dsim0$replicate)) == nrow(d) # OK
# names(dsim0)
# unique(aggregate(burned_sim ~ replicate, dsim0, length)[, "burned_sim"])
# # OK
# dsim <- matrix(dsim0$burned_sim, ncol = max(dsim0$replicate))
# saveRDS(dsim, "data_and_files/data_randomized_fires_mindist_800m_500_matrix.rds")
dsim <- readRDS("data_and_files/data_randomized_fires_mindist_800m_500_matrix.rds")
dsim <- dsim[!out_ids, ]
dim(dsim) # 8228 data points, not 8278

# veg_type distribution and burned proportion

# merge dry forest A with subalpine
burned_veg <- burned_veg0[burned_veg0$vegetation_code != 4, ]
burned_veg[burned_veg$vegetation_class == "Subalpine forest",
           c("area_available_ha", "area_burned_ha")] <-
  colSums(
    burned_veg0[burned_veg0$vegetation_class %in% c("Subalpine forest", "Dry forest A"),
                c("area_available_ha", "area_burned_ha")]
  )
burned_veg$vegetation_class[burned_veg$vegetation_class == "Dry forest B"] <- "Dry forest"
burned_veg$vegetation_class[burned_veg$vegetation_class == "Steppe and grassland"] <- "Grassland"
burned_veg$vegetation_class[burned_veg$vegetation_class == "Anthropogenic prairie and shrubland"] <- "Anthropogenic prairie"


# compute burned proportion
data_vd <- burned_veg[-1, ]
data_vd$vegetation_class <- factor(data_vd$vegetation_class, levels = veg_levels,
                                   labels = veg_labels)
data_vd$prop <- data_vd$area_burned_ha / data_vd$area_available_ha
data_vd$veg_dist <- normalize(data_vd$area_available_ha)



# Spatial correlation models ----------------------------------------------


# no correlation
m0 <- gam(burned ~ s(elevation, k = 3, bs = "cr"), method = "REML",
          family = "binomial", data = d)
plot(m0)

m1 <- gam(burned ~
            s(elevation, k = 3, bs = "cr") +
            s(neighs_bprop, k = 3, bs = "cr"),
          method = "REML", family = "binomial", data = d)
plot(m1)


m2 <- gam(burned ~
            s(elevation, k = 3, bs = "cr") +
            neighs_burned,
          method = "REML", family = "binomial", data = d)
plot(m2)

m3 <- gam(burned ~
            s(elevation, k = 3, bs = "cs") +
            s(slope, k = 3, bs = "cs") +
            s(TPI2k, k = 3, bs = "cs") +
            s(ndvi_mean, k = 3, bs = "cs") +
            s(pp, k = 3, bs = "cs") +
            s(dist_human, k = 3, bs = "cs") +
            neighs_bprop,#neighs_bin,
          method = "REML", family = "binomial", data = d)
m4 <- gam(burned ~
            s(elevation, k = 3, bs = "cs") +
            s(slope, k = 3, bs = "cs") +
            s(TPI2k, k = 3, bs = "cs") +
            s(ndvi_mean, k = 3, bs = "cs") +
            s(pp, k = 3, bs = "cs") +
            s(dist_human, k = 3, bs = "cs"),
          method = "REML", family = "binomial", data = d)
par(mfrow = c(3, 2))
plot(m3)
plot(m4)
par(mfrow = c(1, 1))

# no cambian tanto los resultados cuando condicionamos en los vecinos.
# Se achican los efectos, pero tampoco tanto.

mean(d$neighs_burned)
mean(fitted(m3))

# corr espacial?
library(DHARMa)

rm4 <- simulateResiduals(m4, integerResponse = T, n = 1500)
rm3 <- simulateResiduals(m3, integerResponse = T, n = 1500)

d2 <- d[!is.na(d$neighs_bprop), ]

d2$res3 <- rm3$scaledResiduals
d$res4 <- rm4$scaledResiduals


d2dist <- dist(d2[, c("lat", "long")])
ddist <- dist(d[, c("lat", "long")])

d2rdist <- dist(d2$res3)
drdist <- dist(d$res4)

disdat3 <- data.frame(resdist = d2rdist[lower.tri(d2rdist)],
                      spdist = d2dist[lower.tri(d2dist)])

disdat4 <- data.frame(resdist = drdist[lower.tri(drdist)],
                      spdist = ddist[lower.tri(ddist)])

disdat3 <- disdat3[disdat3$spdist < 100000, ]
disdat4 <- disdat4[disdat4$spdist < 100000, ]

# order and aggregate
disdat3 <- disdat3[order(disdat3$spdist), ]
disdat3$g <- rep(1:500, each = floor(nrow(disdat3 / 500)))
disdat3 <- disdat3[disdat3$g < 500, ]
dagg3 <- aggregate(cbind(resdist, spdist) ~ g, disdat3, mean)

disdat4 <- disdat4[order(disdat4$spdist), ]
disdat4$g <- rep(1:500, each = ceiling(nrow(disdat4 / 500)))
disdat4 <- disdat4[disdat4$g < 500, ]
dagg4 <- aggregate(cbind(resdist, spdist) ~ g, didat4, mean)

ggplot(dagg3, aes(spdist, resdist)) + geom_smooth() + geom_point()


# Falta ver esto. Le cuesta.