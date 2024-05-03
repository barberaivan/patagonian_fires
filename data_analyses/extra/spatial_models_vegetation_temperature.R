# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(terra)
library(mgcv)
library(HDInterval) # HDI

library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr
library(ggh4x)    # varying strip theme for veg_types and all together

library(circular)
library(vegan)    # pca
library(ggdensity) # geom_hdr


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

r2bern <- function(p) {
  var_mu <- var(p)
  var_y <- mean(p * (1 - p))
  return(var_mu / (var_mu + var_y))
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
# function to make new data varying only one predictor.
make_newdata <- function(varying = "temp", data, l = 150, HDI = TRUE,
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

# Data --------------------------------------------------------------------

v0 <- vect("data_and_files/data_spatial_variables_mindist_800m.shp")
pp <- rast("data_and_files/pp_atlas-climatico/geonode_precipitacion_anual.tif")
tt <- rast("data_and_files/temp_atlas-climatico/13_tempera_anual.tif")

v2 <- extract(pp, v0)
v0$pp <- v2$geonode_precipitacion_anual
t2 <- extract(tt, v0)
v0$temp <- t2$`13_tempera_anual`

v <- project(v0, "EPSG:5343")
v <- v[!(is.na(v$pp) | is.na(v$temp)), ]

# Tidy data ---------------------------------------------------------------

d <- as.data.frame(values(v))
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

## Do not transform plantation and prairie into other stuff for the separate
## analyses.

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
  "Anthropogenic prairie",
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


data$vegetation_class <- factor(data$vegetation_class, levels = veg_levels,
                                labels = veg_labels)

flat_coords <- crds(v)


predictors <- c("temp", "slope", "aspect", "TPI2k",
                "pp", "dist_human", "dist_roads") ## ndvi removed!!!

names_frame <- data.frame(variable = predictors,
                          var_name = c("(A) Temperature (°C)",
                                       "(B) Slope (°)",
                                       "(C) Aspect",
                                       "(D) Topographic\nposition",
                                       "(E) Precipitation\n(mm / year)",
                                       "(F) Distance from\nhuman settlements (km)",
                                       "(G) Distance from\nroads (km)"),
                          var_name_simple = c(
                            "Temperature",
                            "Slope",
                            "Aspect",
                            "Topographic position",
                            "Precipitation",
                            "Distance from human settlements",
                            "Distance from roads"
                          ))


# data for vegetation model (removes low-abundance classes)
data_veg <- data[data$vegetation_class %in% veg_labels_sub, ]
data_veg$vegetation_class <- factor(as.character(data_veg$vegetation_class),
                                    levels = veg_labels_sub)
data_veg$veg_num <- as.numeric(data_veg$vegetation_class) - 1
# unique(data_veg$veg_num)
# numeric categories start at zero because mgcv likes it this way.

k_gam <- 4

GGally::ggpairs(data[, c("temp", "elevation", "slope", "TPI2k")])

# Vegetation multivariate model -------------------------------------------

# wet forest is set as the reference.
k <- length(veg_levels_sub) - 1

# we need k (4) linear predictors, which will be for
veg_levels_sub[-1]
# we can simplify the smooths for dry forest.

vegmod <- gam(
  list(veg_num ~ s(temp, k = k_gam, bs = "cr") + s(slope, k = k_gam, bs = "cr") +
                 s(aspect, k = k_gam, bs = "cc") + s(TPI2k, k = k_gam, bs = "cr") +
                 s(pp, k = k_gam, bs = "cr") +
                 s(dist_human, k = k_gam, bs = "cr") + s(dist_roads, k = k_gam, bs = "cr"),

               ~ s(temp, k = k_gam, bs = "cr") + s(slope, k = k_gam, bs = "cr") + # dry forest
                 s(aspect, k = k_gam, bs = "cc") + s(TPI2k, k = k_gam, bs = "cr") +
                 s(pp, k = k_gam, bs = "cr") +
                 s(dist_human, k = k_gam, bs = "cr") + s(dist_roads, k = k_gam, bs = "cr"),

               ~ s(temp, k = k_gam, bs = "cr") + s(slope, k = k_gam, bs = "cr") +
                 s(aspect, k = k_gam, bs = "cc") + s(TPI2k, k = k_gam, bs = "cr") +
                 s(pp, k = k_gam, bs = "cr") +
                 s(dist_human, k = k_gam, bs = "cr") + s(dist_roads, k = k_gam, bs = "cr"),

               ~ s(temp, k = k_gam, bs = "cr") + s(slope, k = k_gam, bs = "cr") +
                 s(aspect, k = k_gam, bs = "cc") + s(TPI2k, k = k_gam, bs = "cr") +
                 s(pp, k = k_gam, bs = "cr") +
                 s(dist_human, k = k_gam, bs = "cr") + s(dist_roads, k = k_gam, bs = "cr")),
  knots = list(aspect = c(0, 360)),
  family = multinom(K = k), data = data_veg
)
saveRDS(vegmod, "data_and_files/vegetation_model_gam_temperature.rds")
## se toma su tiempo
vegmod <- readRDS("data_and_files/vegetation_model_gam_temperature.rds")

# predictions

nsim <- 10000
K <- length(veg_labels_sub)

# new data_veg for predictions.

ndveg0 <- do.call("rbind", lapply(predictors, function(var) {
  make_newdata(var, data = data_veg[, predictors])
}))
ndveg <- ndveg0
# make distances co-vary
ndveg$dist_roads[ndveg0$varying_var == "dist_human"] <-
  ndveg$dist_roads[ndveg0$varying_var == "dist_roads"]
ndveg$dist_human[ndveg0$varying_var == "dist_roads"] <-
  ndveg$dist_human[ndveg0$varying_var == "dist_human"]

# make linear predictor matrix:
lpmat <- predict(vegmod, ndveg, type = "lpmatrix")
lpmat <- lpmat[, 1:(ncol(lpmat) / k)] # all matrices are the same

np <- nrow(ndveg)

# create an array for simulated predicted probabilites.
pred_arr <- array(NA, dim = c(np, K, nsim),   # it's an array with nsim tpm
                  dimnames = list(
                    case = 1:np,
                    veg = 1:K,
                    sim = 1:nsim
                  ))

# the model parameter estimates, under hipothetical infinite replicates of the
# study, are assumed to follow a multivariate normal distribution, with the
# MLE as their mean and a certain variance-covariance matrix.
# We simulate parameter vectors from this distribution:
set.seed(123)
coef_samples <- rmvn(nsim, coef(vegmod),
                     V = vcov(vegmod, unconditional = TRUE, freq = F)) %>% t

# each linear predictor has the same number of coefficients, so
coef_ids <- rep(1:k, each = nrow(coef_samples) / k)

for(i in 1:nsim) {
  # i = 1
  # order the coefficients vector in matrix form. It enters by column, as the
  # default for matrix()
  # i = 1 # just to test the loop
  coef_temp <- coef_samples[, i]
  # turn potential Inf into something reasonable (exp(700) is finite, but a bit
  # above return Inf and breaks the loop)
  coef_temp[coef_temp > 700] <- 700

  coef_mat <- matrix(coef_temp, ncol = k)
  # add reference linear predictor
  linpred <- cbind(rep(0, np), lpmat %*% coef_mat)
  pred_arr[, , i] <- apply(linpred, 1, softmax) %>% t
}

# compute ci and longanize array:
pred_long_ci <- as.data.frame.table(apply(pred_arr, 1:2, ci))
head(pred_long_ci)

# get the MLE in a similar way, but before match the names.
prob_hat <- predict(vegmod, ndveg, type = "response")
dimnames(prob_hat) <- dimnames(pred_arr)[1:2]
pred_long_mle <- as.data.frame.table(prob_hat)
pred_long_mle <- cbind(data.frame("Var1" = rep("p_mle", np)),
                       pred_long_mle)

# merge with ci data_veg:
pred_long <- rbind(pred_long_ci, pred_long_mle)
# a bit wider for ggplot:
pred <- pivot_wider(pred_long, names_from = "Var1", values_from = "Freq")

# merge with predictor variables
ndveg$case <- 1:np
pred$case <- as.numeric(as.character(pred$case))
pred_veg <- left_join(pred, ndveg, by = "case")

pred_veg$veg_cat <- factor(veg_levels_sub[pred_veg$veg],
                           levels = veg_levels_sub,
                           labels = veg_labels_sub)

ids_dist <- grep("dist", pred_veg$varying_var)
pred_veg$varying_val[ids_dist] <- pred_veg$varying_val[ids_dist] / 1000

# plot
P <- length(predictors)
vegplot_list <- vector("list", P)
for(i in 1:P) {
  # i = 1
  vv <- predictors[i]

  dd <- pred_veg[pred_veg$varying_var == vv, ]

  p <-
  ggplot(dd, aes(y = p_mle, ymin = p_lower, ymax = p_upper,
                 x = varying_val, group = veg_cat, color = veg_cat,
                 fill = veg_cat)) +
    geom_ribbon(alpha = 0.25, color = NA) +
    geom_line() +
    scale_color_viridis(option = "A", discrete = TRUE, end = 0.9) +
    scale_fill_viridis(option = "A", discrete = TRUE, end = 0.9) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8)) +
    xlab(names_frame$var_name[i]) +
    # ylim(0, 100) +
    ylab("Expected cover (%)") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1))

  if(!(i %in% c(1, 5))) {p <- p + theme(axis.title.y = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank())}

  if(vv == "aspect") {
    p <- p + scale_x_continuous(breaks = c(0, 90, 180, 270, 360),
                                labels = c("N", "E", "S", "W", "N"))
  }

  vegplot_list[[i]] <- p
  # print(p)
}
# son hermosos

vegplot <- egg::ggarrange(plots = vegplot_list, nrow = 2, ncol = 4)

leg <- ggpubr::get_legend(vegplot_list[[1]] + theme(legend.position = "bottom",
                                                    legend.title = element_blank(),
                                                    legend.text = element_text(size = 8),))

vegplot_leg <- grid.arrange(vegplot, leg, nrow = 2, heights = c(20, 2))

# ggsave("figures/spatial patterns - veg as a function of fire drivers_temperature.png",
#        plot = vegplot_leg,
#        width = 16, height = 12, units = "cm")


# corr entre params?
cc <- cov2cor(vcov(vegmod, unconditional = T))
cc[lower.tri(cc)] %>% range
# cc[lower.tri(cc)] %>% hist(breaks = 30)


# r2 for conditional model
fff <- apply(cbind(rep(0, nrow(data_veg)), fitted(vegmod)), 1, softmax) %>% t
class_wise_r2_full <- apply(fff, 2, r2bern)


# Vegetation univariate models --------------------------------------------

# wet forest is set as the reference.
K <- length(veg_levels_sub)
k <- K-1
# we need k (4) linear predictors, which will be for
nsim <- 10000 # I'd recommend using 10000 for a real application
nr <- 150 # number of points to predict in sequence

## Loop over predictors until we get a nice plot. save model and predictions.
P <- length(predictors)
vegplot_list <- vector("list", P)
vegplot_list_dens <- vector("list", P)

# r2 will be in %
r2table <- names_frame
r2table$r2_weighted <- NA
r2table$r2_plain <- NA
mm <- matrix(NA, P, K)
colnames(mm) <- veg_levels_sub
r2table$r2_veg <- mm

# coordinates for r2
ranges <- apply(as.matrix(data_veg[, predictors]), 2, range)
ranges[, 7:P] <- ranges[, 7:P]
widths <- apply(ranges, 2, diff)
r2table$x <- ranges[1, ] + widths * 0.2
r2table$y <- 0.9
veg_props <- as.numeric(table(data_veg$vegetation_class) / nrow(data_veg))

names(vegplot_list) <- names(vegplot_list_dens) <- predictors

for(p in predictors) {
  print(p)
  # reduce dataset, and name predictor as "x"
  dp <- data_veg[, c(p, "veg_num")]
  names(dp) <- c("x", "veg_num")

  # model settings
  kn = k_gam; bs = "cr"; knots = NULL
  if(p == "aspect") {bs = "cc"; knots = list(x = c(0, 360)) }

  # fit model
  vm <- gam(
    list(veg_num ~ s(x, k = kn, bs = bs),
                 ~ s(x, k = kn, bs = bs),
                 ~ s(x, k = kn, bs = bs),
                 ~ s(x, k = kn, bs = bs)),
    knots = knots,
    family = multinom(K = k), data = dp
  )
  mname <- paste("data_and_files/vegetation_model_gam", "_", p, "_temperature.rds", sep = "")
  saveRDS(vm, mname)
  vm <- readRDS(mname)

  # r2 for the model: bayesian r2 for every class, weigthed by class abundance.

  fff <- apply(cbind(rep(0, nrow(dp)), fitted(vm)), 1, softmax) %>% t
  class_wise_r2 <- apply(fff, 2, r2bern)

  rptex <- paste(format(round(mean(class_wise_r2) * 100, 2), nsmall = 2), "%")
  rwtex <- paste(format(round(sum(class_wise_r2 * veg_props) * 100, 2), nsmall = 2), "%")

  r2table$r2_plain[r2table$variable == p] <- rptex
  r2table$r2_weighted[r2table$variable == p] <- rwtex
  r2table$r2_veg[r2table$variable == p, ] <- class_wise_r2 * 100
  # predictions

  # new data
  ndveg <- data.frame(x = seq(min(dp$x), max(dp$x), length.out = nr))

  #  linear predictor matrix:
  lpmat <- predict(vm, ndveg, type = "lpmatrix")
  lpmat <- lpmat[, 1:(ncol(lpmat) / k)] # all matrices are the same
  np <- nrow(ndveg)

  # create an array for simulated predicted probabilites.
  pred_arr <- array(NA, dim = c(np, K, nsim),
                    dimnames = list(
                      case = 1:np,
                      veg = 1:K,
                      sim = 1:nsim
                    ))
  # Simulate parameter vectors
  set.seed(123)
  coef_samples <- rmvn(nsim, coef(vm),
                       V = vcov(vm, unconditional = TRUE, freq = F)) %>% t
  for(i in 1:nsim) {
    # i = 1
    # order the coefficients vector in matrix form. It enters by column, as the
    # default for matrix()
    # i = 1 # just to test the loop
    coef_temp <- coef_samples[, i]
    # turn potential Inf into something reasonable (exp(700) is finite, but a bit
    # above return Inf and breaks the loop)
    coef_temp[coef_temp > 700] <- 700

    coef_mat <- matrix(coef_temp, ncol = k)
    # add reference linear predictor
    linpred <- cbind(rep(0, np), lpmat %*% coef_mat)
    pred_arr[, , i] <- apply(linpred, 1, softmax) %>% t
  }
  arrname <- paste("data_and_files/vegetation_model_gam_predictions_array_marg", "_", p, "_temperature.rds", sep = "")
  saveRDS(pred_arr, arrname)
  pred_arr <- readRDS(arrname)

  # compute ci and longanize array:
  pred_long_ci <- as.data.frame.table(apply(pred_arr, 1:2, ci))

  # get the MLE in a similar way, but before match the names.
  prob_hat <- predict(vm, ndveg, type = "response")
  dimnames(prob_hat) <- dimnames(pred_arr)[1:2]
  pred_long_mle <- as.data.frame.table(prob_hat)
  pred_long_mle <- cbind(data.frame("Var1" = rep("p_mle", np)),
                         pred_long_mle)

  # merge with ci data:
  pred_long <- rbind(pred_long_ci, pred_long_mle)
  # a bit wider for ggplot:
  pred <- pivot_wider(pred_long, names_from = "Var1", values_from = "Freq")

  # merge with predictor variables
  ndveg$case <- 1:np
  pred$case <- as.numeric(as.character(pred$case))
  pred_veg <- left_join(pred, ndveg, by = "case")

  pred_veg$veg_cat <- factor(veg_levels_sub[pred_veg$veg],
                             levels = veg_levels_sub,
                             labels = veg_labels_sub)

  # plot
  plotcito <-
    ggplot(pred_veg, aes(y = p_mle, ymin = p_lower, ymax = p_upper,
                         x = x, group = veg_cat, color = veg_cat,
                         fill = veg_cat)) +
    geom_ribbon(alpha = 0.25, color = NA) +
    geom_line() +
    geom_text(data = r2table[r2table$variable == p, , drop = F],
              mapping = aes(x, y, label = r2_weighted), size = 2.8,
              inherit.aes = F) +
    scale_color_viridis(option = "A", discrete = TRUE, end = 0.85) +
    scale_fill_viridis(option = "A", discrete = TRUE, end = 0.85) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.margin = margin(t = 2, r = 2, b = 4, l = 2 ,unit = "mm")) +
    xlab(names_frame$var_name[names_frame$variable == p]) +
    ylab("Expected cover (%)") +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, 1), expand = c(0.01, 0))
  # plotcito
  if(!(p %in% names_frame$variable[c(1, 5)])) {
    plotcito <- plotcito + theme(axis.title.y = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank())
  }

  if(p == "aspect") {
    plotcito <- plotcito +
      scale_x_continuous(breaks = c(0, 90, 180, 270, 360),
                         labels = c("N", "E", "S", "W", "N"))
  }

  vegplot_list[[p]] <- plotcito

  # minimal density plot
  dens <- density(dp$x, from = min(dp$x), to = max(dp$x), n = 2^10)
  dens_approx <- approx(dens$x, dens$y, xout = ndveg$x)
  ndveg$dens <- dens_approx$y

  if(p == "aspect") {
    xc <- circular(dp$x, type = "angles", units = "degrees", template = "geographic")
    dens <- density.circular(xc, bw = 10, n = 2 ^ 10)

    # circular hace cosas raras: su seq, dens$x empieza en 90 y va decreciendo
    # (gira a contrarreloj, donde el último número, cercano a 90, es -270 )
    # Entonces, convertimos los números negativos a lo que son + 360.
    densx_cor <- dens$x
    densx_cor[dens$x < 0] <- dens$x[dens$x < 0] + 360

    dens_approx <- approx(densx_cor, dens$y, xout = ndveg$x)
    ndveg$dens <- dens_approx$y
  }

  dplot <- ggplot(ndveg, aes(x = x, y = dens, ymax = dens, ymin = 0)) +
    geom_ribbon(color = NA, alpha = 0.05) +
    geom_line(linewidth = 0.25, alpha = 0.7) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(1, 0, -4, 0, "mm"))
  # dplot
  vegplot_list_dens[[p]] <- dplot
}

# ggarrange(plots = vegplot_list, ncol = 4)

# listas con densities

# make blank density to put above legend, because a blank plot did not work well.
blank <- ggplot(data_veg, aes(x = ndvi_mean)) +
  geom_density(color = "white", alpha = 0, fill = "white") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(1, 0, -4, 0, "mm"))

# get legend and cast as ggplot object
leg <- ggpubr::get_legend(vegplot_list[[1]] + theme(legend.position = "right",
                                                    legend.title = element_blank(),
                                                    legend.text = element_text(size = 8),
                                                    legend.box.margin = margin(t = 2, r = 2, b = 4, l = 2 ,unit = "mm")))
leg <- ggpubr::as_ggplot(leg)

# make a list for ggarrange
aaa <- c(vegplot_list_dens[1:4], vegplot_list[1:4],
         vegplot_list_dens[5:7], list(blank),
         vegplot_list[5:7], list(leg))

# predictions
vegplot_marg_dens <- ggarrange(plots = aaa, nrow = 4,
                               heights = c(0.35, 1, 0.35, 1))

ggsave("figures/S01) vegetation and drivers marginal_temperature.png",
       plot = vegplot_marg_dens,
       width = 16, height = 14, units = "cm")


# Veg models R2 table ----------------------------------------------------

# View(r2table)
r2table_veg <- cbind(r2table$var_name_simple, as.data.frame(round(r2table$r2_veg, 2)))
names(r2table_veg)[1] <- "Variable"

class_wise_r2_full_m <- matrix(round(class_wise_r2_full * 100, 2), nrow = 1)
colnames(class_wise_r2_full_m) <- veg_levels_sub

fffull <- cbind(data.frame(Variable = "Multiple regression"),
                as.data.frame(class_wise_r2_full_m))

r2_export <- rbind(r2table_veg, fffull)
write.csv(r2_export, "data_and_files/vegetation_models_r2.csv")


# Conditional predictions using focal mean --------------------------------

# Using the conditional model to make partial predictions has a problem: when
# the non-plotted covariates are fixed at their mean, they move towards zero the
# probabilities of vegetation types that belong to the extreme environmental
# conditions. To avoid this we can fix the non-plotted covariates at the mean
# of every vegetation type. It's like saying, how do you vary if you are in your
# range? This allows to see how different drivers affect each class, but not
# how a single driver controls all classes.

# copy code from predictions above.

nsim <- 10000
K <- length(veg_labels_sub)

# new data_veg for predictions.

ndveg0 <- do.call("rbind", lapply(predictors, function(var) {
  do.call("rbind", lapply(veg_labels_sub, function(veg) {
    # make new data filtering vegetation and predictor
    ##  var = "temp"; veg = "Shrubland" # test

    if(var != "aspect") {
      nd <- make_newdata(var,
                         data_veg[data_veg$vegetation_class == veg, predictors],
                         HDI = TRUE, prob = 0.98)
    } else {
      nd <- make_newdata(var,
                         data_veg[data_veg$vegetation_class == veg, predictors],
                         HDI = FALSE)
    }

    nd$focal_veg <- factor(veg, levels = veg_labels_sub)
    return(nd)
  }))
}))
# nrow(ndveg0) # 150 * 8 * 5

# change the reference temp for subalpine forest, so its probabilities
# can vary instead of being flat at 1.
filt <- ndveg0$varying_var != "temp" &
        ndveg0$focal_veg == "Subalpine\nforest"
ndveg0[filt, "temp"] <- 1150 # it was 1322

ndveg <- ndveg0
# make distances co-vary
ndveg$dist_roads[ndveg0$varying_var == "dist_human"] <-
  ndveg$dist_roads[ndveg0$varying_var == "dist_roads"]
ndveg$dist_human[ndveg0$varying_var == "dist_roads"] <-
  ndveg$dist_human[ndveg0$varying_var == "dist_human"]

# make linear predictor matrix:
lpmat <- predict(vegmod, ndveg, type = "lpmatrix")
lpmat <- lpmat[, 1:(ncol(lpmat) / k)] # all matrices are the same

np <- nrow(ndveg)

# create an array for simulated predicted probabilites.
pred_arr <- array(NA, dim = c(np, K, nsim),   # it's an array with nsim tpm
                  dimnames = list(
                    case = 1:np,
                    veg = 1:K,
                    sim = 1:nsim
                  ))

set.seed(123)
coef_samples <- rmvn(nsim, coef(vegmod),
                     V = vcov(vegmod, unconditional = TRUE, freq = F)) %>% t

for(i in 1:nsim) {
  print(i)
  # i = 1
  # order the coefficients vector in matrix form. It enters by column, as the
  # default for matrix()
  # i = 1 # just to test the loop
  coef_temp <- coef_samples[, i]
  # turn potential Inf into something reasonable (exp(700) is finite, but a bit
  # above return Inf and breaks the loop)
  coef_temp[coef_temp > 700] <- 700

  coef_mat <- matrix(coef_temp, ncol = k)
  # add reference linear predictor
  linpred <- cbind(rep(0, np), lpmat %*% coef_mat)
  pred_arr[, , i] <- apply(linpred, 1, softmax) %>% t
}

# compute ci and longanize array:
pred_ci <- apply(pred_arr, 1:2, ci)
names(dimnames(pred_ci))[1] <- "summary"
pred_ci <- aperm(pred_ci, c(2, 3, 1))
saveRDS(pred_ci, "data_and_files/vegetation_pred_ci_conditional_focal_temperature.rds")
pred_ci <- readRDS("data_and_files/vegetation_pred_ci_conditional_focal_temperature.rds")

# compute mle prediction
prob_hat <- predict(vegmod, ndveg, type = "response")

# Select only the columns corresponding to the focal vegetation
pred_focal <- matrix(NA, nrow(pred_ci), 3)
colnames(pred_focal) <- c("p_lower", "p_upper", "p_mle")

for(v in 1:K) {
  # v = 2
  rows <- ndveg$focal_veg == veg_labels_sub[v]
  pred_focal[rows, 1:2] <- pred_ci[rows, v, ]
  pred_focal[rows, 3] <- prob_hat[rows, v]
}

# data for prediction
pred_veg <- cbind(ndveg, as.data.frame(pred_focal))

# x breaks for each predictor
xbreaks <- list(
  "temp" = c(400, 1000, 1600),#c(500, 1000, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = seq(0, 1, by = 0.5),#c(0, 0.5, 1),
  "pp" = c(800, 1400, 2000),# seq(500, 2000, by = 500)#c(1000, 2000)
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")


# plot
P <- length(predictors)
pl_cond <- vector("list", P)
pl_dens <- vector("list", P)

for(i in 1:P) {
  # i = 1
  vv <- predictors[i]

  dd <- pred_veg[pred_veg$varying_var == vv, ]

  pp <-
    ggplot(dd, aes(y = p_mle, ymin = p_lower, ymax = p_upper,
                   x = varying_val, group = focal_veg, color = focal_veg,
                   fill = focal_veg)) +
    geom_ribbon(alpha = 0.25, color = NA) +
    geom_line() +
    scale_color_viridis(option = "A", discrete = TRUE, end = 0.9) +
    scale_fill_viridis(option = "A", discrete = TRUE, end = 0.9) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8)) +
    xlab(names_frame$var_name[i]) +
    # ylim(0, 100) +
    ylab("Expected cover (%)") +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, 1), expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(min(data_veg[, vv]),
                                  max(data_veg[, vv])))

  if(!(i %in% c(1, 5))) {pp <- pp + theme(axis.title.y = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank())}

  if(vv == "aspect") {
    pp <- pp + scale_x_continuous(breaks = c(0, 90, 180, 270, 360),
                                  limits = c(0, 360),
                                  labels = c("N", "E", "S", "W", "N"))
  }

  pl_cond[[i]] <- pp
}

# densities
for(i in 1:P) {
  # print(i)
  # i = 1
  vv <- predictors[i]
  dp <- data.frame(x = data_veg[, vv])

  dens <- density(dp$x, from = min(dp$x), to = max(dp$x), n = 2^10)
  xseq <- seq(min(dp$x), max(dp$x), length.out = 150)
  dens_approx <- approx(dens$x, dens$y, xout = xseq)
  datadens <- data.frame(dens = dens_approx$y,
                         x = xseq)

  if(vv == "aspect") {
    xc <- circular(dp$x, type = "angles", units = "degrees", template = "geographic")
    dens <- density.circular(xc, bw = 10, n = 2 ^ 10)

    # circular hace cosas raras: su seq, dens$x empieza en 90 y va decreciendo
    # (gira a contrarreloj, donde el último número, cercano a 90, es -270 )
    # Entonces, convertimos los números negativos a lo que son + 360.
    densx_cor <- dens$x
    densx_cor[dens$x < 0] <- dens$x[dens$x < 0] + 360
    xseq <- seq(0, 360, length.out = 150)
    dens_approx <- approx(densx_cor, dens$y, xout = xseq)
    datadens <- data.frame(dens = dens_approx$y,
                           x = xseq)
  }

  dplot <- ggplot(datadens, aes(x = x, y = dens, ymax = dens, ymin = 0)) +
    geom_ribbon(color = NA, alpha = 0.05) +
    geom_line(linewidth = 0.25, alpha = 0.7) +
    theme_minimal() +
    scale_x_continuous(limits = c(min(xseq), max(xseq))) +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(1, 0, -4, 0, "mm"))
  # print(dplot)
  pl_dens[[i]] <- dplot
}

vegplot_cond <- egg::ggarrange(plots = pl_cond, nrow = 2, ncol = 4)
vegdens <- c(pl_dens[1:4], pl_cond[1:4], pl_dens[5:8], pl_cond[5:8])
# vegplot_cond_dens <- egg::ggarrange(plots = vegdens,
#                                     heights = c(0.35, 1, 0.35, 1),
#                                     nrow = 4, ncol = 4, draw = F)
# leg <- ggpubr::get_legend(pl_cond[[1]] + theme(legend.position = "bottom",
#                                                legend.title = element_blank(),
#                                                legend.text = element_text(size = 8),))
#
# p2 <- grid.arrange(vegplot_cond_dens, leg, nrow = 2,
#                                       heights = c(40, 2))

# ggsave("figures/spatial patterns - veg as a function of fire drivers - cond focal with density_temperature.png",
#        plot = p2,
#        width = 16, height = 14, units = "cm")


## Similar plot using facet_wrap

pl_facet <- vector("list", P)

for(i in 1:P) {
  # i = 1
  vv <- predictors[i]
  dd <- pred_veg[pred_veg$varying_var == vv, ]

  # compute means by veg type to in the plot
  mmeans <- aggregate(data_veg[, vv] ~ vegetation_class, data_veg, mean)
  names(mmeans) <- c("focal_veg", "x")
  mmeans$veg_focal <- factor(mmeans$focal_veg, levels = veg_labels_sub)
  # change reference value for temp at subalpine forest
  if(vv == "temp") {
    mmeans[mmeans$veg_focal == "Subalpine\nforest", "x"] <- 1150
  }


  pp <-
    ggplot(dd, aes(y = p_mle, ymin = p_lower, ymax = p_upper,
                   x = varying_val, group = focal_veg, color = focal_veg,
                   fill = focal_veg)) +
    # mean line
    geom_vline(data = mmeans, aes(xintercept = x, color = focal_veg),
               linewidth = 0.5, linetype = 3) +

    # predictions
    geom_ribbon(alpha = 0.4, color = NA) +
    geom_line() +
    scale_color_viridis(option = "A", discrete = TRUE, end = 0.8) +
    scale_fill_viridis(option = "A", discrete = TRUE, end = 0.8) +
    facet_wrap(vars(focal_veg), ncol = 1, strip.position = "right") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          strip.background = element_rect(fill = "white", color = "white")) +
    xlab(names_frame$var_name[i]) +
    # ylim(0, 100) +
    ylab("Expected cover (%)") +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, 1), expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(min(data_veg[, vv]),
                                  max(data_veg[, vv])),
                       breaks = xbreaks[[vv]])
  # pp

  # remove y axis in >1 column
  if(i != 1) {
    pp <- pp + theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
  }
  if(i < P) {
    pp <- pp + theme(strip.background = element_blank(),
                     strip.text = element_blank())
  }

  # pp
  if(vv == "aspect") {
    pp <- pp + scale_x_continuous(breaks = c(0, 90, 180, 270, 360),
                                  limits = c(0, 360),
                                  labels = c("N", "E", "S", "W", "N"))
  }

  pl_facet[[i]] <- pp
}

vegsep <- egg::ggarrange(plots = pl_facet, ncol = P)
vegsep_dens <- egg::ggarrange(plots = c(pl_dens, pl_facet), ncol = P, nrow = 2,
                              heights = c(1, 7))

ggsave("figures/S02) vegetation drivers conditional_temperature.png",
       plot = vegsep_dens,
       width = 24, height = 14, units = "cm")

# Conditional predictions using focal mean - co-varying predictor ----------

# Similar as above, but co-varying correlated predictors, based on the
# topographic, distance and ndvi models fitted in the fire analyses script.
# However, in the case of subalpine forests we use the global model for
# topographic conditions. Otherwise, temp would be predicted too high,
# and we could not observe the effect of other predictors.

# copy code from predictions above.
nsim <- 10000
K <- length(veg_labels_sub)

# new data_veg for predictions.

ndveg0 <- do.call("rbind", lapply(predictors, function(var) {
  do.call("rbind", lapply(veg_labels_sub, function(veg) {
    # make new data filtering vegetation and predictor
    ##  var = "temp"; veg = "Shrubland" # test

    if(var != "aspect") {
      nd <- make_newdata(var,
                         data_veg[data_veg$vegetation_class == veg, predictors],
                         HDI = TRUE, prob = 0.98)
    } else {
      nd <- make_newdata(var,
                         data_veg[data_veg$vegetation_class == veg, predictors],
                         HDI = FALSE)
    }

    nd$focal_veg <- factor(veg, levels = veg_labels_sub)
    return(nd)
  }))
}))
# nrow(ndveg0) # 150 * 8 * 5

# change the reference temp for subalpine forest, so its probabilities
# can vary instead of being flat at 1.
filt <- ndveg0$varying_var != "temp" &
  ndveg0$focal_veg == "Subalpine\nforest"
ndveg0[filt, "temp"] <- 1150 # it was 1322

ndveg <- ndveg0


# Edit newdata considering correlations.

# Load required models for predictors
mtopo <- readRDS("data_and_files/models_topography_temperature.rds")
mdist <- readRDS("data_and_files/models_distances_temperature.rds")

topo <- c("temp", "slope", "TPI2k")
distv <- predictors[grep("dist", predictors)]

# combinations of topographic variables
combs_topo <- expand.grid(response = topo, predictor = topo)
combs_topo <- combs_topo[combs_topo$response != combs_topo$predictor, ]

# models are fitted with vegetation class, not "focal veg"
ndveg$vegetation_class <- ndveg$focal_veg

for(veg in veg_labels_sub) {
  # veg = "Shrubland"

  # define model type for topography, using the global for subalpine.
  model_type <- ifelse(veg == "Subalpine\nforest", "all", "veg")

  # topography
  for(predictor in topo) {
    # predictor <- "temp"
    rr <- ndveg$varying_var == predictor & ndveg$focal_veg == veg
    responses <- combs_topo$response[combs_topo$predictor == predictor]
    for(response in responses) {
      # response = "slope"
      prediction <- predict(mtopo[[predictor]][[response]][[model_type]],
                            newdata = ndveg[rr, ],
                            type = "response")
      ndveg[rr, response] <- prediction
    }
  }

  # distances
  for(predictor in distv) {
    # predictor <- "dist_human"
    rr <- ndveg$varying_var == predictor & ndveg$focal_veg == veg
    response <- distv[distv != predictor]
    prediction <- predict(mdist[[predictor]][[response]][["veg"]],
                          newdata = ndveg[rr, ],
                          type = "response")
    ndveg[rr, response] <- prediction
  }
}

# check
# View(ndveg)

# make linear predictor matrix:
lpmat <- predict(vegmod, ndveg, type = "lpmatrix")
lpmat <- lpmat[, 1:(ncol(lpmat) / k)] # all matrices are the same

np <- nrow(ndveg)

# create an array for simulated predicted probabilites.
pred_arr <- array(NA, dim = c(np, K, nsim),   # it's an array with nsim tpm
                  dimnames = list(
                    case = 1:np,
                    veg = 1:K,
                    sim = 1:nsim
                  ))

set.seed(123)
coef_samples <- rmvn(nsim, coef(vegmod),
                     V = vcov(vegmod, unconditional = TRUE, freq = F)) %>% t

for(i in 1:nsim) {
  if(i %% 50 == 0) print(i)
  # i = 1
  # order the coefficients vector in matrix form. It enters by column, as the
  # default for matrix()
  # i = 1 # just to test the loop
  coef_temp <- coef_samples[, i]
  # turn potential Inf into something reasonable (exp(700) is finite, but a bit
  # above return Inf and breaks the loop)
  coef_temp[coef_temp > 700] <- 700

  coef_mat <- matrix(coef_temp, ncol = k)
  # add reference linear predictor
  linpred <- cbind(rep(0, np), lpmat %*% coef_mat)
  pred_arr[, , i] <- apply(linpred, 1, softmax) %>% t
}

# compute ci and longanize array:
pred_ci <- apply(pred_arr, 1:2, ci)
names(dimnames(pred_ci))[1] <- "summary"
pred_ci <- aperm(pred_ci, c(2, 3, 1))
saveRDS(pred_ci, "data_and_files/vegetation_pred_ci_conditional_focal_cov_temperature.rds")
pred_ci <- readRDS("data_and_files/vegetation_pred_ci_conditional_focal_cov_temperature.rds")

# compute mle prediction
prob_hat <- predict(vegmod, ndveg, type = "response")

# Select only the columns corresponding to the focal vegetation
pred_focal <- matrix(NA, nrow(pred_ci), 3)
colnames(pred_focal) <- c("p_lower", "p_upper", "p_mle")

for(v in 1:K) {
  # v = 2
  rows <- ndveg$focal_veg == veg_labels_sub[v]
  pred_focal[rows, 1:2] <- pred_ci[rows, v, ]
  pred_focal[rows, 3] <- prob_hat[rows, v]
}

# data for prediction
pred_veg <- cbind(ndveg, as.data.frame(pred_focal))

# x breaks for each predictor
xbreaks <- list(
  "temp" = c(400, 1000, 1600),#c(500, 1000, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = seq(0, 1, by = 0.5),#c(0, 0.5, 1),
  "pp" = c(800, 1400, 2000),# seq(500, 2000, by = 500)#c(1000, 2000)
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")


# plot
P <- length(predictors)
pl_dens <- vector("list", P)

# densities
for(i in 1:P) {
  # print(i)
  # i = 1
  vv <- predictors[i]
  dp <- data.frame(x = data_veg[, vv])

  dens <- density(dp$x, from = min(dp$x), to = max(dp$x), n = 2^10)
  xseq <- seq(min(dp$x), max(dp$x), length.out = 150)
  dens_approx <- approx(dens$x, dens$y, xout = xseq)
  datadens <- data.frame(dens = dens_approx$y,
                         x = xseq)

  if(vv == "aspect") {
    xc <- circular(dp$x, type = "angles", units = "degrees", template = "geographic")
    dens <- density.circular(xc, bw = 10, n = 2 ^ 10)

    # circular hace cosas raras: su seq, dens$x empieza en 90 y va decreciendo
    # (gira a contrarreloj, donde el último número, cercano a 90, es -270 )
    # Entonces, convertimos los números negativos a lo que son + 360.
    densx_cor <- dens$x
    densx_cor[dens$x < 0] <- dens$x[dens$x < 0] + 360
    xseq <- seq(0, 360, length.out = 150)
    dens_approx <- approx(densx_cor, dens$y, xout = xseq)
    datadens <- data.frame(dens = dens_approx$y,
                           x = xseq)
  }

  dplot <- ggplot(datadens, aes(x = x, y = dens, ymax = dens, ymin = 0)) +
    geom_ribbon(color = NA, alpha = 0.05) +
    geom_line(linewidth = 0.25, alpha = 0.7) +
    theme_minimal() +
    scale_x_continuous(limits = c(min(xseq), max(xseq))) +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(1, 0, -4, 0, "mm"))
  # print(dplot)
  pl_dens[[i]] <- dplot
}

pl_facet <- vector("list", P)

for(i in 1:P) {
  # i = 1
  vv <- predictors[i]
  dd <- pred_veg[pred_veg$varying_var == vv, ]

  # compute means by veg type to in the plot
  mmeans <- aggregate(data_veg[, vv] ~ vegetation_class, data_veg, mean)
  names(mmeans) <- c("focal_veg", "x")
  mmeans$veg_focal <- factor(mmeans$focal_veg, levels = veg_labels_sub)
  # change reference value for temp at subalpine forest
  if(vv == "temp") {
    mmeans[mmeans$veg_focal == "Subalpine\nforest", "x"] <- 1150
  }


  pp <-
    ggplot(dd, aes(y = p_mle, ymin = p_lower, ymax = p_upper,
                   x = varying_val, group = focal_veg, color = focal_veg,
                   fill = focal_veg)) +
    # mean line
    geom_vline(data = mmeans, aes(xintercept = x, color = focal_veg),
               linewidth = 0.5, linetype = 3) +

    # predictions
    geom_ribbon(alpha = 0.4, color = NA) +
    geom_line() +
    scale_color_viridis(option = "A", discrete = TRUE, end = 0.8) +
    scale_fill_viridis(option = "A", discrete = TRUE, end = 0.8) +
    facet_wrap(vars(focal_veg), ncol = 1, strip.position = "right") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          strip.background = element_rect(fill = "white", color = "white")) +
    xlab(names_frame$var_name[i]) +
    # ylim(0, 100) +
    ylab("Expected cover (%)") +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, 1), expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(min(data_veg[, vv]),
                                  max(data_veg[, vv])),
                       breaks = xbreaks[[vv]])
  # pp

  # remove y axis in >1 column
  if(i != 1) {
    pp <- pp + theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
  }
  if(i < P) {
    pp <- pp + theme(strip.background = element_blank(),
                     strip.text = element_blank())
  }

  # pp
  if(vv == "aspect") {
    pp <- pp + scale_x_continuous(breaks = c(0, 90, 180, 270, 360),
                                  limits = c(0, 360),
                                  labels = c("N", "E", "S", "W", "N"))
  }

  pl_facet[[i]] <- pp
}

vegsep <- egg::ggarrange(plots = pl_facet, ncol = P)
vegsep_dens <- egg::ggarrange(plots = c(pl_dens, pl_facet), ncol = P, nrow = 2,
                              heights = c(1, 7))

ggsave("figures/S02) vegetation drivers conditional_cov_temperature.png",
       plot = vegsep_dens,
       width = 24, height = 14, units = "cm")

# Ordination ---------------------------------------------------------------

dord <- data_veg[, predictors]
dord$aspect <- cos(data_veg$aspect * pi / 180)
plot(dord$aspect ~ data_veg$aspect)
pc <- prcomp(dord,
             center = TRUE,
             scale. = TRUE)
summary(pc)
print(pc)

library(ggbiplot)
g <- ggbiplot(pc,
              alpha = 0.0,
              obs.scale = 1,
              var.scale = 1,
              groups = data_veg$vegetation_class,
              ellipse = TRUE, ellipse.level = 0.5, ellipse.alpha = 0.1,
              circle = F,
              ellipse.prob = 0.95)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'vertical',
               legend.position = 'right')
print(g)



# Relationships between drivers --------------------------------------------

pp <- ggplot(data_veg, aes(x = pp, y = temp)) +
  geom_hdr(probs = c(0.99, 0.95, 0.8, 0.5, 0.2, 0.05)) +
  geom_point(alpha = 0) +
  theme(legend.position = "bottom")
pp
ggExtra::ggMarginal(pp)

ggplot(data_veg, aes(x = pp, y = temp)) +
  geom_hdr(probs = c(0.99, 0.95, 0.8, 0.5, 0.2, 0.05)) +
  facet_wrap(vars(vegetation_class)) +
  theme(legend.position = c(0.85, 0.25))

ggplot(data_veg, aes(x = pp, y = temp, color = vegetation_class,
                     linetype = vegetation_class)) +
  geom_hdr_lines(probs = c(0.90), alpha = 1, linewidth = 0.7) +
  scale_color_viridis(discrete = T, end = 0.6, option = "A") +
  ylim(100, 2000) +
  xlim(250, 2000) +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("90 % probability highest density regions")

ggplot(data_veg, aes(x = pp, y = temp, color = vegetation_class)) +
  geom_hdr_lines(probs = c(0.90), alpha = 1, linewidth = 0.7) +
  scale_color_viridis(discrete = T, end = 0.6, option = "A") +
  ylim(100, 2000) +
  xlim(250, 2000) +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("90 % probability highest density regions")

ggplot(data_veg, aes(x = pp, y = temp, color = vegetation_class)) +
  geom_hdr_lines(probs = c(0.75), alpha = 1, linewidth = 0.7) +
  scale_color_viridis(discrete = T, end = 0.6, option = "A") +
  ylim(100, 2000) +
  xlim(250, 2000) +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("75 % probability highest density regions")
