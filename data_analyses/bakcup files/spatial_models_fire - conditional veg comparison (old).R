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

v0 <- vect("data/data_spatial_variables_mindist_800m.shp")
pp <- rast("data/pp_atlas-climatico/geonode_precipitacion_anual.tif")
tt <- rast("data/temp_atlas-climatico/13_tempera_anual.tif")
wb <- rast("data/pp_atlas-climatico/geonode_p_ep_balanceh.tif")

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
burned_veg0 <- read.csv("data/data_burned_and_available_area_by_vegetation_dryforest2.csv")


# Tidy data ---------------------------------------------------------------

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
burned_veg <- read.csv("data/data_burned_and_available_area_by_vegetation_dryforest2.csv")
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

# ff <- list.files("data")
# ff <- ff[grep("data_randomized_fires_mindist_800m_", ff)]
# dsim0 <- do.call("rbind", lapply(ff, function(x) {
#   read.csv(paste("data", x, sep = "/"))
# }))
# nrow(dsim0) / length(unique(dsim0$replicate)) == nrow(d) # OK
# names(dsim0)
# unique(aggregate(burned_sim ~ replicate, dsim0, length)[, "burned_sim"])
# # OK
# dsim <- matrix(dsim0$burned_sim, ncol = max(dsim0$replicate))
# saveRDS(dsim, "exports/data_randomized_fires_mindist_800m_500_matrix.rds")
dsim <- readRDS("exports/data_randomized_fires_mindist_800m_500_matrix.rds")
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


# Average burned proportion in real and simulated datasets ----------------

hist(colMeans(dsim))
abline(v = mean(data$burned))


# Univariate fire models --------------------------------------------------

# All loops will run across veg types and then across predictors. This is so
# for convenience to fit a conditional model.

# objects to save
ob <- c("model", "pred", "pred_sim", "r2", "p", "plot")
O <- length(ob)

# pred is the data frame to plot everything
# pred_sim is the matrix of randomized predictions

# create nested lists to save all
firelist <- lapply(veg_model, function(v) {
  # list for variables
  vlist <- lapply(predictors, function(p) {
    # list of objects
    olist <- vector("list", O)
    names(olist) <- ob
    return(olist)
  })
  names(vlist) <- predictors
  return(vlist)
})
names(firelist) <- veg_model

# huge loop
for(v in 1:V) {
  # v = 1
  veg <- veg_model[v]
  print(veg)
  for(p in 1:P) {
    # p = 5
    var <- predictors[p]
    print(var)

    # create local dataset
    if(v == 1) { # all veg types
      rows <- 1:nrow(data)
      dlocal <- data[, c("burned", var)]
      dsim_local <- dsim
    } else {
      rows <- data$vegetation_class == veg
      dlocal <- data[rows, c("burned", var)]
      dsim_local <- dsim[rows, ]
    }
    names(dlocal) <- c("burned", "x")

    # Fit model to observed data

    kgam <- 4; bs <- "cr"; knots <- NULL
    if(var == "aspect") {
      bs <- "cc"; knots <- list(x = c(0, 360))
    }

    mod <- bam(burned ~ s(x, k = kgam, bs = bs), #, pc = mean(dlocal$x)),
               family = "binomial", knots = knots,
               data = dlocal, method = "fREML", discrete = T, nthreads = 8)
    firelist[[v]][[p]]$model <- mod

    # compute prediction (with ci)
    ext <- hdi(dlocal$x, 0.98)
    xseq <- seq(ext["lower"], ext["upper"], length.out = 150)
    if(v == 1) {
      xseq <- seq(min(dlocal$x), max(dlocal$x), length.out = 150)
    }
    if(var == "aspect") {
      xseq <- seq(0, 360, length.out = 150)
    }

    xpred <- data.frame(x = xseq)
    ppp <- predict(mod, xpred, se.fit = T)
    pred_obs <- data.frame(p_mle = plogis(ppp$fit),
                           p_lower = plogis(ppp$fit - qnorm(0.975) * ppp$se.fit),
                           p_upper = plogis(ppp$fit + qnorm(0.975) * ppp$se.fit),
                           varying_val = xseq,
                           varying_var = var)

    # compute observed r2
    firelist[[v]][[p]]$r2 <- r2bern(fitted(mod))

    # Get predicted linear predictor for randomized fires, without intercept
    print("Fitting models to randomized fires")
    sims_lp <- sapply(1:ncol(dsim_local), function(j) {
      if(j %% 50 == 0) print(j)
      # j = 1
      data_sim <- dlocal
      data_sim$burned <- dsim_local[, j]
      mran <-  bam(burned ~ s(x, k = kgam, bs = bs),#, pc = mean(dlocal$x)),
                   family = "binomial", knots = knots,
                   data = data_sim, method = "fREML", discrete = T, nthreads = 8)
      linpred_ran <- predict(mran, xpred, type = "linear", se.fit = F) - coef(mran)[1]
      return(linpred_ran)
    })

    # add the observed intercept
    pred_sim <- plogis(sims_lp + coef(mod)[1])
    firelist[[v]][[p]]$pred_sim <- pred_sim

    # compute quantiles
    pred_sim_q <- apply(pred_sim, 1, quantile,
                        probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
                        method = 8) %>% t %>% as.data.frame
    colnames(pred_sim_q) <- c("lower95", "lower90", "lower80", "median",
                              "upper80", "upper90", "upper95")

    # merge with observed prediction
    pred <- cbind(pred_obs, pred_sim_q)
    firelist[[v]][[p]]$pred <- pred

    # p-values over the sequence
    pvalseq <- sapply(1:nrow(pred_sim), function(j) {
      distrib <- ecdf(pred_sim[j, ] %>% as.numeric())
      perc <- distrib(pred_obs$p_mle[j])
      p <- ifelse(perc >= 0.5, 2 * (1 - perc), 2 * perc)
      return(p)
    })
    firelist[[v]][[p]]$p <- pvalseq

    # Compute density and add to pred
    dens <- density(dlocal$x, from = min(dlocal$x), to = max(dlocal$x), n = 2^10)
    xseq_full <- seq(min(dlocal$x), max(dlocal$x), length.out = 150)
    if(var == "aspect") {
      xseq_full <- xseq
    }
    dens_approx <- approx(dens$x, dens$y, xout = xseq_full)

    # scale maximum density to maximum p
    maxp <- as.numeric(as.matrix(pred[, grep("upper", names(pred))])) %>% max
    dens_factor <- max(dens_approx$y) / maxp
    dens_data <- data.frame(
      density = dens_approx$y / dens_factor,
      varying_val = xseq_full,
      varying_var = var,
      vegetation_class = veg
    )

    if(var == "aspect") {
      xc <- circular(dlocal$x, type = "angles", units = "degrees", template = "geographic")
      dens <- density.circular(xc, bw = 10, n = 2 ^ 10)

      # circular hace cosas raras: su seq, dens$x empieza en 90 y va decreciendo
      # (gira a contrarreloj, donde el último número, cercano a 90, es -270 )
      # Entonces, convertimos los números negativos a lo que son + 360.
      densx_cor <- dens$x
      densx_cor[dens$x < 0] <- dens$x[dens$x < 0] + 360

      dens_approx <- approx(densx_cor, dens$y, xout = xseq_full)
      # no puede interpolar el extremo
      dens_approx$y[is.na(dens_approx$y)] <- mean(dens$y[which.min(densx_cor)],
                                                  dens$y[which.max(densx_cor)])

      dens_factor <- max(dens_approx$y) / maxp

      dens_data <- data.frame(
        density = dens_approx$y / dens_factor,
        varying_val = xseq_full,
        varying_var = var,
        vegetation_class = veg
      )
    }
    firelist[[v]][[p]]$dens <- dens_data

    # plot!

    plotcito <-
    ggplot() +
      # density
      geom_ribbon(
        data = dens_data, mapping = aes(x = varying_val, ymax = density, ymin = 0),
        fill = viridis(1, alpha = 0.2, option = "A", begin = 0.5),
        color = NA
      ) +
      # randomized prediction
      geom_ribbon(
        data = pred, mapping = aes(x = varying_val,
                                   ymax = upper95, ymin = lower95),
        color = NA, alpha = 0.2
      ) +
      geom_ribbon(
        data = pred, mapping = aes(x = varying_val,
                                   ymax = upper80, ymin = lower80),
        color = NA, alpha = 0.2
      ) +
      geom_line(
        data = pred, mapping = aes(x = varying_val, y = median),
        linetype = 2, linewidth = 0.4
      ) +
      # observed prediction
      geom_line(
        data = pred, mapping = aes(x = varying_val, y = p_mle),
        color = viridis(1, option = "A", begin = 0.2)
      ) +
      xlab(names_frame$var_name[p]) +
      ylab("Burn probability") +
      theme(panel.grid.minor = element_blank())

    firelist[[v]][[p]]$plot <- plotcito
  }
}

saveRDS(firelist, "exports/fire_models_univariate.rds")
firelist <- readRDS("exports/fire_models_univariate.rds")

# Univariate fire plots ---------------------------------------------------

# Two figures: topo and (pp, ndvi, dist). Use a verticla facet_wrap for each
# variable.

# settings for varying strip across veg_types (to differentiate the
# marginal row)
a <- element_text(face = "bold", colour = "white", angle = 270, size = 8)
b <- element_text(colour = "black", angle = 270, size = 9)
texts <- list(a, b, b, b, b, b)
c <- element_rect(fill = "gray10", color = "gray10")
d <- element_rect(fill = "white", color = "white")
backgrounds <- list(c, d, d, d, d, d)

# x breaks for each predictor
xbreaks <- list(
  "elevation" = c(500, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = c(0, 0.5, 1),
  "ndvi_mean" = c(0, 0.5, 1),
  "pp" = c(1000, 2000),
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

# factor to multiply r2 y position
r2yfac <- c(
  "elevation" = 0.85,
  "slope" = 0.85,
  "aspect" = 0.15,
  "TPI2k" = 0.15,
  "ndvi_mean" = 0.85,
  "pp" = 0.85,
  "dist_human" = 0.85,
  "dist_roads" = 0.85
)

r2xfac <- c(
  "elevation" = 0.85,
  "slope" = 0.85,
  "aspect" = 0.85,
  "TPI2k" = 0.5,
  "ndvi_mean" = 0.4,
  "pp" = 0.85,
  "dist_human" = 0.85,
  "dist_roads" = 0.85
)



plist_uni <- vector("list", P)

for(p in 1:P) {
  # p = 1
  var <- predictors[p]

  # extract data from all vegs
  pred <- do.call("rbind", lapply(1:V, function(v) {
    d <- firelist[[v]][[p]]$pred
    d$vegetation_class <- veg_model[v]
    return(d)
  }))
  dens_data <- do.call("rbind", lapply(1:V, function(v) {
    d <- firelist[[v]][[p]]$dens
    d$vegetation_class <- veg_model[v]
    return(d)
  }))
  r2 <- sapply(1:V, function(v) {
    firelist[[v]][[p]]$r2
  })

  pred$vegetation_class <- factor(pred$vegetation_class,
                                  levels = veg_model,
                                  labels = veg_model_num)
  dens_data$vegetation_class <- factor(dens_data$vegetation_class,
                                       levels = veg_model,
                                       labels = veg_model_num)

  # set position for r2. Different computation because y starts alway at zero,
  # but it varies in max between vet types
  dxrange <- range(dens_data$varying_val)
  xr2 <- dxrange[1] + r2xfac[p] * (diff(dxrange))

  max_y <- sapply(veg_model_num, function(v) {
    yy <- c(pred$p_mle[pred$vegetation_class == v],
            pred$upper95[pred$vegetation_class == v],
            dens_data$density[dens_data$vegetation_class == v])
    return(max(yy))
  })
  r2data <- data.frame(
    r2 = paste(format(round(r2 * 100, 2), nsmall = 2), "%"),
    vegetation_class = veg_model_num,
    x = xr2,
    y = max_y * r2yfac[p]
  )

  plotcito <-
    ggplot() +
      # density
      geom_ribbon(
        data = dens_data, mapping = aes(x = varying_val, ymax = density, ymin = 0),
        fill = viridis(1, alpha = 0.2, option = "A", begin = 0.5),
        color = NA
      ) +
      # randomized prediction
      geom_ribbon(
        data = pred, mapping = aes(x = varying_val,
                                   ymax = upper95, ymin = lower95),
        color = NA, alpha = 0.2
      ) +
      geom_ribbon(
        data = pred, mapping = aes(x = varying_val,
                                   ymax = upper80, ymin = lower80),
        color = NA, alpha = 0.2
      ) +
      geom_line(
        data = pred, mapping = aes(x = varying_val, y = median),
        linetype = 2, linewidth = 0.4
      ) +
      # observed prediction
      geom_line(
        data = pred, mapping = aes(x = varying_val, y = p_mle),
        color = viridis(1, option = "A", begin = 0.2)
      ) +
      # r2
      geom_text(
        data = r2data,
        mapping = aes(x, y, label = r2), size = 2.5, color = "gray10"
      ) +

      facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right",
                scales = "free_y") +
      xlab(names_frame$var_name[p]) +
      ylab("Burn probability (%)") +
      scale_x_continuous(breaks = xbreaks[[p]],
                         labels = xlabels[[p]]) +
      scale_y_continuous(labels = scales::label_percent(suffix = "")) +
      theme(panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 9),
            axis.title.y = element_text(size = 10))
  # plotcito

  # remove y-axis title for most predictors except elevation and pp
  if(!(var %in% c("elevation", "pp"))) {
    plotcito <- plotcito + theme(axis.title.y = element_blank())
  }

  # remove strip for most predictors except pos topo and dist_roads
  if(var %in% c("TPI2k", "dist_roads")) {
    plotcito <- plotcito +
      theme(strip.background = element_rect(),
            strip.text = element_text()) +
      facet_wrap2(vars(vegetation_class), ncol = 1, strip.position = "right",
                  scales = "free_y",
                  strip = strip_themed(text_y = texts,
                                       background_y = backgrounds))
  }

  plist_uni[[p]] <- plotcito
}

fire_uni_topo <- ggarrange(plots = plist_uni[1:4], ncol = 4)
ggsave("figures/XX spatial patterns - fire by veg and topo - univariate.png",
       plot = fire_uni_topo,
       width = 16, height = 16.5, units = "cm")

fire_uni_pp <- ggarrange(plots = plist_uni[c(6, 5, 7, 8)], ncol = 4)
ggsave("figures/XX spatial patterns - fire by veg and pp-ndvi-dist - univariate.png",
       plot = fire_uni_pp,
       width = 16, height = 16.5, units = "cm")


# Univariate fire plots, fixed y ------------------------------------------

# Two figures: topo and (pp, ndvi, dist). Use a verticla facet_wrap for each
# variable.

# settings for varying strip across veg_types (to differentiate the
# marginal row)
a <- element_text(face = "bold", colour = "white", angle = 270, size = 8)
b <- element_text(colour = "black", angle = 270, size = 9)
texts <- list(a, b, b, b, b, b)
c <- element_rect(fill = "gray10", color = "gray10")
d <- element_rect(fill = "white", color = "white")
backgrounds <- list(c, d, d, d, d, d)

# x breaks for each predictor
xbreaks <- list(
  "elevation" = c(500, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = c(0, 0.5, 1),
  "ndvi_mean" = c(0, 0.5, 1),
  "pp" = c(1000, 2000),
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

# factor to multiply r2 y position
r2yfac <- c(
  "elevation" = 0.8,
  "slope" = 0.8,
  "aspect" = 0.8,
  "TPI2k" = 0.8,
  "ndvi_mean" = 0.8,
  "pp" = 0.8,
  "dist_human" = 0.8,
  "dist_roads" = 0.8
)

r2xfac <- c(
  "elevation" = 0.85,
  "slope" = 0.85,
  "aspect" = 0.85,
  "TPI2k" = 0.85,
  "ndvi_mean" = 0.4,
  "pp" = 0.85,
  "dist_human" = 0.85,
  "dist_roads" = 0.85
)

# get maximum y value at mle prediction with observed data
(max_p <- sapply(1:V, function(v) {
  ff <- firelist[[v]]
  mm <- sapply(1:P, function(p) max(ff[[p]]$pred$p_mle))
  return(max(mm))
}) %>% max)
roof <- ceiling(max_p * 100) / 100
max_p <- 0.4; roof <- 0.4

plist_uni <- vector("list", P)

for(p in 1:P) {
  # p = 1
  var <- predictors[p]

  # extract data from all vegs
  pred <- do.call("rbind", lapply(1:V, function(v) {
    d <- firelist[[v]][[p]]$pred
    d$vegetation_class <- veg_model[v]
    return(d)
  }))

  dens_data <- do.call("rbind", lapply(1:V, function(v) {
    d <- firelist[[v]][[p]]$dens

    # rescale the density to the roof
    dens_factor <- max(d$density) / (roof * 0.98)
    d$density <- d$density / dens_factor

    d$vegetation_class <- veg_model[v]
    return(d)
  }))

  r2 <- sapply(1:V, function(v) {
    firelist[[v]][[p]]$r2
  })

  pred$vegetation_class <- factor(pred$vegetation_class,
                                  levels = veg_model,
                                  labels = veg_model_num)
  dens_data$vegetation_class <- factor(dens_data$vegetation_class,
                                       levels = veg_model,
                                       labels = veg_model_num)

  # set position for r2. Different computation because y starts alway at zero,
  # but it varies in max between vet types
  dxrange <- range(dens_data$varying_val)
  xr2 <- dxrange[1] + r2xfac[p] * (diff(dxrange))

  r2data <- data.frame(
    r2 = paste(format(round(r2 * 100, 2), nsmall = 2), "%"),
    vegetation_class = veg_model_num,
    x = xr2,
    y = roof * r2yfac[p]
  )

  # truncate the ribbon at the roof
  vars_trunc <- c("upper95", "lower95", "upper80", "lower80")
  for(vt in vars_trunc) {
    pred[pred[, vt] > roof, vt] <- roof
  }

  plotcito <-
    ggplot() +
    # density
    geom_ribbon(
      data = dens_data, mapping = aes(x = varying_val, ymax = density, ymin = 0),
      fill = viridis(1, alpha = 0.2, option = "A", begin = 0.5),
      color = NA
    ) +
    # randomized prediction
    geom_ribbon(
      data = pred, mapping = aes(x = varying_val,
                                 ymax = upper95, ymin = lower95),
      color = NA, alpha = 0.2
    ) +
    geom_ribbon(
      data = pred, mapping = aes(x = varying_val,
                                 ymax = upper80, ymin = lower80),
      color = NA, alpha = 0.2
    ) +
    geom_line(
      data = pred, mapping = aes(x = varying_val, y = median),
      linetype = 2, linewidth = 0.4
    ) +
    # observed prediction
    geom_line(
      data = pred, mapping = aes(x = varying_val, y = p_mle),
      color = viridis(1, option = "A", begin = 0.2)
    ) +
    # r2
    geom_text(
      data = r2data,
      mapping = aes(x, y, label = r2), size = 2.5, color = "gray10"
    ) +

    facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right",
               scales = "free_y") +
    xlab(names_frame$var_name[p]) +
    ylab("Burn probability (%)") +
    scale_x_continuous(breaks = xbreaks[[p]],
                       labels = xlabels[[p]]) +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, roof),
                       breaks = seq(0, roof, by = 0.1),
                       expand = c(0.0025, 0.0025)) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 10))
  # plotcito

  # remove y-axis title for most predictors except elevation and pp
  if(!(var %in% c("elevation", "pp"))) {
    plotcito <- plotcito + theme(axis.title.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.text.y = element_blank())
  }

  # remove strip for most predictors except pos topo and dist_roads
  if(var %in% c("TPI2k", "dist_roads")) {
    plotcito <- plotcito +
      theme(strip.background = element_rect(),
            strip.text = element_text()) +
      facet_wrap2(vars(vegetation_class), ncol = 1, strip.position = "right",
                  scales = "free_y",
                  strip = strip_themed(text_y = texts,
                                       background_y = backgrounds))
  }

  plist_uni[[p]] <- plotcito
}

fire_uni_topo <- ggarrange(plots = plist_uni[1:4], ncol = 4)
ggsave("figures/XX spatial patterns - fire by veg and topo - univariate.png",
       plot = fire_uni_topo,
       width = 16, height = 16.5, units = "cm")

fire_uni_pp <- ggarrange(plots = plist_uni[c(6, 5, 7, 8)], ncol = 4)
ggsave("figures/XX spatial patterns - fire by veg and pp-ndvi-dist - univariate.png",
       plot = fire_uni_pp,
       width = 16, height = 16.5, units = "cm")



# NDVI models --------------------------------------------------------------

# used to predict NDVI at different values of topographic variables
# and precipitation. In this way, the multiple model including NDVI
# can show partial effects of covariates but retaining the NDVI effect.
predictors_env <- predictors[predictors != "ndvi_mean"]
rows_comp <- data$vegetation_class %in% veg_labels_sub
data_comp <- data[rows_comp, ]
nv <- veg_labels_sub %>% length

ndvi_models <- list(
  "all" = bam(
    ndvi_mean ~
      s(elevation, k = 3, bs = "cr") +
      s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") +
      s(TPI2k, k = 4, bs = "cr") +
      s(pp, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = betar(), discrete = T, nthreads = 8,
    data = data
  ),
  "veg" = bam(
    ndvi_mean ~
      vegetation_class +
      s(elevation, by = vegetation_class, k = 3, bs = "cr", id = 1) +
      s(slope, by = vegetation_class, k = 4, bs = "cr", id = 2) +
      s(aspect, by = vegetation_class, k = 4, bs = "cc", id = 3) +
      s(TPI2k, by = vegetation_class, k = 4, bs = "cr", id = 4) +
      s(pp, by = vegetation_class, k = 4, bs = "cr", id = 5),
    knots = list(aspect = c(0, 360)),
    family = betar(), discrete = T, nthreads = 8,
    data = data_comp
  )
)

predictors_ndvi <- predictors_env[grep("dist", predictors_env, invert = T)]

# saveRDS(ndvi_models, "exports/models_ndvi.rds")

# NDVI models predictions --------------------------------------------------

# settings for varying strip across veg_types (to differentiate the
# marginal row)
a <- element_text(face = "bold", colour = "white", angle = 270, size = 8)
b <- element_text(colour = "black", angle = 270, size = 9)
texts <- list(a, b, b, b, b, b)
c <- element_rect(fill = "gray10", color = "gray10")
d <- element_rect(fill = "white", color = "white")
backgrounds <- list(c, d, d, d, d, d)

# x breaks for each predictor
xbreaks <- list(
  "elevation" = c(400, 1000, 1600),#c(500, 1000, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = seq(0, 1, by = 0.25),#c(0, 0.5, 1),
  "pp" = c(800, 1400, 2000)# seq(500, 2000, by = 500)#c(1000, 2000)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

pred_ndvi_names <- names_frame$var_name[names_frame$variable %in% predictors_ndvi]
pred_ndvi_names[5] <- "(E) Precipitation\n(mm / year)"

P <- length(predictors_ndvi)
pndvi_list <- vector("list", P)

for(p in 1:P) {
  # p = 1
  var <- predictors_ndvi[p]

  # compute predictions by veg type
  ndvi_preds <- do.call("rbind", lapply(veg_model, function(veg) {
    # veg = "Shrubland"
    if(veg == "All vegetation\ntypes") {
      dlocal <- data[, predictors_ndvi]
    } else {
      dlocal <- data[data$vegetation_class == veg, predictors_ndvi]
    }

    if(var != "aspect") {
      ndvi_pred <- make_newdata(var, dlocal, HDI = TRUE, prob = 0.98)
    } else {
      ndvi_pred <- make_newdata(var, dlocal, HDI = FALSE)
    }
    ndvi_pred$vegetation_class <- veg

    if(veg == "All vegetation\ntypes") {
      m_use <- ndvi_models[["all"]]
    } else {
      m_use <- ndvi_models[["veg"]]
    }
    ndvi_logit <- predict(m_use, ndvi_pred, se.fit = T)

    ndvi_pred$pfit <- plogis(ndvi_logit$fit) %>% as.numeric
    ndvi_pred$plower <- plogis(ndvi_logit$fit - qnorm(0.975) * ndvi_logit$se.fit) %>% as.numeric
    ndvi_pred$pupper <- plogis(ndvi_logit$fit + qnorm(0.975) * ndvi_logit$se.fit) %>% as.numeric

    return(ndvi_pred)
  }))

  # rename vegs
  ndvi_preds$vegetation_class <- factor(ndvi_preds$vegetation_class,
                                        levels = veg_model,
                                        labels = veg_model_num)

  plotcito <-
    ggplot(data = ndvi_preds,
           mapping = aes(x = varying_val, y = pfit,
                         ymin = plower, ymax = pupper)) +
    geom_ribbon(alpha = 0.5, color = NA) +
    geom_line() +
    facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right") +
    xlab(pred_ndvi_names[p]) +
    ylab("NDVI") +
    scale_x_continuous(breaks = xbreaks[[p]],
                       labels = xlabels[[p]]) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 10))
  plotcito

  # remove y-axis title for most predictors except elevation
  if(var != "elevation") {
    plotcito <- plotcito + theme(axis.title.y = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank())
  }

  # remove strip for most predictors except pp
  if(var == "pp") {
    plotcito <- plotcito +
      theme(strip.background = element_rect(),
            strip.text = element_text()) +
      facet_wrap2(vars(vegetation_class), ncol = 1, strip.position = "right",
                  strip = strip_themed(text_y = texts,
                                       background_y = backgrounds))
  }

  pndvi_list[[p]] <- plotcito
}

ndvi_plot_all <- egg::ggarrange(plots = pndvi_list, ncol = length(pndvi_list))

ggsave("figures/S06) NDVI model predictions.png", plot = ndvi_plot_all,
       width = 16, height = 16.5, units = "cm")


# Topographic models ------------------------------------------------------

# to take into account their correlation when making predictions.
topo <- predictors[c(1, 2, 4)]
PT <- length(topo)

# list of formulas
combs_topo <- expand.grid(response = topo, predictor = topo)
combs_topo <- combs_topo[combs_topo$response != combs_topo$predictor, ]
combs_topo$formula_all <- paste(
  combs_topo$response, " ~ ", "s(", combs_topo$predictor, ", k = 4, bs = \"cr\")",
  sep = ""
)
combs_topo$formula_veg <- paste(
  combs_topo$response, " ~ ", "s(", combs_topo$predictor, ", by = vegetation_class, k = 4, bs = \"cr\", id = 1) + vegetation_class",
  sep = ""
)

family_resp <- list("elevation" = gaussian(),
                    "slope" = Gamma(link = "log"),
                    "TPI2k" = mgcv::betar())

# list of models, ordered as
# predictor > response > vegetation (all-veg)
mtopo <- vector("list", PT)
names(mtopo) <- topo

# duplicate dataset to modify slope to use Gamma distribution:
data_topo <- data
data_topo$slope[data$slope == 0] <- 0.0001

for(pred in topo) {
  # pred <- "slope"
  mtopo[[pred]] <- vector("list", 2) # predict two variables
  responses <- combs_topo$response[combs_topo$predictor == pred]
  names(mtopo[[pred]]) <- responses

  # loop by response, to make a list of models by vegetation type
  for(resp in responses) {
    # resp = "elevation"
    mtopo[[pred]][[resp]] <- vector("list", 2) # vegetation types + all
    names(mtopo[[pred]][[resp]]) <- c("all", "veg")

    # get formula position
    formula_row <- combs_topo$response == resp & combs_topo$predictor == pred

    # fit models
    mtopo[[pred]][[resp]][["all"]] <- gam(
      as.formula(combs_topo$formula_all[formula_row]), data = data_topo,
      family = family_resp[[resp]], method = "REML"
    )

    mtopo[[pred]][[resp]][["veg"]] <- gam(
      as.formula(combs_topo$formula_veg[formula_row]),
      data = data_topo[data_topo$vegetation_class %in% veg_labels_sub, ],
      family = family_resp[[resp]], method = "REML"
    )
  }
}

# mtopo[["elevation"]][["slope"]][["all"]] %>% plot
# mtopo[["elevation"]][["TPI2k"]][["veg"]] %>% plot

# saveRDS(mtopo, "exports/models_topography.rds")

# Distance models ----------------------------------------------------------

# to take into account their correlation when making predictions.
distv <- predictors[grep("dist", predictors)]
PD <- length(distv)

# list of formulas
combs_dist <- expand.grid(response = distv, predictor = distv)
combs_dist <- combs_dist[combs_dist$response != combs_dist$predictor, ]
combs_dist$formula_all <- paste(
  combs_dist$response, " ~ ", combs_dist$predictor,
  sep = ""
)
combs_dist$formula_veg <- paste(
  combs_dist$response, " ~ ", combs_dist$predictor, " * vegetation_class",
  sep = ""
)

# list of models, ordered as
# predictor > response > vegetation (all-veg)
mdist <- vector("list", PD)
names(mdist) <- distv

# duplicate dataset to modify slope to use Gamma distribution:
data_dist <- data
data_dist$dist_human[data$dist_human == 0] <- 0.0001
data_dist$dist_roads[data$dist_roads == 0] <- 0.0001

for(pred in distv) {
  # pred <- "slope"
  mdist[[pred]] <- vector("list", 2) # predict two variables
  responses <- combs_dist$response[combs_dist$predictor == pred]
  names(mdist[[pred]]) <- responses

  # loop by response, to make a list of models by vegetation type
  for(resp in responses) {
    # resp = "elevation"
    mdist[[pred]][[resp]] <- vector("list", 2) # vegetation types + all
    names(mdist[[pred]][[resp]]) <- c("all", "veg")

    # get formula position
    formula_row <- combs_dist$response == resp & combs_dist$predictor == pred

    # fit models
    mdist[[pred]][[resp]][["all"]] <- glm(
      as.formula(combs_dist$formula_all[formula_row]), data = data_dist,
      family = Gamma(link = "log")
    )

    mdist[[pred]][[resp]][["veg"]] <- glm(
      as.formula(combs_dist$formula_veg[formula_row]),
      data = data_dist[data_dist$vegetation_class %in% veg_labels_sub, ],
      family = Gamma(link = "log")
    )
  }
}
# mdist[["dist_human"]][["dist_roads"]][["all"]] %>% summary

# saveRDS(mdist, "exports/models_distances.rds")

# Multivariate fire models -----------------------------------------------

# Fit models to observed data
# use the same smoothing parameter across vegetation types, but varying across
# predictors. In elevation and pp, reduce k from 4 to 3 to avoid fitting noisy
# curves.
fire_models_obs <- list(
  "all" = bam(
    burned ~
      s(elevation, k = 3, bs = "cr") +
      s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") +
      s(TPI2k, k = 4, bs = "cr") +
      s(pp, k = 3, bs = "cr") +
      s(ndvi_mean, k = 4, bs = "cr") +
      s(dist_human, k = 4, bs = "cr") +
      s(dist_roads, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data, discrete = T, nthreads = 8
  ),
  "veg" = bam(
    burned ~
      vegetation_class +
      s(elevation, by = vegetation_class, k = 3, bs = "cr", id = 1) +
      s(slope, by = vegetation_class, k = 4, bs = "cr", id = 2) +
      s(aspect, by = vegetation_class, k = 4, bs = "cc", id = 3) +
      s(TPI2k, by = vegetation_class, k = 4, bs = "cr", id = 4) +
      s(pp, by = vegetation_class, k = 3, bs = "cr", id = 5) +
      s(ndvi_mean, by = vegetation_class, k = 4, bs = "cr", id = 6) +
      s(dist_human, by = vegetation_class, k = 4, bs = "cr", id = 7) +
      s(dist_roads, by = vegetation_class, k = 4, bs = "cr", id = 8),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data_veg, discrete = T, nthreads = 8
  )
)

# fit the same models to randomized fires
data_ran <- data # to overwrite burned column
nsim <- ncol(dsim)

fire_models_sim <- list(
  "all" = vector("list", nsim),
  "veg" = vector("list", nsim)
)

for(i in 1:nsim) {
  print(i)
  data_ran$burned <- dsim[, i]

  fire_models_sim[["all"]][[i]] <- bam(
    burned ~
      s(elevation, k = 3, bs = "cr") +
      s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") +
      s(TPI2k, k = 4, bs = "cr") +
      s(pp, k = 3, bs = "cr") +
      s(ndvi_mean, k = 4, bs = "cr") +
      s(dist_human, k = 4, bs = "cr") +
      s(dist_roads, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data_ran, discrete = T, nthreads = 8
  )

  fire_models_sim[["veg"]][[i]] <- bam(
    burned ~
      vegetation_class +
      s(elevation, by = vegetation_class, k = 3, bs = "cr", id = 1) +
      s(slope, by = vegetation_class, k = 4, bs = "cr", id = 2) +
      s(aspect, by = vegetation_class, k = 4, bs = "cc", id = 3) +
      s(TPI2k, by = vegetation_class, k = 4, bs = "cr", id = 4) +
      s(pp, by = vegetation_class, k = 3, bs = "cr", id = 5) +
      s(ndvi_mean, by = vegetation_class, k = 4, bs = "cr", id = 6) +
      s(dist_human, by = vegetation_class, k = 4, bs = "cr", id = 7) +
      s(dist_roads, by = vegetation_class, k = 4, bs = "cr", id = 8),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data_ran[data_ran$vegetation_class %in% veg_labels_sub, ],
    discrete = T, nthreads = 8
  )
}


# saveRDS(fire_models_sim, "exports/fire_models_multivariate.rds")
fire_models_sim <- readRDS("exports/fire_models_multivariate.rds")
# object.size(fire_models_sim) / 1e6 # 2gb

# Multivariate fire predictions --------------------------------------------

# covary topographic and ndvi variables to make predictions?
cov_topo <- F
cov_ndvi <- F
cov_dist <- F

nsim <- ncol(dsim)
# make newdata for predictions

newdata <- list(
  "all" = do.call("rbind", lapply(predictors, function(p) {
    if(p != "aspect") {
      rr <- make_newdata(p, data[, predictors], HDI = TRUE, prob = 0.98)
    } else {
      rr <- make_newdata(p, data[, predictors], HDI = FALSE)
    }
    return(rr)
  })),

  "veg" = do.call("rbind", lapply(veg_labels_sub, function(v) {
    dataloc <- data[data$vegetation_class == v, predictors]
    rr2 <- do.call("rbind", lapply(predictors, function(p) {
      if(p != "aspect") {
        rr <- make_newdata(p, dataloc, HDI = TRUE, prob = 0.98)
      } else {
        rr <- make_newdata(p, dataloc, HDI = FALSE)
      }
      return(rr)
    }))
    rr2$vegetation_class <- v
    return(rr2)
  }))
)
newdata$all$vegetation_class <- veg_model[1]

# alter newdata by making predictors co-vary
if(cov_topo) {
  # in the case of topographic variables, consider their correlation
  # (elevation, slope and aspect.)
  for(predictor in topo) {
    # all vegetation types
    # predictor <- "elevation"
    rows_varying <- newdata$all$varying_var == predictor
    responses <- combs_topo$response[combs_topo$predictor == predictor]
    for(response in responses) {
      # response = "slope"
      prediction <- predict(mtopo[[predictor]][[response]][["all"]],
                            newdata = newdata$all[rows_varying, predictor, drop = F],
                            type = "response")
      newdata$all[rows_varying, response] <- prediction
    }

    # separate vegetation types
    rows_varying <- newdata$veg$varying_var == predictor
    responses <- combs_topo$response[combs_topo$predictor == predictor]
    for(response in responses) {
      prediction <- predict(mtopo[[predictor]][[response]][["veg"]],
                            newdata = newdata$veg[rows_varying, ],
                            type = "response")
      newdata$veg[rows_varying, response] <- prediction
    }
  }
}
if(cov_dist) {
  for(predictor in distv) {
    # all vegetation types
    rows_varying <- newdata$all$varying_var == predictor
    responses <- combs_dist$response[combs_dist$predictor == predictor]
    for(response in responses) {
      prediction <- predict(mdist[[predictor]][[response]][["all"]],
                            newdata = newdata$all[rows_varying, predictor, drop = F],
                            type = "response")
      newdata$all[rows_varying, response] <- prediction
    }

    # separate vegetation types
    rows_varying <- newdata$veg$varying_var == predictor
    responses <- combs_dist$response[combs_dist$predictor == predictor]
    for(response in responses) {
      prediction <- predict(mdist[[predictor]][[response]][["veg"]],
                            newdata = newdata$veg[rows_varying, ],
                            type = "response")
      newdata$veg[rows_varying, response] <- prediction
    }
  }
}
if(cov_ndvi) {
  # all vegetation types
  rows_varying <- newdata$all$varying_var %in% predictors_ndvi
  prediction <- predict(ndvi_models[["all"]],
                        newdata = newdata$all[rows_varying, , drop = F],
                        type = "response")
  newdata$all[rows_varying, "ndvi_mean"] <- prediction

  # all vegetation types
  rows_varying <- newdata$veg$varying_var %in% predictors_ndvi
  prediction <- predict(ndvi_models[["veg"]],
                        newdata = newdata$veg[rows_varying, , drop = F],
                        type = "response")
  newdata$veg[rows_varying, "ndvi_mean"] <- prediction
}


# compute predictions
pred_obs_logit <- list(
  "all" = predict(fire_models_obs$all, newdata$all, se.fit = T),
  "veg" = predict(fire_models_obs$veg, newdata$veg, se.fit = T)
)

# compute mean by group at logit scale to center randomized predictions
pred_obs_logit_means <- list(
  "all" = aggregate(pred_obs_logit$all$fit ~ newdata$all$varying_var,
                    FUN = mean),
  "veg" = aggregate(pred_obs_logit$veg$fit ~
                      newdata$veg$varying_var +
                      newdata$veg$vegetation_class,
                    FUN = mean)
)
colnames(pred_obs_logit_means$all) <- c("predictor", "logit_mean")
colnames(pred_obs_logit_means$veg) <- c("predictor", "vegetation_class", "logit_mean")

pred_obs_p <- list(
  "all" = data.frame(
    p_mle = plogis(pred_obs_logit$all$fit),
    p_lower = plogis(pred_obs_logit$all$fit - qnorm(0.975) * pred_obs_logit$all$se.fit),
    p_upper = plogis(pred_obs_logit$all$fit + qnorm(0.975) * pred_obs_logit$all$se.fit)
  ),
  "veg" = data.frame(
    p_mle = plogis(pred_obs_logit$veg$fit),
    p_lower = plogis(pred_obs_logit$veg$fit - qnorm(0.975) * pred_obs_logit$veg$se.fit),
    p_upper = plogis(pred_obs_logit$veg$fit + qnorm(0.975) * pred_obs_logit$veg$se.fit)
  )
)

# from randomized fires, will be filled later
pred_sim <- list(
  "logit" = list(
    "all" = matrix(NA, nrow(newdata$all), nsim),
    "veg" = matrix(NA, nrow(newdata$veg), nsim)
  ),
  "prob" = list(
    "all" = matrix(NA, nrow(newdata$all), nsim),
    "veg" = matrix(NA, nrow(newdata$veg), nsim)
  ),
  "prob_summ" = list(
    "all" = NULL,
    "veg" = NULL
  )
)

# compute randomized predictions and fire effects with p-values. First,
# for all veg types, and then, separately
for(i in 1:nsim) {
  pred_sim$logit$all[, i] <- predict(fire_models_sim$all[[i]],
                                     newdata$all, type = "link")
  pred_sim$logit$veg[, i] <- predict(fire_models_sim$veg[[i]],
                                     newdata$veg, type = "link")
}

# center predictions at logit scale and transform to probability
for(p in predictors) {
  # p = "elevation"
  # all
  rr <- newdata$all$varying_var == p
  meanfit <- mean(pred_obs_logit$all$fit[rr])
  meanran <- colMeans(pred_sim$logit$all[rr, ])
  mdiff <- meanfit - meanran
  mdiff_mat <- outer(rep(1, sum(rr)), mdiff)
  pred_sim$prob$all[rr, ] <- plogis(pred_sim$logit$all[rr, ] + mdiff_mat)

  # veg
  for(v in veg_labels_sub) {
    # v = "Shrubland"
    rr <- newdata$veg$varying_var == p & newdata$veg$vegetation_class == v
    meanfit <- mean(pred_obs_logit$veg$fit[rr])
    meanran <- colMeans(pred_sim$logit$veg[rr, ])
    mdiff <- meanfit - meanran
    mdiff_mat <- outer(rep(1, sum(rr)), mdiff)
    pred_sim$prob$veg[rr, ] <- plogis(pred_sim$logit$veg[rr, ] + mdiff_mat)
  }
}

# compute quantiles
pred_sim$prob_summ$all <- apply(
  pred_sim$prob$all, 1, quantile,
  probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
  method = 8) %>% t %>% as.data.frame
pred_sim$prob_summ$veg <- apply(
  pred_sim$prob$veg, 1, quantile,
  probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
  method = 8) %>% t %>% as.data.frame

colnames(pred_sim$prob_summ$all) <-
  colnames(pred_sim$prob_summ$veg) <-
    c("lower95", "lower90", "lower80", "median", "upper80", "upper90", "upper95")


# put all predictions in a single df
pred <- cbind(
  rbind(newdata$all, newdata$veg),
  rbind(pred_obs_p$all, pred_obs_p$veg),
  rbind(pred_sim$prob_summ$all, pred_sim$prob_summ$veg)
)

# compute effects
pred_sim_both <- rbind(pred_sim$prob$all, pred_sim$prob$veg)
effects_data <- expand.grid(vegetation_class = veg_model,
                            predictor = predictors,
                            delta = NA, pval = NA, text = NA)

for(p in predictors) {
  for(v in veg_model) {
    # p = "elevation"; v = "Shrubland"

    # get rows to filter predictions
    rr <- pred$varying_var == p & pred$vegetation_class == v

    # get density to compute weights
    dens_data <- firelist[[v]][[p]]$dens
    weights <- approx(dens_data$varying_val,
                      dens_data$density,
                      xout = pred$varying_val[rr], rule = 2)$y %>% normalize

    # pointwise mean in randomized fires
    pred_both <- cbind(pred$p_mle[rr], pred_sim_both[rr, ])
    p_means <- colSums(pred_both * weights)

    # for each curve, get the avg absolute difference with the mean
    delta_mat <- abs(pred_both - p_means) * 100
    delta_means <- colSums(delta_mat * weights)

    # metric and p value
    delta <- delta_means[1]
    pval <- 1 - ecdf(delta_means[-1])(delta)

    # make texts
    text_delta <- paste(format(round(delta, 2), nsmall = 2), "%")
    text_pval <- paste("(", format(round(pval, 3), nsmall = 3), ")", sep = "")
    if(text_pval == "(0.000)") text_pval <- "(< 0.002)"
    text <- paste(text_delta, "\n", text_pval, sep = "")

    # fill table
    rout <- effects_data$predictor == p & effects_data$vegetation_class == v
    effects_data$delta[rout] <- delta
    effects_data$pval[rout] <- pval
    effects_data$text[rout] <- text
  }
}

# save predictions and effects
result <- list("pred" = pred, "effects_data" = effects_data,
               covs = c("topo" = cov_topo, "ndvi" = cov_ndvi, "dist" = cov_dist))
# saveRDS(result, "exports/fire_models_multreg_result_cov-ndvi-topo-dist.rds")
# saveRDS(result, "exports/fire_models_multreg_result_cov0.rds")
result <- readRDS("exports/fire_models_multreg_result_cov-ndvi-topo-dist.rds")

pred <- result$pred
effects_data <- result$effects_data

# Multivariate fire plots, fixed y ----------------------------------------

# Two figures: topo and (pp, ndvi, dist). Use a vertical facet_wrap for each
# variable.

# settings for varying strip across veg_types (to differentiate the
# marginal row)
a <- element_text(face = "bold", colour = "white", angle = 270, size = 8)
b <- element_text(colour = "black", angle = 270, size = 9)
texts <- list(a, b, b, b, b, b)
c <- element_rect(fill = "gray10", color = "gray10")
d <- element_rect(fill = "white", color = "white")
backgrounds <- list(c, d, d, d, d, d)

# x breaks for each predictor
xbreaks <- list(
  "elevation" = c(400, 1000, 1600),#c(500, 1000, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = seq(0, 1, by = 0.25),#c(0, 0.5, 1),
  "ndvi_mean" = c(0, 0.5, 1),
  "pp" = c(800, 1400, 2000),# seq(500, 2000, by = 500)#c(1000, 2000)
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

# factor to multiply derivatives y position
delta_yfac <- c(
  "elevation" = 0.8,
  "slope" = 0.8,
  "aspect" = 0.8,
  "TPI2k" = 0.8,
  "ndvi_mean" = 0.8,
  "pp" = 0.8,
  "dist_human" = 0.8,
  "dist_roads" = 0.8
)

delta_xfac <- c(
  "elevation" = 0.85,
  "slope" = 0.85,
  "aspect" = 0.85,
  "TPI2k" = 0.85,
  "ndvi_mean" = 0.2,
  "pp" = 0.85,
  "dist_human" = 0.85,
  "dist_roads" = 0.85
)

P <- length(predictors)
plist_multi <- vector("list", P)

# # to compare effects later, save derivatives
# delta_list <- vector("list", P)
# ## instead of derivatives, compute the weighted difference between the curve
# ## fitted with real vs simulated data.
# ## The p value may be for this or for every point in the curve.

# get maximum y value at mle prediction with observed data
max_p <- max(pred$p_mle)
roof <- ceiling(max_p * 100) / 100
max_p <- roof <- 0.4

for(p in 1:P) {
  # p = 1
  var <- predictors[p]

  # extract data from focal predictor
  pred_local <- pred[pred$varying_var == var, ]

  # effects
  effects_data_local <- effects_data[effects_data$predictor == var, ]

  # extract density from the univariate firelist
  dens_data <- do.call("rbind", lapply(1:V, function(v) {
    # v = 1; p = 1
    d <- firelist[[v]][[p]]$dens
    d$vegetation_class <- veg_model[v]

    # rescale the density to the roof
    dens_factor <- max(d$density) / (roof * 0.98)
    d$density <- d$density / dens_factor

    return(d)
  }))

  pred_local$vegetation_class <- factor(pred_local$vegetation_class,
                                  levels = veg_model,
                                  labels = veg_model_num)
  dens_data$vegetation_class <- factor(dens_data$vegetation_class,
                                       levels = veg_model,
                                       labels = veg_model_num)

  # truncate the ribbon at the roof
  vars_trunc <- c("upper95", "lower95", "upper80", "lower80")
  for(vt in vars_trunc) {
    pred_local[pred_local[, vt] > roof, vt] <- roof
  }

  # set position for delta. Different computation because y starts always at zero,
  # but it varies in max between veg types
  dxrange <- range(dens_data$varying_val)
  x_delta <- dxrange[1] + delta_xfac[p] * (diff(dxrange))

  effects_data_local$x <- x_delta
  effects_data_local$y <- roof * delta_yfac[p]

  effects_data_local$vegetation_class <- factor(
    effects_data_local$vegetation_class, levels = veg_model,
    labels = veg_model_num
  )

  plotcito <-
    ggplot() +
    # density
    geom_ribbon(
      data = dens_data, mapping = aes(x = varying_val, ymax = density, ymin = 0),
      fill = viridis(1, alpha = 0.2, option = "A", begin = 0.5),
      #fill = viridis(1, alpha = 0.45, option = "D", begin = 0.7), # no me gustó
      #fill = viridis(1, alpha = 0.2, option = "D", begin = 0.45),
      color = NA
    ) +
    # randomized prediction
    geom_ribbon(
      data = pred_local, mapping = aes(x = varying_val,
                                 ymax = upper95, ymin = lower95),
      color = NA, alpha = 0.2
    ) +
    geom_ribbon(
      data = pred_local, mapping = aes(x = varying_val, y = p_mle,
                                 ymax = upper80, ymin = lower80),
      color = NA, alpha = 0.2, fill = "black"
    ) +
    geom_line(
      data = pred_local, mapping = aes(x = varying_val, y = median),
      linetype = 2, linewidth = 0.4
    ) +
    # observed prediction
    geom_line(
      data = pred_local, mapping = aes(x = varying_val, y = p_mle),
      color = viridis(1, option = "A", begin = 0.2)
    ) +
    # delta
    geom_text(
      data = effects_data_local,
      mapping = aes(x, y, label = text), size = 2.2, color = "gray10"
    ) +

    facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right",
               scales = "free_y") +
    xlab(names_frame$var_name[p]) +
    ylab("Burn probability (%)") +
    scale_x_continuous(breaks = xbreaks[[p]],
                       labels = xlabels[[p]]) +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, roof),
                       breaks = seq(0, roof, by = 0.1),
                       expand = c(0.0025, 0.0025)) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 10))
  plotcito

  # remove y-axis title for most predictors except elevation and pp
  if(!(var %in% c("elevation", "pp"))) {
    plotcito <- plotcito + theme(axis.title.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.text.y = element_blank())
  }

  # if(var == "pp") {
  #   plotcito <- plotcito + xlab(bquote("Precipitation\n(mm"~yr^{-1}~")"))
  # }

  # remove strip for most predictors except pos topo and dist_roads
  if(var %in% c("TPI2k", "dist_roads")) {
    plotcito <- plotcito +
      theme(strip.background = element_rect(),
            strip.text = element_text()) +
      facet_wrap2(vars(vegetation_class), ncol = 1, strip.position = "right",
                  scales = "free_y",
                  strip = strip_themed(text_y = texts,
                                       background_y = backgrounds))
  }

  plist_multi[[p]] <- plotcito
}

fire_multi_topo <- ggarrange(plots = plist_multi[1:4], ncol = 4)
ggsave("figures/04) spatial patterns - environment topo.png",
       plot = fire_multi_topo,
       width = 16, height = 16.5, units = "cm")

fire_multi_pp <- ggarrange(plots = plist_multi[c(6, 5, 7, 8)], ncol = 4)
ggsave("figures/05) spatial patterns - environment pp.png",
       plot = fire_multi_pp,
       width = 16, height = 16.5, units = "cm")

# Multivariate fire plots - vegetation types side by side -----------------

roof <- 0.4

# x breaks for each predictor
xbreaks <- list(
  "elevation" = c(400, 1000, 1600),#c(500, 1000, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = seq(0, 1, by = 0.25),#c(0, 0.5, 1),
  "ndvi_mean" = c(0, 0.5, 1),
  "pp" = c(800, 1400, 2000),# seq(500, 2000, by = 500)#c(1000, 2000)
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

plist_side <- vector("list", P)

for(p in 1:P) {
  # p = 1
  var <- predictors[p]

  # extract data from all vegs
  pred_local <- pred[pred$varying_var == var, ]

  pred_local$vegetation_class <- factor(pred_local$vegetation_class,
                                        levels = veg_labels_sub)

  plotcito <-
    ggplot() +

    # observed prediction global
    geom_line(
      data = pred_local[is.na(pred_local$vegetation_class), ],
      mapping = aes(x = varying_val, y = p_mle),
      color = "black", alpha = 0.2, linewidth = 1.2
    ) +

    # observed prediction vegs
    geom_line(
      data = pred_local[!is.na(pred_local$vegetation_class), ],
      mapping = aes(x = varying_val, y = p_mle,
                    color = vegetation_class, group = vegetation_class)
    ) +
    scale_color_viridis(option = "A", discrete = TRUE, end = 0.9) +    xlab(names_frame$var_name[p]) +

    xlab(names_frame$name2[p]) +
    ylab("Burn probability (%)") +
    scale_x_continuous(breaks = xbreaks[[p]],
                       labels = xlabels[[p]]) +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, roof),
                       breaks = seq(0, roof, by = 0.1),
                       expand = c(0.0025, 0.0025)) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 10),
          legend.position = "none")
  plotcito

  # remove y-axis title for most predictors except elevation and pp
  if(!(var %in% c("elevation", "pp"))) {
    plotcito <- plotcito + theme(axis.title.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.text.y = element_blank())
  }

  plist_side[[p]] <- plotcito
}

plist_side <- plist_side[c(1:4, 6, 5, 7, 8)]

pans <- ggarrange(plots = plist_side, nrow = 2)
pans <- ggpubr::as_ggplot(pans)

# get legend and cast as ggplot object
leg <- ggpubr::get_legend(plist_side[[1]] + theme(legend.position = "bottom",
                                                  legend.title = element_blank(),
                                                  legend.text = element_text(size = 8),
                                                  legend.box.margin = margin(t = 2, r = 2, b = 4, l = 2 ,unit = "mm")))
leg <- ggpubr::as_ggplot(leg)

vegs_side <- ggarrange(pans, leg, nrow = 2,
                               heights = c(1, 0.1))

ggsave("figures/07) XX curves side by side_cov-topo-ndvi-dist.png",
       plot = vegs_side,
       width = 16, height = 11, units = "cm")

# Spatial effects summary --------------------------------------------------

effects_plot <- effects_data

# remove distances, for their absurd effect in dry forest.
inn <- grep("dist", effects_plot$predictor, invert = T)
effects_plot <- effects_plot[inn, ]
q <- max(effects_plot$delta) / max(effects_plot$pval)
q <- 10
effects_plot$pval_scaled <- effects_plot$pval * q

# effects_plot$vegetation_class <- factor(effects_plot$vegetation_class,
#                                levels = veg_model_num,
#                                labels = veg_model2)
effects_plot$predictor2 <- factor(effects_plot$predictor,
                        levels = names_frame$variable,
                        labels = names_frame$name3)
effects_plot$predictor3 <- factor(effects_plot$predictor,
                                levels = names_frame$variable,
                                labels = names_frame$name4)

effects_plot$delta_label <- "difference"
effects_plot$p_label <- "p-value"

ggplot(effects_plot) +
  geom_bar(aes(x = vegetation_class, y = delta, color = delta_label,
               fill = delta_label),
           stat = "identity", alpha = 0.4, linewidth = 0.4) +
  scale_fill_viridis(option = "A", discrete = TRUE, name = NULL) +
  scale_color_viridis(option = "A", discrete = TRUE, name = NULL) +

  ggnewscale::new_scale_colour() +
  ggnewscale::new_scale_fill() +

  geom_point(aes(x = vegetation_class, y = pval_scaled, color = p_label,
                 fill = p_label)) +
  scale_fill_viridis(option = "A", discrete = TRUE, name = NULL, begin = 0.5) +
  scale_color_viridis(option = "A", discrete = TRUE, name = NULL, begin = 0.5) +

  facet_wrap(vars(predictor2), nrow = 2) +
  scale_y_continuous(
    name = "Burn probability difference (%)",
    limits = c(0, q * 1.05), breaks = seq(0, 10, by = 2.5), expand = c(0.01, 0.01),
    sec.axis = sec_axis(~ . / q,
                        name = "p-value",
                        breaks = round(seq(0, 1, by = 1/4), 2))
  ) +
  xlab("Vegetation type") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))

# ggsave("figures/06) XX spatial effects summary - predictors_cov0.png",
#        width = 16, height = 14, units = "cm")
ggsave("figures/06) XX spatial effects summary - predictors_cov-topo-ndvi-dist.png",
       width = 16, height = 14, units = "cm")


ggplot(effects_plot) +
  geom_bar(aes(x = predictor3, y = delta, color = delta_label,
               fill = delta_label),
           stat = "identity", alpha = 0.4, linewidth = 0.4) +
  scale_fill_viridis(option = "A", discrete = TRUE, name = NULL) +
  scale_color_viridis(option = "A", discrete = TRUE, name = NULL) +

  ggnewscale::new_scale_colour() +
  ggnewscale::new_scale_fill() +

  geom_point(aes(x = predictor3, y = pval_scaled, color = p_label,
                 fill = p_label)) +
  scale_fill_viridis(option = "A", discrete = TRUE, name = NULL, begin = 0.5) +
  scale_color_viridis(option = "A", discrete = TRUE, name = NULL, begin = 0.5) +

  facet_wrap(vars(vegetation_class), nrow = 2) +
  scale_y_continuous(
    name = "Burn probability difference (%)",
    limits = c(0, q*1.05), breaks = seq(0, 10, by = 2.5), expand = c(0.01, 0.01),
    sec.axis = sec_axis(~ . / q,
                        name = "p-value",
                        breaks = round(seq(0, 1, by = 1/4), 2))
  ) +
  xlab("Predictor") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))

# ggsave("figures/06) XX spatial effects summary - veg_cov0.png",
#        width = 16, height = 14, units = "cm")
ggsave("figures/06) XX spatial effects summary - veg_cov-topo-ndvi-dist.png",
       width = 16, height = 14, units = "cm")



# Multivariate plot, co-varying predictors or not --------------------------
# (To compare easily)

cov0 <- readRDS("exports/fire_models_multreg_result_cov0.rds")$pred
covall <- readRDS("exports/fire_models_multreg_result_cov-ndvi-topo-dist.rds")$pred
cov0$cov <- "fixed"
covall$cov <- "co-varying"

datacomp <- rbind(cov0, covall)
datacomp$cov <- factor(datacomp$cov)

# settings for varying strip across veg_types (to differentiate the
# marginal row)
a <- element_text(face = "bold", colour = "white", angle = 270, size = 8)
b <- element_text(colour = "black", angle = 270, size = 9)
texts <- list(a, b, b, b, b, b)
c <- element_rect(fill = "gray10", color = "gray10")
d <- element_rect(fill = "white", color = "white")
backgrounds <- list(c, d, d, d, d, d)

# x breaks for each predictor
xbreaks <- list(
  "elevation" = c(500, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = c(0, 0.5, 1),
  "ndvi_mean" = c(0, 0.5, 1),
  "pp" = c(1000, 2000),
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

max_p <- 0.4; roof <- 0.4

predictors_cov <- predictors[predictors != "ndvi_mean"]
Pcov <- length(predictors_cov)
plist_cov <- vector("list", Pcov)

names_frame_cov <- names_frame[names_frame$variable != "ndvi_mean", ]
names_frame_cov$var_name_cov <- paste("(", LETTERS[1:Pcov], ") ", names_frame_cov$name2, sep = "")

for(p in 1:Pcov) {
  # p = 1
  var <- predictors_cov[p]

  # extract data from all vegs
  pred_local <- datacomp[datacomp$varying_var == var, ]

  pred_local$vegetation_class <- factor(pred_local$vegetation_class,
                                       levels = veg_model,
                                       labels = veg_model_num)

  plotcito <-
    ggplot() +

    # observed prediction
    geom_line(
      data = pred_local, mapping = aes(x = varying_val, y = p_mle,
                                      color = cov)
    ) +
    scale_color_viridis(discrete = T, begin = 0.2, end = 0.5, option = "A") +
    facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right",
               scales = "free_y") +
    xlab(names_frame_cov$var_name_cov[p]) +
    ylab("Burn probability (%)") +
    scale_x_continuous(breaks = xbreaks[predictors_cov][[p]],
                       labels = xlabels[predictors_cov][[p]]) +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, roof),
                       breaks = seq(0, roof, by = 0.1),
                       expand = c(0.0025, 0.0025)) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 10),
          legend.position = "none")
  plotcito

  # remove y-axis title for most predictors except elevation and pp
  if(var != "elevation") {
    plotcito <- plotcito + theme(axis.title.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.text.y = element_blank())
  }

  # remove strip for most predictors except pos topo and dist_roads
  if(p == Pcov) {
    plotcito <- plotcito +
      theme(strip.background = element_rect(),
            strip.text = element_text(),
            legend.position = "right",
            legend.title = element_blank()) +
      facet_wrap2(vars(vegetation_class), ncol = 1, strip.position = "right",
                  scales = "free_y",
                  strip = strip_themed(text_y = texts,
                                       background_y = backgrounds))
  }

  plist_cov[[p]] <- plotcito
}


fire_multi_ndvi_varfix <- ggarrange(plots = plist_cov, ncol = Pcov)
ggsave("figures/S07) burn probability varying predictors or not.png",
       plot = fire_multi_ndvi_varfix,
       width = 24, height = 16.5, units = "cm")

# Burn prob ~ vegetation, marginal ---------------------------------------

# plot with vegetation distribution and burn probability, both randomized
# and estimated. (no need to fit model, just average!)

sumlen <- function(x) return(c("sum" = sum(x), "len" = length(x)))
dveg <- do.call("data.frame", aggregate(burned ~ vegetation_class, data, sumlen))
pobs_fit <- glm(cbind(burned.sum, burned.len - burned.sum) ~ vegetation_class - 1,
                family = "binomial", data = dveg)

props_data <- data.frame(veg = dveg$vegetation_class,
                         pobs = dveg$burned.sum / dveg$burned.len)
# r2
r2bern(props_data$pobs) * 100

# compute the same for simulated fires
sim_sums <- as.matrix(aggregate(dsim ~ vegetation_class, data, sum)[, -1])
sim_lens <- as.matrix(aggregate(dsim ~ vegetation_class, data, length)[, -1])
sim_prob <- sim_sums / sim_lens

# compute quantiles
sim_prob_q <- apply(sim_prob, 1, quantile,
                    probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
                    method = 8) %>% t %>% as.data.frame
colnames(sim_prob_q) <- c("lower95", "lower90", "lower80", "median",
                          "upper80", "upper90", "upper95")
sim_prob_q$max <- apply(sim_prob, 1, max)

# p-values
# min is 2 * 1 / ncol(dsim) = 0.004
props_data$p <- sapply(1:nrow(sim_prob), function(j) {
  distrib <- ecdf(sim_prob[j, ] %>% as.numeric())
  perc <- distrib(props_data$pobs[j])
  p <- ifelse(perc >= 0.5, 2 * (1 - perc), 2 * perc)
  return(p)
}) %>% round(digits = 3) %>% format(nsmall = 3)
props_data$p[props_data$p == "0.000"] <- "< 0.004"
props_data <- cbind(props_data, sim_prob_q)

# show distributions
props_dens <- as.data.frame(cbind(vegetation_class = dveg$vegetation_class,
                                  as.data.frame(sim_prob)))
# names(props_dens)
props_dens_long <- pivot_longer(props_dens, all_of(2:ncol(props_dens)),
                                values_to = "p", names_to = "id")

# Compute the categorical effect and p-value
vegeff_marg <- categorical_effect2(y = props_data$pobs * 100,
                                   w = data_vd$veg_dist)
vegeff_marg_sim <- sapply(1:ncol(dsim), function(j) {
  categorical_effect2(y = sim_prob[, j] * 100,
                      w = data_vd$veg_dist)
})
vegeff_marg_pval <- 1 - ecdf(vegeff_marg_sim)(vegeff_marg)

# put into a df
effdata <- data.frame(
  veg = "Shrubland",
  y = 0.35,
  text = paste(
    format(round(vegeff_marg, 2), nsmall = 2), " % ",
    "(", format(round(vegeff_marg_pval, 3), nsmall = 3), ")", sep = ""
  )
)

# two in one:

props_dens_long$label_violin <- "Burn probability in\nrandomized fires"
data_vd$label_vd <- "Relative abundance of\nvegetation types"
props_data$label_point <- "Burn probability\nin observed fires"

global_burn_prob <- mean(data$burned)

veg_marg <-
  ggplot(props_dens_long, aes(x = vegetation_class, y = p)) +

  # veg distribution
  geom_bar(data = data_vd,
           mapping = aes(x = vegetation_class, y = veg_dist,
                         fill = label_vd),
           stat = "identity",
           position = position_dodge2(width = 2, padding = 0.05),
           width = 0.7) +
  # scale_fill_viridis(discrete = T, alpha = 0.3, option = "D", begin = 0.6) +
  scale_fill_viridis(discrete = T, alpha = 0.2, option = "A", begin = 0.5) +
  # scale_fill_viridis(discrete = T, alpha = 0.45, option = "D", begin = 0.7) +

  # mean burn probability in pixels
  geom_hline(yintercept = mean(data$burned),
             linetype = "dashed", linewidth = 0.2, alpha = 0.8) +

  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +

  geom_violin(linewidth = 0.3, alpha = 0.2,
              mapping = aes(fill = label_violin, color = label_violin)) +
  scale_fill_manual(values = "black") +
  scale_color_manual(values = "black") +

  ggnewscale::new_scale_color() +

  geom_point(data = props_data, mapping = aes(x = veg, pobs, color = label_point),
             inherit.aes = F, size = 2.5) +
  # scale_color_viridis(discrete = T, option = "A", begin = 0.2,
  #                     guide = guide_legend(order = 1)) +
  scale_color_viridis(#discrete = T, option = "D", begin = 0,
                      discrete = T, option = "A", begin = 0.15,
                      guide = guide_legend(order = 1)) +

  # individual p-values
  geom_text(data = props_data, mapping = aes(x = veg, y = max + 0.025,
                                             label = p),
            color = "gray10",
            inherit.aes = F, size = 2.8) +

  # global effect
  geom_text(data = effdata, mapping = aes(x = veg, y = y, label = text),
            color = "gray10",
            inherit.aes = F, size = 2.8) +
  ylab("Burn probability and\nrelative abundance (%)") +
  xlab("Vegetation type") +
  scale_y_continuous(labels = scales::label_percent(suffix = ""),
                     expand = c(0, 0), limits = c(0, 0.45),
                     breaks = seq(0, 0.4, by = 0.1)) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 9),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 11, hjust = 0.5),
        axis.title.x = element_text(size = 10, vjust = 2),
        plot.title = element_text(hjust = -0.05),
        legend.title = element_blank()) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2))
veg_marg

effdata2 <- effdata
effdata2$veg <- 1
effdata2$y <- 0.175
effdata2$text <- "1.83 %\n(0.196)"

fill_colors <- c(
  viridis(5, option = "A", end = 0.85), # the most abundant
  viridis(2, option = "A", begin = 0.96, end = 1)
)

# Vegetation marginal -----------------------------------------------------

props_data2 <- props_data
names(props_data)[1] <- "vegetation_class"

veg_marg2 <-
  ggplot() +
  # mean burn probability in pixels
  geom_hline(yintercept = mean(data$burned),
             linetype = "dashed", linewidth = 0.3, alpha = 0.8) +

  geom_bar(data = props_data,
           mapping = aes(vegetation_class, pobs, fill = vegetation_class),
           stat = "identity", alpha = 0.7, color = "black", linewidth = 0.26,
           width = 0.9) +
  scale_fill_manual(values = fill_colors) +

  # individual p-values
  geom_text(data = props_data, mapping = aes(x = vegetation_class,
                                             y = 0.1425,
                                             label = p),
            color = "gray10",
            inherit.aes = F, alpha = 0.85, size = 2.8) +

  # global effect
  geom_text(data = effdata2, mapping = aes(x = veg, y = y, label = text),
            color = "gray10", alpha = 0.85, size = 2.8, #3
            inherit.aes = F) +
  ylab("Burn probability (%)") +
  xlab("Vegetation type") +
  scale_y_continuous(labels = scales::label_percent(suffix = ""),
                     expand = c(0, 0), limits = c(0, 0.2),
                     breaks = seq(0, 0.2, by = 0.05)) +
  labs(tag = "A") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 9),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 11, hjust = 0.5),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.05),
        legend.position = "none",
        plot.margin = margin(2, 2, 4, 2, "mm")) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2))
veg_marg2



# Burn prob ~ vegetation, conditional -------------------------------------

predictors_env <- predictors[predictors != "ndvi_mean"]

#summary(ndvi_model)

# get means to compare. As pp and elevation have high effects on vegetation types,
# a common value is unrealistic, so we make 2 conditions: low and high.
# The constrain for each common value is that they must not be outside the
# [0.15, 0.85] percentiles in the remaining communities
# (after removing the extreme).
qq <- aggregate(cbind(pp, elevation) ~ vegetation_class, data_comp, quantile,
                probs = c(0.15, 0.85), method = 8)

concensus <- function(x) { # get consensus when col 1 is the lower and col 2 is the upper
  # x <- qq$pp
  x <- x[, grep("15|85", colnames(x))]
  colnames(x) <- c("min", "max")
  # find communities from the extremes
  id_lower <- which.min(x[, "max"])
  id_upper <- which.max(x[, "min"])

  # find who may live with the extremes.
  neigh_lower <- which(x[, "min"] <= x[id_lower, "max"])
  neigh_upper <- which(x[, "max"] >= x[id_upper, "min"])

  low_val <- mean(c(min(x[neigh_lower, "max"]), max(x[neigh_lower, "min"])))
  high_val <- mean(c(min(x[neigh_upper, "max"]), max(x[neigh_upper, "min"])))

  return(c(low_val, high_val))
}

# table to show means
veg_means <- aggregate(as.matrix(data_comp[, predictors_env]) ~ vegetation_class,
                       data_comp, mean)
veg_means$aspect <- aggregate(aspect ~ vegetation_class,
                              data_comp, mean_circular_deg)[, "aspect"]
global_means <- lapply(veg_means[, -1], function(x) mean(x))
global_means$aspect <- mean_circular_deg(veg_means$aspect)

# make values to predict (replace global_means for pp and elevation)
ref_values <- global_means
ref_values$pp <- concensus(qq$pp)
ref_values$elevation <- concensus(qq$elevation)

# for dist_roads, use the median, so it's not so extreme for dry forests
ref_values$dist_roads <- median(veg_means$dist_roads, method = 8)

# compare this values with slope and TPI being predicted by elevation
ref_values2 <- ref_values

ref_values2$slope <- predict(mtopo[["elevation"]][["slope"]][["all"]],
                             data.frame(elevation = ref_values$elevation),
                             type = "response") %>% as.numeric
ref_values2$TPI2k <- predict(mtopo[["elevation"]][["TPI2k"]][["all"]],
                             data.frame(elevation = ref_values$elevation),
                             type = "response") %>% as.numeric

names1 <- rep(names(ref_values), sapply(ref_values, length))
names2 <- rep(names(ref_values2), sapply(ref_values2, length))

# check visually:
ref_df <- data.frame(val = c(do.call("c", ref_values),
                             do.call("c", ref_values2)),
                     var = c(names1, names2),
                     group = rep(c("fixed", "cov"),
                                 c(length(names1), length(names2))))

data_comp_long <- pivot_longer(
  data_comp, all_of(which(names(data_comp) %in% predictors_env)),
  names_to = "var", values_to = "val"
)
data_comp_long$vari <- factor(data_comp_long$var,
                              levels = names_frame$variable[-5],
                              labels = names_frame$name2[-5])
ref_df$vari <- factor(ref_df$var,
                      levels = names_frame$variable[-5],
                      labels = names_frame$name2[-5])

# ggplot(data_comp_long, aes(x = val, fill = vegetation_class, color = vegetation_class)) +
#   geom_density(alpha = 0.25) +
#   scale_color_viridis(option = "A", discrete = TRUE, end = 0.95) +
#   scale_fill_viridis(option = "A", discrete = TRUE, end = 0.95) +
#
#   # ggnewscale::new_scale_color() +
#   geom_vline(data = ref_df, mapping = aes(xintercept = val,
#                                           linetype = group, group = group),
#              linewidth = 0.3,
#              show.legend = F) +
#   scale_linetype_manual(values = c("dotted", "dashed")) +
#
#   facet_wrap(vars(vari), scales = "free", nrow = 2,
#              strip.position = "bottom") +
#   theme(legend.position = c(0.88, 0.23),
#         legend.title = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.x = element_blank(),
#         strip.placement = "outside",
#         strip.background = element_rect(color = "white", fill = "white"),
#         strip.text = element_text(vjust = 1),
#         panel.spacing.y = unit(4, "mm"),
#         panel.grid.minor = element_blank())
#
# ggsave("figures/XX predictors densities by veg type_cov-topo.png",
#        width = 16, height = 13, units = "cm")

# newdata to predict
vars_unique <- ref_values2 # choose: 2 co-varies topography, ref_values, not.
vars_unique$vegetation_class <- veg_labels_sub
pd <- expand.grid(vars_unique)

# remove combinations that should not exist
out <- (
         pd$elevation == ref_values2$elevation[1] &
         (
           pd$slope == ref_values2$slope[2] |
           pd$TPI2k == ref_values2$TPI2k[2]
         )
       ) |
       (
         pd$elevation == ref_values2$elevation[2] &
         (
           pd$slope == ref_values2$slope[1] |
           pd$TPI2k == ref_values2$TPI2k[1]
         )
       )
pd <- pd[!out, ]

# # write values:
# write.csv(pd[pd$vegetation_class == "Shrubland", predictors_env],
#           "exports/common_environment_values_used_cov-topo.csv",
#           row.names = F)

# Here, use the NDVI model to predict NDVI at these environments
pd$ndvi_mean <- predict(ndvi_models$veg, pd, type = "response")

pd$elevation_class <- factor(
  ifelse(pd$elevation == min(pd$elevation), "Low", "High"),
  levels = c("High", "Low")
)

pd$pp_class <- factor(
  ifelse(pd$pp == min(pd$pp), "Low", "High"),
  levels = c("Low", "High")
)

pd$show <- factor(
  ifelse(
    (pd$vegetation_class == "Grassland" & pd$pp_class == "High") |
    (pd$vegetation_class == "Wet forest" & pd$pp_class == "Low") |
    (pd$vegetation_class == "Subalpine\nforest" & pd$elevation_class == "Low") |
    (pd$vegetation_class == "Dry forest" & pd$elevation_class == "High"),
    "no", "yes"),
  levels = c("yes", "no")
)

pd$pfit <- predict(fire_models_obs$veg, pd, "response")

# compute predictions from randomized fires, to compare all predictions with
# chance
nsim <- ncol(dsim)
pmat <- matrix(NA, nrow(pd), nsim)
for(i in 1:nsim) {
  pmat[, i] <- predict(fire_models_sim$veg[[i]], pd, "response")
}

pd$pval <- sapply(1:nrow(pd), function(i) {
  pcum <- ecdf(pmat[i, ])(pd$pfit[i])
  pval <- ifelse(pcum < 0.5, pcum * 2, (1 - pcum) * 2)
})
pd$pval <- format(round(pd$pval, 3), nsmall = 3)

# make titles as nested facets
pd$elevation_tit <- "Elevation"
pd$pp_tit <- "Precipitation"

# bring the observed points
names(props_data)[1] <- "vegetation_class"
pdobs <- left_join(pd, props_data[, 1:2], by = "vegetation_class")

# get max probability to label pvalue
pdobs$pmax <- apply(as.matrix(pdobs[, c("pfit", "pobs")]), 1, max)

# put labels in each panel
panlabs <- pdobs[!duplicated(pdobs[, c("elevation_class", "pp_class")]),
                 c("elevation_class", "pp_class",
                   "elevation_tit", "pp_tit")]
panlabs$vegetation_class[panlabs$pp_class == "Low"] <- "Wet forest"
panlabs$vegetation_class[panlabs$pp_class == "High"] <- "Grassland"

pdobs$pfit %>% range # check useful limit
panlabs$pfit <- 18
panlabs$label <- c("(3)", "(1)", "(4)", "(2)")#c("C", "A", "D", "B")

# Compute categorical effect at each condition
vegeffs <- panlabs
vegeffs$pfit <- 14.75
vegeffs$eff <- NA
vegeffs$pval <- NA
vegeffs$text <- NA
# import conditional vegetation model to estimate vegetation weights
vegmod <- readRDS("exports/vegetation_model_gam.rds")
pd_uni <- pdobs[, c(predictors_env,
                    "elevation_class", "pp_class", "elevation_tit", "pp_tit")]
pd_uni <- pd_uni[!duplicated(pd_uni), ]
vpred <- predict(vegmod, pd_uni, type = "response")

# compute effects
for(cond in 1:nrow(vegeffs)) {
  # cond = 1
  # get vegetation types at each condition
  eee <- vegeffs$elevation_class[cond]
  ppp <- vegeffs$pp_class[cond]
  rows_focal <- pdobs$elevation_class == eee &
                pdobs$pp_class == ppp &
                pdobs$show == "yes"

  pdsub <- pdobs[rows_focal, ]
  pmat_sub <- pmat[rows_focal, ]

  veg_ids <- which(veg_labels_sub %in% pdsub$vegetation_class)

  w <- vpred[cond, veg_ids]

  eff_obs <- categorical_effect2(y = pdsub$pfit, w = w) * 100

  # compute the effect across all simulations, for p-val
  eff_sim <- sapply(1:ncol(pmat_sub), function(j) {
    categorical_effect2(y = pmat_sub[, j], w = w) * 100
  })
  pval <- 1 - ecdf(eff_sim)(eff_obs)

  vegeffs$eff[cond] <- eff_obs
  vegeffs$pval[cond] <- pval
  vegeffs$text[cond] <- paste(
    format(round(eff_obs, 2), nsmall = 2), " %\n",
    "(", format(round(pval, 3), nsmall = 3), ")", sep = ""
  )

}


# Vegetation conditional --------------------------------------------------

# plot
veg_cond <-
ggplot(pdobs[pdobs$show == "yes", ],
       aes(x = vegetation_class, y = pfit * 100,
           fill = vegetation_class)) +

  geom_hline(yintercept = mean(data_comp$burned) * 100,
             linetype = "dashed", linewidth = 0.3, alpha = 0.8) +

  geom_bar(stat = "identity", width = 0.9, color = "black",
           alpha = 0.7, linewidth = 0.3) +

  # p values
  geom_text(aes(x = vegetation_class, y = 17,#(pmax + 0.02) * 100,
                label = pval), alpha = 0.85,
            data = pdobs[pdobs$show == "yes", ],
            size = 2.8, inherit.aes = F) +

  #scale_fill_viridis(option = "A", discrete = TRUE, end = 0.85) +
  scale_fill_manual(values = fill_colors[1:5]) +

  # labels for panels
  geom_text(aes(x = vegetation_class, y = pfit, label = label),
            data = panlabs,
            size = 4.2, inherit.aes = F) +

  # effects
  geom_text(aes(x = vegetation_class, y = pfit, label = text),
            data = vegeffs, alpha = 0.85,
            size = 2.8, inherit.aes = F) +

  facet_nested(rows = vars(elevation_tit, elevation_class),
               cols = vars(pp_tit, pp_class)) +

  xlab("Vegetation type") +
  ylab("Burn probability (%)") +
  scale_y_continuous(expand = c(0, 0, 0, 0),
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 5)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  # ggtitle("B") +
  labs(tag = "C") +
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11,
                                  margin = margin(0.01, 0.01, 1, 1, "mm")),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = -0.05))

veg_cond


veg_dist <-
  ggplot(data_vd, aes(x = "", y = veg_dist, fill = vegetation_class)) +
  scale_fill_manual(values = fill_colors) +
  geom_bar(stat = "identity", width = 1,
           alpha = 0.7, linewidth = 0.3, color = "black") +
  coord_polar("y", start = 0) +
  labs(tag = "B",
       title = "Relative abundance\nin the landscape") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title.position = "panel",
        plot.title = element_text(vjust = 0, hjust = 0.5, size = 9),
        panel.border = element_blank())
veg_dist

veg_all <- grid.arrange(
  grobs = list(veg_marg2,
               veg_dist + theme(plot.margin = margin(1, 1, 20, 1, "mm")),
               veg_cond),
  widths = c(2.1, 1),
  heights = c(1.3, 2),
  layout_matrix = matrix(c(1, 2, 3, 3), 2, 2, byrow = T)
)

ggsave("figures/06) veg_all.png", plot = veg_all,
       width = 16, height = 21, units = "cm")


# veg cond not deleting veg types -----------------------------------------

# plot
veg_cond_all <-
  ggplot(pdobs,
         aes(x = vegetation_class, y = pfit * 100,
             fill = vegetation_class,
             alpha = show,
             color = show)) +
  scale_alpha_manual(values = c(0.85, 0.3)) +
  scale_color_manual(values = c("black", "gray")) +

  geom_hline(yintercept = mean(data_comp$burned) * 100,
             linetype = "dashed", linewidth = 0.3, alpha = 0.8) +

  geom_bar(stat = "identity", width = 0.7,
           linewidth = 0.3) +

  # p values
  geom_text(aes(x = vegetation_class, y = 20,
                label = pval, alpha = show),
            data = pdobs,
            size = 2.8, inherit.aes = F) +

  # scale_fill_viridis(discrete = T, end = 0.9, option = "A") +
  scale_fill_viridis(option = "A", discrete = TRUE, end = 0.85) +
  # scale_fill_viridis(option = "C", discrete = TRUE, end = 0.85) +
  # scale_fill_viridis(discrete = T, end = 0.9, option = "A") +

  # observed proportions
  ggnewscale::new_scale_fill() +
  geom_point(data = pdobs,
             mapping = aes(x = vegetation_class, pobs * 100, alpha = show),
             inherit.aes = F, size = 2.5,
             #color = viridis(1, option = "A", begin = 0.15)) +
             color = viridis(1, option = "A", begin = 0.15)) +

  # # labels for panels
  # geom_text(aes(x = vegetation_class, y = 27, label = label),
  #           data = panlabs,
  #           size = 4.2, inherit.aes = F) +

  # # effects
  # geom_text(aes(x = vegetation_class, y = 25, label = text),
  #           data = vegeffs, alpha = 0.85,
  #           size = 2.8, inherit.aes = F) +

  facet_nested(rows = vars(elevation_tit, elevation_class),
               cols = vars(pp_tit, pp_class)) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.9) +

  xlab("Vegetation type") +
  ylab("Burn probability (%)") +
  scale_y_continuous(expand = c(0, 0, 0, 0),
                     limits = c(0, 22),
                     breaks = seq(0, 20, by = 5)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11,
                                  margin = margin(0.01, 0.01, 1, 1, "mm")),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
veg_cond_all

ggsave("figures/06) vegetation - conditional all types.png",
       width = 16, height = 14, units = "cm")

# Write reference values --------------------------------------------------

# save reference values and means by veg:
r1 <- sapply(ref_values2, function(x) {x[1]})
r2 <- sapply(ref_values2, function(x) {
  if(length(x) == 2) return(x[2]) else return(NA)
})
names(r1) <- names(r2) <- predictors_env
ref_mat <- rbind(r1, r2)

vegm_mat <- as.matrix(veg_means[, -1])
rownames(vegm_mat) <- veg_levels_sub
veg_text <- vegm_mat
veg_text[,] <- NA

# compute percentiles of the reference values in the corresponding distribution
for(p in 1:length(predictors_env)) {
  for(v in 1:nv) {
    # p = 2; v = 1
    mtext <- format(round(vegm_mat[v, p], 2), nsmall = 2)

    ee <- ecdf(data_comp[data_comp$vegetation_class == veg_labels_sub[v],
                         predictors_env[p]])
    p1 <- ee(ref_mat[1, p])

    # no percentiles for aspect
    # if(p != which(predictors_env == "aspect")) {
      p1text <- format(round(p1 * 100, 2), nsmall = 2)
    # } else {
      # p1text <- ""
    # }

    if(!is.na(ref_mat[2, p])) {
      p2 <- ee(ref_mat[2, p])
      p2text <- format(round(p2 * 100, 2), nsmall = 2)
      separ <- ", "
      percent2 <- " %"
    } else {
      separ <- p2text <- percent2 <- ""
    }

    veg_text[v, p] <- paste(mtext, "\n(", p1text, " %", separ,
                            p2text, percent2, ")", sep = "")

    if(p == which(predictors_env == "aspect")) {
      veg_text[v, p] <- mtext
    }
  }
}

# format ref values
rr <- apply(ref_mat, 2, function(x) format(round(x, 2), nsmall = 2))
export_table <- rbind(veg_text, rr) %>% t
rownames(export_table) <- names_frame$name2[names_frame$variable %in% predictors_env]
# View(export_table)
write.csv(export_table,
          "exports/common_environment_reference_and_vegetation_means_cov-topo.csv")

## NDVI reference values
ndvi_ref_long <- pd[, c("ndvi_mean", "pp_class", "elevation_class",
                        "vegetation_class")]
ndvi_ref <- pivot_wider(ndvi_ref_long, names_from = "vegetation_class",
                        values_from = "ndvi_mean")

# add percentile in the NDVI distribution
ndvi_ref2 <- ndvi_ref

for(v in 1:nv) {
  # v = 1
  veg <- veg_labels_sub[v]
  ee <- ecdf(data_comp[data_comp$vegetation_class == veg, "ndvi_mean"])
  ndvi_vals <- ndvi_ref[, veg, drop = T]
  percs <- ee(ndvi_vals)

  texts <- character(nrow(ndvi_ref))
  for(r in 1:nrow(ndvi_ref)) {
    # r = 1
    texts[r] <- paste(
      format(round(ndvi_vals[r], 3), nsmall = 3),
      "\n(",
      format(round(percs[r] * 100, 2), nsmall = 2),
      " %)",
      sep = ""
    )
  }

  ndvi_ref2[, veg] <- texts
}
# View(ndvi_ref2)

write.csv(ndvi_ref2, "exports/common_environments_ndvi_ref_cov-topo.csv",
          row.names = F)


# Correlation between fire drivers ----------------------------------------

predictors2 <- predictors
predictors2[predictors == "ndvi_mean"] <- "pp"
predictors2[predictors == "pp"] <- "ndvi_mean"
newlabs <- c("elev", "slope", "aspect", "top. pos.",
             "pp", "ndvi", "dist settl.", "dist roads")

# function to add smooth to ggpairs, from
# https://stackoverflow.com/questions/37889222/change-colors-in-ggpairs-now-that-params-is-deprecated
smooth_gg <- function(data, mapping, points = list(), smooth = list(), ...) {
  ggplot(data = data, mapping = mapping, ...) +
    do.call(geom_point, points) +
    do.call(geom_smooth, smooth)
}

cc <- viridis(1, option = "A", begin = 0.5)

# Plot
for(v in 1:V) {
  # v = 1
  if(v == 1) {
    rows <- 1:nrow(data)
  } else {
    rows <- data$vegetation_class == veg_model[v]
  }
  datav <- data[rows, predictors2]
  pa1 <-
    ggpairs(
      datav,
      lower = list(continuous = wrap(smooth_gg,
                                     points = list(size = 2, colour = "black", alpha = 0.05),
                                     smooth = list(method = "gam", se = F, linewidth = 0.7, colour = cc,
                                                   formula = y ~ s(x, k = 4, bs = "cr")))),
      upper = list(continuous = wrap("cor", size = 2.5)),
      columnLabels = newlabs,
      axisLabels = "none"
    ) +
    theme(panel.grid.minor = element_blank()) +
    ggtitle(veg_model2[v])

  nn <- paste("figures/", "pairs_by_veg_", v, ".png", sep = "")
  ggsave(filename = nn, plot = pa1,
         width = 16, height = 16, units = "cm")
}


# Correlations with temperature -water -------------------------------------

GGally::ggpairs(data[, c("pp", "temp", "wb", "elevation")])


cor.test(data$elevation, data$temp_wc) # -0.77, WordlClim
cor.test(data$elevation, data$temp) # -0.67 , Atlas climatico digital (Bianchi and Cravero)
cor.test(data$pp, data$temp) # -0.72

# temperature ~ elevation
ggplot(data, aes(elevation, temp)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam")

# temperature ~ pp
ggplot(data, aes(pp, temp)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam")

# pp ~ lat
ggplot(data, aes(lat, pp)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam")

# pp ~ lat
ggplot(data, aes(long, pp)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam")

GGally::ggpairs(data[, c("pp", "temp", "elevation", "slope", "TPI2k")])

# pp:temp = -0.717
# elev:temp = -0.673
# elev:pp = 0.183

# no usar temp porque se asocia muy fuertemente a dos cosas.
mean_ci(data$temp, name = "temp_") # 6.198543 [2.171532, 9.738004]
range(data$temp)
mean_ci(data$pp, name = "pp")      # 1256.412 [580.000, 2508.525]
range(data$pp)

# Lo bueno de la temp es que también considera el gradiente latitudinal,
# cuyo eff es 1/3 del efecto de la elev en magnitud, pero la corr con pp
# es tan alta que trae problemas.

m_elev <- summary(lm(temp ~ elevation, data = data))
m_pp <- summary(lm(temp ~ pp, data = data))
m_ppelev <- summary(lm(temp ~ pp + elevation, data = data))

m_elev$r.squared
m_pp$r.squared
m_ppelev$r.squared # high R2!

m_wbelev <- summary(lm(temp ~ wb + elevation, data = data))
m_wbelev$r.squared # 0.835



# Ordination - environment and vegetation ---------------------------------

data_veg |> names()

data_veg$northing <- cos(data_veg$aspect * pi / 180)
hist(data_veg$northing)

data_ord <- data_veg[, c("elevation", "slope", "northing", "TPI2k",
                         "pp", "dist_human", "dist_roads")]
data_ord_z <- scale(data_ord)

ord <- vegan::rda(data_ord_z)
plot(ord)

biplot(ord, display = c("sites",
                        "species"),
       type = c("text",
                "points"))

ord$Ybar |> str()

cc <- vegan::scores(ord)
cc$sites

sumord <- summary(ord)
imp1 <- round(sumord$cont$importance["Proportion Explained", "PC1"] * 100, 2)
imp2 <- round(sumord$cont$importance["Proportion Explained", "PC2"] * 100, 2)

lab1 <- paste("Axis 1 (", imp1, " %)", sep = "")
lab2 <- paste("Axis 2 (", imp2, " %)", sep = "")

cc$sites

ordtab <- data.frame(
  x1 = cc$sites[, 1],
  x2 = cc$sites[, 2],
  veg = data_veg$vegetation_class
)



ggplot(ordtab, aes(x1, x2, fill = veg, color = veg, shape = veg)) +
  # geom_point() +
  ggdensity::geom_hdr(probs = c(0.90), alpha = 0.2) +
  scale_fill_viridis(option = "B", end = 0.8, discrete = T, direction = -1,
                     name = "Vegetation\ntype") +
  scale_color_viridis(option = "B", end = 0.8, discrete = T, direction = -1,
                      name = "Vegetation\ntype") +
  xlab(lab1) +
  ylab(lab2) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/XX environmantal variables ordintaion with veg type.png",
       width = 14, height = 10, units = "cm")


# Refit models to compute R2 ----------------------------------------------

# vegetation model
vegmod <- glm(burned ~ vegetation_class, family = "binomial", data = data_veg)
r2bern(fitted(vegmod))
0.01507833

# environmental model
envmod <- bam(
    burned ~
      s(elevation, k = 3, bs = "cr") +
      s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") +
      s(TPI2k, k = 4, bs = "cr") +
      s(pp, k = 3, bs = "cr") +
      s(ndvi_mean, k = 4, bs = "cr") +
      s(dist_human, k = 4, bs = "cr") +
      s(dist_roads, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data_veg, discrete = T, nthreads = 8
  )
r2bern(fitted(envmod))
0.08354483

jointmod <- bam(
    burned ~
      vegetation_class +
      s(elevation, by = vegetation_class, k = 3, bs = "cr", id = 1) +
      s(slope, by = vegetation_class, k = 4, bs = "cr", id = 2) +
      s(aspect, by = vegetation_class, k = 4, bs = "cc", id = 3) +
      s(TPI2k, by = vegetation_class, k = 4, bs = "cr", id = 4) +
      s(pp, by = vegetation_class, k = 3, bs = "cr", id = 5) +
      s(ndvi_mean, by = vegetation_class, k = 4, bs = "cr", id = 6) +
      s(dist_human, by = vegetation_class, k = 4, bs = "cr", id = 7) +
      s(dist_roads, by = vegetation_class, k = 4, bs = "cr", id = 8),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data_veg, discrete = T, nthreads = 8
)
r2bern(fitted(jointmod))
0.1272726

# Likelihood ration test

lmtest::lrtest(envmod, jointmod)

# #     Df  LogLik     Df  Chisq Pr(>Chisq)
# 1 15.949 -1612.6
# 2 72.532 -1488.1 56.583  249.01  < 2.2e-16 ***

lmtest::lrtest(vegmod, jointmod)

# #     Df  LogLik     Df  Chisq Pr(>Chisq)
# 1  5.000 -1826.9
# 2 72.532 -1488.1 67.532 677.65  < 2.2e-16 ***


# Comparing vegetation types while controlling environment -----------------

# As vegetation types interact with physical variables, comparing them in 
# a single point may no be very representative of their differences. 
# Here I compute a Partial Prediction Plot for the vegetation effect, which
# marginalises over the remaining predictors. However, to avoid the problem
# of the correlation among vegetation and environment, I subset the data in 
# four combinations of high and low elevation and precipitation. Those four 
# subsets are tanken as the distribution of the environment, and the 
# predictions of burn probability will be marginal to those 4 distributions.
# This can be thought of as a partitioned partial dependence plot.

data_veg$elev_class <- NA
data_veg$pp_class <- NA

filter_elev_low <- data_veg$elevation >= 500 & data_veg$elevation < 900
filter_elev_high <- data_veg$elevation >= 900 & data_veg$elevation < 1300
data_veg$elev_class[filter_elev_low] <- "elev_low"#"Low: [500, 900) m a.s.l."
data_veg$elev_class[filter_elev_high] <- "elev_high"#High: [900, 1300) m a.s.l."

filter_pp_low <- data_veg$pp >= 750 & data_veg$pp < 1000
filter_pp_high <- data_veg$pp >= 1000 & data_veg$pp < 1500
data_veg$pp_class[filter_pp_low] <- "pp_low"#"Low: [750, 1000) mm"
data_veg$pp_class[filter_pp_high] <- "pp_high"#"High: [1000, 1500) mm"

# subset data in those ranges

datasub <- data_veg[!is.na(data_veg$elev_class) &
                    !is.na(data_veg$pp_class), ]

table(datasub$elev_class, datasub$pp_class, datasub$vegetation_class)
# OK

datasub$envir <- paste(datasub$elev_class, datasub$pp_class, sep = "/")
envirs <- unique(datasub$envir)
datasub$envir <- factor(datasub$envir, levels = envirs)

# Communities in each environment
envir_table <- data.frame(
  "elev_high/pp_high" = c("Wet forest", "Subalpine\nforest", "Shrubland"),
  "elev_low/pp_low" = c("Dry forest", "Shrubland", "Grassland"),  
  "elev_high/pp_low" = c("Subalpine\nforest", "Shrubland", "Grassland"),
  "elev_low/pp_high"  = c("Wet forest", "Dry forest", "Shrubland")
)

# proportion of each vegetation type within each environment, 
# to compute weighted average and effect size
envir_prop <- envir_table
for(e in 1:4) {
  # e = 1
  tabb <- table(datasub$vegetation_class[datasub$envir == envirs[e]])
  props <- normalize(tabb[envir_table[, e]])
  envir_prop[, e] <- props
}

pdata_env <- lapply(1:4, function(e) {
  datalocal0 <- datasub[datasub$envir == envirs[e], ]
  datalocal <- datalocal0[datalocal0$vegetation_class %in% envir_table[, e], ]
  return(datalocal[, c("vegetation_class", predictors)])
})

names(pdata_env) <- envirs
# lapply(pdata_env, function(x) unique(x$vegetation_class))

# Create a prediction list with the 4 envirs. Each environment will have a 
# prediction table, with the predicted probs, observed marginal probs, 
# and p-values. 
# In addition, they will have the effect size and its p-value, in a separate
# (1-row) df.

props_data2 <- props_data
names(props_data)[1] <- "vegetation_class"

predlist <- lapprops_data2predlist <- lapply(1:4, function(e) {
  # e = 1
  print(e)
  ptable <- data.frame(
    vegetation_class = factor(envir_table[, e], levels = veg_labels_sub),
    pfit = NA,
    pval = NA,
    envir = envirs[e]
  )
  
  # bring marginal burn prob
  ptable <- left_join(ptable, props_data, by = "vegetation_class")  
  
  # table to fill with pfit distribution under null models
  pfit_sim <- matrix(NA, np, nsim)
  
  # Expand prediction data
  pdata_local <- do.call("rbind", lapply(1:np, function(v) {
    ppp <- pdata_env[[e]]
    ppp$vegetation_class <- factor(ptable$vegetation_class[v], 
                                   levels = veg_labels_sub)
    return(ppp)
  }))
  
  # predict and aggregate with observed model
  pfit_all <- predict(fire_models_obs$veg, pdata_local, type = "response")
  pfit_means <- tapply(pfit_all, pdata_local$vegetation_class, mean)
  ptable$pfit <- pfit_means[!is.na(pfit_means)] 
  
  # Predictions from null models
  psim_all <- sapply(1:nsim, function(i) {
    predict(fire_models_sim$veg[[i]], pdata_local, type = "response")
  })
  
  psim_means <- aggregate(psim_all ~ vegetation_class, pdata_local, mean)[, -1] |> as.matrix()
  
  # The p-value for each vegetation type is the p-value for the difference between
  # the weighted average of burn probs across veg types and its value. 
  fitted_mean <- as.numeric(t(ptable$pfit) %*% envir_prop[, e])
  diff_obs <- abs(ptable$pfit - fitted_mean)
  
  diff_sim <- abs(t(psim_means) - as.numeric(t(psim_means) %*% envir_prop[, e]))
  
  ptable$pval <- sapply(1:3, function(v) {
    distrib <- ecdf(diff_sim[, v])
    return(1 - distrib(diff_obs[v]))
  })
  
  # Effect size with its pvalue
  efftable <- data.frame(
    effect = NA, 
    pval = NA,
    envir = envirs[e],
    pmean = fitted_mean * 100
  )
  
  # Compute the categorical effect and p-value
  efftable$effect <- categorical_effect2(y = ptable$pfit * 100,
                                         w = envir_prop[, e])
  eff_distrib <- sapply(1:nsim, function(j) {
    categorical_effect2(y = psim_means[, j] * 100,
                        w = envir_prop[, e])
  })
  efftable$pval <- 1 - ecdf(eff_distrib)(efftable$effect)
  
  # turn probabilities into percentages:
  ptable$pfit <- ptable$pfit * 100
  ptable$pobs <- ptable$pobs * 100
  
  return(
    list(
      ptable = ptable,
      efftable = efftable
    )
  )
  
})

ddd <- datasub[, c("envir", "elev_class", "pp_class")]
ddd <- ddd[!duplicated(ddd), ]

# Extract predicted proababilities
pred_table <- do.call("rbind", lapply(predlist, function(x) x[["ptable"]]))
pred_table <- left_join(pred_table, ddd, by = "envir")

eff_table <- do.call("rbind", lapply(predlist, function(x) x[["efftable"]]))
eff_table <- left_join(eff_table, ddd, by = "envir")


# Vegetation partitioned PDP ---------------------------------------------

fill_colors <- c(
  viridis(5, option = "A", end = 0.85), # the most abundant
  viridis(2, option = "A", begin = 0.96, end = 1)
)

# Replace labels for elevation and pp class
pred_table$elev_class2 <- factor(
  pred_table$elev_class, 
  levels = c("elev_high", "elev_low"),
  labels = c("High\n[900; 1300) m a.s.l.", "Low\n[500; 900) m a.s.l.")
) 
eff_table$elev_class2 <- factor(
  eff_table$elev_class, 
  levels = c("elev_high", "elev_low"),
  labels = c("High\n[900; 1300) m a.s.l.", "Low\n[500; 900) m a.s.l.")
) 
pred_table$elev_title <- "Elevation"
eff_table$elev_title <- "Elevation"

## pp
pred_table$pp_class2 <- factor(
  pred_table$pp_class, 
  levels = c("pp_low", "pp_high"),
  labels = c("Low\n[750; 1000) mm", "High\n[1000; 1500) mm")
) 
eff_table$pp_class2 <- factor(
  eff_table$pp_class, 
  levels = c("pp_low", "pp_high"),
  labels = c("Low\n[750; 1000) mm", "High\n[1000; 1500) mm")
) 
pred_table$pp_title <- "Precipitation"
eff_table$pp_title <- "Precipitation"

# put labels in each panel
eff_table$label <- c("B", "C", "A", "D")
eff_table$vegetation_class <- factor(
  c("Grassland", "Wet forest", "Wet forest", "Grassland"),
  levels = veg_labels_sub
)
eff_table$y_label <- 17
eff_table$y_effect <- 2.8

eff_table$eff_label <- paste(
  round(eff_table$effect, 2), " %\n",
  "(", round(eff_table$pval, 3), ")", sep = ""
)

# plot
veg_ppdp <-
  
  # fitted probs data
  ggplot(pred_table,
         aes(x = vegetation_class, y = pfit, fill = vegetation_class)) +
  
  # average prob by environment
  geom_hline(data = eff_table, mapping = aes(yintercept = pmean),
             linetype = "dashed", linewidth = 0.3, alpha = 0.8) +
  
  # fitted probs bar
  geom_bar(stat = "identity", width = 0.9, color = "black",
           alpha = 0.7, linewidth = 0.3) +
  
  # p values
  geom_text(aes(x = vegetation_class, y = 17,#(pmax + 0.02) * 100,
                label = pval), alpha = 0.85,
            size = 2.8, inherit.aes = F) +
  
  #scale_fill_viridis(option = "A", discrete = TRUE, end = 0.85) +
  scale_fill_manual(values = fill_colors[1:5]) +
  
  # labels for panels
  geom_text(aes(x = vegetation_class, y = y_label, label = label),
            data = eff_table,
            size = 4.2, inherit.aes = F) +
  
  # effects
  geom_text(aes(x = vegetation_class, y = y_effect, label = eff_label),
            data = eff_table, alpha = 0.85,
            size = 2.8, inherit.aes = F) +
  
  facet_nested(rows = vars(elev_title, elev_class2),
               cols = vars(pp_title, pp_class2)) +
  
  xlab("Vegetation type") +
  ylab("Burn probability (%)") +
  scale_y_continuous(expand = c(0, 0, 0, 0),
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 5)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  # ggtitle("B") +
  labs(tag = "C") +
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11,
                                  margin = margin(0.01, 0.01, 1, 1, "mm")),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = -0.05))

veg_ppdp


veg_dist <-
  ggplot(data_vd, aes(x = "", y = veg_dist, fill = vegetation_class)) +
  scale_fill_manual(values = fill_colors) +
  geom_bar(stat = "identity", width = 1,
           alpha = 0.7, linewidth = 0.3, color = "black") +
  coord_polar("y", start = 0) +
  labs(tag = "B",
       title = "Relative abundance\nin the landscape") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title.position = "panel",
        plot.title = element_text(vjust = 0, hjust = 0.5, size = 9),
        panel.border = element_blank())
veg_dist

veg_all <- grid.arrange(
  grobs = list(veg_marg2,
               veg_dist + theme(plot.margin = margin(1, 1, 20, 1, "mm")),
               veg_ppdp),
  widths = c(2.1, 1),
  heights = c(1.3, 2),
  layout_matrix = matrix(c(1, 2, 3, 3), 2, 2, byrow = T)
)

ggsave("figures/06) veg_all.png", plot = veg_all,
       width = 16, height = 21, units = "cm")



veg_all <- grid.arrange(
  grobs = list(veg_marg2,
               veg_dist + theme(plot.margin = margin(1, 1, 20, 1, "mm")),
               veg_ppdp),
  widths = c(2.1, 1),
  heights = c(1.3, 2),
  layout_matrix = matrix(c(1, 2, 3, 3), 2, 2, byrow = T)
)
ggsave("figures/06) veg_all_ppdp.png", plot = veg_all,
       width = 16, height = 21, units = "cm")
