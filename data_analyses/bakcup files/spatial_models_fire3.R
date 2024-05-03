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
v2 <- extract(pp, v0)
v0$pp <- v2$geonode_precipitacion_anual
v <- project(v0, "EPSG:5343")
out_ids <- is.na(v$pp)
v <- v[!out_ids, ]
# veg_type distribution
burned_veg0 <- read.csv("data_and_files/data_burned_and_available_area_by_vegetation_dryforest2.csv")


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


# correlation elev~temp
cor.test(data$elevation, data$temp) # -0.77

# Fire ~ vegetation -------------------------------------------------------

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
names(props_dens)
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

veg_burn <-
ggplot(props_dens_long, aes(x = vegetation_class, y = p)) +
  geom_violin(linewidth = 0.3, alpha = 0.2, fill = "black") +
  geom_point(data = props_data, mapping = aes(x = veg, pobs),
             inherit.aes = F, size = 2,
             color = viridis(1, option = "A", begin = 0.2)) +
  geom_label(data = props_data, mapping = aes(x = veg, y = upper95 * 1.6,
                                              label = p),
             color = "gray10",
             inherit.aes = F, size = 2.5) +
  ylab("Burned percentage (%)") +
  scale_y_continuous(labels = scales::label_percent(suffix = ""),
                     expand = c(0, 0), limits = c(0, 0.45),
                     breaks = seq(0, 0.4, by = 0.1)) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(hjust = -0.05)) +
  ggtitle("A")
veg_burn


# Vegetation type distribution

veg_dist <-
  ggplot(data_vd, aes(x = vegetation_class, y = veg_dist)) +
  geom_bar(stat = "identity",
           position = position_dodge2(width = 2, padding = 0.05),
           width = 0.7,
           fill = viridis(1, alpha = 0.7, option = "A", begin = 0.5)) +
  # scale_color_viridis(discrete = TRUE, option = "A", end = 0.5) +
  # scale_fill_viridis(discrete = TRUE, option = "A", end = 0.5) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 9),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10, vjust = 2),
        plot.title = element_text(hjust = -0.05)) +
  ylab("Relative abundance\nin the landscape (%)") +
  xlab("Vegetation type") +
  ggtitle("B") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45),
                     breaks = seq(0, 0.4, by = 0.1),
                     labels = scales::label_percent(suffix = ""))
veg_dist

fireveg <- ggarrange(veg_burn, veg_dist, ncol = 1)


# two in one:

props_dens_long$label_violin <- "Burn probability in\nrandomized fires"
data_vd$label_vd <- "Relative abundance of\nvegetation types"
props_data$label_point <- "Burn probability\nin observed fires"

global_burn_prob <- mean(data$burned)

veg_burn2 <-
  ggplot(props_dens_long, aes(x = vegetation_class, y = p)) +

  # veg distribution
  geom_bar(data = data_vd,
           mapping = aes(x = vegetation_class, y = veg_dist,
                         fill = label_vd),
           stat = "identity",
           position = position_dodge2(width = 2, padding = 0.05),
           width = 0.7) +
  scale_fill_viridis(discrete = T, alpha = 0.3, option = "A", begin = 0.5) +

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
             inherit.aes = F, size = 3) +
  scale_color_viridis(discrete = T, option = "A", begin = 0.2,
                      guide = guide_legend(order = 1)) +

  # individual p-values
  geom_text(data = props_data, mapping = aes(x = veg, y = max + 0.025,
                                             label = p),
             color = "gray10",
             inherit.aes = F, size = 2.5) +

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
        axis.title.y = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10, vjust = 2),
        plot.title = element_text(hjust = -0.05),
        legend.title = element_blank()) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2))
veg_burn2

ggsave("figures/04) spatial patterns - fire by veg type.png", plot = veg_burn2,
       width = 16, height = 9, units = "cm")

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

saveRDS(firelist, "data_and_files/fire_models_univariate.rds")
firelist <- readRDS("data_and_files/fire_models_univariate.rds")

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



# NDVI model --------------------------------------------------------------

# used to predict NDVI at different values of topographic variables
# and precipitation. In this way, the multiple model including NDVI
# can show partial effects of covariates but retaining the NDVI effect.

predictors_env <- predictors[predictors != "ndvi_mean"]
rows_comp <- data$vegetation_class %in% veg_labels_sub
data_comp <- data[rows_comp, ]
nv <- veg_labels_sub %>% length

ndvi_models <- lapply(1:nv, function(v) {
  # v = 1
  rows_local <- data_comp$vegetation_class == veg_labels_sub[v]
  data_comp_local <- data_comp[rows_local, ]
  model <- bam(
    ndvi_mean ~
      s(elevation, k = 4, bs = "cr") +
      s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") +
      s(TPI2k, k = 4, bs = "cr") +
      s(pp, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = betar(), discrete = T, nthreads = 8,
    data = data_comp_local
  )
  return(model)
})
names(ndvi_models) <- veg_labels_sub

# add model for all vegetation types
ndvi_all <- bam(
  ndvi_mean ~
    s(elevation, k = 4, bs = "cr") +
    s(slope, k = 4, bs = "cr") +
    s(aspect, k = 4, bs = "cc") +
    s(TPI2k, k = 4, bs = "cr") +
    s(pp, k = 4, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  family = betar(), discrete = T, nthreads = 8,
  data = data
)

ndvi_models_6 <- c(list(ndvi_all), ndvi_models)
names(ndvi_models_6)[[1]] <- veg_model[1]

predictors_ndvi <- predictors_env[grep("dist", predictors_env, invert = T)]

# # checking smooth centering
# m1 <- ndvi_models_6[[1]]
#
# tt <- predict(m1, type = "terms")
# intercept <- coef(m1)[1]
# colMeans(tt) # Means are close to zero
# colSums(tt) # Sums are far from zero
# hist(rowSums(tt))
# rowSums(tt) %>% mean
#
# tt_int <- intercept + tt
# p_terms <- plogis(tt_int)
# p_terms %>% colMeans
# p_terms %>% colMeans %>% mean
# mean(fitted(m1))
# mean(data$ndvi_mean)

# NDVI model predictions --------------------------------------------------

ndvi_pred <- do.call("rbind", lapply(veg_labels_sub, function(veg) {
  # veg = "Shrubland"
  dlocal <- data_comp[data_comp$vegetation_class == veg, predictors_ndvi]
  nnn <- do.call("rbind", lapply(predictors_ndvi, function(p) {
    # p = "elevation"
    if(p != "aspect") {
      rr <- make_newdata(p, dlocal, HDI = TRUE, prob = 0.98)
    } else {
      rr <- make_newdata(p, dlocal, HDI = FALSE)
    }
    return(rr)
  }))
  nnn$vegetation_class <- veg
  return(nnn)
}))
# nrow(newdata_comp) # 150 * length(predictors_ndvi) * 5 # OK

ndvi_pred$pfit <- NA
ndvi_pred$plower <- NA
ndvi_pred$pupper <- NA

for(v in 1:nv) {
  rows_local <- ndvi_pred$vegetation_class == veg_labels_sub[v]
  ndfit <- predict(ndvi_models[[v]], ndvi_pred[rows_local, ], se.fit = T)
  ndvi_pred$pfit[rows_local] <- plogis(ndfit$fit)
  ndvi_pred$plower[rows_local] <- plogis(ndfit$fit - qnorm(0.975) * ndfit$se.fit)
  ndvi_pred$pupper[rows_local] <- plogis(ndfit$fit + qnorm(0.975) * ndfit$se.fit)
}


pred_ndvi_names <- names_frame$var_name[names_frame$variable %in% predictors_ndvi]
pred_ndvi_names[5] <- "(E) Precipitation\n(mm / year)"

ndvi_pred$varying_var <- factor(ndvi_pred$varying_var,
                                levels = predictors_ndvi,
                                labels = pred_ndvi_names)
ndvi_pred$vegetation_class <- factor(ndvi_pred$vegetation_class,
                                     levels = veg_labels_sub)


# x breaks for each predictor
xbreaks <- list(
  "elevation" = c(500, 1000, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = seq(0, 1, by = 0.25),#c(0, 0.5, 1),
  "pp" = seq(500, 2000, by = 500)#c(1000, 2000)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

pndvi_list <- vector("list", length(predictors_ndvi) + 1)
for(i in 1:length(predictors_ndvi)) {
  # i = 3
  ddd <- ndvi_pred[ndvi_pred$varying_var == pred_ndvi_names[i], ]

  pp <- ggplot(ddd, aes(x = varying_val, y = pfit,
                        ymin = plower, ymax = pupper,
                        color = vegetation_class, fill = vegetation_class)) +
    geom_ribbon(color = NA, alpha = 0.4) +
    geom_line() +
    scale_color_viridis(discrete = T, end = 0.85, option = "A") +
    scale_fill_viridis(discrete = T, end = 0.85, option = "A") +
    xlab(pred_ndvi_names[i]) +
    ylim(min(ndvi_pred$plower), max(ndvi_pred$pupper)) +
    scale_x_continuous(breaks = xbreaks[[i]],
                       labels = xlabels[[i]]) +
    # limits = c(0, 360),
    # expand = c(0, 0)) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10))
  pndvi_list[[i]] <- pp
}

legleg <- ggpubr::get_legend(pndvi_list[[1]] + theme(legend.position = "right")) %>% ggpubr::as_ggplot()
pndvi_list[[length(predictors_ndvi) + 1]] <- legleg

lablab <- textGrob("NDVI", gp = gpar(fontsize = 12), hjust = -0.2, rot = 90)
ndviplot <- egg::ggarrange(plots = pndvi_list, nrow = 2, left = lablab)

ggsave("figures/S06) NDVI model predictions only vegs.png", plot = ndviplot,
       width = 16, height = 12, units = "cm")



# NDVI model predictions with all veg types -------------------------------

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
  "elevation" = c(500, 1000, 1500),
  "slope" = c(0, 20, 40, 60),
  "aspect" = c(0, 90, 180, 270, 360),
  "TPI2k" = seq(0, 1, by = 0.25),#c(0, 0.5, 1),
  "pp" = seq(500, 2000, by = 500)#c(1000, 2000)
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

    lpmat_full <- predict(ndvi_models_6[[veg]], ndvi_pred, type = "lpmatrix")
    use_cols <- c(1, grep(var, colnames(lpmat_full)))

    fit <- lpmat_full[, use_cols] %*% coef(ndvi_models_6[[veg]])[use_cols]
    Sigma <- vcov(ndvi_models_6[[veg]], unconditional = T, freq = F)[use_cols, use_cols]
    se.fit <- apply(lpmat_full[, use_cols], MARGIN = 1,
                    FUN = function(x) sqrt(t(x) %*% Sigma %*% x))

    ndvi_pred$pfit <- plogis(fit) %>% as.numeric
    ndvi_pred$plower <- plogis(fit - qnorm(0.975) * se.fit) %>% as.numeric
    ndvi_pred$pupper <- plogis(fit + qnorm(0.975) * se.fit) %>% as.numeric

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

# Multivariate fire models --------------------------------------------------

# here we do not save models and we don't save densities. The same densities
# computed in univariate models will be used.
ob_cond <- c("model", "pred", "pred_sim", "r2", "p")
O_cond <- length(ob_cond)

# pred is the data frame to plot everything
# pred_sim is the matrix of randomized predictions

# create nested lists to save all
firelist_cond <- lapply(veg_model, function(v) {
  olist <- vector("list", O_cond)
  names(olist) <- ob_cond
  return(olist)
})
names(firelist_cond) <- veg_model

P <- length(predictors)

for(v in 1:V) {
  # v = 4
  veg <- veg_model[v]
  print(veg)

  # create local dataset
  if(v == 1) { # all veg types
    rows <- 1:nrow(data)
    dlocal <- data[, c("burned", predictors)]
    dsim_local <- dsim
  } else {
    rows <- data$vegetation_class == veg
    dlocal <- data[rows, c("burned", predictors)]
    dsim_local <- dsim[rows, ]
  }
  names(dlocal) <- c("burned", predictors)

  # Fit model to observed data
  mod <- bam(
    burned ~
      s(elevation, k = 4, bs = "cr") +
      s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") +
      s(TPI2k, k = 4, bs = "cr") +
      s(pp, k = 4, bs = "cr") +
      s(ndvi_mean, k = 4, bs = "cr") +
      s(dist_human, k = 4, bs = "cr") +
      s(dist_roads, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = dlocal, discrete = T, nthreads = 8
  )
  firelist_cond[[v]]$model <- mod

  # compute prediction (with ci)

  # new data
  newdata <- do.call("rbind", lapply(predictors, function(p) {
    if(p != "aspect") {
      rr <- make_newdata(p, dlocal, HDI = TRUE, prob = 0.98)
    } else {
      rr <- make_newdata(p, dlocal, HDI = FALSE)
    }
    return(rr)
  }))

  # in all cases where NDVI is not the varying variable, predict NDVI from
  # the NDVI model, but using only the basis functions of the varying variable
  # (average effect over the remaining predictors).
  # This does not apply to predictions for distances, which do not affect NDVI
  rows_ndvi_predict <- !(newdata$varying_var %in% c("ndvi_mean", "dist_human",
                                                    "dist_roads"))
  lp_ndvi_full <- predict(ndvi_models_6[[v]], newdata, type = "lpmatrix")
  lp_ndvi <- lp_ndvi_full
  lp_ndvi[,] <- 0
  lp_ndvi[, 1] <- 1

  for(pn in predictors_ndvi) {
    # pn = "elevation"
    rows <- newdata$varying_var == pn
    cols <- grep(pn, colnames(lp_ndvi))
    lp_ndvi[rows, cols] <- lp_ndvi_full[rows, cols]
  }

  newdata$ndvi_mean[rows_ndvi_predict] <-
    plogis(lp_ndvi[rows_ndvi_predict, ] %*% coef(ndvi_models_6[[v]]))

  # The same computation for burn prob: instead of fixing the covariates at their
  # means, we simply remove from the design matrix (lpmatrix) the non-varying
  # predictors.
  lpmat_full <- predict(mod, newdata, type = "lpmatrix")
  lpmat <- lpmat_full
  lpmat[,] <- 0
  lpmat[, 1] <- 1 # intercept
  for(p in 1:P) {
    #p = 1
    focal_var_bin <- as.numeric(newdata$varying_var == predictors[p])
    focal_cols <- grep(predictors[p], colnames(lpmat))
    lpmat_local <- lpmat_full[, focal_cols]
    # make zero the non-focal
    lpmat_local <- lpmat_local * focal_var_bin
    lpmat[, focal_cols] <- lpmat_local
  }

  # But we keep the NDVI effect in the case of physical variables, as their
  # effect is mediated by NDVI, but we are interested in the whole effect.
  ndvi_cols <- grep("ndvi_mean", colnames(lpmat))
  lpmat[rows_ndvi_predict, ndvi_cols] <- lpmat_full[rows_ndvi_predict, ndvi_cols]

  # compute predictions
  fit <- lpmat %*% coef(mod)
  Sigma <- vcov(mod, unconditional = T, freq = F)
  se.fit <- apply(lpmat, MARGIN = 1, FUN = function(x) sqrt(t(x) %*% Sigma %*% x))

  pred_ppp <- data.frame(p_mle = plogis(fit),
                         p_lower = plogis(fit - qnorm(0.975) * se.fit),
                         p_upper = plogis(fit + qnorm(0.975) * se.fit))

  # now do the same but fixing NDVI at its mean effect (no co-variation, fixed="f")
  newdata_f <- do.call("rbind", lapply(predictors, function(p) {
    if(p != "aspect") {
      rr <- make_newdata(p, dlocal, HDI = TRUE, prob = 0.98)
    } else {
      rr <- make_newdata(p, dlocal, HDI = FALSE)
    }
    return(rr)
  }))

  lpmat_full_f <- predict(mod, newdata_f, type = "lpmatrix")
  lpmat_f <- lpmat_full_f
  lpmat_f[,] <- 0
  lpmat_f[, 1] <- 1 # intercept

  for(p in 1:P) {
    #p = 1
    focal_var_bin <- as.numeric(newdata_f$varying_var == predictors[p])
    focal_cols <- grep(predictors[p], colnames(lpmat_f))
    lpmat_local <- lpmat_full_f[, focal_cols]
    # make zero the non-focal
    lpmat_local <- lpmat_local * focal_var_bin
    lpmat_f[, focal_cols] <- lpmat_local
  }

  # compute predictions
  fit_f <- lpmat_f %*% coef(mod)
  se.fit_f <- apply(lpmat_f, MARGIN = 1, FUN = function(x) sqrt(t(x) %*% Sigma %*% x))

  pred_ppp_f <- data.frame(p_mle_f = plogis(fit_f),
                           p_lower_f = plogis(fit_f - qnorm(0.975) * se.fit_f),
                           p_upper_f = plogis(fit_f + qnorm(0.975) * se.fit_f))

  pred_obs <- cbind(newdata, pred_ppp, pred_ppp_f) # _f has fixed NDVI
  # check: # plot(p_mle ~ p_mle_f, pred_obs)

  # compute observed r2
  firelist_cond[[v]]$r2 <- r2bern(fitted(mod))

  # Get predicted linear predictor for randomized fires
  sims_lp <- sapply(1:ncol(dsim_local), function(j) {
    if(j %% 50 == 0) print(j)
    # j = 1
    data_sim <- dlocal
    data_sim$burned <- dsim_local[, j]
    mran <-  bam(
      burned ~
        s(elevation, k = 4, bs = "cr") +
        s(slope, k = 4, bs = "cr") +
        s(aspect, k = 4, bs = "cc") +
        s(TPI2k, k = 4, bs = "cr") +
        s(pp, k = 4, bs = "cr") +
        s(ndvi_mean, k = 4, bs = "cr") +
        s(dist_human, k = 4, bs = "cr") +
        s(dist_roads, k = 4, bs = "cr"),
      knots = list(aspect = c(0, 360)),
      family = "binomial", data = data_sim, discrete = T, nthreads = 8
    )
    linpred_ran <- (lpmat %*% coef(mran)) #  - coef(mran)[1] # do not subtract
    return(linpred_ran)
  })

  # # add the observed intercept (did not work when NDVI was in the way)
  # pred_sim <- plogis(sims_lp + coef(mod)[1])

  # For every prediction (varying_var), force equal mean along the curve for
  # observed and randomized, at logit scale.

  sims_lp_back <- sims_lp # just backup
  for(p in predictors) {
    # p = "elevation"
    rr <- newdata$varying_var == p
    meanfit <- mean(fit[rr])
    meanran <- colMeans(sims_lp_back[rr, ])
    mdiff <- meanfit - meanran
    mdiff_mat <- outer(rep(1, sum(rr)), mdiff)
    sims_lp[rr, ] <- sims_lp_back[rr, ] + mdiff_mat
  }

  # # check centering
  # rrr <- newdata$varying_var == "aspect"
  # plot(plogis(fit[rrr]) ~ newdata$aspect[rrr], type = "l")
  # for(i in 1:ncol(sims_lp)) {
  #   lines(plogis(sims_lp[rrr, i]) ~ newdata$aspect[rrr], col = rgb(1, 0, 0, 0.05))
  # }

  pred_sim <- plogis(sims_lp)
  firelist_cond[[v]]$pred_sim <- pred_sim

  # compute quantiles
  pred_sim_q <- apply(pred_sim, 1, quantile,
                      probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
                      method = 8) %>% t %>% as.data.frame
  colnames(pred_sim_q) <- c("lower95", "lower90", "lower80", "median",
                            "upper80", "upper90", "upper95")

  # merge with observed prediction
  pred <- cbind(pred_obs, pred_sim_q)
  firelist_cond[[v]]$pred <- pred

  # p-values over the sequence
  pvalseq <- sapply(1:nrow(pred_sim), function(j) {
    distrib <- ecdf(pred_sim[j, ] %>% as.numeric())
    perc <- distrib(pred_obs$p_mle[j])
    p <- ifelse(perc >= 0.5, 2 * (1 - perc), 2 * perc)
    return(p)
  })
  firelist_cond[[v]]$p <- pvalseq
}
saveRDS(firelist_cond, "data_and_files/fire_models_multivariate_varying-NDVI-at-prediction.rds")
firelist_cond <- readRDS("data_and_files/fire_models_multivariate_varying-NDVI-at-prediction.rds")


# Fitted values in pooled model -------------------------------------------

mm <- firelist_cond[[1]]$model
tt <- predict(mm, type = "terms")
colMeans(tt)
colMeans(tt) %>% sum

# Multivariate fire plots ---------------------------------------------------

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

# factor to multiply derivatives y position
delta_yfac <- c(
  "elevation" = 0.85,
  "slope" = 0.85,
  "aspect" = 0.15,
  "TPI2k" = 0.85,
  "ndvi_mean" = 0.85,
  "pp" = 0.85,
  "dist_human" = 0.85,
  "dist_roads" = 0.85
)

delta_xfac <- c(
  "elevation" = 0.85,
  "slope" = 0.85,
  "aspect" = 0.85,
  "TPI2k" = 0.35,
  "ndvi_mean" = 0.2,
  "pp" = 0.85,
  "dist_human" = 0.85,
  "dist_roads" = 0.85
)

# get sd by variable to standardize derivatives (mad is harder to interpret,
# and the circular equivalent is also harder.)
# mads <- apply(as.matrix(data[, predictors]), 2, mad)
# mads1 <- apply(as.matrix(data[, predictors]), 2, mad, constant = 1)
sds <- apply(as.matrix(data[, predictors]), 2, sd)
# aspect_circ <- as.circular(data$aspect, type = "angles", template = "geographic",
#                            units = "degrees", rotation = "clock")
aspect_circ <- as.circular(data$aspect * pi / 180, type = "angles", template = "geographic",
                           units = "radians", rotation = "clock")
# in radians
sds["aspect"] <- sd(aspect_circ) * 180 / pi # sd always in radians,

P <- length(predictors)
plist_multi <- vector("list", P)

# to compare effects later, save derivatives
delta_list <- vector("list", P)

## instead of derivatives, compute the weighted difference between the curve
## fitted with real vs simulated data.
## The p value may be for this or for every point in the curve.


for(p in 1:P) {
  # p = 6
  var <- predictors[p]

  # extract data from all vegs
  pred <- do.call("rbind", lapply(1:V, function(v) {
    # observed predictions
    d <- firelist_cond[[v]]$pred
    ds <- firelist_cond[[v]]$pred_sim

    rows_use <- d$varying_var == var

    d <- d[rows_use, ]
    ds <- ds[rows_use, ]

    d$vegetation_class <- veg_model[v]
    colnames(ds) <- paste("p_mle_sim_", 1:ncol(ds), sep = "")

    return(cbind(d, ds))
  }))

  # extract density from the univariate firelist
  dens_data <- do.call("rbind", lapply(1:V, function(v) {
    # v = 1; p = 1
    d <- firelist[[v]][[p]]$dens
    d$vegetation_class <- veg_model[v]

    # rescale the density to the new y_max
    maxp <- max(pred[pred$vegetation_class == veg_model[v], c("p_mle", "upper95")])
    dens_factor <- max(d$density) / maxp
    d$density <- d$density / dens_factor

    return(d)
  }))

  pred$vegetation_class <- factor(pred$vegetation_class,
                                  levels = veg_model,
                                  labels = veg_model_num)
  dens_data$vegetation_class <- factor(dens_data$vegetation_class,
                                       levels = veg_model,
                                       labels = veg_model_num)

  # compute deltas with p-values
  delta_data <- do.call("rbind", lapply(1:V, function(i) {
    # i = 6
    cols_use <- grep("p_mle", colnames(pred))
    remove <- grep("p_mle_f", colnames(pred))
    cols_use <- cols_use[cols_use != remove]

    dens_local <- dens_data[dens_data$vegetation_class == veg_model_num[i], ]
    pred_local <- pred[pred$vegetation_class == veg_model_num[i], ]
    pred_prob <- as.matrix(pred_local[, cols_use])

    weights <- approx(dens_local$varying_val,
                      dens_local$density,
                      xout = pred_local$varying_val, rule = 2)$y %>% normalize

    # pointwise mean in randomized fires
    p_means <- colSums(pred_prob * weights)

    # for each curve, get the avg absolute difference with the mean
    delta_mat <- abs(pred_prob - p_means) * 100
    delta_means <- colSums(delta_mat * weights)

    # metric and p value
    delta <- delta_means[1]
    pval <- 1 - ecdf(delta_means[-1])(delta)

    # make texts
    text_delta <- paste(format(round(delta, 2), nsmall = 2), "%")
    text_pval <- paste("(", format(round(pval, 3), nsmall = 3), ")", sep = "")
    if(text_pval == "(0.000)") text_pval <- "(< 0.002)"
    text <- paste(text_delta, "\n", text_pval, sep = "")

    res <- data.frame(delta = delta,
                      pval = pval,
                      text = text,
                      vegetation_class = veg_model_num[i])
    return(res)
  }))

  delta_data$predictor <- predictors[p]
  delta_list[[p]] <- delta_data

  # set position for delta. Different computation because y starts always at zero,
  # but it varies in max between veg types
  dxrange <- range(dens_data$varying_val)
  x_delta <- dxrange[1] + delta_xfac[p] * (diff(dxrange))

  max_y <- sapply(veg_model_num, function(v) {
    yy <- c(pred$p_mle[pred$vegetation_class == v],
            pred$upper95[pred$vegetation_class == v],
            dens_data$density[dens_data$vegetation_class == v])
    return(max(yy))
  })

  delta_data$x <- x_delta
  delta_data$y <- max_y * delta_yfac[p]

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
      data = pred, mapping = aes(x = varying_val, y = p_mle,
                                 ymax = upper80, ymin = lower80),
      color = NA, alpha = 0.2, fill = "black"
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
    # delta
    geom_text(
      data = delta_data,
      mapping = aes(x, y, label = text), size = 2.2, color = "gray10"
    ) +
    ## choose which pvalue: text1 or text2

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
  plotcito

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

  plist_multi[[p]] <- plotcito
}

fire_multi_topo <- ggarrange(plots = plist_multi[1:4], ncol = 4)
ggsave("figures/XX 05) spatial patterns - fire by veg and topo - mult reg - varying NDVI.png",
       plot = fire_multi_topo,
       width = 16, height = 16.5, units = "cm")

fire_multi_pp <- ggarrange(plots = plist_multi[c(6, 5, 7, 8)], ncol = 4)
ggsave("figures/XX 06) spatial patterns - fire by veg and pp-ndvi-dist - mult reg - varying NDVI.png",
       plot = fire_multi_pp,
       width = 16, height = 16.5, units = "cm")


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

# to compare effects later, save derivatives
delta_list <- vector("list", P)
## instead of derivatives, compute the weighted difference between the curve
## fitted with real vs simulated data.
## The p value may be for this or for every point in the curve.

# get maximum y value at mle prediction with observed data
max_p <- sapply(1:V, function(v) firelist_cond[[v]]$pred$p_mle) %>% max
roof <- ceiling(max_p * 100) / 100
max_p <- 0.4; roof <- 0.4

for(p in 1:P) {
  # p = 1
  var <- predictors[p]

  # extract data from all vegs
  pred <- do.call("rbind", lapply(1:V, function(v) {
    # observed predictions
    d <- firelist_cond[[v]]$pred
    ds <- firelist_cond[[v]]$pred_sim

    rows_use <- d$varying_var == var

    d <- d[rows_use, ]
    ds <- ds[rows_use, ]

    d$vegetation_class <- veg_model[v]
    colnames(ds) <- paste("p_mle_sim_", 1:ncol(ds), sep = "")

    return(cbind(d, ds))
  }))

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

  pred$vegetation_class <- factor(pred$vegetation_class,
                                  levels = veg_model,
                                  labels = veg_model_num)
  dens_data$vegetation_class <- factor(dens_data$vegetation_class,
                                       levels = veg_model,
                                       labels = veg_model_num)

  # compute deltas with p-values
  delta_data <- do.call("rbind", lapply(1:V, function(i) {
    # i = 6
    cols_use <- grep("p_mle", colnames(pred))
    remove <- grep("p_mle_f", colnames(pred))
    cols_use <- cols_use[cols_use != remove]

    dens_local <- dens_data[dens_data$vegetation_class == veg_model_num[i], ]
    pred_local <- pred[pred$vegetation_class == veg_model_num[i], ]
    pred_prob <- as.matrix(pred_local[, cols_use])

    weights <- approx(dens_local$varying_val,
                      dens_local$density,
                      xout = pred_local$varying_val, rule = 2)$y %>% normalize

    # pointwise mean in randomized fires
    p_means <- colSums(pred_prob * weights)

    # for each curve, get the avg absolute difference with the mean
    delta_mat <- abs(pred_prob - p_means) * 100
    delta_means <- colSums(delta_mat * weights)

    # metric and p value
    delta <- delta_means[1]
    pval <- 1 - ecdf(delta_means[-1])(delta)

    # make texts
    text_delta <- paste(format(round(delta, 2), nsmall = 2), "%")
    text_pval <- paste("(", format(round(pval, 3), nsmall = 3), ")", sep = "")
    if(text_pval == "(0.000)") text_pval <- "(< 0.002)"
    text <- paste(text_delta, "\n", text_pval, sep = "")

    res <- data.frame(delta = delta,
                      pval = pval,
                      text = text,
                      vegetation_class = veg_model_num[i])
    return(res)
  }))

  delta_data$predictor <- predictors[p]
  delta_list[[p]] <- delta_data

  # set position for delta. Different computation because y starts always at zero,
  # but it varies in max between veg types
  dxrange <- range(dens_data$varying_val)
  x_delta <- dxrange[1] + delta_xfac[p] * (diff(dxrange))

  delta_data$x <- x_delta
  delta_data$y <- roof * delta_yfac[p]

  # truncate the ribbon at the roof
  vars_trunc <- c("upper95", "lower95", "upper80", "lower80")
  for(vt in vars_trunc) {
    pred[pred[, vt] > roof, vt] <- roof
  }
  # max(as.matrix(pred[, vars_trunc]))

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
      data = pred, mapping = aes(x = varying_val, y = p_mle,
                                 ymax = upper80, ymin = lower80),
      color = NA, alpha = 0.2, fill = "black"
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
    # delta
    geom_text(
      data = delta_data,
      mapping = aes(x, y, label = text), size = 2.2, color = "gray10"
    ) +
    ## choose which pvalue: text1 or text2

    facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right",
               scales = "free_y") +
    xlab(names_frame$var_name[p]) +
    ylab("Burn probability (%)") +
    scale_x_continuous(breaks = xbreaks[[p]],
                       labels = xlabels[[p]]) +
    scale_y_continuous(labels = scales::label_percent(suffix = ""),
                       limits = c(0, roof),
                       breaks = seq(0, roof, by = 0.1),
                       expand = c(0.01, 0.01)) +
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
ggsave("figures/05) spatial patterns - fire by veg and topo - mult reg - varying NDVI_fixedy.png",
       plot = fire_multi_topo,
       width = 16, height = 16.5, units = "cm")

fire_multi_pp <- ggarrange(plots = plist_multi[c(6, 5, 7, 8)], ncol = 4)
ggsave("figures/06) spatial patterns - fire by veg and pp-ndvi-dist - mult reg - varying NDVI_fixedy.png",
       plot = fire_multi_pp,
       width = 16, height = 16.5, units = "cm")


# Spatial effects summary --------------------------------------------------

delta_plot <- do.call("rbind", delta_list)
write.csv(delta_plot,
          "data_and_files/fire_models_multivariate_effects_summary.csv",
          row.names = F)

# remove distances, for their absurd effect in dry forest.
inn <- grep("dist", delta_plot$predictor, invert = T)
delta_plot <- delta_plot[inn, ]
q <- max(delta_plot$delta) / max(delta_plot$pval)
q <- 10
delta_plot$pval_scaled <- delta_plot$pval * q

delta_plot$vegetation_class <- factor(delta_plot$vegetation_class,
                               levels = veg_model_num,
                               labels = veg_model2)
delta_plot$predictor2 <- factor(delta_plot$predictor,
                        levels = names_frame$variable,
                        labels = names_frame$name3)
delta_plot$predictor3 <- factor(delta_plot$predictor,
                                levels = names_frame$variable,
                                labels = names_frame$name4)

delta_plot$delta_label <- "difference"
delta_plot$p_label <- "p-value"

ggplot(delta_plot) +
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
    sec.axis = sec_axis(~ . / q,
                        name = "p-value",
                        breaks = round(seq(0, 1, by = 1/4), 2))
  ) +
  xlab("Vegetation type") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figures/XX spatial effects summary - predictors.png",
       width = 16, height = 14, units = "cm")


ggplot(delta_plot) +
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
    sec.axis = sec_axis(~ . / q,
                        name = "p-value",
                        breaks = round(seq(0, 1, by = 1/4), 2))
  ) +
  xlab("Predictor") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figures/XX spatial effects summary - veg.png",
       width = 16, height = 14, units = "cm")



# Multivariate plot, varying NDVI and not ---------------------------------
# (To compare easily)

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
  "pp" = c(500, 1000, 1500),#c(1000, 2000),
  "dist_human" = c(0, 10, 20, 30),
  "dist_roads" = c(0, 10, 20, 30)
)
xlabels <- xbreaks
xlabels$aspect <- c("N", "E", "S", "W", "N")

Pndvi <- length(predictors_ndvi)
plist_ndvi <- vector("list", Pndvi)

for(p in 1:Pndvi) {
  # p = 1
  var <- predictors_ndvi[p]

  # extract data from all vegs
  pred <- do.call("rbind", lapply(1:V, function(v) {
    d <- firelist_cond[[v]]$pred
    d <- d[d$varying_var == var, ]
    d$vegetation_class <- veg_model[v]
    return(d)
  }))

  # force the same mean as when NDVI is varied (across the curve) for the
  # prediction where with fixed NDVI, at logit scale
  for(veg in unique(pred$vegetation_class)) {
    # get logit-mean
    # veg = "Wet forest"
    rr <- pred$vegetation_class == veg
    logit_vary <- qlogis(pred$p_mle[rr])
    logit_fixed <- qlogis(pred$p_mle_f[rr])
    mdif <- mean(logit_vary) - mean(logit_fixed)

    logit_fixed <- logit_fixed + mdif

    # replace centred curve
    pred$p_mle_f[rr] <- plogis(logit_fixed)
  }

  # longanize
  ccc <- which(names(pred) %in% c("p_mle", "p_mle_f"))
  pred_long <- pivot_longer(pred, all_of(ccc), names_to = "NDVI",
                            values_to = "pburn")

  pred_long$NDVI <- factor(pred_long$NDVI, levels = c("p_mle", "p_mle_f"),
                           labels = c("varying", "fixed"))

  pred_long$vegetation_class <- factor(pred_long$vegetation_class,
                                       levels = veg_model,
                                       labels = veg_model_num)

  plotcito <-
    ggplot() +

    # observed prediction
    geom_line(
      data = pred_long, mapping = aes(x = varying_val, y = pburn,
                                      color = NDVI)
    ) +
    scale_color_viridis(discrete = T, begin = 0.2, end = 0.5, option = "A") +
    facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right",
               scales = "free_y") +
    xlab(names_frame$var_name[p]) +
    ylab("Burn probability (%)") +
    scale_x_continuous(breaks = xbreaks[predictors_ndvi][[p]],
                       labels = xlabels[predictors_ndvi][[p]]) +
    scale_y_continuous(labels = scales::label_percent(suffix = "")) +
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
    plotcito <- plotcito + theme(axis.title.y = element_blank())
  }

  # remove strip for most predictors except pos topo and dist_roads
  if(var == "pp") {
    plotcito <- plotcito +
      theme(strip.background = element_rect(),
            strip.text = element_text(),
            legend.position = "right") +
      facet_wrap2(vars(vegetation_class), ncol = 1, strip.position = "right",
                  scales = "free_y",
                  strip = strip_themed(text_y = texts,
                                       background_y = backgrounds)) +
      xlab("(E) Precipitation\n(mm / year)")
  }

  plist_ndvi[[p]] <- plotcito
}


fire_multi_ndvi_varfix <- ggarrange(plots = plist_ndvi, ncol = 5)
ggsave("figures/S07) burn probability varying NDVI or not.png",
       plot = fire_multi_ndvi_varfix,
       width = 24, height = 16.5, units = "cm")


# Burn prob ~ veg at the same environmental conditions --------------------

predictors_env <- predictors[predictors != "ndvi_mean"]
rows_comp <- data$vegetation_class %in% veg_labels_sub
data_comp <- data[rows_comp, ]
nv <- veg_labels_sub %>% length

modcomp <- bam(
  burned ~
    vegetation_class +
    s(elevation, by = vegetation_class, k = 4, bs = "cr") +
    s(slope, by = vegetation_class, k = 4, bs = "cr") +
    s(aspect, by = vegetation_class, k = 4, bs = "cc") +
    s(TPI2k, by = vegetation_class, k = 4, bs = "cr") +
    s(pp, by = vegetation_class, k = 4, bs = "cr") +
    s(ndvi_mean, by = vegetation_class, k = 4, bs = "cr") +
    s(dist_human, by = vegetation_class, k = 4, bs = "cr") +
    s(dist_roads, by = vegetation_class, k = 4, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  family = "binomial", discrete = T, nthreads = 8,
  data = data_comp
)

ndvi_models <- lapply(1:nv, function(v) {
  # v = 1
  rows_local <- data_comp$vegetation_class == veg_labels_sub[v]
  data_comp_local <- data_comp[rows_local, ]
  model <- bam(
    ndvi_mean ~
      s(elevation, k = 4, bs = "cr") +
      s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") +
      s(TPI2k, k = 4, bs = "cr") +
      s(pp, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = betar(), discrete = T, nthreads = 8,
    data = data_comp_local
  )
  return(model)
})
names(ndvi_models) <- veg_labels_sub

#summary(ndvi_model)

# get means to compare. As pp and elevation have high effects on vegetation types,
# a common value is unrealistic, so we make 2 conditions: low and high.
# The constrain for each common value is that they must not be outside the
# [0.15, 0.85] percentiles in the remaining communities
# (after removing the extreme).
#
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

  return(c("low" = low_val, "high" = high_val))
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

# check visually:
vals_loose <- do.call("c", ref_values)
ref_df <- data.frame(val = vals_loose,
                     var = names(vals_loose))
ref_df$var[grep("elevation", ref_df$var)] <- "elevation"
ref_df$var[grep("pp", ref_df$var)] <- "pp"

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

ggplot(data_comp_long, aes(x = val, fill = vegetation_class, color = vegetation_class)) +
  geom_density(alpha = 0.25) +
  geom_vline(data = ref_df, mapping = aes(xintercept = val),
             linetype = "dashed", linewidth = 0.3) +
  facet_wrap(vars(vari), scales = "free", nrow = 2,
             strip.position = "bottom") +
  # scale_fill_viridis(option = "H", discrete = TRUE, end = 0.8) +
  # scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  scale_color_viridis(option = "A", discrete = TRUE, end = 0.95) +
  scale_fill_viridis(option = "A", discrete = TRUE, end = 0.95) +
  theme(legend.position = c(0.88, 0.23),
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(vjust = 1),
        panel.spacing.y = unit(4, "mm"),
        panel.grid.minor = element_blank())

ggsave("figures/XX predictors densities by veg type.png",
       width = 16, height = 13, units = "cm")

# save values:
write.csv(expand.grid(ref_values),
          "data_and_files/common_environment_values_used.csv",
          row.names = F)

# newdata to predict
vars_unique <- ref_values
vars_unique$vegetation_class <- veg_labels_sub
pd <- expand.grid(vars_unique)

# Here, use the NDVI model to predict NDVI at these environments
pd$ndvi_mean <- NA
for(v in 1:nv) {
  rows_local <- pd$vegetation_class == veg_labels_sub[v]
  pd$ndvi_mean[rows_local] <- predict(ndvi_models[[v]],
                                      pd[rows_local, ], "response")
}

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

pd$pfit <- predict(modcomp, pd, "response")

# compute predictions from randomized fires, to compare all predictions with
# chance
nsim <- ncol(dsim)
pmat <- matrix(NA, nrow(pd), nsim)

# fit model and get conditional probabilities
modcomp_sim <- lapply(1:nsim, function(i) {

  print(i)
  dlocal <- data_comp
  dlocal$burned <- dsim[rows_comp, i]

  mmod <- bam(
    burned ~
      vegetation_class +
      s(elevation, by = vegetation_class, k = 4, bs = "cr") +
      s(slope, by = vegetation_class, k = 4, bs = "cr") +
      s(aspect, by = vegetation_class, k = 4, bs = "cc") +
      s(TPI2k, by = vegetation_class, k = 4, bs = "cr") +
      s(pp, by = vegetation_class, k = 4, bs = "cr") +
      s(ndvi_mean, by = vegetation_class, k = 4, bs = "cr") +
      s(dist_human, by = vegetation_class, k = 4, bs = "cr") +
      s(dist_roads, by = vegetation_class, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = "binomial", discrete = T, nthreads = 8,
    data = dlocal
  )

  return(mmod)
})
saveRDS(modcomp_sim, "data_and_files/counterfactual_fire_models_sim_with_ndvi.rds")
modcomp_sim <- readRDS("data_and_files/counterfactual_fire_models_sim_with_ndvi.rds")

for(i in 1:nsim) {
  pmat[, i] <- predict(modcomp_sim[[i]], pd, "response")
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
panlabs$pfit <- 17
panlabs$label <- c("C", "A", "D", "B")

# Compute categorical effect at each condition
vegeffs <- panlabs
vegeffs$pfit <- 2.7
vegeffs$eff <- NA
vegeffs$pval <- NA
vegeffs$text <- NA
# import conditional vegetation model to estimate vegetation weights
vegmod <- readRDS("data_and_files/vegetation_model_gam.rds")
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


# plot
ggplot(pdobs[pdobs$show == "yes", ],
       aes(x = vegetation_class, y = pfit * 100,
           fill = vegetation_class)) +

  geom_hline(yintercept = mean(data_comp$burned) * 100,
             linetype = "dashed", linewidth = 0.3, alpha = 0.8) +

  geom_bar(stat = "identity", width = 0.7, color = "black",
           alpha = 0.7, linewidth = 0.3) +

  # p values
  geom_text(aes(x = vegetation_class, y = 17,#(pmax + 0.02) * 100,
                label = pval), alpha = 0.85,
            data = pdobs[pdobs$show == "yes", ],
            size = 2.8, inherit.aes = F) +

  # scale_fill_viridis(discrete = T, end = 0.9, option = "A") +
  scale_fill_viridis(option = "A", discrete = TRUE, end = 0.85) +
  # scale_fill_viridis(option = "C", discrete = TRUE, end = 0.85) +
  # scale_fill_viridis(discrete = T, end = 0.9, option = "A") +

  # observed proportions
  ggnewscale::new_scale_fill() +
  geom_point(data = pdobs[pdobs$show == "yes", ],
             mapping = aes(x = vegetation_class, pobs * 100),
             inherit.aes = F, size = 2,
             color = viridis(1, option = "A", begin = 0.15)) +

  # labels for panels
  geom_label(aes(x = vegetation_class, y = pfit, label = label),
            data = panlabs,
            size = 4.2, inherit.aes = F) +

  # effects
  geom_text(aes(x = vegetation_class, y = pfit, label = text),
            data = vegeffs, alpha = 0.85,
            size = 2.8, inherit.aes = F) +

  facet_nested(rows = vars(elevation_tit, elevation_class),
               cols = vars(pp_tit, pp_class)) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.9) +

  xlab("Vegetation type") +
  ylab("Burn probability (%)") +
  scale_y_continuous(expand = c(0, 0, 0, 0),
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 5)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave("figures/07) burn probability at common environments.png",
       width = 16, height = 13, units = "cm")
write.csv(pd, "data_and_files/common_environments_table.csv",
          row.names = F)


# save reference values and means by veg:
r1 <- sapply(ref_values, function(x) {x[1]})
r2 <- sapply(ref_values, function(x) {
  if(length(x) == 2) return(x[2]) else return(NA)
})
names(r1) <- names(r2) <- predictors_env
ref_mat <- rbind(r1, r2)
ref_mat[1, "dist_roads"] <- median(veg_means$dist_roads, method = 8)

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
    p1text <- format(round(p1 * 100, 2), nsmall = 2)

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

    if(p == "aspect") {
      veg_text[v, p] <- mtext
    }
  }
}

# format ref values
rr <- apply(ref_mat, 2, function(x) format(round(x, 2), nsmall = 2))
export_table <- rbind(veg_text, rr) %>% t
rownames(export_table) <- names_frame$name2[names_frame$variable %in% predictors_env]

write.csv(export_table,
          "data_and_files/common_environment_reference_and_vegetation_means.csv")


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

write.csv(ndvi_ref2, "data_and_files/common_environments_ndvi_ref.csv",
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



# Spatially explicit analysis -----------------------------------------------

# Make a few maps to show
# 1) Multiple model including the spatial smoothing term,
# 2) A model with only a spatial term

# 4) Maps with probability > and < the average (0.0577) significantly different
#    from chance (p <= 0.1) for the models with spatial terms.

# study area to clip images
sa <- vect("/home/ivan/Insync/patagonian_fires/study_area/study_area.shp")
sa <- project(sa, "EPSG:5343")
fires_vec <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
fires_vec <- project(fires_vec, "EPSG:5343")
# randomized fires to fit models many times
dsim <- readRDS("data_and_files/data_randomized_fires_mindist_800m_500_matrix.rds")

# image at 500 m to predict
data_im <- rast("data_and_files/study_area_5343_with_predictors.tif")
crds_im <- crds(data_im) %>% as.data.frame()

# remap class (dry forest A = araucaria) as subalpine
data_im$vegetation_valdivian[data_im$vegetation_valdivian == 4] <- 2
data_im$vegetation_valdivian <- as.factor(data_im$vegetation_valdivian)
crds_im <- crds(data_im) %>% as.data.frame()

# subset and rename columns
data_im <- subset(data_im, c("vegetation_valdivian",
                             "elevation", "slope", "aspect", "TPI2k",
                             "pp", "ndvi_mean_max", "dist_humans", "dist_roads"))
names(data_im) <- c("vegetation",
                    "elevation", "slope", "aspect", "tpi",
                    "pp", "ndvi", "dist_humans", "dist_roads")

# data to fit model
data_vec_points <- project(v0, "EPSG:5343")

data_vec_points$vegetation[data_vec_points$vegetation == 4] <- 2
names(data_vec_points)
names(data_im)

data_vec_points <- data_vec_points[, c("burned", "vegetation",
                                       "elevation", "slope", "aspect", "TPI2k",
                                       "pp", "ndvi_mean_", "dist_human", "dist_roads")]
names(data_vec_points) <- c("burned", "vegetation",
                            "elevation", "slope", "aspect", "tpi",
                            "pp", "ndvi", "dist_humans", "dist_roads")

data_vec_points$vegetation <- factor(data_vec_points$vegetation)

# add coordinates
crds_vec <- crds(data_vec_points) %>% as.data.frame()
data_vec <- cbind(as.data.frame(data_vec_points), crds_vec)

# turn image to data.frame (careful, it returns many pixels with NA)
im_values <- cbind(values(data_im, na.rm = F), crds(data_im, na.rm = F)) %>% as.data.frame()
inside_pixels <- !apply(im_values, 1, anyNA)
ok_veg <- im_values$vegetation > 1
ok_veg[is.na(ok_veg)] <- FALSE
valid_pixels <- inside_pixels & ok_veg
data_im_df <- im_values[valid_pixels, ]
data_im_df$vegetation <- factor(data_im_df$vegetation)

# settings for GAMs
scale_q <- diff(range(data_vec$x)) / diff(range(data_vec$y))
k_x <- 4
k_y <- ceiling(k_x / scale_q)

# Fit models to observed data
m_pred <- bam(
  burned ~
    vegetation - 1 +
    s(elevation, k = 4, bs = "cr") + s(slope, k = 4, bs = "cr") +
    s(aspect, k = 4, bs = "cc") + s(tpi, k = 4, bs = "cr") +
    s(pp, k = 4, bs = "cr") + s(ndvi, k = 4, bs = "cr") +
    s(dist_humans, k = 4, bs = "cr") + s(dist_roads, k = 4, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  family = "binomial", data = data_vec, discrete = T, nthreads = 8
)

m_pred_sp <- bam(
  burned ~
    vegetation - 1 +
    s(elevation, k = 4, bs = "cr") + s(slope, k = 4, bs = "cr") +
    s(aspect, k = 4, bs = "cc") + s(tpi, k = 4, bs = "cr") +
    s(pp, k = 4, bs = "cr") + s(ndvi, k = 4, bs = "cr") +
    s(dist_humans, k = 4, bs = "cr") + s(dist_roads, k = 4, bs = "cr") +
    # spatial term
    te(x, y, bs = "cr", k = c(k_x, k_y)),
  knots = list(aspect = c(0, 360)),
  family = "binomial", data = data_vec, discrete = T, nthreads = 8
)

m_sp <- bam(
  burned ~
    # just spatial term
    te(x, y, bs = "cr", k = c(k_x, k_y)),
  family = "binomial", data = data_vec, discrete = T, nthreads = 8
)

# predictions
preds_df <- data.frame(m_pred = numeric(nrow(data_im_df)),
                       m_pred_sp = numeric(nrow(data_im_df)),
                       m_sp = numeric(nrow(data_im_df)))
preds_df$m_pred <- predict(m_pred, data_im_df, "response")
preds_df$m_pred_sp <- predict(m_pred_sp, data_im_df, "response")
preds_df$m_sp <- predict(m_sp, data_im_df, "response")

# add the spatial effect removing the explained effect
cc <- coef(m_pred_sp)
coef_veg <- cc[grep("vegetation", names(cc))]
coef_sp <- cc[grep("te", names(cc))]
lpm <- predict(m_pred_sp, data_im_df, "lpmatrix")
lpm_veg <- lpm[, grep("vegetation", names(cc))]
lpm_sp <- lpm[, grep("te", names(cc))]

preds_df$m_pred_sp_only <- plogis(mean(lpm_veg %*% coef_veg) +
                                  lpm_sp %*% coef_sp)


## Fit models to observed and predicted data

# Function to fit 3 models and compute 4 types of predictions,
# varying the burned points data (binary vector)
data_fit <- data_vec
dsim_ext <- cbind(data_vec$burned, dsim)
colnames(dsim_ext) <- c("observed", paste("sim", 1:ncol(dsim), sep = "_"))

spat_fit <- function(burned) {

  # assign new layer
  data_fit$burned <- burned

  # Fit models to observed data
  m_pred <- bam(
    burned ~
      vegetation - 1 +
      s(elevation, k = 4, bs = "cr") + s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") + s(tpi, k = 4, bs = "cr") +
      s(pp, k = 4, bs = "cr") + s(ndvi, k = 4, bs = "cr") +
      s(dist_humans, k = 4, bs = "cr") + s(dist_roads, k = 4, bs = "cr"),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data_fit, discrete = T, nthreads = 8
  )

  m_pred_sp <- bam(
    burned ~
      vegetation - 1 +
      s(elevation, k = 4, bs = "cr") + s(slope, k = 4, bs = "cr") +
      s(aspect, k = 4, bs = "cc") + s(tpi, k = 4, bs = "cr") +
      s(pp, k = 4, bs = "cr") + s(ndvi, k = 4, bs = "cr") +
      s(dist_humans, k = 4, bs = "cr") + s(dist_roads, k = 4, bs = "cr") +
      # spatial term
      te(x, y, bs = "cr", k = c(k_x, k_y)),
    knots = list(aspect = c(0, 360)),
    family = "binomial", data = data_fit, discrete = T, nthreads = 8
  )

  m_sp <- bam(
    burned ~
      # just spatial term
      te(x, y, bs = "cr", k = c(k_x, k_y)),
    family = "binomial", data = data_fit, discrete = T, nthreads = 8
  )

  # predictions
  preds_mat <- matrix(NA, nrow(data_im_df), 4)
  colnames(preds_mat) <- c("m_pred", "m_pred_sp", "m_pred_sp_only", "m_sp")

  preds_mat[, "m_pred"] <- predict(m_pred, data_im_df, "response")
  preds_mat[, "m_pred_sp"] <- predict(m_pred_sp, data_im_df, "response")
  preds_mat[, "m_sp"] <- predict(m_sp, data_im_df, "response")

  # add the spatial effect removing the predictors' effects
  cc <- coef(m_pred_sp)
  coef_veg <- cc[grep("vegetation", names(cc))]
  coef_sp <- cc[grep("te", names(cc))]
  lpm <- predict(m_pred_sp, data_im_df, "lpmatrix")
  lpm_veg <- lpm[, grep("vegetation", names(cc))]
  lpm_sp <- lpm[, grep("te", names(cc))]

  preds_mat[, "m_pred_sp_only"] <- plogis(mean(lpm_veg %*% coef_veg) +
                                          lpm_sp %*% coef_sp)

  return(preds_mat)
}

# spat_fit(data_vec$burned) # test
# preds_list <- lapply(1:ncol(dsim_ext), function(i) {
#   print(i)
#   spat_fit(dsim_ext[, i])
# })
# names(preds_list) <- colnames(dsim_ext)
# preds_arr <- abind::abind(preds_list, along = 3)
# saveRDS(preds_arr, "data_and_files/spatial_models_predictions.rds")
preds_arr <- readRDS("data_and_files/spatial_models_predictions.rds")

pred_means <- apply(preds_arr, 1:2, mean)
pred_perc <- apply(preds_arr, 1:2, function(x) ecdf(x[2:length(x)])(x[1]))
pred_sign <- apply(pred_perc, 2, function(x) as.numeric(x > 0.5) + 1)
pred_p <- (pred_perc >= 0.5) * (1 - pred_perc) * 2 +
          (pred_perc < 0.5) * pred_perc * 2
preds_prob <- preds_arr[, , 1]

pred_bin <- pred_sign
# pred_bin[pred_sign == 1] <- "more than expected"
# pred_bin[pred_sign == -1] <- "less than expected"
pred_bin[pred_p > 0.1] <- NA

# map predictions and p-values
data_im_pred <- subset(data_im, (1:ncol(preds_prob)) + 1) # to avoid the factor
vv <- values(data_im_pred)
vv[] <- NA
vv[valid_pixels, ] <- as.matrix(preds_prob)
values(data_im_pred) <- vv
names(data_im_pred) <- colnames(preds_prob)
data_im_pred_crop <- crop(data_im_pred, sa, mask = T)
names(data_im_pred_crop) <- c("Covariates",
                              "Covariates + spatial",
                              "Spatial - covariates",
                              "Spatial only")
probs_map <- ggplot() +
  geom_spatraster(data = data_im_pred_crop, maxcell = 1e6) +
  scale_fill_viridis(na.value = "transparent", option = "F",
                     direction = 1, begin = 0, end = 1,
                     guide = "colourbar",
                     name = "burn\nprobability",
                     limits = c(0, 1)) +
  geom_spatvector(data = sa, linewidth = 0.3, color = "black", fill = NA) +
  facet_wrap(~ lyr, nrow = 1) +
  theme(panel.border = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        strip.background = element_rect(fill = "white", color = "white")) +
  scale_y_continuous(expand = c(0.02, 0.05)) #+
probs_map
# ggsave("figures/XX spatial predictions.png",
#        width = 20, height = 15, units = "cm")

# p values
data_im_p <- data_im_pred
vv <- values(data_im_p)
vv[] <- NA
vv[valid_pixels, ] <- as.matrix(pred_p)
values(data_im_p) <- vv
names(data_im_p) <- colnames(pred_p)
data_im_p_crop <- crop(data_im_p, sa, mask = T)
names(data_im_p_crop) <- names(data_im_pred_crop)

p_map <- ggplot() +
  geom_spatraster(data = data_im_p_crop, maxcell = 1e6) +
  scale_fill_viridis(na.value = "transparent", option = "A",
                     direction = -1, begin = 0, end = 0.6,
                     guide = "colourbar",
                     name = "p value",
                     limits = c(0, 0.15)) +
  geom_spatvector(data = sa, linewidth = 0.3, color = "black", fill = NA) +
  facet_wrap(~ lyr, nrow = 1) +
  theme(panel.border = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        strip.background = element_rect(fill = "white", color = "white")) +
  scale_y_continuous(expand = c(0.02, 0.05))
p_map

# regions with high and low burn probability
"more than expected"
pred_bin[pred_sign == -1] <- "less than expected"

data_im_pbin <- data_im_pred
vv <- values(data_im_pbin)
vv[] <- NA
vv[valid_pixels, ] <- pred_bin
values(data_im_pbin) <- vv
for(l in 1:4) {
  levels(data_im_pbin[[l]]) <- data.frame(ids = 1:2,
                                          class = c("low",
                                                    "high"))
}
data_im_pbin_crop <- crop(data_im_pbin, sa, mask = T)
names(data_im_pbin_crop) <- names(data_im_pred_crop)

pbin_map <- ggplot() +
  geom_spatraster(data = data_im_pbin_crop, maxcell = 1e6) +
  scale_fill_viridis(discrete = TRUE,
                     na.value = "transparent", option = "A",
                     na.translate = F, name = "burn\nprobability",
                     direction = 1, begin = 0, end = 0.6) +
  geom_spatvector(data = sa, linewidth = 0.3, color = "black", fill = NA) +
  facet_wrap(~ lyr, nrow = 1) +
  theme(panel.border = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        # legend.title = element_blank(),
        strip.background = element_rect(fill = "white", color = "white")) +
  scale_y_continuous(expand = c(0.02, 0.05))
pbin_map

# study area with fires
aa <- ggplot() +
  geom_spatvector(data = sa, linewidth = 0.3, color = "black", fill = NA) +
  geom_spatvector(data = fires_vec, linewidth = 0.3, color = "red", fill = "red") +
  theme(panel.border = element_blank(),
        legend.position = "right",
        legend.direction = "vertical") +
  scale_y_continuous(expand = c(0.02, 0.05))

map_prob_fires <- egg::ggarrange(aa, probs_map, nrow = 1, widths = c(1, 4))
ggsave("figures/XX spatial predictions burn prob.png", plot = map_prob_fires,
       width = 22, height = 15, units = "cm")

map_p_fires <- egg::ggarrange(aa, p_map, nrow = 1, widths = c(1, 4))
ggsave("figures/XX spatial predictions p value.png", plot = map_p_fires,
       width = 22, height = 15, units = "cm")

map_pbin_fires <- egg::ggarrange(aa, pbin_map, nrow = 1, widths = c(1, 4))
ggsave("figures/XX spatial predictions p value binary.png", plot = map_pbin_fires,
       width = 22, height = 15, units = "cm")

