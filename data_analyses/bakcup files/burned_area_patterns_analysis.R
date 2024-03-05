# In this code the burned_area_patterns_analysis is repeated but using
# data separating the dry forest in A and B (araucaria and cypress)

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr
library(ggh4x)    # varying strip theme for veg_types and all together

library(scales)   # log scale
library(circular) # density.circular, for aspect
library(lubridate)

library(rgdal)       # readOGR
library(rgeos)       # gIntersection, gArea...
library(maptools)    # UnionSpatialPolygons

library(mgcv)
library(scam)
library(DHARMa)

# functions ---------------------------------------------------------------

# Normalize
normalize <- function(x) x / sum(x)

# Functions to summaryze
mean_ci <- function(x, ci = 0.95, name = "mu") {
  p_low <- (1 - ci) / 2
  p_high <- ci + p_low

  qq <- quantile(x, probs = c(p_low, p_high), method = 8)
  result <- c(qq[1], mean(x), qq[2])
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

median_ci <- function(x, ci = 0.95, name = "mu") {
  p_low <- (1 - ci) / 2
  p_high <- ci + p_low

  qq <- quantile(x, probs = c(p_low, p_high), method = 8)
  result <- c(qq[1], median(x), qq[2])
  names(result) <- paste(rep(name, 3), c("lower", "median", "upper"), sep = "_")
  result
}

# To share legend (not used I think)
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# custom ggplot theme -----------------------------------------------------

# from https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_mine <- function() {
  font <- "Arial"   #assign font family up front
  marg <- 2 # figure margin in mm

  theme_bw() %+replace%    #replace elements we want to change

    theme(

      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 16,                #set font size
        #face = 'bold',            #bold typeface
        hjust = -0.1,                #left align
        vjust = 1),

      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size

      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 12),

      # para separar el eje y de los nros
      axis.title.y = element_text(
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),

      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size

      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9, family = font),

      strip.text = element_text(size = 12, family = font, color = "white"),
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),

      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

theme_set(theme_mine())



# Data for spatial analyses -----------------------------------------------

# veg_type distribution
burned_veg0 <- read.csv("data_burned_and_available_area_by_vegetation_dryforest2.csv")

# merge dry forest A with subalpine
burned_veg <- burned_veg0[burned_veg0$vegetation_code != 4, ]
burned_veg[burned_veg$vegetation_class == "Subalpine forest",
           c("area_available_ha", "area_burned_ha")] <-
  colSums(
    burned_veg0[burned_veg0$vegetation_class %in% c("Subalpine forest", "Dry forest A"),
     c("area_available_ha", "area_burned_ha")]
  )
burned_veg$vegetation_class[burned_veg$vegetation_class == "Dry forest B"] <- "Dry forest"


# spatial variables marginal to veg_type
spatial_vars <- read.csv("data_spatial_variables.csv")

# sample 20000 from each class
s1 <- sample(which(spatial_vars$class == "Burned"), size = 20000)
s2 <- sample(which(spatial_vars$class == "Burnable"), size = 20000)
spatial_vars <- spatial_vars[c(s1, s2), ]
nrow(spatial_vars)

# rename vegetation
spatial_vars$vegetation_valdivian <- 99 # altogether

# spatial variables conditional to veg_type
spatial_vars_veg0 <- read.csv("data_spatial_variables_by-veg.csv")

# spatial variables in randomized fires
drandom1 <- read.csv("data_spatial_variables_burned_simulations.csv")
drandom2 <- read.csv("data_spatial_variables_burned_simulations_2.csv")
drandom <- rbind(drandom1, drandom2)
rm("drandom1", "drandom2")

# order and filter columns
spatial_vars_veg <- spatial_vars_veg0[, c("class", "vegetation_valdivian",
                                          "elevation", "aspect", "slope", "TPI2k",
                                          "ndvi", "pp",
                                          "dist_humans", "dist_roads")]

drandom <- drandom[, c("vegetation_valdivian",
                       "elevation", "aspect", "slope", "TPI2k",
                       "ndvi", "pp",
                       "dist_humans", "dist_roads",
                       "simulation")]

spatial_vars_veg <- spatial_vars_veg[spatial_vars_veg$vegetation_valdivian > 1, ]
drandom <- drandom[drandom$vegetation_valdivian > 1, ]

# select ~ 10000 pixels by veg and class
spatial_vars_veg <- do.call(
  "rbind",
  lapply(unique(spatial_vars_veg$vegetation_valdivian), function(v) {

    # v = 2

    f1 <- which(spatial_vars_veg$vegetation_valdivian == v &
                spatial_vars_veg$class == "Burned")

    f2 <- which(spatial_vars_veg$vegetation_valdivian == v &
                spatial_vars_veg$class == "Burnable")

    n1 <- min(length(f1), 10000); n2 <- min(length(f2), 10000)

    s1 <- sample(f1, n1); s2 <- sample(f2, n2)

    return(spatial_vars_veg[c(s1, s2), ])
}))

# merge both data sets (marginal and conditional)
spatdata0 <- rbind(spatial_vars_veg,
                   spatial_vars[, names(spatial_vars_veg)])
names(spatdata0)[2] <- "vegetation_code"
drandom <- rename(drandom, vegetation_code = vegetation_valdivian)

# get name of vegetation type
spatdata <- left_join(spatdata0,
                      burned_veg0[, c("vegetation_code", "vegetation_class")],
                      by = "vegetation_code")
spatdata$vegetation_class[spatdata$vegetation_code == 99] <- "All vegetation types"

drandom <- left_join(drandom,
                     burned_veg0[, c("vegetation_code", "vegetation_class")],
                     by = "vegetation_code")

# Merge subalpine and dry forest A by sampling in accordance to their relative
# abundances

# Take 10000 pixels from Subalpine and Dry forest A in a proportion equal to
# their abundance in the landscape.
prop_sub_dry_burn <- normalize(burned_veg0$area_burned_ha[burned_veg0$vegetation_class %in% c("Subalpine forest", "Dry forest A")])
prop_sub_dry_marg <- normalize(burned_veg0$area_available_ha[burned_veg0$vegetation_class %in% c("Subalpine forest", "Dry forest A")])

sub_burned <- which(spatdata$vegetation_class == "Subalpine forest" &
                    spatdata$class == "Burned")
sub_marg <- which(spatdata$vegetation_class == "Subalpine forest" &
                  spatdata$class == "Burnable")
dry_burned <- which(spatdata$vegetation_class == "Dry forest A" &
                      spatdata$class == "Burned")
dry_marg <- which(spatdata$vegetation_class == "Dry forest A" &
                    spatdata$class == "Burnable")

nn <- length(sub_burned)

set.seed(493876)

bb <- c(
  sample(sub_burned, size = round(nn * prop_sub_dry_burn[1])),
  sample(dry_burned, size = round(nn * prop_sub_dry_burn[2]))
)

mm <- c(
  sample(sub_marg, size = round(nn * prop_sub_dry_marg[1])),
  sample(dry_marg, size = round(nn * prop_sub_dry_marg[2]))
)

spatdata <- rbind(
  spatdata[!spatdata$vegetation_class %in% c("Subalpine forest", "Dry forest A"), ],
  spatdata[c(bb, mm), ]
)

spatdata$vegetation_class[spatdata$vegetation_class == "Dry forest B"] <- "Dry forest"
spatdata$vegetation_class[spatdata$vegetation_class == "Dry forest A"] <- "Subalpine forest"

# rename anthrop type
burned_veg$vegetation_class[burned_veg$vegetation_class == "Anthropogenic prairie and shrubland"] <-
  "Anthropogenic prairie"
spatdata$vegetation_class[spatdata$vegetation_class == "Anthropogenic prairie and shrubland"] <-
  "Anthropogenic prairie"


# rename in simulated dataset
drandom$vegetation_class[drandom$vegetation_class == "Dry forest B"] <- "Dry forest"
drandom$vegetation_class[drandom$vegetation_class == "Dry forest A"] <- "Subalpine forest"
drandom$vegetation_class[drandom$vegetation_class == "Anthropogenic prairie and shrubland"] <-
  "Anthropogenic prairie"

# make levels
veg_levels <- c("Wet forest",
                "Subalpine forest",
                "Plantation",
                "Dry forest",
                "Shrubland",
                "Anthropogenic prairie",
                "Steppe and grassland",
                # "Non burnable",
                "All vegetation types")

veg_levels_ord <- c("All vegetation types",
                    "Wet forest",
                    "Subalpine forest",
                    "Plantation",
                    "Dry forest",
                    "Shrubland",
                    "Anthropogenic prairie",
                    "Steppe and grassland")

veg_levels2 <- c("Wet forest",
                 "Subalpine\nforest",
                 "Plantation",
                 "Dry forest",
                 "Shrubland",
                 "Anthropogenic\nprairie",
                 "Steppe and\ngrassland",
                 # "Non burnable",
                 "All vegetation\ntypes")
plant_anth_id <- c(3, 7)

veg_levels3 <- c("All vegetation\ntypes", # this one first
                 "Wet forest",
                 "Subalpine\nforest",
                 "Plantation",
                 "Dry forest",
                 "Shrubland",
                 "Anthropogenic\nprairie",
                 "Steppe and\ngrassland")

veg_levels3_num <- c("(1) All vegetation\ntypes", # this one first
                     "(2) Wet forest",
                     "(3) Subalpine\nforest",
                     "(3) Plantation",
                     "(4) Dry forest",
                     "(5) Shrubland",
                     "(5) Anthropogenic\nprairie",
                     "(6) Steppe and\ngrassland")

veg_levels0 <- c("Wet forest",
                "Subalpine forest",
                "Plantation",
                "Dry forest",
                "Shrubland",
                "Anthropogenic prairie",
                "Steppe and grassland")

veg_levels0_sub <- c("Wet forest",
                     "Subalpine forest",
                     "Dry forest",
                     "Shrubland",
                     "Steppe and grassland")


# Spatial patterns: Distribution of spatial variables in burned and available area ----------

# Vegetation type

burned_veg$prob_av <- c(NA, normalize(burned_veg$area_available_ha[-1]))
burned_veg$prob_b <- c(NA, normalize(burned_veg$area_burned_ha[-1]))
veg_dist_data <- data.frame(
  vegetation_class = factor(
    rep(burned_veg$vegetation_class[-1], 2),
    levels = veg_levels, labels = veg_levels2
  ),
  prob = c(burned_veg$prob_av[-1], burned_veg$prob_b[-1]),
  class = factor(
    rep(c("Burnable", "Burned"), each = nrow(burned_veg) - 1),
    levels = c("Burnable", "Burned")
    ),
  variable = "Vegetation type"
)

# compute overlap
ol_num <- 1 - sum(na.omit(abs(burned_veg$prob_av - burned_veg$prob_b))) / 2

# compute p-value for overlap from simulations.
drandom_veg <- drandom
drandom_veg$vegetation_class <- factor(drandom_veg$vegetation_class,
                                       levels = burned_veg$vegetation_class[-1])
nsim <- length(unique(drandom$simulation))
veg_ov_sim <- sapply(unique(drandom$simulation), function(s) {
  # s = 0
  tt <- table(drandom_veg$vegetation_class[drandom_veg$simulation == s])
  prop <- tt / sum(tt)
  ov <- 1 - sum(na.omit(abs(burned_veg$prob_av[-1] - prop))) / 2
  return(ov)
})

hist(veg_ov_sim, xlim = c(0.5, 1)); abline(v = ol_num)
p_val_veg <- sum(veg_ov_sim <= ol_num) / nsim
data_ol <- data.frame(
  prob = 0.35, vegetation_class = 6.5,
  lab = paste(round(ol_num * 100, 2), " % ", "(",
              round(p_val_veg, 3), ")", sep = "")
)

# plot
veg_dist <-
ggplot(veg_dist_data,
       aes(x = vegetation_class, y = prob,
           colour = class, fill = class)) +
  geom_bar(stat = "identity",
           # position = "identity",
           position = position_dodge2(width = 2, padding = 0.05),
           alpha = 0.3, width = 0.7) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  #facet_wrap(vars(variable)) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", #c(0.85, 0.8),
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.67, size = 8),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = 7),
        strip.text.x = element_text(margin = margin(3,0,3,0, "mm")),
        plot.title = element_text(hjust = -0.05)) +
  ylab("Proportion of burnable or\nburned area") +
  xlab("Vegetation type") +
  geom_text(data = data_ol, inherit.aes = FALSE, size = 3,
            mapping = aes(x = vegetation_class, y = prob, label = lab),
            colour = "black") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45),
                     breaks = seq(0, 0.4, by = 0.1))
veg_dist


# proportion burned by vegetation
burned_veg$burned_prop <- burned_veg$area_burned_ha / burned_veg$area_available_ha
burned_veg$vegetation_class2 <- factor(burned_veg$vegetation_class,
                                      levels = veg_levels, labels = veg_levels2)
veg_prop <-
  ggplot(burned_veg[burned_veg$vegetation_code > 1, ],
         aes(x = vegetation_class2, y = burned_prop)) +
  geom_bar(stat = "identity",
           # position = "identity",
           # position = position_dodge2(width = 2, padding = 0.05),
           color = "black", width = 0.4, size = 0.4) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", #c(0.85, 0.8),
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.67, size = 8),
        axis.title.x = element_blank(),# element_text(vjust = 7),
        axis.title.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(margin = margin(3,0,3,0, "mm")),
        plot.title = element_text(hjust = -0.05)) +
  ylab("Burned proportion\n(burned area / burnable area)") +
  xlab("Vegetation type") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.155),
                     breaks = c(0, 0.05, 0.10, 0.15)) +
  ggtitle("A")
veg_prop


veg_plot <- ggarrange(veg_prop,
                      veg_dist + ggtitle("B"),
                      nrow = 2)
ggsave("figures/05) vegetation distribution and proportion.jpeg",
       plot = veg_plot,
       width = 16, height = 13, units = "cm")


# Fig 6 and 7: -------------------------------------------------------------------
# Densities and overlap by veg type -----------------------------------------

var_names <- c("elevation", "slope", "aspect", "TPI2k",
               "ndvi", "pp",
               "dist_humans_km", "dist_roads_km")

n_var <- length(var_names)
n_veg <- length(veg_levels) # includes overall

# list to save densities
dlist <- vector(mode = "list", length = n_var)
names(dlist) <- var_names

# list to save overlaps
overlap_data_raw <- expand.grid(vegetation_class = veg_levels,
                            variable = names(dlist),
                            overlap_num = NA,
                            p_value = NA,
                            overlap_char = NA,
                            value = NA,   # variables to plot alongside dplots
                            density = NA)

# rescale distances to km
spatdata$dist_humans_km <- spatdata$dist_humans / 1000
spatdata$dist_roads_km <- spatdata$dist_roads / 1000

drandom$dist_humans_km <- drandom$dist_humans / 1000
drandom$dist_roads_km <- drandom$dist_roads / 1000

# get ranges to evaluate variables
ranges <- data.frame(vari = var_names,
                     mindens = rep(0, n_var),
                     maxdens = c(max(spatdata$elevation) * 1.1, 90, 360, # topo
                                 1, # tpi
                                 1, 2300, # ndvi and pp
                                 40, 30)#, # dist (km)
                     )

# Compute densities and overlap
for(i in 1:n_var) {

  # i = 6
  v <- names(dlist)[i]
  print(v)

  veg_list <- vector(mode = "list", length = n_veg)
  names(veg_list) <- veg_levels

  # loop over veg_types
  for(veg in veg_levels) {

    # veg <- "Wet forest"
    print(veg)

    dat <- spatdata %>% filter(vegetation_class == veg)

    d_burned <- density(dat[dat$class == "Burned", v],
                        from = ranges$mindens[ranges$vari == v],
                        to = ranges$maxdens[ranges$vari == v])
    d_av <- density(dat[dat$class == "Burnable", v],
                    from = ranges$mindens[ranges$vari == v],
                    to = ranges$maxdens[ranges$vari == v])

    if(v == "aspect") {
      filter_b <- dat$class == "Burned"
      filter_av <- dat$class == "Burnable"

      aspect_burned <- circular(dat$aspect[filter_b],
                                type = "angles", units = "degrees",
                                template = "geographic")
      aspect_av <- circular(dat$aspect[filter_av],
                            type = "angles", units = "degrees",
                            template = "geographic")

      bw_circ <- 20 # smaller = smoother

      d_burned <- density.circular(aspect_burned, bw = bw_circ)
      d_av <- density.circular(aspect_av, bw = bw_circ)

    }

    # densities df
    veg_list[[veg]] <- data.frame(variable = v,
                             vegetation_class = veg,
                             value = c(d_burned$x, d_av$x),
                             density = c(d_burned$y, d_av$y),
                             class = factor(rep(c("Burned", "Burnable"),
                                                each = length(d_av$x)),
                                            levels = c("Burned", "Burnable")))

    # Compute delta and overlap

    # approximate densities on a common x sequence
    xseq <- seq(ranges$mindens[i], ranges$maxdens[i], length.out = 200)
    if(v == "aspect") {
      xseq <- seq(max(d_av$x), min(d_av$x), length.out = 200)
    }
    diff_size <- abs(unique(diff(xseq))[1])

    d_burned_pred <- approx(d_burned$x, d_burned$y, xseq, method = "linear",
                            yleft = 0, yright = 0)
    d_av_pred <- approx(d_av$x, d_av$y, xseq, method = "linear",
                        yleft = 0, yright = 0)

    # circular densities are not normalized, so we get the normalizing factor.
    # (And anyway, this ensures all approximate densities are exactly
    # normalize.)
    height_b <- rowMeans(cbind(d_burned_pred$y[-1], d_burned_pred$y[-length(xseq)]))
    height_av <- rowMeans(cbind(d_av_pred$y[-1], d_av_pred$y[-length(xseq)]))

    auc_b <- sum(height_b * diff_size)
    auc_av <- sum(height_av * diff_size)

    norm_factor <- 1 / mean(auc_b, auc_av)

    den_diff <- height_b - height_av
    delta <-  ((abs(den_diff) * diff_size) %>% sum) / 2 * norm_factor
    overlap <- 1 - delta

    # fill overlap df
    fff <- overlap_data_raw$variable == v & overlap_data_raw$vegetation_class == veg
    overlap_data_raw$overlap_num[fff] <- overlap
    overlap_data_raw$overlap_char[fff] <- paste(round(overlap * 100, 2), "%")
    overlap_data_raw$value[fff] <- quantile(xseq, 0.8)
    overlap_data_raw$density[fff] <- max(c(d_burned$y, d_av$y)) * 0.8

  } # end loop across veg_types

  dlist[[v]] <- do.call("rbind", veg_list)

} # end loop across variables

# merge in one df
densdata0 <- do.call("rbind", dlist)

## Compute overlap distribution from randomized fires to get the p-value

# Make an overlap matrix where rows correspond to overlap_data_raw
ov_matrix <- matrix(NA, nrow(overlap_data_raw), nsim)

# subset available densities once
dens_av <- densdata0[densdata0$class == "Burnable", ]

# loop over replicates
for(i in 1:nsim) { # :nsim
  #i = 1
  print(i)
  data_sim0 <- drandom[drandom$simulation == (i-1), ]

  # loop over vegetation classes
  for(veg in veg_levels) {
    # veg = "All vegetation types"

    if(veg != "All vegetation types") {
      data_sim1 <- data_sim0[data_sim0$vegetation_class == veg, ]
      dens_ref1 <- dens_av[dens_av$vegetation_class == veg, ]
    } else {
      data_sim1 <- data_sim0
      dens_ref1 <- dens_av[dens_av$vegetation_class == "All vegetation types", ]
    }

    # are there enough with this veg type?
    if(nrow(data_sim1) < 100) next

    # loop over variables
    for(var in var_names) {
      # var = "aspect"
      data_sim <- data_sim1[, var]
      dens_ref <- dens_ref1[dens_ref1$variable == var, ]

      # compute density, overlap, and save.
      if(var != "aspect") {
        den_sim <- density(data_sim,
                           from = ranges$mindens[ranges$vari == var],
                           to = ranges$maxdens[ranges$vari == var])
      }

      if(var == "aspect") {
        aspect_sim <- circular(data_sim,
                               type = "angles", units = "degrees",
                               template = "geographic")

        bw_circ <- 20 # smaller = smoother

        den_sim <- density.circular(aspect_sim, bw = bw_circ)
      }

      # approximate density to the same sequence as reference density
      xseq <- dens_ref$value
      diff_size <- abs(unique(diff(xseq))[1])

      d_sim <- approx(den_sim$x, den_sim$y, xseq, method = "linear",
                      yleft = 0, yright = 0)

      # circular densities are not normalized, so we get the normalizing factor.
      # (And anyway, this ensures all approximate densities are exactly
      # normalize.)
      height_sim <- rowMeans(cbind(d_sim$y[-1], d_sim$y[-length(xseq)]))
      height_av <- rowMeans(cbind(dens_ref$density[-1], dens_ref$density[-length(xseq)]))

      auc_sim <- sum(height_sim * diff_size)
      auc_av <- sum(height_av * diff_size)

      norm_factor <- 1 / mean(auc_sim, auc_av)

      den_diff <- height_sim - height_av
      delta <-  ((abs(den_diff) * diff_size) %>% sum) / 2 * norm_factor
      overlap <- 1 - delta

      ## check overlap makes sense:
      # plot(d_sim$y ~ xseq, type = "l",
      #      ylim = c(0, max(c(height_sim, height_av)) * 1.05))
      # lines(dens_ref$density ~ dens_ref$value, col = 2)

      # fill overlap matrix
      fff <- overlap_data_raw$variable == var &
             overlap_data_raw$vegetation_class == veg
      ov_matrix[fff, i] <- overlap
    }
  }
}

# saveRDS(ov_matrix, "overlap_matrix_randomizations.rds")
naes <- apply(ov_matrix, 1, function(x) sum(is.na(x)))
aggregate(naes ~ vegetation_class, overlap_data_raw, range)
# OK, dry forests have enough data (na = 16)

# now, the p-values
overlap_data_raw$p_value <- sapply(1:nrow(ov_matrix), function(r) {
  x <- na.omit(ov_matrix[r, ])
  if(length(x) > 100) {
    p <- sum(x <= overlap_data_raw$overlap_num[r]) / length(x)
    return(p)
  } else return(NA)
})

# minimum p-value
1/400
1/(400-16)
# pero en los plots luego lo redondea a 0.002,
# así que pondré < 0.002

p3 <- round(overlap_data_raw$p_value, 3)
p3[p3 == 0] <- "<0.002"

overlap_data_raw$p_char <- paste("(", p3, ")", sep = "")
overlap_data_raw$overlap_char <- paste(
  round(overlap_data_raw$overlap_num * 100, 2), "%"
)

overlap_data_raw$ov_p_char <- paste(
  overlap_data_raw$overlap_char, "\n",
  overlap_data_raw$p_char,
  sep = ""
)

# Details for plotting
# good variable names:
names_frame <- data.frame(variable = var_names,
                          var_name = c("(A) Elevation (m a.s.l.)",
                                       "(B) Slope (°)",
                                       "(C) Aspect",
                                       "(D) Topographic position",
                                       "(A) NDVI max",
                                       "(B) Precipitation\n(mm / year)",
                                       "(C) Distance from\nhuman settlements (km)",
                                       "(D) Distance from\nroads (km)"))

# rename variabes in ov data
overlap_data <- left_join(overlap_data_raw, names_frame, by = "variable")
densdata <- left_join(densdata0, names_frame, by = "variable")

# move overlap values to avoid overlapping the densities
overlap_data$value[overlap_data$variable == "ndvi"] <- 0.5
overlap_data$value[overlap_data$variable == "ndvi" &
                   overlap_data$vegetation_class == "Steppe and grassland"] <- 0.85
overlap_data$density[overlap_data$variable == "ndvi"] <- overlap_data$density[overlap_data$variable == "ndvi"] * 1.1
overlap_data$value[overlap_data$variable == "elevation"] <- 1950
overlap_data$density[overlap_data$variable == "elevation"] <- overlap_data$density[overlap_data$variable == "elevation"] * 1.1
overlap_data$density[overlap_data$variable == "aspect"] <-
  max(overlap_data$density[overlap_data$variable == "aspect"])
overlap_data$value[overlap_data$variable == "aspect"] <- -145
overlap_data$value[overlap_data$variable == "pp"] <- 1900

# put p below overlap:
overlap_data$density_p <- overlap_data$density * 0.85
# different for aspect
overlap_data$density_p[overlap_data$variable == "aspect"] <-
  overlap_data$density[overlap_data$variable == "aspect"]# * 1.15

# aspect needs the x-value for p different from overlap.
overlap_data$value_p <- overlap_data$value
overlap_data$value_p[overlap_data$variable == "aspect"] <- -145 - 35 - 35 # mirror around -180

# Make plots

# rename vegetations
densdata$vegetation_class <- factor(densdata$vegetation_class,
                                    levels = veg_levels_ord,
                                    labels = veg_levels3_num)
overlap_data$vegetation_class <- factor(overlap_data$vegetation_class,
                                        levels = veg_levels_ord,
                                        labels = veg_levels3_num)

plant_anth_id <- c(4, 7) # id for plantation and anthrop in veg_levels3

plist <- vector(mode = "list", length = n_var)
names(plist) <- var_names

for(i in 1:n_var) {
  # i = 3

  ## OJO, ACÁ DEFINIR SI EL PLOT VA CON O SIN PLANTATION Y ANTHROP

  dd <- densdata[densdata$var_name == names_frame$var_name[i] &
                 !(densdata$vegetation_class %in% veg_levels3_num[plant_anth_id]), ]
  dd$vegetation_class %>% levels
  ov_d <- overlap_data[overlap_data$var_name == names_frame$var_name[i] &
                       !(overlap_data$vegetation_class %in% veg_levels3_num[plant_anth_id]),]


  plist[[i]] <- ggplot(dd,
                       aes(x = value, y = density, ymin = 0, ymax = density,
                           colour = class, fill = class)) +
    geom_line(size = 0.3) +
    geom_ribbon(alpha = 0.2, colour = NA) +

    facet_grid(rows = vars(vegetation_class), cols = vars(var_name),
               scales = "free") +

    scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
    scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
    geom_text(data = ov_d, hjust = 0.5,
              mapping = aes(x = value, y = density, label = overlap_char),
              size = 2.1, inherit.aes = FALSE) +
    geom_text(data = ov_d, hjust = 0.5,
              mapping = aes(x = value_p, y = density_p, label = p_char),
              size = 2.1, inherit.aes = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),

          strip.background.y = element_blank(),
          strip.text.y = element_blank(),

          strip.background.x = element_rect(color = "white", fill = "white"),
          strip.text.x = element_text(size = 7, color = "black"),

          legend.title = element_blank(),
          legend.position = "none",

          axis.title.y = element_blank(),
          axis.title.x = element_blank(),

          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),

          axis.text.x = element_text(size = 7),

          plot.margin = unit(c(1, 1, 1, 1), "mm"))
  # print(plist[[i]])

  if(var_names[i] == "aspect") {
    # i = 5
    # i = 3
    plist[[i]] <-
      ggplot(dd,
                         aes(x = value, y = density, ymin = 0, ymax = density,
                             colour = class, fill = class)) +
      geom_line(size = 0.3) +
      geom_ribbon(alpha = 0.2, colour = NA) +

      facet_grid(rows = vars(vegetation_class), cols = vars(var_name)) +
      # facet_wrap(vars(vegetation_class)) +

      scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
      scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
      geom_text(data = ov_d, hjust = 0.5,
                mapping = aes(x = value, y = density, label = overlap_char),
                size = 2.1, inherit.aes = FALSE) +
      geom_text(data = ov_d, hjust = 0.5,
                mapping = aes(x = value_p, y = density_p, label = p_char),
                size = 2.1, inherit.aes = FALSE) +
      # scale_x_continuous(breaks = c(-270, -180, -90, 0, 90),
      #                    labels = c("E", "S", "W", "N", "E")) +
      scale_x_continuous(breaks = c(-270, -180, -90, 0),
                         labels = c("E", "S", "W", "N")) +
      scale_y_continuous(limits = c(0, max(dd$density) * 1)) +

      coord_polar(start = pi / 2) +

      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),

            strip.background.y = element_blank(),
            strip.text.y = element_blank(),

            strip.background.x = element_rect(color = "white", fill = "white"),
            strip.text.x = element_text(size = 7, color = "black"),

            legend.title = element_blank(),
            legend.position = "none",

            axis.title.y = element_blank(),
            axis.title.x = element_blank(),

            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),

            axis.text.x = element_text(size = 4.5),

            plot.margin = unit(c(1, 1, 1, 1), "mm"))
    # plist[[i]]

  }
}


# add veg label in rightmost plots (dist_roads and TPI)



## ATENTION this list must have length = length(veg_levels) - 2
texts <- list(a, b, b, b, b, b, b)
c <- element_rect(fill = "gray10", color = "gray10")
d <- element_rect(fill = "white", color = "white")
backgrounds <- list(c, d, d, d, d, d, d)

plist[["dist_roads_km"]] <- plist[["dist_roads_km"]] +
  theme(strip.background.y = element_rect(),
        strip.text.y = element_text()) +
  facet_grid2(rows = vars(vegetation_class), cols = vars(var_name),
              scales = "free",
              strip = strip_themed(text_y = texts,
                                   background_y = backgrounds))

plist[["TPI2k"]] <- plist[["TPI2k"]] +
  theme(strip.background.y = element_rect(),
        strip.text.y = element_text()) +
  facet_grid2(rows = vars(vegetation_class), cols = vars(var_name),
              scales = "free",
              strip = strip_themed(text_y = texts,
                                   background_y = backgrounds)) +
  scale_x_continuous(breaks = c(0, 0.5, 1))


# add y axis name at leftmost plot (elevation and ndvi)
plist[["elevation"]] <- plist[["elevation"]] +
  theme(axis.title.y = element_text(size = 10)) +
  ylab("Probability density")

plist[["ndvi"]] <- plist[["ndvi"]] +
  theme(axis.title.y = element_text(size = 10)) +
  ylab("Probability density")

# Edit NDVI labels
plist[["ndvi"]] <- plist[["ndvi"]] +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# Edit pp labels
plist[["pp"]] <- plist[["pp"]] +
  scale_x_continuous(breaks = c(0, 1000, 2000))

# Edit Elevation labels
plist[["elevation"]] <- plist[["elevation"]] +
  scale_x_continuous(breaks = c(0, 1000, 2000))


# plot vegetation dist to get legend
veg_dist3 <-
  ggplot(veg_dist_data,
         aes(x = vegetation_class, y = prob,
             colour = class, fill = class)) +
  geom_bar(stat = "identity",
           # position = "identity",
           position = position_dodge2(width = 2, padding = 0.05),
           alpha = 0.3,
           size = 0.3) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5) +
  #facet_wrap(vars(variable)) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", #c(0.85, 0.8),
        axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.67),
        axis.title.x = element_text(vjust = 7),
        axis.title.y = element_blank(),
        strip.text.x = element_text(margin = margin(3,0,3,0, "mm"))) +
  ylab("Probability mass") +
  xlab("Vegetation type") +
  geom_text(data = data_ol, inherit.aes = FALSE,
            mapping = aes(x = vegetation_class, y = prob, label = lab),
            colour = "black")

veg_dist3 <- veg_dist3 +
  theme(legend.position = "right",
        legend.key.size = unit(5, 'mm'),
        legend.spacing.y = unit(1, 'mm'),
        legend.text = element_text(size = 8,
                                   margin = margin(t = 1, b = 1,
                                                   r = 2,
                                                   unit = "mm"))) +
  guides(fill = guide_legend(byrow = TRUE), color = guide_legend(byrow = TRUE))

leg <- ggpubr::get_legend(veg_dist3)

# merge core plots
dens_topo_0 <- egg::ggarrange(plots = plist[c("elevation", "slope", "aspect", "TPI2k")],
                              ncol = n_var / 2)
dens_prod_0 <- egg::ggarrange(plots = plist[c("ndvi", "pp", "dist_humans_km", "dist_roads_km")],
                              ncol = n_var / 2)

# add legends
# dens_topo_1 <- grid.arrange(dens_topo_0, leg, nrow = 2, heights = c(20, 1))
# dens_prod_1 <- grid.arrange(dens_prod_0, leg, nrow = 2, heights = c(20, 1))

dens_topo_1 <- grid.arrange(dens_topo_0, leg, nrow = 1, widths = c(20, 3))
dens_prod_1 <- grid.arrange(dens_prod_0, leg, nrow = 1, widths = c(20, 3))


ggsave("figures/06) spatial patterns - topo by veg type.png",
       plot = dens_topo_1,
       width = 16, height = 15, units = "cm") # width was 17
ggsave("figures/07) spatial patterns - prod and dist by veg type.png",
       plot = dens_prod_1,
       width = 16, height = 15, units = "cm") # width was 17

# corregir posición de los números.

# 07) Burned proportion as a function of environmental variables across veg types --------


dd <- spatdata[spatdata$class == "Burnable" &
               spatdata$vegetation_class != "All vegetation types", ]

dd$northing <- cos(dd$aspect * (pi / 180))
nn <- c(names_frame$variable, "northing")
nn <- nn[nn != "aspect"]

# longanize
ddl <- pivot_longer(dd, which(names(dd) %in% nn),
                    values_to = "x", names_to = "var")
# summaries by veg
veg_data0 <- aggregate(x ~ var + vegetation_class, ddl, mean_ci)
veg_data0 <- do.call(data.frame, veg_data0)
colnames(veg_data0)[grep("x", colnames(veg_data0))] <- c("lower", "mean", "upper")

veg_data <- left_join(
  veg_data0,
  burned_veg[, c("vegetation_class", "burned_prop", "area_available_ha")],
  by = "vegetation_class"
)

veg_data$vegetation_class <- factor(veg_data$vegetation_class,
                                    levels = veg_levels)


nf2 <- names_frame
nf2[names_frame$variable == "aspect", c("variable", "var_name")] <- c("northing", "Northing")

veg_data$var_name <- factor(veg_data$var,
                            levels = nf2$variable,
                            labels = nf2$var_name)

dddd <- veg_data[!veg_data$vegetation_class %in% levels(veg_data$vegetation_class)[c(3, 6, 8)], ]

(veg_means <- ggplot(dddd, aes(x = mean, xmin = lower, xmax = upper,
                                   y = burned_prop,
                                   fill = vegetation_class)) +
  geom_linerange(data = dddd,
                 mapping = aes(x = mean, xmin = lower, xmax = upper,
                               y = burned_prop,
                               colour = vegetation_class)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_viridis(discrete = T, option = "B", end = 0.9) +
  scale_colour_viridis(discrete = T, option = "A", end = 0.9) +

  facet_wrap(vars(var_name), scales = "free_x",
             strip.position = "bottom",
             ncol = 2) +
  ylab("Burned proportion") +
  ylim(0, 0.15) +
  theme(legend.position = "right",
        # legend.position = "bottom",
        # legend.box = "horizontal",
        legend.spacing.y = unit(0.1, 'mm'),
        legend.text = element_text(margin = margin(l = -2, r = 2, unit = "mm")),
        axis.title.x = element_blank(),
        strip.text.x = element_text(color = "black", size = 11),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.placement = "outside") +
  # guides(fill = guide_legend(ncol = 4, nrow = 2, byrow = TRUE)) +
  scale_x_continuous(n.breaks = 4)
)

ggsave("figures/07) spatial patterns - means by veg type_c.png",
       plot = veg_means,
       width = 16, height = 16, units = "cm") # width was 17


(veg_means <- ggplot(veg_data_long, aes(x = x, y = burned_prop,
                                        fill = vegetation_class)) +
    geom_point(size = 3, shape = 21) +
    scale_fill_viridis(discrete = T, option = "B", end = 1) +
    facet_wrap(vars(var_name), scales = "free_x",
               strip.position = "bottom",
               ncol = 2) +
    ylab("Burned proportion") +
    ylim(0, 0.15) +
    theme(legend.position = "right",
          # legend.position = "bottom",
          # legend.box = "horizontal",
          legend.spacing.y = unit(0.1, 'mm'),
          legend.text = element_text(margin = margin(l = -2, r = 2, unit = "mm")),
          axis.title.x = element_blank(),
          strip.text.x = element_text(color = "black", size = 11),
          strip.background = element_rect(color = "white", fill = "white"),
          strip.placement = "outside") +
    # guides(fill = guide_legend(ncol = 4, nrow = 2, byrow = TRUE)) +
    scale_x_continuous(n.breaks = 4)
)

ggsave("figures/07) spatial patterns - means by veg type_b.png",
       plot = veg_means,
       width = 16, height = 16, units = "cm") # width was 17


# boxplots by veg
ddl <- pivot_longer(dd, which(names(dd) %in% nn),
                    values_to = "x", names_to = "var")
ddl$vegetation_class <- factor(ddl$vegetation_class, levels = veg_levels)
ddl$var_name <- factor(ddl$var, levels = nf2$variable, labels = nf2$var_name)

(veg_bb <- ggplot(ddl, aes(x = vegetation_class, y = x,
                                     fill = vegetation_class)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = T, option = "B", end = 1) +
    facet_wrap(vars(var_name), scales = "free_y",
               strip.position = "top",
               ncol = 2) +
    theme(legend.position = "right",
          # legend.position = "bottom",
          # legend.box = "horizontal",
          legend.spacing.y = unit(0.1, 'mm'),
          legend.text = element_text(margin = margin(l = -2, r = 2, unit = "mm")),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.x = element_text(color = "black", size = 11),
          strip.background = element_rect(color = "white", fill = "white"),
          strip.placement = "outside") #+
    # guides(fill = guide_legend(ncol = 4, nrow = 2, byrow = TRUE)) +
    # scale_x_continuous(n.breaks = 4)
)

ggsave("figures/07) spatial patterns - boxplot by veg type_b.png",
       plot = veg_bb,
       width = 16, height = 16, units = "cm") # width was 17


# Interannual temporal patterns -------------------------------------------


## Prepare data

burned_annual <- read.csv("data_and_files/data_burned_area_by_year.csv")
burned_veg <- read.csv("data_and_files/data_burned_and_available_area_by_vegetation_dryforest2.csv")

# bring fwi data
fwi_data <- read.csv("data_and_files/data_climate_interannual_fwi.csv")
fwi_data$date <- as.Date(fwi_data$date, format = "%Y-%m-%d")
fwi_data$year <- format(fwi_data$date, format = "%Y") %>% as.numeric
fwi_data$month <- format(fwi_data$date, format = "%m") %>% as.numeric
fwi_data$fseason <- NA
fwi_data$fseason <- fwi_data$year
fwi_data$fseason[fwi_data$month == 12] <- fwi_data$year[fwi_data$month == 12] + 1

fwi_agg <- aggregate(fwi ~ fseason, fwi_data[fwi_data$month %in% c(12, 1:3), ],
                     FUN = mean)
names(fwi_agg) <- c("year", "fwi")

# remaining climatic variables
climate_long <- read.csv("data_and_files/data_climate_interannual.csv")
climate <- pivot_wider(climate_long[, c("variable", "value", "year")],
                       names_from = "variable",
                       values_from = "value")

# merge with fwi
climate <- cbind(climate[climate$year > 1998, ], fwi = fwi_agg$fwi)

# merge with fire data
burned_annual_clim <- left_join(
  burned_annual, climate, by = "year"
)
burned_annual_clim_long <- pivot_longer(
  burned_annual_clim, which(names(burned_annual_clim) %in% c("pp", "temp",
                                                             "vpd", "wind", "fwi")),
  names_to = "clim_var", values_to = "clim_value"
)

# better names for variables
burned_annual_clim_long$clim_var2 <- plyr::revalue(
  burned_annual_clim_long$clim_var,
  replace = c(
    "fwi" = "Fire Weather\nIndex",
    "pp" = "Precipitation (mm)",
    "temp" = "Temperature (°C)",
    "vpd" = "Vapour pressure\ndeficit (kPa)"
  )
)
burned_annual_clim_long$clim_var2 <- factor(
  burned_annual_clim_long$clim_var2,
  levels = c(
    "Temperature (°C)",
    "Precipitation (mm)",
    "Vapour pressure\ndeficit (kPa)",
    "Fire Weather\nIndex",
    "wind"
  )
)

# turn area into proportion
burned_annual_clim_long$area_perc <- burned_annual_clim_long$area_ha / 2194363 * 100
burned_annual_clim$area_perc <- burned_annual_clim$area_ha / 2194363 * 100

burned_annual_clim_long$area_perc_log <- log10(burned_annual_clim_long$area_perc)
burned_annual_clim$area_perc_log <- log10(burned_annual_clim$area_perc)


## Burned proportion as a function of climate  (monotonic GAMs)

# get coordinates in x for r2
ranges <- apply(burned_annual_clim[, c("pp", "temp", "vpd", "wind", "fwi")], 2, range)
widths <- apply(ranges, 2, diff)

# r2 for burned area
r2_area <- data.frame(clim_var2 = unique(burned_annual_clim_long$clim_var2),
                      clim_var = unique(burned_annual_clim_long$clim_var),
                      r2 = NA,
                      r2ori = NA,
                      clim_value = ranges[1, ] + widths / 2,
                      y = 3.6)#max(burned_annual_clim_long$area_perc) * 1.1)

for(v in 1:nrow(r2_area)) {
  # v = 1
  var_name <- r2_area$clim_var[v]
  d <- burned_annual_clim[, c("area_perc_log", "area_perc", var_name)]
  names(d)[ncol(d)] <- "clim_var"
  basis <- ifelse(var_name == "pp", "mpd", "mpi")
  m <- scam(area_perc_log ~ s(clim_var, bs = basis, k = 4), data = d) # y ~ s(x, bs = "cs", k = 4)
  mu <- fitted(m)
  var_y <- sigma(m) ^ 2
  rr <- round(var(mu) / (var(mu) + var_y) * 100, 2) %>% format(nsmall = 2)
  r2_area$r2[v] <- paste(rr, "%")

  # r2 at the original scale
  m_ori <- scam(area_perc ~ s(clim_var, bs = basis, k = 4), data = d,
                family = Gamma(link = "log")) # y ~ s(x, bs = "cs", k = 4)
  mu_ori <- fitted(m_ori)
  var_y_ori <- m_ori$family$variance(mu_ori) %>% mean # just mu ^ 2,https://pj.freefaculty.org/guides/stat/Regression-GLM/Gamma/GammaGLM-01.pdf
  rr_ori <- round(var(mu_ori) / (var(mu_ori) + var_y_ori) * 100, 2) %>% format(nsmall = 2)
  r2_area$r2ori[v] <- paste(rr_ori, "%")

  # check autocorr
  r <- simulateResiduals(m)
  plot(acf(r$scaledResiduals), main = v)
}

# clim vars in columns:
clim_area <-
  ggplot(burned_annual_clim_long[burned_annual_clim_long$clim_var != "wind", ],
         aes(x = clim_value, y = area_perc)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              # method.args = list(family = Gamma(link = "log")),
              color = viridis(1, begin = 0.2, option = "B"),
              fill = viridis(1, begin = 0.2, option = "B"),
              alpha = 0.2, linewidth = 0.8) +
  geom_point(shape = 19, alpha = 0.7, size = 2.5) +
  geom_text(aes(x = clim_value, y = y, label = r2), size = 3.0,
            data = r2_area[r2_area$clim_var != "wind", ]) +
  facet_wrap(vars(clim_var2), scales = "free_x", nrow = 1) +
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 11)) +
  xlab("Climatic variable") +
  ylab("Burned proportion (%)") +
  scale_x_continuous(n.breaks = 4) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c(0.001, 0.01, 0.1, 1)) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
  #               labels = trans_format("log10", math_format(10 ^ .x))) +
  ggtitle("B")
clim_area


## Number of fires as a function of climate (monotonic GAMs)

# models and r2 for number of fires. As scam() does not estimate theta for
# the neg bin, we used negbin(). This function is to be used with known theta.
# So we first take an initial value from a gam(nb()) fit and optimize the
# likelihood.
theta_optimizer <- function(th, data, basis) {
  mod <- scam(fires ~ s(clim_var, bs = basis, k = 4), data = data,
              family = negbin(theta = th))
  return(-as.numeric(logLik.scam(mod)))
}

r2_fires <- r2_area
r2_fires$y <- 38

# save predictions from scams in a list
preds_list <- vector("list", nrow(r2_fires))

for(v in 1:nrow(r2_fires)) {
  # v = 5
  var_name <- r2_area$clim_var[v]
  print(var_name)
  d <- burned_annual_clim[, c("fires", var_name)]
  names(d)[ncol(d)] <- "clim_var"
  basis <- ifelse(var_name == "pp", "mpd", "mpi")

  # unconstrained gam to get initial guess for theta.
  m0 <- gam(fires ~ s(clim_var, bs = "cr", k = 4), data = d,
            family = nb(link = "log"))
  theta_init <- m0$family$getTheta(TRUE)

  # optimize theta fitting a scam
  opthet <- optim(theta_init, theta_optimizer, data = d, basis = basis,
                  method = "Brent",
                  lower = theta_init * 0.1,
                  upper = min(theta_init * 10, 30))

  # fit final scam with optimal theta.
  m <- scam(fires ~ s(clim_var, bs = basis, k = 4), data = d,
            family = negbin(theta = opthet$par))

  mu <- fitted(m)
  theta <- opthet$par
  var_y <- mean(mu + mu ^ 2 / theta) # neg bin variance
  rr <- round(var(mu) / (var(mu) + var_y) * 100, 2) %>% format(nsmall = 2)
  r2_fires$r2[v] <- paste(rr, "%")

  # compute predictions
  nd <- data.frame(clim_var = seq(min(d$clim_var), max(d$clim_var),
                                  length.out = 150))
  pp <- predict(m, nd, se.fit = T)
  nd$mu <- exp(pp$fit)
  nd$mu_lower <- exp(pp$fit - qnorm(0.975) * pp$se.fit)
  nd$mu_upper <- exp(pp$fit + qnorm(0.975) * pp$se.fit)
  nd$clim_value <- nd$clim_var
  nd$clim_var <- var_name
  nd$clim_var2 <- r2_fires$clim_var2[v]

  preds_list[[v]] <- nd

  # check autocorr (make residuals by hand)
  ysim <- matrix(rnbinom(nrow(d) * 3000, size = theta, mu = mu),
                 nrow(d), 3000)
  r <- createDHARMa(ysim, d$fires, mu)
  plot(r)
  plot(acf(r$scaledResiduals))
}

preds_fires <- do.call("rbind", preds_list)
r2_fires$y <- 41

clim_fires <-
  ggplot(burned_annual_clim_long[burned_annual_clim_long$clim_var != "wind", ],
         aes(x = clim_value, y = fires)) +
  # geom_smooth(
  #             # method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
  #             # method.args = list(family = mgcv::nb(link = "log")),
  #             # method = scam::scam, formula = y ~ s(x, bs = "mpi", k = 4),
  #             # method.args = list(family = "poisson"),#mgcv::negbin(theta = c(-10), link = "log")),
  #             color = viridis(1, begin = 0.2, option = "B"),
  #             fill = viridis(1, begin = 0.2, option = "B"),
  #             alpha = 0.2, linewidth = 0.8) +
  # smooth by hand, with scam
  geom_ribbon(data = preds_fires[preds_fires$clim_var != "wind", ],
              mapping = aes(x = clim_value, ymin = mu_lower, ymax = mu_upper),
              fill = viridis(1, begin = 0.2, option = "B"),
              alpha = 0.2, color = NA, inherit.aes = F) +
  geom_line(data = preds_fires[preds_fires$clim_var != "wind", ],
            mapping = aes(x = clim_value, y = mu),
            color = viridis(1, begin = 0.2, option = "B"),
            inherit.aes = F) +

  geom_point(shape = 19, alpha = 0.7, size = 2.5) +
  geom_text(aes(x = clim_value, y = y, label = r2), size = 3.0,
            data = r2_fires[r2_fires$clim_var != "wind", ]) +
  facet_wrap(vars(clim_var2), scales = "free_x", nrow = 1,
             strip.position = "bottom") +
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(color = "black", size = 11),
        strip.background = element_rect(color = "white", fill = "white"),
        axis.title.y = element_text(hjust = 0.5, vjust = 5.2, size = 11),
        strip.placement = "outside",
        axis.text.x = element_text(size = 9),
        axis.title.x = element_blank()) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(breaks = seq(0, 30, 10)) +
  ylab("Number of fires")
clim_fires

## Burned area time series

burned_annual$area_perc <- burned_annual$area_ha / 2194363 * 100

mintemp <- min(climate$temp)
scaletemp <- max(climate$temp - mintemp)
climate2 <- climate
lowertemp <- 0.5
climate2$temp <- (climate$temp - mintemp) / scaletemp + lowertemp
# scale temperature between 1 and 2

# y = (x - mintemp) / scaletemp + 1
# y - 1 = (x - mintemp) / scaletemp
# (y - 1) * scaletemp = (x - mintemp)
# (y - 1) * scaletemp + mintemp = x -

# make temp scale vary between 18 and 22 when perc varies from 0.5 to 2
data_scales <- data.frame(perc = c(0.5, 2),
                          temp = c(18, 21))
cc <- coef(lm(perc ~ temp, data = data_scales))
climate2 <- climate
climate2$temp <- cc[1] + climate$temp * cc[2]
# invert:
# perc = cc[1] + temp * cc[2]
# temp = (perc - cc[1]) / cc[2]

barplot(c(1, 1) ~ c(1, 2), col = viridis(2, option = "A", end = 0.3))
temp_color <- viridis(2, option = "A", end = 0.4)[2]

ts_area_temp <-
  ggplot(burned_annual, aes(x = year, y = area_perc)) +

  # Burned area and number of fires
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_smooth(mapping = aes(x = year, y = area_perc),
              linetype = "dashed", method = "glm",
              method.args = list(family = Gamma(link = "log")),
              se = FALSE, linewidth = 0.3,
              color = "black") +

  # temp
  geom_line(data = climate2, mapping = aes(x = year, y = temp),
            color = temp_color, linewidth = 0.45) +
  geom_smooth(data = climate2, mapping = aes(x = year, y = temp),
              color = temp_color, linetype = "dotted", method = "lm",
              se = FALSE, linewidth = 0.5) +
  geom_point(data = climate2, mapping = aes(x = year, y = temp),
             fill = temp_color, size = 2.2, shape = 21) +
  xlab("Year") +
  ggtitle("A") +
  scale_y_continuous(
    name = "Burned proportion (%)",
    sec.axis = sec_axis(~ (. - cc[1]) / cc[2],
                        name = "Temperature (°C)",
                        breaks = seq(18, 21, 1))
  ) +
  theme(plot.title = element_text(hjust = 0),
        axis.text.y.right = element_text(color = temp_color),
        axis.title.y.right = element_text(color = temp_color),
        axis.ticks.y.right = element_line(color = temp_color),
        axis.title = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
ts_area_temp

## number of fires time series (with FWI)

minfwi <- min(climate$fwi)
scalefwi <- max(climate$fwi - minfwi)
climate2 <- climate
climate2$fwi <- (climate$fwi - minfwi) / scalefwi + 1

fwi_color <- viridis(2, option = "A", end = 0.4)[2]

ts_fires_fwi <-
  ggplot(burned_annual, aes(x = year, y = fires)) +

  # Burned area and number of fires
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_smooth(mapping = aes(x = year, y = fires),
              linetype = "dashed", method = "gam",
              formula = y ~ x,
              method.args = list(family = nb(link = "log")),
              se = FALSE, linewidth = 0.3,
              color = "black") +

  # fwi
  geom_line(data = climate, mapping = aes(x = year, y = fwi),
            color = fwi_color, linewidth = 0.45) +
  geom_smooth(data = climate, mapping = aes(x = year, y = fwi),
              color = fwi_color, linetype = "dotted", method = "lm",
              se = FALSE, linewidth = 0.5) +
  geom_point(data = climate, mapping = aes(x = year, y = fwi),
             fill = fwi_color, size = 2.2, shape = 21) +
  xlab("Year") +
  scale_y_continuous(
    name = "Number of fires",
    sec.axis = sec_axis(~ .,
                        name = "Fire Weather Index",
                        breaks = c(10, 20))#seq(5, 20, 5))
  ) +
  theme(plot.title = element_text(hjust = 0),
        axis.text.y.right = element_text(color = fwi_color),
        axis.title.y.right = element_text(color = fwi_color),
        axis.ticks.y.right = element_line(color = fwi_color),
        axis.title = element_text(size = 11),
        axis.title.y = element_text(vjust = 1.7))
ts_fires_fwi


# figure with all climatic results
clim_fig <- ggarrange(ts_area_temp, ts_fires_fwi, clim_area, clim_fires,
                      ncol = 1)
ggsave("figures/03) temporal patterns_monotonic.jpeg",
       clim_fig,
       width = 16, height = 20, units = "cm")

# Data for text

av_area <- burned_veg$area_available_ha[-1] %>% sum

(burned_annual$area_ha %>% mean) / av_area
(burned_annual$area_ha %>% mean)
(burned_annual$area_ha %>% summary)

burned_annual$area_ha[order(burned_annual$area_ha)]

prop_ord <- burned_annual$area_ha[order(burned_annual$area_ha)] / av_area * 100
plot(ecdf(prop_ord))
min(prop_ord)
max(prop_ord)

15000 / av_area

burned_annual$fires %>% mean
burned_annual$fires %>% range

# Temporal trends tests
acf(burned_annual$area_ha)
acf(burned_annual$fires)
acf(climate$fwi)
acf(climate$pp)
acf(climate$temp)
acf(climate$vpd)


# no temporal correlation, use raw models
m_perc <- glm(area_perc ~ year, data = burned_annual, family = Gamma(link = "log"))
m_fires <- MASS::glm.nb(fires ~ year, data = burned_annual)
m_fwi <- lm(fwi ~ year, data = climate)
m_pp <- lm(pp ~ year, data = climate)
m_temp <- lm(temp ~ year, data = climate)
m_vpd <- lm(vpd ~ year, data = climate)

# check dharmas
mlist <- list(m_area, m_fires, m_fwi, m_pp, m_temp, m_vpd)
lapply(mlist, function(m) {
  r <- simulateResiduals(m)
  plot(r)
  plot(acf(r$scaledResiduals))
}) # all ok

coef_table <- rbind(
  summary(m_perc)$coefficients["year", ],
  summary(m_fires)$coefficients["year", ],
  summary(m_fwi)$coefficients["year", ],
  summary(m_pp)$coefficients["year", ],
  summary(m_temp)$coefficients["year", ],
  summary(m_vpd)$coefficients["year", ]
)

# compute bayesian r2 for all models
r2trends <- numeric(6)

# number of fires (neg bin)
mu <- fitted(m_fires)
theta <- m_fires$theta
var_y <- mean(mu + mu ^ 2 / theta) # neg bin variance
rr <- round(var(mu) / (var(mu) + var_y) * 100, 3) %>% format(nsmall = 3)
r2trends[2] <- rr

# burned proportion (Gamma)
mu <- fitted(m_perc)
var_y <- m_perc$family$variance(mu) %>% mean
rr <- round(var(mu) / (var(mu) + var_y) * 100, 3) %>% format(nsmall = 3)
r2trends[1] <- rr

# r2 for normal models (climate)
r2norm <- function(model) {
  mu <- fitted(model)
  sigma <- sigma(model)
  rr <- var(mu) / (var(mu) + sigma ^ 2)
  return(rr)
}
r2trends[3:6] <- do.call("c", lapply(list(m_fwi, m_pp, m_temp, m_vpd), function(x) {
  rr <- r2norm(x)
  rr <- round(rr * 100, 3) %>% format(nsmall = 3)
  return(rr)
}))


coef_export <- cbind(
  variable = c("Burned percentage", "Number of fires",
    "Fire weather index", "Precipitation", "Temperature", "Vapour pressure deficit"),
  as.data.frame(format(round(coef_table, digits = 3), nsmall = 3)),
  r2 = r2trends
)

# rownames(coef_export) <- coef_export$variable
# coef_export[c("Burned area", "Number of fires"), ]$Estimate <- round(exp(coef_export[c("Burned area", "Number of fires"), ]$Estimate), 4)
write.csv(coef_export, "trends.csv", row.names = F)


# Theil-Sen regression (Thomas pide)
library(RobustLinearReg)
theil_sen_regression()

burned_annual$area_log <- log(burned_annual$area_ha)
mts_area <- theil_sen_regression(area_log ~ year, data = burned_annual)
mts_fwi <- theil_sen_regression(fwi ~ year, data = climate)
mts_pp <- theil_sen_regression(pp ~ year, data = climate)
mts_temp <- theil_sen_regression(temp ~ year, data = climate)
mts_vpd <- theil_sen_regression(vpd ~ year, data = climate)

coef_table_ts <- rbind(
  summary(mts_area)$coefficients["year", ],
  summary(mts_fwi)$coefficients["year", ],
  summary(mts_pp)$coefficients["year", ],
  summary(mts_temp)$coefficients["year", ],
  summary(mts_vpd)$coefficients["year", ]
)

# da muy parecido.

# Fig 3: ------------------------------------------------------------------
# Intraannual fire activity -----------------------------------------------

# get clipped fires database

# fires_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/patagonian_fires/patagonian_fires.shp")
# study_area_wgs <- readOGR("/home/ivan/Insync/Burned area mapping/patagonian_fires/study_area/study_area.shp")
# proj_posgar <- "+proj=tmerc +lat_0=-90 +lon_0=-72 +k=1 +x_0=1500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
# fires_posgar <- spTransform(fires_wgs, proj_posgar)
# sa_posgar <- spTransform(study_area_wgs, proj_posgar)
# fires_clipped <- do.call("rbind", lapply(1:nrow(fires_posgar), function(i) {
#   print(i)
#   fire_clip <- gIntersection(polygons(fires_posgar[i, ]), polygons(sa_posgar))
#   row.names(fire_clip) <- row.names(fires_posgar[i, ])
#   fire_df <- SpatialPolygonsDataFrame(fire_clip,
#                                       data = fires_posgar@data[i, , drop = F])
#   # recompute area
#   fire_df$area_ha <- gArea(polygons(fire_clip)) * 0.0001 # m2 to ha
#
#   return(fire_df)
# }))
# fires_clipped_data <- fires_clipped@data
# write.csv(fires_clipped_data, "data_fires_clipped.csv")

fires_clipped_data <- read.csv("data_and_files/data_fires_clipped.csv")[, -1]
remove <- which(fires_clipped_data$obs == "uncertain_year")
dmonth <- fires_clipped_data[-remove, ]
dmonth$month_num <- format(as.Date(dmonth$date, format = "%Y-%m-%d"),
                           format = "%m") %>% as.numeric
dmonth$month <- factor(dmonth$month_num, levels = c(7:12, 1:6),
                       labels = month.name[c(7:12, 1:6)])
nrow(dmonth)
#### Fires count: fill with zeroes and then analize and aggregate.

# aggregate length by month and year
data_my <- aggregate(area_ha ~ month_num + month + year, dmonth, length)
colnames(data_my)[4] <- "fires"
data_my$id <- paste(data_my$month_num, data_my$year, sep = "_")

# add zeroes
ym <- expand.grid(month_num = 1:12,
                  year = 1999:2022)
ym$month <- factor(ym$month_num,
                   levels = c(7:12, 1:6))
ym$id <- paste(ym$month_num, ym$year, sep = "_")


data_my <- left_join(ym,
                     data_my[, c("id", "fires")],
                     by = "id")
data_my$fires[is.na(data_my$fires)] <- 0

# model
m_count <- gam(fires ~ s(month_num, bs = "cc",
                         k = length(unique(data_my$month_num))),
               data = data_my, family = nb(link = "log"),
               method = "REML", gamma = 1,
               knots = list(month_num = c(0, 12)))

# aggregate and predict by month
count_data <- aggregate(fires ~ month_num + month, data_my, mean)

# merge datasets by date
start <- as.Date("2015-07-01"); end <- as.Date("2016-06-01")
date_seq <- seq(start, end, 1)
date_seq_focal <- date_seq[day(date_seq) == 1]#date_seq[day(date_seq) == 15]
count_data$date <- date_seq_focal
# focal has only 12 dates; the other is a dayly sequence for GAM

# make continuous value for month
month_pred <- data.frame(month_num = c(seq(7, 18,
                                           length.out = length(date_seq))))
pp <- predict(m_count, month_pred, se.fit = T)
month_pred$count_mle <- exp(pp$fit)
month_pred$count_lower <- exp(pp$fit - qnorm(0.975) * pp$se.fit)
month_pred$count_upper <- exp(pp$fit + qnorm(0.975) * pp$se.fit)
month_pred$date <- date_seq

ggplot(month_pred, aes(x = date, y = count_mle,
                       ymin = count_lower, ymax = count_upper)) +
  geom_ribbon(alpha = 0.2, col = NA) + geom_line() +
  geom_point(data = count_data, mapping = aes(x = date, y = fires),
             inherit.aes = F) +
  scale_x_date(labels = date_format("%b"),
               breaks = "1 month",
               minor_breaks = "1 month")


#### Fires size

# aggregate size by month
size_data <- aggregate(area_ha ~ month_num + month, dmonth, mean)
colnames(size_data)[3] <- "size"
size_data$date <- date_seq_focal

# model fit to raw data
m_size <- gam(area_ha ~ s(month_num, bs = "cc",
                         k = length(unique(dmonth$month_num))),
               data = dmonth, family = Gamma(link = "log"),
               method = "REML", gamma = 1,
               knots = list(month_num = c(0, 12)))

# predict over the same values as for count
pp <- predict(m_size, month_pred, se.fit = T)
month_pred$size_mle <- exp(pp$fit)
month_pred$size_lower <- exp(pp$fit - qnorm(0.975) * pp$se.fit)
month_pred$size_upper <- exp(pp$fit + qnorm(0.975) * pp$se.fit)

ggplot(month_pred, aes(x = date, y = size_mle,
                       ymin = size_lower, ymax = size_upper)) +
  geom_ribbon(alpha = 0.2, col = NA) + geom_line() +
  geom_point(data = size_data, mapping = aes(x = date, y = size),
             inherit.aes = F) +
  scale_x_date(labels = date_format("%b"),
               breaks = "1 month",
               minor_breaks = "1 month")


# Merge both datasets and predictions, for count and size.
# Count will be scaled, so that they can live in the same space.
q <- max(month_pred$size_upper) / max(month_pred$count_upper)

month_pred_scaled <- month_pred
count_data_scaled <- count_data

month_pred_scaled$count_mle <- month_pred$count_mle * q
month_pred_scaled$count_lower <- month_pred$count_lower * q
month_pred_scaled$count_upper <- month_pred$count_upper * q

count_data_scaled$count <- count_data$fires * q

# longanize
ddwide <- cbind(count_data_scaled, size = size_data$size)
ddlong <- pivot_longer(ddwide, which(names(ddwide) %in% c("count", "size")),
                       names_to = "var", values_to = "y")

# preds
nn <- c("date", "mle", "lower", "upper", "var")
ss <- month_pred[, c("date", "size_mle", "size_lower", "size_upper")]
ss$var <- "size"
cc <- month_pred_scaled[, c("date", "count_mle", "count_lower", "count_upper")]
cc$var <- "count"
names(ss) <- names(cc) <- nn
pplong <- rbind(ss, cc)

ddlong$var <- factor(ddlong$var, levels = c("size", "count"),
                     labels = c("size", "number of fires"))
pplong$var <- factor(pplong$var, levels = c("size", "count"),
                     labels = c("size", "number of fires"))


intra_fire <-
  ggplot() +
  geom_ribbon(data = pplong, mapping = aes(x = date, y = mle,
                                           ymin = lower, ymax = upper,
                                           fill = var),
              alpha = 0.1, color = NA) +
  geom_line(data = pplong, mapping = aes(x = date, y = mle,
                                         color = var),
            alpha = 0.5) +

  geom_point(data = ddlong,
             mapping = aes(x = date, y = y, colour = var, shape = var,
                           group = var), size = 3, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "B", end = 0.5, direction = -1) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60, vjust = 0.5)) +
  scale_y_continuous(
    name = "Mean fire size (ha)",
    sec.axis = sec_axis(~ . / q, name = "Mean number of fires")
  )  +
  scale_x_date(labels = date_format("%b"),
               breaks = "1 month",
               minor_breaks = "1 month")+
  xlab("Month")
intra_fire

# climate variables

date_table <- data.frame(
  date = date_seq_focal,
  month = month(date_seq_focal)
)

# mean over the study area
clim_intra <- read.csv("data_and_files/data_climate_monthly.csv")
clim_intra <- left_join(
  clim_intra,
  date_table,
  by = "month"
)
clim_intra$variable <- plyr::revalue(
  clim_intra$variable,
  c("prec" = "pp", "tavg" = "temp")
)

# percentiles data
clim_intra_perc <- read.csv("data_and_files/data_climate_monthly_percentiles.csv")
clim_intra_perc <- separate(clim_intra_perc, "var",
                            into = c("month", "variable", "percentile"),
                            sep = "_")
clim_intra_perc$variable <- plyr::revalue(
  clim_intra_perc$variable,
  c("prec" = "pp", "tavg" = "temp")
)
clim_intra_perc$month <- as.numeric(clim_intra_perc$month)
clim_intra_perc <- left_join(clim_intra_perc, date_table, by = "month")


# scale for plot with 2 y-axes
coef_clim <- 2

clim_intra$y_scaled <- clim_intra$value
clim_intra$y_scaled[clim_intra$variable == "temp"] <-
  clim_intra$value[clim_intra$variable == "temp"] * coef_clim

clim_intra_perc$y_scaled <- clim_intra_perc$value
clim_intra_perc$y_scaled[clim_intra_perc$variable == "temp"] <-
  clim_intra_perc$value[clim_intra_perc$variable == "temp"] * coef_clim

# widenize by hand
clim_intra_perc <- clim_intra_perc[order(clim_intra_perc$variable,
                                         clim_intra_perc$date), ]
ids <- !duplicated(clim_intra_perc[, c("variable", "date")])

cperc_wide <- clim_intra_perc[ids, c("variable", "date")]
cperc_wide$lower <- clim_intra_perc[clim_intra_perc$percentile == "p5",
                                    "y_scaled"]
cperc_wide$upper <- clim_intra_perc[clim_intra_perc$percentile == "p95",
                                    "y_scaled"]

maxtemp <- max(clim_intra_perc$value[clim_intra_perc$variable == "temp"])

intra_clim <-
  ggplot() +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = variable),
              data = cperc_wide, color = NA, alpha = 0.3) +
  geom_line(data = clim_intra,
            mapping = aes(x = date, y = y_scaled,
                          colour = variable,
                          group = variable)) +
  geom_point(data = clim_intra,
             mapping = aes(x = date, y = y_scaled,
                           colour = variable, shape = variable,
                           group = variable),
             size = 3, alpha = 0.7) +
  scale_color_viridis(discrete = TRUE, option = "D", end = 0.5, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "D", end = 0.5, direction = -1) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60, vjust = 0.5),
        plot.margin = unit(c(2, 2, 2, 4), "mm")) +
  xlab("Month") +
  scale_y_continuous(
    name = "Precipitation (mm)",
    sec.axis = sec_axis(~ . / coef_clim, name = "Temperature (°C)",
                        breaks = seq(0, 15, by = 15))
  ) +
  scale_x_date(labels = date_format("%b"),
               breaks = "1 month",
               minor_breaks = "1 month")
intra_clim

# Fire size distribution

fires_wgs <- readOGR("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
proj_posgar <- "+proj=tmerc +lat_0=-90 +lon_0=-72 +k=1 +x_0=1500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
fires_posgar <- spTransform(fires_wgs, proj_posgar)

size_data <- data.frame("area_abs" = fires_posgar$area_ha,
                        "year" = fires_posgar$year)

size_data <- size_data[order(size_data$area_abs, decreasing = TRUE), , drop = F]
size_data$area_prop <- cumsum(size_data$area_abs) / sum(size_data$area_abs)
size_data$number_abs <- 1:nrow(size_data)
size_data$number_prop <- size_data$number_abs / max(size_data$number_abs)

size_props <-
  ggplot(size_data, aes(y = area_prop, x = number_prop)) +
  geom_point(
    color = "black",#viridis(1, begin = 0.2, option = "B"),
    size = 2, alpha = 0.6
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  ylab("Cumulative proportion\nof burned area") +
  xlab("Cumulative proportion of fires")
size_props

# in which years occurred the fires burning 90 % of the burned?
years_fires10p <- size_data$year[size_data$number_prop < 0.11] %>% unique
years_fires10p[order(years_fires10p)]
years_fires_half_area <- size_data$year[size_data$area_prop <= 0.5] %>% unique
years_fires_half_area



# log-log plot

# log-size range:
# 1.004864 4.432912
# 10 regular classes in width, not in N, from 1 to 4.5 at log scale.
# every class represented by its mean at the log scale.

size_data$area_abs_log10 <- log10(size_data$area_abs)
k <- 10
size_limits <- seq(1, 4.5, length.out = k+1)

size_data$size_class <- NA

for(i in 1:nrow(size_data)) {
  #i = 3
  cl <- as.numeric(size_data$area_abs_log10[i] > size_limits) %>% sum
  size_data$size_class[i] <- cl
}

# aggregate
mean_length <- function(x) c("mean" = mean(x), "length" = length(x))
size_data_agg <- do.call("data.frame",
                         aggregate(cbind(area_abs_log10, area_abs) ~
                                     size_class,
                                   size_data, mean_length))

size_data_agg$freq <- size_data_agg$area_abs.length / sum(size_data_agg$area_abs.length)

size_freq <-
  ggplot(size_data_agg, aes(x = area_abs.mean, y = freq)) +
  # geom_smooth(method = "glm", formula = y ~ log10(x),
  #             method.args	= list(family = Gamma(link = "log"))) +
  geom_smooth(method = "lm", se = TRUE,
              color = viridis(1, begin = 0.2, option = "B"),
              fill = viridis(1, begin = 0.2, option = "B"),
              alpha = 0.2, size = 0.8) +
  geom_point(size = 2.5, alpha = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  ylab("Relative frequency") +
  xlab("Size class (ha)")


# fire regime plots

size_season <-
ggarrange(size_props + ggtitle("A"),
          size_freq + ggtitle("B"),
          intra_fire + theme(legend.position = "bottom") +
            ggtitle("C") + xlab("Month"),
          intra_clim + theme(legend.position = "bottom") + ggtitle("D"),
          nrow = 2)

ggsave("figures/03) size distribution and seasonality.jpeg",
       size_season,
       width = 16, height = 14, units = "cm", dpi = 300)

# Burned area global results ----------------------------------------------

dburn_freq <- read.csv("data_and_files/data_reburn_frequency.csv")[, -1]
dburn_veg <- read.csv("data_and_files/data_burned_and_available_area_by_vegetation.csv")

# Total burned area
(burned_total_ha <- dburn_veg$area_burned_ha[-1] %>% sum)
# 126593.8

# Total burnable area
(burnable_ha <- dburn_veg$area_available_ha[-1] %>% sum)
# 2194363 ha

# percentage burned
round(burned_total_ha / burnable_ha * 100, 2)
# 5.77 %
24 / 0.0577
# FRI = 415.9445 = 416 años

# burnable / total
sum(dburn_veg$area_available_ha[-1]) / sum(dburn_veg$area_available_ha) * 100



# reburn:
options(scipen = 999)
dburn_freq
round(dburn_freq$prob * 100, digits = 2)

# add row for freq = 0
data_ex <- data.frame(fire_freq = 0,
                      area_ha = burnable_ha - burned_total_ha,
                      prob = NA)
reburns <- rbind(data_ex, dburn_freq)
reburns$prob <- reburns$area_ha / sum(reburns$area_ha)
# cumsum(reburns$prob[6:1])

fire_num_avg <- sum(reburns$fire_freq * reburns$prob)
fire_annual_prob <- fire_num_avg / length(1999:2022)
(fire_return <- 1 / fire_annual_prob)


# Comparison of annual burn rates -----------------------------------------

# Paritsis et al. 2013, Lagos region, 1984:2010
total_rate <- 440 / 19392 * 100
annual_rate <- total_rate / length(1984:2010)
# 0.084

# Holz et al. 2012,



