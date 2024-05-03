# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr
library(ggh4x)    # varying strip theme for veg_types and all together

library(scales)   # log scale
library(lubridate)

library(terra)
# library(rgdal)    # readOGR
# library(rgeos)    # gIntersection, gArea...
# library(maptools) # UnionSpatialPolygons     ## all replaced by terra

library(mgcv)
library(scam)
library(DHARMa)

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


# Interannual temporal patterns -------------------------------------------


## Prepare data

burned_annual <- read.csv("data/data_burned_area_by_year.csv")
burned_veg <- read.csv("data/data_burned_and_available_area_by_vegetation_dryforest2.csv")

# bring fwi data
fwi_data <- read.csv("data/data_climate_interannual_fwi.csv")
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
climate_long <- read.csv("data/data_climate_interannual.csv")
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
write.csv(coef_export, "exports/trends.csv", row.names = F)


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

# Intraannual fire activity -----------------------------------------------


fires <- vect(file.path("..", "patagonian_fires/patagonian_fires.shp"))
# # get clipped fires database (NOT USED)
# study_area <- vect(file.path("..", "study_area/study_area.shp"))
# fires_clipped <- crop(fires, study_area)
# fires_clipped$area_ha <- expanse(fires_clipped, unit = "ha")

# DO NOT use clipped fires here
fires_data <- as.data.frame(fires)
remove <- which(fires_data$obs == "uncertain_year")
dmonth <- fires_data[-remove, ]
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
data_my$year_factor <- factor(data_my$year)

# model
m_count <- gam(fires ~
                 s(month_num, bs = "cc",
                   k = length(unique(data_my$month_num))) +
                 s(year_factor, bs = "re"),
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
                                           length.out = length(date_seq))),
                         year_factor = "1999")

# compute prediction for an average year
lpmat <- predict(m_count, month_pred, type = "lpmatrix")
b <- coef(m_count)
V <- vcov(m_count, unconditional = T)
sigma_year <- gam.vcomp(m_count)["s(year_factor)", "std.dev"]
pars_use <- c(1, grep("month_num", names(b)))
bsim <- rmvn(10000, b, V)[, pars_use] %>% t
mu_sim <- exp(lpmat[, pars_use] %*% bsim + 0.5 * sigma_year ^ 2)

month_pred$count_mle <- exp(lpmat[, pars_use] %*% b[pars_use] +
                            0.5 * sigma_year ^ 2)
month_pred$count_lower <- apply(mu_sim, 1, quantile, 0.025)
month_pred$count_upper <- apply(mu_sim, 1, quantile, 0.975)
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

dmonth$year_factor <- factor(dmonth$year)

# model fit to raw data
m_size <- gam(area_ha ~
                s(month_num, bs = "cc",
                  k = length(unique(dmonth$month_num))) +
                s(year_factor, bs = "re"),
               data = dmonth, family = Gamma(link = "log"),
               method = "REML", gamma = 1,
               knots = list(month_num = c(0, 12)))

lpmat <- predict(m_size, month_pred, type = "lpmatrix")
b <- coef(m_size)
V <- vcov(m_size, unconditional = T)
sigma_year <- gam.vcomp(m_size)["s(year_factor)", "std.dev"]
pars_use <- c(1, grep("month_num", names(b)))
bsim <- rmvn(10000, b, V)[, pars_use] %>% t
mu_sim <- exp(lpmat[, pars_use] %*% bsim + 0.5 * sigma_year ^ 2)

month_pred$size_mle <- exp(lpmat[, pars_use] %*% b[pars_use] +
                           0.5 * sigma_year ^ 2)
month_pred$size_lower <- apply(mu_sim, 1, quantile, 0.025)
month_pred$size_upper <- apply(mu_sim, 1, quantile, 0.975)

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
                     labels = c("size", "n° fires"))
pplong$var <- factor(pplong$var, levels = c("size", "count"),
                     labels = c("size", "n° fires"))


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
    sec.axis = sec_axis(~ . / q, name = "Mean n° of fires")
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

# mean over the study area, temperature from worldclim
temp_intra <- read.csv("data/data_climate_monthly.csv")
temp_intra <- temp_intra[temp_intra$variable == "tavg", ]
names(temp_intra)[3] <- "mean"
temp_intra$variable <- "temp"

temp_intra_perc <- read.csv("data/data_climate_monthly_percentiles.csv")
temp_intra_perc <- temp_intra_perc[order(temp_intra_perc$var), ]
rows_temp_5 <- grep("tavg_p5$", temp_intra_perc$var)
rows_temp_95 <- grep("tavg_p95$", temp_intra_perc$var)

temp_intra$lower <- temp_intra_perc$value[rows_temp_5]
temp_intra$upper <- temp_intra_perc$value[rows_temp_95]

# bring precipitation from climatic atlas
pp_intra <- read.csv("data/pp_atlas-climatico/data_climate_monthly_pp_with_percentiles_atlas.csv")
pp_intra$variable <- "pp"
pp_intra <- pp_intra[, c("month", "variable", "pp_mean", "pp_lower", "pp_upper")]
names(pp_intra) <- names(temp_intra)

clim_intra_0 <- rbind(temp_intra, pp_intra)

clim_intra <- left_join(
  clim_intra_0,
  date_table,
  by = "month"
)

# scale for plot with 2 y-axes
coef_clim <- 2

clim_intra_sc <- clim_intra
rows_sc <- clim_intra$variable == "temp"
cols_sc <- names(clim_intra) %in% c("mean", "lower", "upper")

clim_intra_sc[rows_sc, cols_sc] <- clim_intra[rows_sc, cols_sc] * coef_clim

intra_clim <-
  ggplot(clim_intra_sc,
         mapping = aes(date, mean, ymin = lower, ymax = upper,
                       fill = variable, color = variable, group = variable)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  geom_point(size = 3, alpha = 0.7) +
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
size_data <- data.frame("area_abs" = fires_data$area_ha,
                        "year" = fires_data$year)

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
# 10% of fires:

row_lim <- which.min(abs(size_data$number_prop - 0.1))
size_sub <- size_data[1:row_lim, ]
# number of fires:
nrow(size_sub)
size_sub$year[order(size_sub$year)] %>% unique # years
size_sub$year[order(size_sub$year)] %>% unique %>% length # 12 years
size_sub$area_prop %>% max # 79 %

# eitght largest
size_data[1:8, ]

years_fires10p[order(years_fires10p)]
years_fires10p[order(years_fires10p)] %>% length # 12 years
# and 50 %?
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
# size_limits <- c(0, seq(1, 4.5, length.out = k+1))

size_data$size_class <- NA

for(i in 1:nrow(size_data)) {
  #i = 3
  cl <- as.numeric(size_data$area_abs_log10[i] > size_limits) %>% sum
  # after clipping, there are fires < 10 ha, so we convert them to that class
  if(cl == 0) cl <- 1
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
size_freq

# fire regime plots

size_season <-
ggarrange(size_props + ggtitle("A"),
          size_freq + ggtitle("B"),
          intra_fire + theme(legend.position = "bottom") +
            ggtitle("C") + xlab("Month"),
          intra_clim + theme(legend.position = "bottom") + ggtitle("D"),
          nrow = 2)

ggsave("figures/02) size distribution and seasonality_b.jpeg",
       size_season,
       width = 16, height = 14, units = "cm", dpi = 300)

# Burned area global results ----------------------------------------------

dburn_freq <- read.csv("data/data_reburn_frequency.csv")[, -1]
dburn_veg <- read.csv("data/data_burned_and_available_area_by_vegetation.csv")

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



