library(tidyverse)
library(terra)

v0 <- vect("data_and_files/data_spatial_variables_mindist_800m.shp")
pp <- rast("data_and_files/pp_atlas-climatico/geonode_precipitacion_anual.tif")

v2 <- terra::extract(pp, v0)

v0$pp_atlas <- v2$geonode_precipitacion_anual

mm <- lm(pp ~ pp_atlas, data = v0)
plot(v0$pp ~ v0$pp_atlas, col = rgb(0, 0, 0, 0.1), pch = 19)
abline(0, 1, col = "blue", lwd = 2)
abline(coef(mm), col = 2)

ggplot(as.data.frame(v0), aes(pp_atlas, pp)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linewidth = 0.6) +
  geom_smooth(method = "lm", color = "red", linewidth = 0.6, se = F) +
  coord_fixed() +
  theme_bw() +
  ylab("pp worldclim") +
  xlab("pp atlas") +
  theme(panel.grid.minor = element_blank())



