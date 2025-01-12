library(ggplot2)
library(viridis)
normalize <- function(x) x / sum(x)



N <- 100
a1 <- qlogis(0.2)
a2 <- qlogis(0.45)
b <- - 1 / 100

colcol <- viridis(2, end = 0.5)
colcol_shade <- viridis(2, end = 0.5, alpha = 0.05)


# Experiment 1 ------------------------------------------------------------

# vegetation type distribution functions
f1 <- function(x) ifelse(x < 1000, 1, 0)
f2 <- function(x) ifelse(x >= 1000, 1, 0)

png("figures/XX simulation correlation.png", width = 12, height = 17,
    units = "cm", res = 300)

par(mfrow = c(2, 1))

# Vegetation type ~ elevation
curve(plogis(a1 + b * (x-1000)), from = 500, to = 1500, ylim = c(0, 1),
      col = rgb(0, 0, 0, 0),
      ylab = "Relative cover", xlab = NA)
curve(f1(x), add = T, col = colcol[1], lwd = 1.5)
curve(f2(x), add = T, col = colcol[2], lwd = 1.5)
text(750, 0.8, "Type A", col = colcol[1])
text(1250, 0.8, "Type B", col = colcol[2])

# Burn probability ~ elevation
curve(plogis(a1 + b * (x-1000)), from = 500, to = 1500, ylim = c(0, 1),
      col = rgb(0, 0, 0, 0),
      ylab = "Burn probability", xlab = "Elevation (m a.s.l.)")
curve(plogis(a1 + b * (x-1000)), from = 500, to = 1000,
      add = T, col = colcol[1], lwd = 1.5)
curve(plogis(a2 + b * (x-1000)), from = 1000, to = 1500,
      add = T, col = colcol[2], lwd = 1.5)

for(i in 1:200) {
  x1 <- runif(N, 500, 1000)
  x2 <- runif(N, 1000, 1500)

  p1 <- plogis(a1 + b * (x1 - 1000))
  p2 <- plogis(a2 + b * (x2 - 1000))

  y1 <- rbinom(N, 1, p1)
  y2 <- rbinom(N, 1, p2)

  m1 <- glm(y1 ~ x1, family = "binomial")
  m2 <- glm(y2 ~ x2, family = "binomial")

  cc1 <- coef(m1)
  cc2 <- coef(m2)

  curve(plogis(cc1[1] + cc1[2] * x), from = 500, to = 1000,
        add = T, col = colcol_shade[1], lwd = 1)
  curve(plogis(cc2[1] + cc2[2] * x), from = 1000, to = 1500,
        add = T, col = colcol_shade[2], lwd = 1)
}


par(mfrow = c(1, 1))

dev.off()

# Experiment 2.1 ------------------------------------------------------------

cuts1 <- seq(500, 1000, by = 100)
cuts2 <- seq(1000, 1500, by = 100)

xmid1 <- seq(550, 950, 100)
xmid2 <- seq(1050, 1450, 100)

png("figures/XX simulation coarse environment_01.png", width = 12, height = 9,
    units = "cm", res = 300)

# Burn probability ~ elevation
curve(plogis(a1 + b * (x-1000)), from = 500, to = 1500, ylim = c(0, 1),
      col = rgb(0, 0, 0, 0),
      ylab = "Burn probability", xlab = "Elevation (m a.s.l.)")

abline(v = xmid1, col = colcol[1], lty = 2, lwd = 0.7)
abline(v = xmid2, col = colcol[2], lty = 2, lwd = 0.7)

curve(plogis(a1 + b * (x-1000)), from = 500, to = 1000,
      add = T, col = colcol[1], lwd = 1.5)
curve(plogis(a2 + b * (x-1000)), from = 1000, to = 1500,
      add = T, col = colcol[2], lwd = 1.5)

for(i in 1:200) {
  x1 <- runif(N, 500, 1000)
  x2 <- runif(N, 1000, 1500)

  x1cat <- cut(x1, breaks = cuts1)
  x2cat <- cut(x2, breaks = cuts2)

  x1fit <- as.numeric(x1cat) * 100 + 450
  x2fit <- as.numeric(x2cat) * 100 + 950

  p1 <- plogis(a1 + b * (x1 - 1000))
  p2 <- plogis(a2 + b * (x2 - 1000))

  y1 <- rbinom(N, 1, p1)
  y2 <- rbinom(N, 1, p2)

  m1 <- glm(y1 ~ x1fit, family = "binomial")
  m2 <- glm(y2 ~ x2fit, family = "binomial")

  cc1 <- coef(m1)
  cc2 <- coef(m2)

  curve(plogis(cc1[1] + cc1[2] * x), from = 500, to = 1000,
        add = T, col = colcol_shade[1], lwd = 1)
  curve(plogis(cc2[1] + cc2[2] * x), from = 1000, to = 1500,
        add = T, col = colcol_shade[2], lwd = 1)
}
dev.off()

# Experiment 2.2 ------------------------------------------------------------

N <- 100
K <- 4
a <- qlogis(runif(K, 0.2, 0.45))
b <- - 1 / 600

colcol <- viridis(K, end = 0.8)
colcol_shade <- viridis(K, end = 0.8, alpha = 0.05)

xwidth <- 500
xmin <- 0
xmax <- 2000

cuts <- sapply(1:K, function(k) {
  seq(xmin + xwidth * (k-1), xwidth * k, by = xwidth / 2)
})

xmid <- cuts[1:2, ] + xwidth / 4


## Plot
png("figures/XX simulation coarse environment_01.png", width = 12, height = 9,
    units = "cm", res = 300)

# Burn probability ~ elevation
curve(plogis(a[1] + b * (x-1000)), from = xmin, to = xmax, ylim = c(0, 1),
      col = rgb(0, 0, 0, 0),
      ylab = "Burn probability", xlab = "Elevation (m a.s.l.)")

for(k in 1:K) {
  # x data
  abline(v = xmid[, k], col = colcol[k], lty = 2, lwd = 0.7)

  # curves from simulated data
  for(i in 1:200) {
    x <- runif(N, cuts[1, k], cuts[3, k])
    xcat <- cut(x, breaks = cuts[, k])
    xfit <- as.numeric(as.character(factor(xcat, labels = as.character(xmid[, k]))))

    p <- plogis(a[k] + b * (x - 1000))
    y <- rbinom(N, 1, p)

    m <- glm(y ~ xfit, family = "binomial")
    cc <- coef(m)

    curve(plogis(cc[1] + cc[2] * x),  from = cuts[1, k], to = cuts[3, k],
          add = T, col = colcol_shade[k], lwd = 1)
  }

  # true function
  curve(plogis(a[k] + b * (x-1000)), from = cuts[1, k], to = cuts[3, k],
        add = T, col = colcol[k], lwd = 1.5)
}
dev.off()


# Batch 2 -----------------------------------------------------------------

N <- 100
a1 <- qlogis(0.2)
a2 <- qlogis(0.45)
b <- - 1 / 100

colcol <- viridis(2, end = 0.5)
colcol_shade <- viridis(2, end = 0.5, alpha = 0.05)




# Graph 1 ------------------------------------------------------------


a1 <- qlogis(0.15)
a2 <- qlogis(0.45)
b <- - 1 / 100
bveg <- - 1 / 75#132

# vegetation type distribution functions
png("figures/XX correlation_01.png", width = 12, height = 17,
    units = "cm", res = 300)

par(mfrow = c(2, 1))

# Vegetation type ~ elevation
curve(plogis(bveg * (x-1000)), from = 500, to = 1500, ylim = c(0, 1),
      col = rgb(0, 0, 0, 0),
      ylab = "Relative cover", xlab = NA)
curve(plogis(bveg * (x - 1000)), add = T, col = colcol[1], lwd = 1.5)
curve(1 - plogis(bveg * (x - 1000)), add = T, col = colcol[2], lwd = 1.5)

text(600, 0.8, "Type A", col = colcol[1])
text(1400, 0.8, "Type B", col = colcol[2])

# Mean elevation:
xlarge <- runif(1e7, 500, 1500)
w1 <- plogis(bveg * (xlarge-1000)) |> normalize()
xmean1 <- sum(xlarge * w1)
w2 <- (1 - plogis(bveg * (xlarge-1000))) |> normalize()
xmean2 <- sum(xlarge * w2)
abline(v = c(xmean1, xmean2), col = colcol, lwd = 0.8, lty = 2)
abline(v = 1000, col = "black", lwd = 1.2, lty = 2)


# Burn probability ~ elevation
curve(plogis(a1 + b * (x-1000)), from = 500, to = 1500, ylim = c(0, 1),
      col = rgb(0, 0, 0, 0),
      ylab = "Burn probability", xlab = "Elevation (m a.s.l.)")
curve(plogis(a1 + b * (x-1000)), from = 500, to = 1500,
      add = T, col = colcol[1], lwd = 1.5)
curve(plogis(a2 + b * (x-1000)), from = 500, to = 1500,
      add = T, col = colcol[2], lwd = 1.5)

abline(v = c(xmean1, xmean2), col = colcol, lwd = 0.8, lty = 2)
abline(v = 1000, col = "black", lwd = 1.2, lty = 2)

par(mfrow = c(1, 1))

dev.off()


# Burn prob means marginal to elevation:

x1 <- sample(xlarge, size = 1e7, replace = T, prob = w1)
x2 <- sample(xlarge, size = 1e7, replace = T, prob = w2)

pm1 <- mean(plogis(a1 + b * (x1 - 1000)))
pm2 <- mean(plogis(a2 + b * (x2 - 1000)))
pm <- c(pm1, pm2)


pcond <- plogis(c(a1, a2))

png("figures/XX correlation_02.png", width = 12, height = 10,
    units = "cm", res = 300)
barplot(pm ~ c("A", "B"), ylab = "Burn probability",
        xlab = NA, col = colcol, ylim = c(0, 1), font.main = 1,
        main = "At the elevation where\neach type occurs (averaged)")
dev.off()

png("figures/XX correlation_03.png", width = 12, height = 10,
    units = "cm", res = 300)
barplot(pcond ~ c("A", "B"), ylab = "Burn probability",
        xlab = "Vegetation type", col = colcol, ylim = c(0, 1), font.main = 1,
        main = "At elevation = 1000 m a.s.l.")
dev.off()




# Batch 3 -----------------------------------------------------------------


colcol <- viridis(2, end = 0.5)
colcol_shade <- viridis(2, end = 0.5, alpha = 0.05)

a1 <- qlogis(0.15)
a2 <- qlogis(0.45)
b <- - 1 / 200

lwd0 <- 0.5
lwd1 <- 1.2

A <- ggplot() +

  # full range curves
  geom_function(fun = function(x) plogis(a1 + b * (x-1000)),
                mapping = aes(colour = "A", linetype = "Absent"),
                linewidth = lwd0,
                xlim = c(500, 1500)) +
  geom_function(fun = function(x) plogis(a2 + b * (x-1000)),
                mapping = aes(colour = "B", linetype = "Absent"),
                linewidth = lwd0,
                xlim = c(500, 1500)) +

  # focal curves
  geom_function(fun = function(x) plogis(a1 + b * (x-1000)),
                mapping = aes(colour = "A", linetype = "Present"),
                linewidth = lwd1,
                xlim = c(500, 1200)) +
  geom_function(fun = function(x) plogis(a2 + b * (x-1000)),
                mapping = aes(colour = "B", linetype = "Present"),
                linewidth = lwd1,
                xlim = c(800, 1500)) +

  scale_color_viridis(discrete = TRUE, end = 0.5,
                      name = "Vegetation\ntype") +

  scale_linetype_manual(values = c("dashed", "solid"), name = "") +

  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  ggtitle("(1) Overlap")

B <- ggplot() +

  # full range curves
  geom_function(fun = function(x) plogis(a1 + b * (x-1000)),
                colour = colcol[1], xlim = c(500, 1500),
                linetype = "dashed", linewidth = lwd0) +
  geom_function(fun = function(x) plogis(a2 + b * (x-1000)),
                colour = colcol[2], xlim = c(500, 1500),
                linetype = "dashed", linewidth = lwd0) +

  # focal curves
  geom_function(fun = function(x) plogis(a1 + b * (x-1000)),
                colour = colcol[1], xlim = c(500, 1000),
                linewidth = lwd1) +
  geom_function(fun = function(x) plogis(a2 + b * (x-1000)),
                colour = colcol[2], xlim = c(1000, 1500),
                linewidth = lwd1) +

  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  ylab("Burn probability") +
  ggtitle("(2) Adjacent")
# B

C <- ggplot() +

  # full range curves
  geom_function(fun = function(x) plogis(a1 + b * (x-1000)),
                colour = colcol[1], xlim = c(500, 1500),
                linetype = "dashed", linewidth = lwd0) +
  geom_function(fun = function(x) plogis(a2 + b * (x-1000)),
                colour = colcol[2], xlim = c(500, 1500),
                linetype = "dashed", linewidth = lwd0) +

  # focal curves
  geom_function(fun = function(x) plogis(a1 + b * (x-1000)),
                colour = colcol[1], xlim = c(500, 800),
                linewidth = lwd1) +
  geom_function(fun = function(x) plogis(a2 + b * (x-1000)),
                colour = colcol[2], xlim = c(1200, 1500),
                linewidth = lwd1) +

  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  xlab("Elevation (m a.s.l.)") +
  ggtitle("(3) Disjoint")
C

library(patchwork)

out <- A / B / C +
  plot_layout(guides = "collect") & # Unifica las leyendas si aplica
  theme(legend.spacing.y = unit(0.2, "cm")) # Ajusta el espacio entre elementos de la leyenda

# Muestra el resultado
ggsave("figures/XX correlation 3 cases_b.png", plot = out,
       width = 10, height = 17, units = "cm")
