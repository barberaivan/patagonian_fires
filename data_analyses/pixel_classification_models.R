
# Script heredado de
# "/home/ivan/Insync/Burned area mapping/Long mapping approach/glms fit nbr_min - nbr_delta - dnbr.R"

# Acá edito, reduzco, traduzco y mejoro los plots para que vayan como material
# suplementario.

# Comments ----------------------------------------------------------------

# En este script ajustaré muchos glms a datos de NBR para predecir prob de quemado.
# En el script de GEE 'dNBR mapping 7 (SNIC p02)', en el que mapeaba solo con
# dNBR, hice la siguiente descripción:

# ----
# Mapeo de fuegos usando solo dNBR con algunos criterios para semillas

# Pasos:
#
#   1) Definir valores de dNBR con los que se estima una p(burn)
# >0.2 y >0.8, para definir semillas.
# bien ese unbral. Podria ser 0.9 tambien.
#
# p(burn) = 0.2 when dNBR = 0.1877488
# p(burn) = 0.5 when dNBR = 0.2917249
# p(burn) = 0.8 when dNBR = 0.3957009
#
# Del GLM que fue ajustado con todos los datos:
#   p(burn) = 1 / (1 + exp(-(-3.89 + 13.33 * dNBR)))
#
# 2) Mascara para burnt: p(burnt) >0.5.
# Tamben mascarear NDVI_max <0.35 para sacar altoandina y estepa.
#
# 3) Mascara para semilla:
#   p(burnt) >0.8 &
#   dNBR = max(dNBR en años vecinos) &
#     NDVI_max > 0.5 (para eliminar pastizales muy secos).
# ----

# Actualmente, la idea es la siguiente:

# Candidatos:
#   nbr_min bajo | dnbr alto
#   NDVI_max (all years) > 0.35

#   para cada índice puedo ajustar un modelo y definir una p(burn) baja
#   como umbral. La idea sería que los candidatos tenga alta sensibilidad,
#   a costa de un alto error de comisión (baja especificidad).

#   NOTA sobre dNBR: este índice trae problemas para detectar el año de quemado
#   cuando el fuego ocurre muy al principio o al final de la temporada.
#   Aún peor es si hay dos fuegos así muy cerca (problema de Alerces).
#   Si con semillas de años focales seguimos teniendo ese problema, podemos
#   considerar como válidos solo los dNBR que tengan una alta señal de quemado en
#   el año focal en un radio de 500 m aprox.
#   Pero primero probaría sin esa restricción.
#   OJO, porque sin esa restricción quizás haya que poner un umbral de dNBR no
#   tan permisivo

# Semillas

#   con dNBR tuve el siguiente criterio:
#   p(dNBR_burned) > 0.8,
#   dNBR = máximo de años vecinos
#   NDVI_max (de todos los años) > 0.5  saca pastizales secos y altoandinas

#   Con NBR:
#
#   p(NBR_min burned) > umbral alto
#   p(NBR_delta burned) > umbral alto
#
#   Estos umbrales podrían ser la menor p que lleva el error de omisión a 0.01



# Datos, paquetes y funciones ---------------------------------------------

library(arm)
library(tidyverse); theme_set(theme_bw())
library(caret)
library(plyr)
library(viridis)
library(ggpubr)
library(gridExtra)

train <- read.csv("data/pixel_classification_training_data.csv")
names_to_change <- which(colnames(train) %in% c("nbr_delta", "nbr_min", "dnbr"))
colnames(train)[names_to_change] # "dnbr"      "nbr_delta" "nbr_min"
colnames(train)[names_to_change] <- c("dnbr_mean", "dnbr_min", "nbr_min")

# elimino un errorcito muy evidente
errorcito <- with(train, which(nbr_min < -0.2 & Burned == 0))
train <- train[-errorcito, ]

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


# Exploration -------------------------------------------------------------

# vars que quiero mirar:
# names(cho)
# antes se llamaban así: c("nbr_delta", "nbr_min", "dnbr")
pred_names <- c("dnbr_min", "nbr_min", "dnbr_mean")

train_long <- pivot_longer(train, which(names(train) %in% pred_names),
                           names_to = "variable",
                           values_to = "value")
train_long$variable <- revalue(train_long$variable,
                               c("nbr_min" = "NBR_min",
                                 "dnbr_min" = "dNBR_min",
                                 "dnbr_mean" = "dNBR_mean"))

ggplot(#train_long[train_long$Fire_ID == "Lolog_2002", ],
       train_long,
       aes(x = value, y = Burned, colour = ndvi_max)) +
  geom_jitter(alpha = 0.35) +
  facet_wrap(vars(variable), nrow = 2, scales = "free") +
  geom_smooth(data = train_long,
              mapping = aes(x = value, y = Burned, colour = "black"),
              color = "black", size = 0.5,
              method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE) +
  scale_color_gradient(low = "blue", high = "red", trans = "log")
# en Lolog 2002 hay varios quemados con NBR muy alto. Eso podría ser por falta
# de imágenes?

ss <- subset(train, subset = Fire_ID == "Lolog_2002" & Burned == 1 & nbr_min > 0.1)
pairs(ss[, pred_names])
# muestra un patrón razonable, no da para sacarlos.

# Modelos de nbr(s) * ndvi ---------------------------------------------------

m_min <- glm(Burned ~ nbr_min * ndvi_max, data = train, family = "binomial")
m_delta <- glm(Burned ~ dnbr_min * ndvi_max, data = train, family = "binomial")
m_dnbr <- glm(Burned ~ dnbr_mean * ndvi_max, data = train, family = "binomial")

nrep = 150
p_min <- expand.grid(
  nbr_min = seq(min(train$nbr_min), max(train$nbr_min), length = nrep),
  ndvi_max = c(0.5, 0.7, 0.9)
)

p_delta <- expand.grid(
  dnbr_min = seq(min(train$dnbr_min), max(train$dnbr_min), length = nrep),
  ndvi_max = c(0.5, 0.7, 0.9)
)

p_dnbr <- expand.grid(
  dnbr_mean = seq(min(train$dnbr_mean), max(train$dnbr_mean), length = nrep),
  ndvi_max = c(0.5, 0.7, 0.9)
)

p_min$p <- predict(m_min, p_min, type = "response")
p_delta$p <- predict(m_delta, p_delta, type = "response")
p_dnbr$p <- predict(m_dnbr, p_dnbr, type = "response")

names(p_min)[1] <- names(p_delta)[1] <- names(p_dnbr)[1] <- "x"

p_all <- rbind(p_min, p_delta, p_dnbr)
p_all$variable <- rep(c("NBR_min", "dNBR_min", "dNBR_mean"), each = nrow(p_min))

plot1 <-
ggplot(train_long, aes(x = value, y = Burned, colour = ndvi_max)) +
  geom_jitter(alpha = 0.08, height = 0.2) +
  #scale_color_gradient(low = "blue", high = "red") +
  geom_line(data = p_all, mapping = aes(x = x, y = p,
                                        group = as.factor(ndvi_max))) +
  facet_wrap(vars(variable), nrow = 1) +
  scale_colour_viridis(direction = -1, option = "B") +
  xlab("NBR-based index") +
  ylab("Burn probability") +
  labs(colour = "NDVI_max") +
  theme(legend.position = "right",
        legend.title = element_text()) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  ggtitle("A")
plot1

# kesigue? ----------------------------------------------------------------

# Evaluar sensibilidad y especificidad en función de umbrales de probabilidad
# para cada modelo. 30 fold CV.
# nrow(train) / 20

k <- 20
set.seed(1)
fold_id <- rep(1:k, length.out = nrow(train)) %>%
  sample(size = nrow(train), replace = FALSE)

# Función para evaluar sens-spec en función del umbral
cutss <- function(cutoff, prediction, y_fac) {
  cutoff <- cutoff
  classif <- (prediction > cutoff) %>% as.numeric() %>% factor(levels = c(0, 1))
  ss <- confusionMatrix(classif, y_fac, mode = "sens_spec", positive = "1") %$%
    byClass %>% '['(1:2)
  return(ss)
}

cut_acc <- function(cutoff, prediction, y_fac) {
  cutoff <- cutoff
  classif <- (prediction > cutoff) %>% as.numeric() %>% factor(levels = c(0, 1))
  ss <- confusionMatrix(classif, y_fac, mode = "sens_spec", positive = "1") %$%
    overall %>% '['(1)
  return(ss)
}


# Tabla con probabilidades predichas a partir de modelos ajustados sin las filas
# de interés

ptable <- data.frame(
  nbr_min = numeric(nrow(train)),
  dnbr_min = numeric(nrow(train)),
  dnbr_mean = numeric(nrow(train))
)

sapply(1:k, function(i) {
  #i = 1
  print(i)

  train_set = train[fold_id != i, ]
  test_set = train[fold_id == i, ]

  m_min_cv <- glm(Burned ~ nbr_min * ndvi_max, data = train_set,
                  family = "binomial")
  m_delta_cv <- glm(Burned ~ dnbr_min * ndvi_max, data = train_set,
                    family = "binomial")
  m_dnbr_cv <- glm(Burned ~ dnbr_mean * ndvi_max, data = train_set,
                   family = "binomial")

  ptable[fold_id == i, 'nbr_min'] <<- predict(m_min_cv, test_set, type = "response")
  ptable[fold_id == i, 'dnbr_min'] <<- predict(m_delta_cv, test_set, type = "response")
  ptable[fold_id == i, 'dnbr_mean'] <<- predict(m_dnbr_cv, test_set, type = "response")

}) %>% invisible

ptable$y_fac <- train$Burned %>% as.character %>% factor(levels = c("0", "1"))
# summary(ptable)


# Curva sens-spec para todos los modelos

cut_seq <- seq(0, 1, length.out = 101)

senspec <- sapply(1:3, function(m) {
  #m = 1
  sstable <- sapply(cut_seq, cutss, prediction = ptable[, m],
                    y_fac = ptable$y_fac, simplify = "array") %>% t()
  return(sstable)
}, simplify = FALSE)

str(senspec)

senspec <- do.call('rbind', senspec) %>% as.data.frame
senspec$model <- rep(names(ptable)[1:3], each = length(cut_seq))
senspec$model <- revalue(senspec$model,
                         c("nbr_min" = "NBR_min",
                           "dnbr_min" = "dNBR_min",
                           "dnbr_mean" = "dNBR_mean"))

senspec$threshold <- rep(cut_seq, 3)

senspec_long <- pivot_longer(senspec, 1:2, names_to = "performance",
                             values_to = "senspec")

senspec_long$error_value <- 1 - senspec_long$senspec
senspec_long$error_name <- revalue(senspec_long$performance, c("Sensitivity" = "Omission",
                                      "Specificity" = "Comission"))

# sens - spec
ggplot(senspec_long, aes(x = threshold, y = senspec, colour = performance)) +
  geom_line(size = 0.7) +
  ylim(0.750, 1) + facet_wrap(vars(model), nrow = 1) +
  ggtitle("Performance de modelos (20-fold CV)") +
  ylab("Sensibilidad o especificidad") +
  xlab("Umbral de clasificación [P(y = 1|x)]") +
  theme(legend.position = "bottom") + labs(colour = NULL)

# com - om
cand_lines <- data.frame(threshold = c(0.20, 0.25, 0.25),
                         model = c("dNBR_mean", "dNBR_min", "NBR_min"))


plot2 <-
ggplot(senspec_long, aes(x = threshold, y = error_value, colour = error_name)) +
  geom_line(size = 0.7) +
  ylim(0, 0.15) +
  facet_wrap(vars(model), nrow = 1) +
  ggtitle("B") +
  ylab("Error rate") +
  xlab("Burn probability threshold") +
  theme(legend.position = "right",
        panel.grid.minor = element_blank()) +
  scale_color_viridis(discrete = TRUE, option = "B", end = 0.5) +
  labs(colour = NULL) +
  geom_vline(data = cand_lines,
             mapping = aes(xintercept = threshold),
             linetype = 3, alpha = 0.8, size = 0.5) +
  geom_vline(xintercept = 0.98, linetype = 2, alpha = 0.8,
             size = 0.5)

plot2

# Save plots --------------------------------------------------------------

marg <- 2
ggarrange(plot1 +
            theme(plot.margin = unit(c(marg, marg + 4.5, marg, marg), "mm"),
                  axis.text = element_text(size = 7)),
          plot2 +
            theme(axis.text = element_text(size = 7)),
          nrow = 2)
ggsave("figures/S03) pixel level classication.jpeg",
       width = 16, height = 11, units = "cm", dpi = 300)





# Umbrales para semillas --------------------------------------------------

# Omisión y comisión en wide
senspec$Comission <- 1 - senspec$Specificity
senspec$Omission <- 1 - senspec$Sensitivity


# buscamos error de comisión muy bajo:
# menor o igual a 0.001 (redondeando)
senspec[senspec$Specificity >= 0.998 & senspec$model == "nbr_min", ]
senspec[senspec$Comission <= 0.0015 & senspec$model == "nbr_min", ]

# nbr_min: su mínima p es 0.84
#   especificidad = 0.9989142
#   sensibilidad = 0.9542125
#   comisión = 0.0010857763
senspec[senspec$Specificity >= 0.998 & senspec$model == "dnbr_min", ]
senspec[senspec$Comission <= 0.0015 & senspec$model == "dnbr_min", ]

# dnbr_min con p >= 0.98 tiene
#   especificidad = 0.9267399
#   sensibilidad = 0.9989142
#   comisión = 0.001085776

# semilla =
#   p_nbr_min >= 0.84 | p_dnbr_min >= 0.98


# seed para dnbr_mean en años en que hay pocas imagenes
# (en otros años no se usará como criterio de semilla)
senspec[round(senspec$Comission, 3) <= 0.001 & senspec$model == "dnbr_mean", ]
# dnbr_mean con p >= 0.91 tiene
#   especificidad = 0.9994571
#   sensibilidad = 0.8876679
#   comisión = 0.0005428882

# Umbrales para candidates ------------------------------------------------

# Buscamos sensibilidad muy alta
# Un criterio puede ser max(omisión) = 0.001
# O max(comision) <= 0.10

senspec[round(senspec$Comission, 2) <= 0.1 & senspec$model == "nbr_min", ]
# nbr_min: su mínima p es 0.06
#   especificidad =
#   sensibilidad =
#   comision = 0.1015200869
#   omisión =  0.007326007           chequear el patrón con
#                                    region growing

senspec[round(senspec$Comission, 2) <= 0.1 & senspec$model == "dnbr_min", ]
# dnbr_min con p >= 0.04 tiene
#   especificidad =
#   sensibilidad =
#   comisión =
#   Omisión =

senspec[round(senspec$Comission, 2) <= 0.1 & senspec$model == "dnbr_mean", ]
# dnbr_mean con p >= 0.18 tiene
#   especificidad = 0.8973941
#   sensibilidad = 0.980464
#   comisión = 0.1026059
#   Omisión = 0.01953602

# Si ponemos este umbral de dNBR para candidato, deberíamos restringirlo a que
# tenga alguna semilla en un radio de 500 m
# pensar más eso


# candidato =

#   {opción 1}
#     p_nbr_min >= 0.06 | p_dnbr_min >= 0.04 | p_dnbr >= 0.18

#   {opción 2}
#     p_nbr_min >= 0.06 | p_dnbr_min >= 0.04 |
#     [p_dnbr >= 0.18 & max(p_nbr_min) >= 0.06 en cierto radio]

# Curva de accuracy para todos los modelos ----------------------------------

accu <- sapply(1:3, function(m) {
  #m = 1
  accu_table <- sapply(cut_seq, cut_acc, prediction = ptable[, m],
                       y_fac = ptable$y_fac, simplify = "array") %>% t()
  return(accu_table)
}, simplify = T) %>% as.data.frame

str(accu)
names(accu) <- names(ptable)[1:3]
accu$threshold <- cut_seq

accu_long <- pivot_longer(accu, 1:3, names_to = "model",
                          values_to = "accuracy")

ggplot(accu_long, aes(x = threshold, y = accuracy, colour = model)) +
  geom_line(size = 0.7) +
  ylim(0.9, 1) +
  ggtitle("Performance de modelos (20-fold CV)") +
  ylab("Accuracy") +
  xlab("Umbral de clasificación [P(y = 1|x)]") +
  theme(legend.position = "right") + labs(colour = NULL)


# Son buenísimos. De hecho, es increíblemente bueno el dnbr_min.
# Pero ambos mucho mejores que el dnbr_mean.
# Aún así, dnbrsubió su accuracy agregando el NDVI (un 1 %)


# Coefs para copiar -------------------------------------------------------

coef(m_min)
coef(m_delta)
coef(m_dnbr)

