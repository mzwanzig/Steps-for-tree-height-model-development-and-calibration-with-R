# R-Script for tree height model development for A. araucana and N. pumilio
# Please open the file with encoding UTF-8
# Author: Xinying Zhou, Edited: Martin Zwanzig
#
# Supplement to the following manuscript submitted to Annals of Forest Science:
# Zhou X, Kutchartt E, Hernández J, Corvalán P, Promis Á, Zwanzig M (2022) Determination of
# optimal tree height models and calibration designs for Araucaria araucana and Nothofagus
# pumilio in mixed stands affected to different levels by anthropogenic disturbance in
# south-central Chile.

# Packages ----
library(ggplot2)
library(ggpubr)
library(cowplot)
library(aplot)
library(grid)
library(gridExtra)
library(lmfor)
library(nlstools)
library(modelr)
library(rcompanion)
library(caret)
library(lemon)
library(ggpubr)
library(psych)
library(leaps)
library(fastDummies)
library(minpack.lm)
library(nlme)
library(dplyr)
library(MASS)
library(factoextra)
library(svglite)
library(readxl)

# load data ----
rawdata <- readxl::read_xlsx("AA_NP_stem_diameter_height.xlsx", sheet = "dbh_height_data") # loads horizontal point sampling measurements 
colnames(rawdata) <- c('Site', 'Stand', 'Plot', 'Species', 'd', 'h')

kmall <- read.csv("kmeanall.csv") # this is a database of stand variabales, independent of species, calculated from the raw data and the stand variables formulae in the supplementary
kmaa <- read.csv("kmeanaa.csv") # this is a database of stand variabales of A. araucana, calculated from the raw data and the stand variables formulae in the supplementary
kmnp <- read.csv("kmeannp.csv") # this is a database of stand variabales of N. pumilio, calculated from the raw data and the stand variables formulae in the supplementary

# data cleaning/selection
rawdata_nond <- rawdata[- which(rawdata$Species == 'nd'), ] # excluding Nothofagus dombeyi
rawdata_d <- rawdata_nond[- which(rowSums(is.na(rawdata_nond)) > 0), ] # considering complete cases only
rawdata_h <- rawdata_d[- which(rawdata_d$h == ''), ] # considering only observations with measured tree heights
database <- unique(rawdata_h) # Remove one duplicate row

# change the data type of the variables
database$h <- as.numeric(database$h)
database$Site <- factor(database$Site, levels = c('La Fusta', 'Conguillio', 'Malalcahuello', 'El Naranjo'))
database$Stand <- factor(database$Stand, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'))

# build database for AA
database_AA <- database[which(database$Species == 'aa'), ]
fitting_AA <- database_AA[database_AA$Stand %in% c('1', '3', '4', '5', '6', '7', '8', '9', '10', '11'), ]
calibrating_AA <- database_AA[database_AA$Stand %in% c('2', '12'), ]

# build database for NP
database_NP <- database[which(database$Species == 'np'), ]
fitting_NP <- database_NP[database_NP$Stand %in% c('1', '3', '4', '5', '6', '7', '8', '9', '10', '11'), ]
calibrating_NP <- database_NP[database_NP$Stand %in% c('2', '12'), ]

# create scatter plots for different databases ----
fitaa <- ggplot(data = fitting_AA, aes(x = d, y = h, colour = Site)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Diameter at the breast height (cm)') + ylab('Height (m)') +
  labs(title = '(b) Fitting_AA') +
  guides(colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'))

caliaa <- ggplot(data = calibrating_AA, aes(x = d, y = h, colour = Site)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Diameter at the breast height (cm)') + ylab('Height (m)') +
  labs(title = '(c) Validating_AA') +
  guides(colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'))

totalaa <- ggplot(data = database_AA, aes(x = d, y = h, colour = Site)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Diameter at the breast height (cm)') + ylab('Height (m)') +
  labs(title = '(a) Database_AA') +
  guides(colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'))

fitnp <- ggplot(data = fitting_NP, aes(x = d, y = h, colour = Site)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Diameter at the breast height (cm)') + ylab('Height (m)') +
  labs(title = '(e) Fitting_NP') +
  guides(colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'))

calinp <- ggplot(data = calibrating_NP, aes(x = d, y = h, colour = Site)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Diameter at the breast height (cm)') + ylab('Height (m)') +
  labs(title = '(f) Validating_NP') +
  guides(colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'))

totalnp <- ggplot(data = database_NP, aes(x = d, y = h, colour = Site)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Diameter at the breast height (cm)') + ylab('Height (m)') +
  labs(title = '(d) Database_NP') +
  guides(colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'))

fig.3_point <- grid_arrange_shared_legend(totalaa, totalnp, fitaa, fitnp, caliaa, calinp, ncol = 2, nrow = 3)

# kmeans ----

# because different stand variables have different scales, all stand variables are first scaled
kmallsca <- scale(kmall) # scale is by default applied column-wise (separately for each variable)
kmaasca <- scale(kmaa)[1:12, ]
kmnpsca <- scale(kmnp)[1:12, ]
kmcom <- cbind(kmallsca, kmaasca, kmnpsca) # forming a database of all stand variables

# kmeans
set.seed(2021)
kmcomre <- kmeans(kmcom, 4, nstart = 25)

fig.2_kmeans <- fviz_cluster(kmcomre, data = kmcom,
             palette = c('#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3'),
             ellipse.type = 'euclid', star.plot = TRUE, repel = TRUE,
             ggtheme = theme_minimal()) +
  ggtitle(label = '') + theme(legend.position = 'none')

ggsave(file = 'fig.2_kmeans.svg', plot = fig.2_kmeans, width = 6, height = 4) # vector figure for Fig.2


# plot of stand variables in each stand that contribute significantly to clustering ----

svclu <- read.csv("figa1.csv") # this a complete database of all stand variables and disturbance levels collated as needed for the supplementary fig.a1
svclu$Störung <- factor(svclu$Störung, levels = c('none', 'low', 'medium', 'high'))

RNNP <- ggplot(svclu, aes(x = factor(Bestand), y = RNNP)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung, scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('RNNP') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

RNAA <- ggplot(svclu, aes(x=factor(Bestand), y = RNAA)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung, scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('RNAA') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

NAA <- ggplot(svclu, aes(x = factor(Bestand), y = NAA)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung, scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('NAA') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

HdNP <- ggplot(svclu, aes(x = factor(Bestand), y = HdNP)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung, scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('HdNP') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

Hdom <- ggplot(svclu, aes(x = factor(Bestand), y = Hdom)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung, scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('Hdom') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

HdomAA <- ggplot(svclu, aes( x = factor(Bestand), y = HdomAA)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung, scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('HdomAA') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

HdomNP <- ggplot(svclu, aes(x = factor(Bestand), y = HdomNP)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung, scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('HdomNP') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

DBHdomNP <- ggplot(svclu, aes(x = factor(Bestand), y = DBHdomNP)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.4) +
  facet_wrap(~ Störung,scales = 'free', nrow = 1) +
  xlab('Stand') + ylab('DBHdomNP') +
  theme(strip.text = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

fig.a1_supsvbar <- plot_grid(RNNP, RNAA, NAA, HdNP, Hdom, HdomAA, HdomNP, DBHdomNP,
          ncol = 2, labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'))


# selection of simple models ----

# A. araucana - Gompertz model - example
startHDgomperz(fitting_AA$d, fitting_AA$h) # get initial values
SM12AA <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)),
              start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
              data = fitting_AA) # model fitting
# different initial values have been tried to reach a global minimum, but it was very likely that some initial values prevent the convergence of the model

overview(SM12AA) # overview model for PRSE calculation

# define a funtion mpc() to evaluate the goodness-of-fit of the model
mpc <- function(model, data){
  RMSE <- modelr::rmse(model = model, data = data)
  Rsquare<-modelr::rsquare(model = model, data = data)
  MAE<-modelr::mae(model = model, data = data)
  return(data.frame(RMSE = RMSE, Rsquare = Rsquare, MAE = MAE))
}

# evaluation the goodness-of-fit of SM12AA
mpc(SM12AA, fitting_AA)

# 10-fold cross-validation repeated 10 times to validate the predictive performance of the model in fitting database using RMSE
set.seed(2021)
folds <- createFolds(fitting_AA$Stand, k = 10)
output <- c()

j <- 1
ncv = 10
while (j <= ncv) {
  for (i in 1:10) {
    fold_train <- fitting_AA[- folds[[i]], ]
    fold_test <- fitting_AA[folds[[i]], ]
    fold_model <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), data = fold_train,
                      start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333))
    RMSE <- modelr::rmse(model = fold_model, data = fold_test)
    output <- rbind(output, cbind(RMSE))
  }
  j <- j + 1
}

output <- data.frame(output)
cat(mean(output$RMSE))

# N. pumilio - Näslund model - example
startHDnaslund(fitting_NP$d, fitting_NP$h)
SM1NP <- nls(h ~ 1.3 + (d ^ 2) / ((a + b * d) ^ 2), 
             data = fitting_NP,
             start = list(a = 3.3212685, b = 0.1995079))
overview(SM1NP)
mpc(SM1NP, fitting_NP)

set.seed(2021)
folds <- createFolds(fitting_NP$Stand, k = 10)
output <- c()
j <- 1
ncv = 10
while(j <= ncv) {
  for(i in 1:10) {
    fold_train <- fitting_NP[-folds[[i]], ]
    fold_test <- fitting_NP[folds[[i]], ]
    fold_model <- nls(h ~ 1.3 + (d ^ 2) / ((a + b * d) ^ 2),
                      data = fold_train,
                      start = list(a = 3.3212685, b = 0.1995079))
    RMSE <- modelr::rmse(model = fold_model, data = fold_test)
    output <- rbind(output, cbind(RMSE))
  }
  j <- j+1
}

output <- data.frame(output)
cat(mean(output$RMSE))


# residual diagnostics of the simple models ----

# A. araucana
plot(SM12AA) # plot showed a clear pattern --> WLS

weight1AA <- 1 / (fitting_AA$d ^ 0.5)
weight2AA <- 1 / (fitting_AA$d ^ 1)
weight3AA <- 1 / (fitting_AA$d ^ 1.5)
weight4AA <- 1 / (fitting_AA$d ^ 2)
weight5AA <- 1 / (fitting_AA$d ^ 2.5)
weight6AA <- 1 / (fitting_AA$d ^ 3)

SM12AA1 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)),
               weights = weight1AA,
               start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
               data = fitting_AA)
SM12AA2 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)),
               weights = weight2AA,
               start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
               data = fitting_AA)
SM12AA3 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)),
               weights = weight3AA,
               start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
               data = fitting_AA)
SM12AA4 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)),
               weights = weight4AA,
               start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
               data = fitting_AA)
SM12AA5 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)),
               weights = weight5AA,
               start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
               data = fitting_AA)
SM12AA6 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)),
               weights = weight6AA,
               start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
               data = fitting_AA)

RDAA0 <- plot(SM12AA, main = '(a) unweighted', col = 'black', pch = 20, cex = 0.7,
                              xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
RDAA1 <- plot(SM12AA1, main = '(b) k = 0.5', col = 'black', pch = 20, cex = 0.7,
              xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
RDAA2 <- plot(SM12AA2, main = '(c) k = 1', col = 'black', pch = 20, cex = 0.7,
              xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
RDAA3 <- plot(SM12AA3, main = '(d) k = 1.5', col = 'black', pch = 20, cex = 0.7,
              xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
RDAA4 <- plot(SM12AA4, main = '(e) k = 2', col = 'black', pch = 20, cex = 0.7,
              xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
RDAA5 <- plot(SM12AA5, main = '(f) k = 2.5', col = 'black', pch = 20, cex = 0.7,
              xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
RDAA6 <- plot(SM12AA6, main = '(g) k = 3', col = 'black', pch = 20, cex = 0.7,
              xlab = 'Predicted height (m)', ylab = 'Studentized residuals')

fig.a2_wls <- plot_grid(NULL, RDAA0, NULL, RDAA1, RDAA2, RDAA3, RDAA4, RDAA5, RDAA6)

plot(SM12AA2, Stand ~ resid(.), abline = 0) # plot showed that the effect of different stands was remained and further consideration was required

# residual diagnostic process for the simple model of N. pumilio was the same as that of A. araucana


# inclusion of stand variables in the simple models ----

# A. araucana

# Fitting a separate nls model for each stand can be done using the function nlsList(). However, since our simple model was weighted, nlsList() did not work here.

# read the parameters of the models fitted for each stand
AA1 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '1', ])
AA3 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '3', ])
AA4 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '4', ])
AA5 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '5', ])
AA6 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '6', ])
AA7 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '7', ])
AA8 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '8', ])
AA9 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '9', ])
AA10 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '10', ])
AA11 <- nls(h ~ 1.3 + a * exp(-b * exp(-c * d)), weights = 1 / (d ^ 1), 
           start = list(a = 46.410000000, b = 2.174185939, c = 0.008847333),
           data = fitting_AA[fitting_AA$Stand == '11', ])

parmsaa <- data.frame(AA1$m$getPars(), AA3$m$getPars(), AA4$m$getPars(),
                      AA5$m$getPars(), AA6$m$getPars(), AA7$m$getPars(),
                      AA8$m$getPars(), AA9$m$getPars(), AA10$m$getPars(),
                      AA11$m$getPars())
print(t(parmsaa))

# input stand variables database
AASV <- read.csv("AASVresults.csv") # database of stand variables, including stand variables not related to species and stand variables related to A. araucana

# check normal distribution of stand variables
# for example: NAA - number of AA trees per hectare
shapiro.test(AASV$NAA) # p-value < 0.05 --> not normal distribution
shapiro.test(log10(AASV$NAA)) # try logarithmic transformation --> p-value > 0.05
ggqqplot(log10(AASV$NAA))

# re-input the database that have been tested for normal distribution and mathematically transformed
AASVres <- read.csv("AASVres.csv") # AASVresults.csv + parameters a, b and c of stands in fitting database

# correlation test

# selection the stand variables to be used in the method 1 (2-phase approach with correlation analysis in the second phase)
AASV1 <- AASVres[, c(3, 4, 7, 8, 11, 12, 15, 17 : 32, 35, 36, 42, 43, 44)]
# correlation test (Pearson)
aasvcor <- corr.test(AASV1)
# get the results of the correlation test (r- and p-value)
aasvcorr <- aasvcor$r[, 26 : 28]
aasvcorp <- aasvcor$p[26 : 28, ]
aasvcorcom <- data.frame(aasvcorr, t(aasvcorp))
colnames(aasvcorcom) <- c('a(r)', 'b(r)', 'c(r)', 'a(p)', 'b(p)', 'c(p)')

# latitude and longitude could not be converted into a normal distribution by simple arithmetic, so the Spearman correlation test should be used
aaalat <- cor.test(AASVres$a, AASVres$BG, method = 'spearman')
aablat <- cor.test(AASVres$b, AASVres$BG, method = 'spearman')
aaclat <- cor.test(AASVres$c, AASVres$BG, method = 'spearman')
aaalon <- cor.test(AASVres$a, AASVres$LG, method = 'spearman')
aablon <- cor.test(AASVres$b, AASVres$LG, method = 'spearman')
aaclon <- cor.test(AASVres$c, AASVres$LG, method = 'spearman')

lat <- c(aaalat$estimate, aablat$estimate, aaclat$estimate, aaalat$p.value, aablat$p.value, aaclat$p.value)
lon <- c(aaalon$estimate, aablon$estimate, aaclon$estimate, aaalon$p.value, aablon$p.value, aaclon$p.value)

# combine all the correlation test results
aasvfin <- rbind(aasvcorcom[1 : 25, ], lat = lat, lon = lon, aasvcorcom[26 : 28, ])

# visualising the results of correlation tests
# load results database with added correlation classes and categories of the stand variables
enaasvres <- read.csv("aasvresen.csv") # results of correlation analysis of stand variables and parameters, including parameters, stand variables, r- and p-value, anthropogenic disturbance levels and group
# give the order of the stand variables and so on
enaasvres$SV <- factor(enaasvres$SV, levels = c('DBHmAA', 'DBHm', 'HdomAA', 'Hdom', 'DBHdomAA',
                                                'DBHdom', 'HDmaxAA', 'RangeHDAA', 'RangeDBHAA',
                                                'DBHmaxAA', 'DBHmax', 'DqAA', 'Dq', 'VardAA',
                                                'GCdAA', 'HdAA', 'SDIAA', 'NlogAA', 'Nlog',
                                                'BALAA', 'BAAA', 'BA', 'lng', 'lat', 'ALT', 
                                                'RNAA', 'RBAAA'))
enaasvres$level <- factor(enaasvres$level, levels = c('hoch', 'mäßig', 'schwach', 'unwesentlich'))
enaasvres$Gruppe <- factor(enaasvres$Gruppe, levels = c('Tree', 'Stand', 'DBH diversity', "Density and competition", "GI", "MS"))
enaasvres <- enaasvres[1 : 81, ]
# plot
aacorplot <- ggplot(enaasvres, aes(x = SV, y = r, fill = level)) +
  geom_bar(stat = 'identity', width = 0.5) +
  facet_grid(Gruppe ~ Parameter, scales = 'free', space = 'free') +
  ylim(-1, 1) +
  coord_flip() +
  theme(legend.position = 'none') +
  ylab('Correlation coefficient') + xlab('Stand variables') +
  scale_fill_manual(values = c("#000000", "#8D8C8B", "#C0C0C0", "#DCDCDC")) +
  theme(strip.text.y = element_text(size = 7, face = 'bold'),
        strip.text.x = element_text(size = 8, face = 'bold'))

# method 1 - inclusion of stand variables
# input database with DBH, height and stand variables
AASVGM <- read.csv("enAASVGM.csv") # integrated original fitting database, stand variables and parameters of A. araucana
GMAA1p7 <- nls(Hoehe ~ 1.3 + (a0 + a1 * HdAA + a2 * HDmaxAA) * exp(-(b0 + b1 * RangeHDAA) * exp(-(c0 + c1 * HdAA) * BHD)),
               start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, a2 = 0, b1 = 0, c1 = 0),
               data = AASVGM)
overview(GMAA1p7)
# parameter a0 and b1 is not significant
# try without b1 first
GMAA1p7_2 <- nls(Hoehe ~ 1.3 + (a0+ a1 * HdAA + a2 * HDmaxAA) * exp(-b0 * exp(-(c0 + c1 * HdAA) * BHD)),
                 start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, a2 = 0, c1 = 0),
                 data = AASVGM)
overview(GMAA1p7_2)
mpc(GMAA1p7_2, AASVGM)

# method 2 - 2-phase approach, with  best subsets regression in the second phase
# selection stand variables that will be used in method 2
AASV2a <- AASVGM[, c(8, 9, 12, 13, 16, 17, 20, 22 : 37, 38, 39, 40, 41, 50, 47)]
AASV2b <- AASVGM[, c(8, 9, 12, 13, 16, 17, 20, 22 : 37, 38, 39, 40, 41, 50, 48)]
AASV2c <- AASVGM[, c(8, 9, 12, 13, 16, 17, 20, 22 : 37, 38, 39, 40, 41, 50, 49)]
AAregsuba <- regsubsets(a ~., data = AASV2a, nbest = 1, nvmax = 4)
AAregsubb <- regsubsets(b ~., data = AASV2b, nbest = 1, nvmax = 4)
AAregsubc <- regsubsets(c ~., data = AASV2c, nbest = 1, nvmax = 4)

# save .eps files for revising the figure again
setEPS()
postscript('aaa.eps', width = 8.0, height = 4.0)
plot(AAregsuba, main = 'Parameter a', scale = 'adjr2')
dev.off()

setEPS()
postscript('aab.eps', width = 8.0, height = 4.0)
plot(AAregsubb, main = 'Parameter b', scale = 'adjr2')
dev.off()

setEPS()
postscript('aac.eps', width = 8.0, height = 4.0)
plot(AAregsubc, main = 'Parameter c', scale = 'adjr2')
dev.off()

# inclusion of stand variable - method 2
GMAA2p6 <- nls(Hoehe ~ 1.3 + (a0 + a1 * HdomAA) * exp(-(b0 + b1 * RangeHDAA + b2 * DBHdomAA) * exp(-(c0) * BHD)),
               start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, b1 = 0, b2 = 0),
               data = AASVGM)
summary(GMAA2p6)
# parameter b0 and b1 is not significant
# try without b1 first
GMAA2p6_2 <- nls(Hoehe ~ 1.3 + (a0 + a1 * HdomAA) * exp(-(b0 + b2 * DBHdomAA) * exp(-(c0) * BHD)),
                 start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, b2 = 0),
                 data = AASVGM)
summary(GMAA2p6_2)
mpc(GMAA2p6_2, AASVGM)

# method 3 - try all combinations by looping through
# selection stand variables that will be used in method 3
AASV3 <- AASVGM[, c(8, 9, 12, 13, 16, 17, 20, 22 : 37, 40, 41, 50)]

# loop
bestresultrmse <- 5 # set an initial value of RMSE
bestresultr2 <- 0.5 # set an initial value of R2
bestresultmae <- 4 # set an initial value of MAE
bestcomb <- c() # the best combination of stand variables will be stored
cout <- 0 # calulate how many times the loop has been run

for (a in 0 : 2) {
  a_index <- combn(16, a)
  for (i in 1 : dim(a_index)[2]) {
    a_var <- list(0, 0, 0, 0, 0)
    if(a != 0) {
      for (idx in 1 : a) {
        a_var[idx] <- AASV3[a_index[, i][idx]]
      }
    }
    
    for (b in 0 : (2 - a)) {
      b_index <- combn(26, b)
      for (j in 1 : dim(b_index)[2]) {
        b_var <- list(0, 0, 0, 0, 0)
        if(b != 0) {
          for (idx in 1 : b) {
            b_var[a + idx] <- AASV3[b_index[, j][idx]]
          }
        }
        
        c <- 2 - a - b
        c_index <- combn(26, c)
        for (k in 1 : dim(c_index)[2]) {
          cout <- cout + 1
          c_var <- list(0, 0, 0, 0, 0)
          if(c != 0) {
            for (idx in 1 : c) {
              c_var[a + b + idx] <- AASV3[c_index[, k][idx]]
            }
          }
          
          modelaa <- try(nls(Hoehe ~ 1.3 + (a0 + a1 * a_var[[1]] + a2 * a_var[[2]]) * exp(-(b0 + a1 * b_var[[1]] + a2 * b_var[[2]]) * exp(-(c0 + a1 * c_var[[1]] + a2 * c_var[[2]]) * BHD)),
                             start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, a2 = 0),
                             data = AASVGM), silent = TRUE)
          
          if('try-error' %in% class(modelaa)) next
          
          res <- mpc(modelaa, AASVGM)
          RMSE <- res$RMSE
          R2 <- res$Rsquare
          MAE <- res$MAE
          
          if((RMSE < bestresultrmse) && (R2 > bestresultr2) && (MAE < bestresultmae)) {
            bestresultrmse <- RMSE
            bestresultr2 <- R2
            bestresultmae <- MAE
            bestcomb <- c(length(a_index[, i]), length(b_index[, j]), length(c_index[, k]),
                          a_index[, i], b_index[, j], c_index[, k])
          }
        }
      }
    }
  }
}

# print result of method 3
cat('bestcomb: ', bestcomb, '\n')

# get stand variables
colnames(AASV3[13])
colnames(AASV3[7])

# inclusion of stand variable - method 3
GMAA3p5 <- nls(Hoehe ~ 1.3 + (a0 + a1 * HdomAA) * exp(-(b0 + b1 * RangeDBHAA) * exp(-(c0) * BHD)),
               start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, b1 = 0),
               data = AASVGM)
summary(GMAA3p5)
mpc(GMAA3p5, AASVGM)

# the process of inclusion stand variables is exactly the same for N. pumilio, so it will not be repeated here, just using A. araucana as an example to show the code


# the effects of the stand variables included ----

# setting the range of DBH
x <- c(0 : 250)

# simulation will based on the maximum and minimum values of HdomAA in different stands, which are assumed to be 18, 23, 28 and 33, RangeDBHAA will use the average value of 168.5 cm
yHdomAA18 <- 1.3 + (17.231019 + 0.330946 * 18) * exp(-(0.886679 + 0.011143 * 168.5) * exp(-(0.026182) * x))
yHdomAA23 <- 1.3 + (17.231019 + 0.330946 * 23) * exp(-(0.886679 + 0.011143 * 168.5) * exp(-(0.026182) * x))
yHdomAA28 <- 1.3 + (17.231019 + 0.330946 * 28) * exp(-(0.886679 + 0.011143 * 168.5) * exp(-(0.026182) * x))
yHdomAA33 <- 1.3 + (17.231019 + 0.330946 * 33) * exp(-(0.886679 + 0.011143 * 168.5) * exp(-(0.026182) * x))

# simulation will based on the maximum and minimum values of RangeDBHAA in different stands, which are assumed to be 122, 137, 152, 167, 182, 197 and 212, HdomAA will use the average value of 25.17 m
yRangeDBHAA122 <- 1.3 + (17.231019 + 0.330946 * 25.17) * exp(-(0.886679 + 0.011143 * 122) * exp(-(0.026182) * x))
yRangeDBHAA137 <- 1.3 + (17.231019 + 0.330946 * 25.17) * exp(-(0.886679 + 0.011143 * 137) * exp(-(0.026182) * x))
yRangeDBHAA152 <- 1.3 + (17.231019 + 0.330946 * 25.17) * exp(-(0.886679 + 0.011143 * 152) * exp(-(0.026182) * x))
yRangeDBHAA167 <- 1.3 + (17.231019 + 0.330946 * 25.17) * exp(-(0.886679 + 0.011143 * 167) * exp(-(0.026182) * x))
yRangeDBHAA182 <- 1.3 + (17.231019 + 0.330946 * 25.17) * exp(-(0.886679 + 0.011143 * 182) * exp(-(0.026182) * x))
yRangeDBHAA197 <- 1.3 + (17.231019 + 0.330946 * 25.17) * exp(-(0.886679 + 0.011143 * 197) * exp(-(0.026182) * x))
yRangeDBHAA212 <- 1.3 + (17.231019 + 0.330946 * 25.17) * exp(-(0.886679 + 0.011143 * 212) * exp(-(0.026182) * x))

dfaasv <- data.frame(x, yHdomAA18, yHdomAA23, yHdomAA28, yHdomAA33,
                     yRangeDBHAA122, yRangeDBHAA137, yRangeDBHAA152, yRangeDBHAA167, yRangeDBHAA182, yRangeDBHAA197, yRangeDBHAA212)

# visualisation
fig.4_HdomAA <- ggplot(dfaasv) +
  geom_line(aes(x = x, y = yHdomAA18, linetype = '18 m')) +
  geom_line(aes(x = x, y = yHdomAA23, linetype = '23 m')) +
  geom_line(aes(x = x, y = yHdomAA28, linetype = '28 m')) +
  geom_line(aes(x = x, y = yHdomAA33, linetype = '33 m')) +
  xlab('DBH (cm)') + ylab('Height (m)') +
  scale_linetype_discrete(name = 'HdomAA') +
  theme(legend.direction = 'horizontal',
        legend.position = c(0.6, 0.1),
        legend.background = element_blank(),
        legend.key = element_blank())

fig.4_RangeDBHAA <- ggplot(dfaasv) +
  geom_line(aes(x = x, y = yRangeDBHAA122, linetype = '122 cm')) +
  geom_line(aes(x = x, y = yRangeDBHAA137, linetype = '137 cm')) +
  geom_line(aes(x = x, y = yRangeDBHAA152, linetype = '152 cm')) +
  geom_line(aes(x = x, y = yRangeDBHAA167, linetype = '167 cm')) +
  geom_line(aes(x = x, y = yRangeDBHAA182, linetype = '182 cm')) +
  geom_line(aes(x = x, y = yRangeDBHAA197, linetype = '197 cm')) +
  geom_line(aes(x = x, y = yRangeDBHAA212, linetype = '212 cm')) +
  xlab('DBH (cm)') + ylab('Height (m)') +
  guides(linetype = guide_legend(ncol = 3)) +
  scale_linetype_discrete(name = 'RangeDBHAA') +
  theme(legend.direction = 'horizontal',
        legend.position = c(0.6, 0.15),
        legend.background = element_blank(),
        legend.key = element_blank())

fig.4_aasveffects <- plot_grid(fig.4_HdomAA, fig.4_RangeDBHAA)

# the process of visualizing the effect of the stand variables included by N. pumilio model on the H-D relationship is also the same


# comparison of H-D models among different sites ----

# also only use A. araucana as an example

# add dummy variables according to location
AAdum <- dummy_cols(AASVGM, select_columns = 'Standort')
# for 4 different sites only 3 dummy variables are needed
AAdum3 <- AAdum[, -52]

# fit the full and reduced models
AAreduce <- nls(Hoehe ~ 1.3 + (a0 + a1 * HdomAA) * exp(-(b0 + b1*RangeDBHAA) * exp(-(c0) * BHD)),
                start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, b1 = 0),
                data = AASVGM)
overview(AAreduce) # get data for F test and Lakkis and Jones test

AAfull <- nls(Hoehe ~ 1.3 + (a0 + a1 * HdomAA + m1 * Standort_Conguillio + m2 * `Standort_La Fusta` + m3 * Standort_Malalcahuello) * exp(-(b0 + b1 * RangeDBHAA + n1 * Standort_Conguillio + n2 * `Standort_La Fusta` + n3 * Standort_Malalcahuello) * exp(-(c0 + p1 * Standort_Conguillio + p2 * `Standort_La Fusta` + p3 * Standort_Malalcahuello) * BHD)),
              start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, b1 = 0, m1 = 0, m2 = 0, m3 = 0, n1 = 0, n2 = 0, n3 = 0, p1 = 0, p2 = 0, p3 = 0),
              data = AAdum3)
overview(AAfull)

# calculation of the test criteria
qf(0.95, 387 - 378, 378) # nonlinear extra sum of squares
qchisq(0.95, df = 387 - 378) # Lakkis-Jones test

# residual diagnosis
plot(AAreduce) # show a clear pattern
# WLS
AAreducew <- nls(Hoehe ~ 1.3 + (a0 + a1 * HdomAA) * exp(-(b0 + b1 * RangeDBHAA) * exp(-(c0) * BHD)),
                 weights = 1 / BHD, data = AASVGM, 
                 start = list(a0 = 46.410000000, b0 = 2.174185939, c0 = 0.008847333, a1 = 0, b1 = 0))
plot(AAreducew) # distribution is more even than when unweighted
# visualisation of results
plotaareduce <- plot(AAreduce, main = 'unweighted', col = 'black', pch = 20, cex = 0.7,
                     xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
plotaareducew <- plot(AAreducew, main = expression( w[i] == 1 / DBH[i]), col = 'black', pch = 20, cex = 0.7,
                      xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
fig.a3_gmresi <- plot_grid(plotaareduce, plotaareducew)



# NLME models ----
# A. araucana - simple NLME model
# try all possibilities of 'random = ' and choose the model with the lowest AIC as the simple NLME model for A. araucana
# calculate how many combinations firsr
choose(3, 1) + choose(3, 2) + choose(3, 3) # total 7 combinations for 3 parameters

# fit 7 NLME models
BMAAmix_1 <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                  start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                  data = fitting_AA, fixed = a + b + c ~ 1, random = a + b + c ~ 1 | Stand)

BMAAmix_2 <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                  start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                  data = fitting_AA, fixed = a + b + c ~ 1, random = a + b ~ 1 | Stand)

BMAAmix_3 <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                  start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                  data = fitting_AA, fixed = a + b + c ~ 1, random = a + c ~ 1 | Stand)

BMAAmix_4 <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                  start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                  data = fitting_AA, fixed = a + b + c ~ 1, random = b + c ~ 1 | Stand)

BMAAmix_5 <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                  start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                  data = fitting_AA, fixed = a + b + c ~ 1, random = a ~ 1 | Stand)

BMAAmix_6 <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                  start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                  data = fitting_AA, fixed = a + b + c ~ 1, random = b ~ 1 | Stand)

BMAAmix_7 <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                  start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                  data = fitting_AA, fixed = a + b + c ~ 1, random = c ~ 1 | Stand)

AIC(BMAAmix_1, BMAAmix_2, BMAAmix_3, BMAAmix_4, BMAAmix_5, BMAAmix_6, BMAAmix_7) # BMAAmix_5 has the lowest AIC among 7 models

# and there is a nested relationship between BMAAmix_5 and BMAAmix_1 / BMAAmix_2 / BMAAmix_3, so the LRT test is required
# BMAAmix_1 and BMAAmix_5 as the example
performance::performance_lrt(BMAAmix_1,BMAAmix_5)

# residual diagnosis
plot(BMAAmix_5) # clear pattern exists
# try different variance functions

BMAAmix_5_power <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                        start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                        data = fitting_AA, fixed = a + b + c ~ 1, random = a ~ 1 | Stand,
                        weights = varPower(form = ~ d))

BMAAmix_5_exp <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                        start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                        data = fitting_AA, fixed = a + b + c ~ 1, random = a ~ 1 | Stand,
                        weights = varExp(form = ~ d))

BMAAmix_5_constpower <- nlme(h ~ 1.3 + a * exp(-b * exp(-c * d)),
                      start = c(a = 46.410000000, b = 2.174185939, c = 0.008847333),
                      data = fitting_AA, fixed = a + b + c ~ 1, random = a ~ 1 | Stand,
                      weights = varConstPower(form = ~ d))

# check AIC and plot
AIC(BMAAmix_5_power, BMAAmix_5_exp, BMAAmix_5_constpower) # BMAAmix_5_power has the lowest AIC

plot(BMAAmix_5_power)
plot(BMAAmix_5_exp)
plot(BMAAmix_5_constpower)
# distribution of the rasiduals is about even across the three plots, so BMAAmix_5_power with the lowest AIC is used

# visualization
plotbmaamix <- plot(BMAAmix_5, main ='unweighted', col = 'black', pch = 20, cex = 0.7,
                    xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
plotbmaamix_power <- plot(BMAAmix_5_power, main ='power variance function', col = 'black', pch = 20, cex = 0.7,
                          xlab = 'Predicted height (m)', ylab = 'Studentized residuals')
fig.a4_nlmebmaaresi <- plot_grid(plotbmaamix, plotbmaamix_power)

# the process of developing a generalized NLME model for A. araucana is the same and therefore will not be repeated
# processes of developing NLME models for N. pumilio will not be repeated either


# calibration and height prediction ----

# see also the supplementary material entitled 'Calibration for Araucaria araucana.pdf'
# calibration simple NLME model of A. araucana in Bestand 2 as the example
calibrating_AA2 <- calibrating_AA[calibrating_AA$Stand == '2', ]

# calibration
.G <- VarCorr(BMAAmix_5_power, rdig = 7) # read variance and correlation components
.estimate <- as.numeric(.G[, 1])
.cov_parm <- c('var_u', 'Residual')
.CEP <- data.frame(.cov_parm, .estimate)
var_e <- .CEP[2, 2] # sigma^2:the scaling factor, given by the value of the residual variance of the model
D <- .CEP[1, 2] # D: the structure of the variance-covariance matrix among stands, since there is only one random parameter in SMAA, D is the variance of this parameter
parms <- data.frame(t(BMAAmix_5_power$coefficients$fixed)) # estimated values for the fixed parameters of SMAA
rho <- coef(BMAAmix_5_power$modelStruct$varStruct, un = FALSE) # the coefficient of the variance function

output <- c()
print <- c()
for(i in 1 : 8){
  for(q in 1 : 1000){
    plot <- sample_n(calibrating_AA2, i, replace=FALSE) # selection of i (i=1~8) trees from stand 2 as sample trees
    plotother <- setdiff(calibrating_AA2, plot) # trees in stand 2 that are not sample trees
    matcor <- diag((plot$d) ^ (2 * rho), ncol = i,nrow = i) # Gi (diagonal matrix describing the non-constant variance)
    r = var_e * matcor # Ri  (variance-covariance matrix within-stand)
    y <- as.matrix(plot$h) # observations of tree heights
    Z <- as.matrix(exp(- parms$b * exp(- parms$c * plot$d))) # determination of the partial derivatives with respect to random parameters
    fxBb <- as.matrix(1.3 + (parms$a) * exp(- parms$b * exp(- parms$c * plot$d))) # determination of the estimated tree height only with fixed effects
    bi <- D %*% t(Z) %*% solve(r + Z %*% D %*% t(Z)) %*% ((y - fxBb)) # determination of the estimated value of a random parameter
    plotother$pre.h <- 1.3 + (parms$a) * exp(-parms$b * exp(-parms$c * plotother$d)) + bi * (exp(- parms$b * exp(- parms$c * plotother$d))) # prediction of tree height for trees that are not sample trees (with mixed effects)
    bias <- sum(plotother$pre.h - plotother$h) / length(plotother$h)
    biasper <- 100 * (sum((plotother$h - plotother$pre.h) / length(plotother$h))) / mean(plotother$h)
    mape <- sum(abs(plotother$h - plotother$pre.h) / plotother$h)
    mpe <- sum((plotother$h - plotother$pre.h) / plotother$h)
    mspe <- sum((plotother$h - plotother$pre.h) ^ 2)
    output <- rbind(output, cbind(i, bias, biasper, mape, mpe, mspe))}
  
  output <- data.frame(output)
  print <- rbind(print, cbind(mean(output$i), mean(output$bias), mean(output$biasper), mean(output$mape), mean(output$mpe), mean(output$mspe)))
  output <- c()
}

print(print)

# 'print' is the results of the prediction performance of the different calibration designs
# the results are not the final results and need to be recalculated based on the results in stand 12
# if you want your result to be reproducible, you can use set.seed() before the loop

# height prediction

# ONLS
startHDgomperz(calibrating_AA2$d, calibrating_AA2$h)
ONLSAAs2 <- nls(h ~ 1.3 + a * exp(- b * exp(- c * d)),
                start = list(a = 55.510000000, b = 1.795581979, c = 0.007193447),
                data = calibrating_AA2)
mpc(ONLSAAs2, calibrating_AA2)

calibrating_AA2$onlspre.h <- ONLSAAs2$m$predict()
calibrating_AA2$fixbmpre.h <- 1.3 + 25.832936 * exp(-2.622668 * exp(-0.025004 * calibrating_AA2$d))

# since a plot is needed to show the example --> seed
set.seed(2021)
ploteg <- sample_n(calibrating_AA2, 5, replace = FALSE)
plotegother <- setdiff(calibrating_AA2, ploteg)

matcor <- diag((ploteg$d) ^ (2 * rho), ncol = 5,nrow = 5)
r <- var_e*matcor
y <- as.matrix(ploteg$h)
Z <- as.matrix(exp(- parms$b * exp(- parms$c * ploteg$d)))
fxBb <- as.matrix(1.3 + parms$a * exp(- parms$b * exp(- parms$c * ploteg$d)))
bi <- D %*% t(Z) %*% solve(r + Z %*% D %*% t(Z)) %*% ((y - fxBb))
plotegother$pre.h <- 1.3 + parms$a * exp(- parms$b * exp(- parms$c * plotegother$d)) + bi * exp(- parms$b * exp(- parms$c * plotegother$d))

# visualisation
fig.5_aas2bm <- ggplot() +
  geom_line(aes(x = d, y = onlspre.h, linetype = 'ONLS'), data = calibrating_AA2) +
  geom_line(aes(x = d, y = fixbmpre.h, linetype = 'Fixed effects'), data = calibrating_AA2) +
  geom_line(aes(x = d, y = pre.h, linetype = 'Calibrated'), data = plotegother) +
  geom_point(aes(x = d, y = h, shape = 'Date'), data = calibrating_AA2) +
  geom_point(aes(x = d, y = h, shape = 'Calibration data'), data = ploteg) +
  scale_shape_manual(name = '', values = c(4, 1)) +
  scale_linetype_discrete(name = '') +
  theme(legend.position = 'bottom') +
  xlab('DBH (cm)') + ylab ('Height (m)') +
  theme(legend.key = element_rect(fill = 'white'),
        plot.title = element_text()) +
  ggtitle('Prediction of height in stand 2 with simple NLME H-D model')

# the rest of the calibration and tree height prediction process is all the same


# THE END ----