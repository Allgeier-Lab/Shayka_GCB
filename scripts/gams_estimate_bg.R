library(arrR)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(MASS)
library(moments)
library(rstatix)
library(plotrix)
library(data.table)
library(mgcv) # GAM, predict.gam

# read in data if not already present in R environment
seafloor_20m <- read.csv("seafloor_20m.csv") # model output for seafloor 20m conditions
seafloor_1m <- read.csv("seafloor_1m.csv") # model output for seafloor 1m conditions
shayka_sg_data <- read.csv("bladeareadata.csv") # empirical data collected for 20m and 1m conditions

# filter model data to relevant distance
seafloor_20m <- filter(seafloor_20m, timestep > 0, distance == 20)
seafloor_1m <- filter(seafloor_1m, timestep > 0, distance == 1)

# create column of predicted values in empirical data
shayka_sg_data$bgprod <- 0

# rename bgweight column to bg_biomass for GAM model syntax (change back later)
shayka_sg_data <- shayka_sg_data %>% rename(bg_biomass=bgweight)

# plot BG/Prod Regression 20m
plot((bg_production)~bg_biomass, data = seafloor_20m)
title(main = "20m BG Biomass ~ Production")

# run 20m models
gam_mod_20m <- gam((bg_production)~s((bg_biomass), k = 10), data = seafloor_20m)
k.check(gam_mod_20m)

# predict over empirical dataframe
shayka_sg_data[shayka_sg_data$Distance == 20, ]$bgprod <- predict.gam(gam_mod_20m,
                                                                   newdata = shayka_sg_data[shayka_sg_data$Distance == 20, ])

# generate fitted  predictions across simulated range
newd <- data.frame(bg_biomass = seq(from = 10,
                                    to = 1600, length.out = 2500))
newd$pred <- predict.gam(gam_mod_20m, newdata = newd)

# plot fitted predictions
points(pred~bg_biomass, data = newd, col = "blue", cex = 0.5)

# plot predicted points onto existing plot
points(bgprod~bg_biomass, data = shayka_sg_data[shayka_sg_data$Distance == 20, ],
       col = "red", cex = 1, pch = 19)

# plot BG/Prod Regression 1m
plot((bg_production)~bg_biomass, data = seafloor_1m)
title(main = "1m BG Biomass ~ Production")

# run 1m models
gam_mod_1m <- gam((bg_production)~s((bg_biomass), k = 10), data = seafloor_1m)
k.check(gam_mod_1m)

# predict over dataframe
shayka_sg_data[shayka_sg_data$Distance == 1, ]$bgprod <- predict.gam(gam_mod_1m,
                                                                      newdata = shayka_sg_data[shayka_sg_data$Distance == 1, ])

# generate fitted  predictions across simulated range
newd <- data.frame(bg_biomass = seq(from = 10,
                                    to = 1600, length.out = 2500))
newd$pred <- predict.gam(gam_mod_1m, newdata = newd)

# plot fitted predictions
points(pred~bg_biomass, data = newd, col = "blue", cex = 0.5)

# plot predicted points onto existing plot
points(bgprod~bg_biomass, data = shayka_sg_data[shayka_sg_data$Distance == 1, ],
       col = "red", cex = 1, pch = 19)

# restore original column names
shayka_sg_data <- shayka_sg_data %>% rename(bgweight = bg_biomass)

# return seagrass data with estimated BG production
write.csv(shayka_sg_data, "bladeareaweight_BG.csv")
