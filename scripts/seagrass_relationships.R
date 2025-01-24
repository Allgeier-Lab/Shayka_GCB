library(arrR)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(MASS)
library(moments)
library(rstatix)
library(plotrix)
library(data.table)
library(mgcv)
library(stringr)

# load nutrient input
FISH_NUTR_INPUT <- read.csv("FISH_NUTR_INPUT_20M.csv")
FISH_NUTR_INPUT <- FISH_NUTR_INPUT[, 2]

#### 1m Seafloor ####

# load parameters
parameters <- arrR::default_parameters
starting_values <- arrR::default_starting

# update parameters according to empirical ranges
shayka_sg_data <- read.csv("bladeareadata.csv")
parameters$ag_biomass_max <- max(shayka_sg_data$agweight)
parameters$bg_biomass_max <- max(shayka_sg_data$bgweight)
parameters$ag_biomass_min <- min(shayka_sg_data$agweight)
parameters$bg_biomass_min <- min(shayka_sg_data$bgweight)

# determine simulation conditions
# one iterations equals 120 minutes
min_per_i <- 120

# run the model for x years
years <- 10
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass once each day
days <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days

# save results at final timestep
days <- 1
save_each <- (24 / (min_per_i / 60)) * days

# create 5 reef cell in center of seafloor
reef_matrix <- matrix(data = c(0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values <- get_req_nutrients(bg_biomass = starting_values$bg_biomass,
                                   ag_biomass = starting_values$ag_biomass,
                                   parameters = parameters)

starting_values$nutrients_pool <- stable_values$nutrients_pool

starting_values$detritus_pool <- stable_values$detritus_pool

# start at minimal seagrass to establish gradient
starting_values$ag_biomass <- parameters$ag_biomass_min
starting_values$bg_biomass <- parameters$bg_biomass_min

# create seafloor, fishpop
# change  starting values and parameters for 1m seagrass conditions
starting_values$pop_n <- 0
parameters$seagrass_slough = 0.75 * arrR::default_parameters$seagrass_slough
parameters$seagrass_slope = arrR::default_parameters$seagrass_slope
parameters$seagrass_thres = -0.175

# create seafloor
input_seafloor <- setup_seafloor(dimensions = c(50, 50), grain = 1,
                                 reef = reef_matrix, starting_values = starting_values)

# create fishpop
input_fishpop <- setup_fishpop(seafloor = input_seafloor, species = rep(1, starting_values$pop_n),
                               starting_values = starting_values,
                               parameters = parameters)

# run model for 1m
result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                               parameters = parameters, movement = "behav",
                               max_i = max_i, min_per_i = min_per_i, nutrients_input = c((FISH_NUTR_INPUT[1:1032] / 3),
                                                                                         (FISH_NUTR_INPUT[1033:1536] / 2.8),
                                                                                         (FISH_NUTR_INPUT[1537:2148] / 2.5),
                                                                                         (FISH_NUTR_INPUT[2149:2796] / 2.25),
                                                                                         (FISH_NUTR_INPUT[2797:max_i] / 1.9)),
                               seagrass_each = seagrass_each, save_each = save_each)

# define reef coords and radius for calculating at reef
reef_x <- 0.5
reef_y <- -0.5
calc_range <- 20

# define near-reef cells (1m distance)
dummy_near_reef_seafloor <- filter(result$seafloor, (
  (x == 0.5 & y == 0.5) | (x == 0.5 & y == -1.5) |
    (x == 1.5 & y == -0.5) | (x == -0.5 & y == -0.5) |
    (x == 1.5 & y == 0.5) | (x == -0.5 & y == -1.5) |
    (x == 1.5 & y == -1.5) | (x == -0.5 & y == 0.5)))
dummy_near_reef_seafloor$distance <- 1
result$seafloor <- dummy_near_reef_seafloor

# create export data frame 'seafloor'
seafloor <- data.frame()
result$seafloor$xy_int <- interaction(result$seafloor$x, result$seafloor$y)
result$seafloor <- result$seafloor[order(result$seafloor$xy_int, result$seafloor$timestep), ]
result$seafloor <- result$seafloor[, c("timestep", "ag_biomass", "bg_biomass", "ag_production",
                                       "bg_production", "distance", "xy_int")]

# calculate difference in production between timesteps
for (xy in unique(result$seafloor$xy_int)){
  temp_result_seafloor <- result$seafloor

  # filter functional result to particular xy pair
  result$seafloor <- filter(temp_result_seafloor, xy_int == xy)

  # calculate diff in production between timesteps for that xy pair
  for (i in seq(max_i, save_each, by = -save_each)) {
    result$seafloor[result$seafloor$timestep ==(i), ]$ag_production =
      result$seafloor[result$seafloor$timestep == i,]$ag_production -
      result$seafloor[result$seafloor$timestep == (i - save_each),]$ag_production

    result$seafloor[result$seafloor$timestep ==(i), ]$bg_production =
      result$seafloor[result$seafloor$timestep == i,]$bg_production -
      result$seafloor[result$seafloor$timestep == (i - save_each),]$bg_production
  }

  # store lagged production result
  seafloor <- rbind(seafloor, result$seafloor)

  # restore original result
  result$seafloor <- temp_result_seafloor
}

seafloor_1m <- seafloor

#### plot model fits AG:BG for 1m ####
# establish graph output and layout parameters
nf <- layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = T))

# plot empirical points
plot(bgweight~agweight, data = filter(shayka_sg_data, Distance == 1),
     ylab = "BG Biomass (g/m2)",
     xlab = "AG Biomass (g/m2)",
     col = "black", cex = 0.5)

# plot model points
points(bg_biomass~ag_biomass, data = filter(seafloor, distance == 1, ag_biomass < 924),
       ylab = "BG Biomass (g/m2)",
       xlab = "AG Biomass (g/m2)",
       col = "red", cex = 0.5)

# generate and plot model of model points
loglm_mod_1m <- lm(bg_biomass~log(ag_biomass), data = filter(seafloor, ag_biomass < 924, distance == 1))
newd <- data.frame(ag_biomass = seq(min(filter(seafloor, distance == 1)$ag_biomass),
                                    max(filter(seafloor, distance == 1)$ag_biomass), by = 1))
Xp <- predict(loglm_mod_1m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(loglm_mod_1m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "pink", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="pink", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="pink", lty=2)

# add legend
ans <- legend(x = grconvertX(0.6, "npc", "user"), y = grconvertY(0.33, "npc", "user"),
              legend = c("Empirical", "Model"),
              pch = c(1, 1), cex = 1, pt.cex = 2,
              y.intersp = 0.6, x.intersp = 0.83, ncol = 2, xpd = NA,
              text.width = c(0, 0), plot = F)
# place legend on plot
r <- ans$rect
legend(x = grconvertX(0.575, "npc", "user"), y = grconvertY(0.4, "npc", "user"),
       legend = c("Empirical", "Model"),
       col = c("black", "red"),
       pch = c(1, 1),
       y.intersp = 0.75, x.intersp = 0.5, ncol = 1, bty = "n", xpd = NA,
       text.width = 0)
rect(r$left, r$top - r$h, r$left + r$w, r$top)
title(main = "1 meter AG/BG Regressions")

# generate and plot model of empirical points
loglm_emp_1m <- lm(bgweight~log(agweight), data = filter(shayka_sg_data, Distance == 1))
newd <- data.frame(agweight = seq(min(filter(shayka_sg_data, Distance == 1)$agweight), max(filter(shayka_sg_data, Distance == 1)$agweight), by = 1))
Xp <- predict(loglm_emp_1m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(loglm_emp_1m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "gray", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="gray", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="gray", lty=2)

# print summaries of model fit
print("EMP LM R2")
print(summary(loglm_emp_1m)$r.squared)
print("MOD LM R2")
print(summary(loglm_mod_1m)$r.squared)

#### plot model fits AG:Prod 1m ####
# plot range of modeled points
plot(ag_production~ag_biomass, data = filter(seafloor,ag_biomass < 924, distance == 1),
     ylab = "AG Production (g/m2/day)",
     xlab = "AG Biomass (g/m2)",
     col = "red", cex = 0.5)

# generate and plot model of modeled points
lm_mod_1m <- lm(ag_production~(ag_biomass), data = filter(seafloor, distance == 1, ag_biomass < 924))
newd <- data.frame(ag_biomass = seq(min(filter(seafloor, distance == 1)$ag_biomass), max(filter(seafloor, distance == 1)$ag_biomass), by = 1))
Xp <- predict(lm_mod_1m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(lm_mod_1m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "pink", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="pink", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="pink", lty=2)

# plot empirical points
points(coreprodweight~agweight, data = filter(shayka_sg_data, Distance == 1),
       ylab = "AG Production (g/m2/day)",
       xlab = "AG Biomass (g/m2)",
       col = "black", cex = 0.5)
title(main = "1 meter AG/Prod Regressions")

# generate and plot model of empirical points
lm_emp_1m <- lm(coreprodweight~(agweight), data = filter(shayka_sg_data, Distance == 1))
newd <- data.frame(agweight = seq(min(filter(shayka_sg_data, Distance == 1)$agweight), max(filter(shayka_sg_data, Distance == 1)$agweight), by = 1))
Xp <- predict(lm_emp_1m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(lm_emp_1m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "gray", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="gray", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="gray", lty=2)

# print summaries of models
print("EMP LM R2")
print(summary(loglm_emp_1m)$r.squared)
print("MOD LM R2")
print(summary(loglm_mod_1m)$r.squared)


#### 20m Seafloor ####

# load nutrient input csv
FISH_NUTR_INPUT <- read.csv("FISH_NUTR_INPUT_20m.csv")
FISH_NUTR_INPUT <- FISH_NUTR_INPUT[, 2]


# load default parameters
parameters <- arrR::default_parameters
starting_values <- arrR::default_starting

# update parameters with experimental data
shayka_sg_data <- read.csv("bladeareadata.csv")
parameters$ag_biomass_max <- max(shayka_sg_data$agweight)
parameters$bg_biomass_max <- max(shayka_sg_data$bgweight)
parameters$ag_biomass_min <- min(shayka_sg_data$agweight)
parameters$bg_biomass_min <- min(shayka_sg_data$bgweight)

# establish parameters for model run
# one iterations equals 120 minutes
min_per_i <- 120

# run the model for x years
years <- 10
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass once each day
days <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days

# save results at final timestep
days <- 1
save_each <- (24 / (min_per_i / 60)) * days

# create 5 reef cell in center of seafloor
reef_matrix <- matrix(data = c(0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values <- get_req_nutrients(bg_biomass = starting_values$bg_biomass,
                                   ag_biomass = starting_values$ag_biomass,
                                   parameters = parameters)

starting_values$nutrients_pool <- stable_values$nutrients_pool

starting_values$detritus_pool <- stable_values$detritus_pool

# start at minimal seagrass to establish gradient
starting_values$ag_biomass <- parameters$ag_biomass_min
starting_values$bg_biomass <- parameters$bg_biomass_min

# modify parameters to model 20m seagrass dynamics
starting_values$pop_n <- 0
parameters$seagrass_slough = 0.05 * arrR::default_parameters$seagrass_slough
parameters$seagrass_slope = arrR::default_parameters$seagrass_slope
parameters$seagrass_thres = -0.4

# create seafloor
input_seafloor <- setup_seafloor(dimensions = c(50, 50), grain = 1,
                                 reef = reef_matrix, starting_values = starting_values)

# create fishpop
input_fishpop <- setup_fishpop(seafloor = input_seafloor, species = rep(1, starting_values$pop_n),
                               starting_values = starting_values,
                               parameters = parameters)
# run the model
result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                               parameters = parameters, movement = "behav",
                               max_i = max_i, min_per_i = min_per_i, nutrients_input = c((FISH_NUTR_INPUT[1:2844] / 12),
                                                                                         (FISH_NUTR_INPUT[2845: 4272] / 11.5),
                                                                                         (FISH_NUTR_INPUT[4273: 4908] / 10.5),
                                                                                         (FISH_NUTR_INPUT[4909: 5220] / 9.75),
                                                                                         (FISH_NUTR_INPUT[5221: 5712] / 9),
                                                                                         (FISH_NUTR_INPUT[5713: 6120] / 8.25),
                                                                                         (FISH_NUTR_INPUT[6121: max_i] / 7.5)),# 6500),
                               seagrass_each = seagrass_each, save_each = save_each)

# define reef coords and radius for calculating at reef
reef_x <- 0.5
reef_y <- -0.5
calc_range <- 20

# calculate primary production at 20m using circle (also observe 1m dynamics
# under 20m conditions)
dummy_result <- result
dummy_result$seafloor$x <- dummy_result$seafloor$x - reef_x
dummy_result$seafloor$y <- dummy_result$seafloor$y - reef_y
dummy_result$seafloor <- dummy_result$seafloor[apply(dummy_result$seafloor[1:2]^2,1,sum) <=
                                                 (calc_range) ^ 2,]
dummy_result$seafloor <- dummy_result$seafloor[apply(dummy_result$seafloor[1:2]^2,1,sum) >
                                                 (calc_range - 1) ^ 2,]
# select near-reef cells
dummy_near_reef_seafloor <- filter(result$seafloor, (
  (x == 0.5 & y == 0.5) | (x == 0.5 & y == -1.5) |
    (x == 1.5 & y == -0.5) | (x == -0.5 & y == -0.5) |
    (x == 1.5 & y == 0.5) | (x == -0.5 & y == -1.5) |
    (x == 1.5 & y == -1.5) | (x == -0.5 & y == 0.5)))
dummy_near_reef_seafloor$distance <- 1
dummy_result$seafloor$distance <- 20
dummy_result$seafloor$xy_int <- interaction(dummy_result$seafloor$x, dummy_result$seafloor$y)
xy_int_list_20m <- as.character(unique(dummy_result$seafloor$xy_int)[1:9])
dummy_result$seafloor <- filter(dummy_result$seafloor, str_detect(xy_int, paste(xy_int_list_20m, collapse="|")))

dummy_near_reef_seafloor$xy_int <- interaction(dummy_near_reef_seafloor$x, dummy_near_reef_seafloor$y)
result$seafloor <- rbind(dummy_near_reef_seafloor, dummy_result$seafloor)

# define export seafloor
seafloor <- data.frame()
result$seafloor <- result$seafloor[order(result$seafloor$xy_int, result$seafloor$timestep), ]
result$seafloor <- result$seafloor[, c("timestep", "ag_biomass", "bg_biomass", "ag_production",
                                       "bg_production", "distance", "xy_int")]

# loop through all cells and calculate daily production
for (xy in unique(result$seafloor$xy_int)) {
  temp_result_seafloor <- result$seafloor

  # filter functional result to particular xy pair
  result$seafloor <- filter(temp_result_seafloor, xy_int == xy)

  # calculate diff in prod between timesteps for that xy pair
  for (i in seq(max_i, save_each, by = -save_each)) {
    result$seafloor[result$seafloor$timestep ==(i), ]$ag_production =
      result$seafloor[result$seafloor$timestep == i,]$ag_production -
      result$seafloor[result$seafloor$timestep == (i - save_each),]$ag_production

    result$seafloor[result$seafloor$timestep ==(i), ]$bg_production =
      result$seafloor[result$seafloor$timestep == i,]$bg_production -
      result$seafloor[result$seafloor$timestep == (i - save_each),]$bg_production
  }

  # store lagged prod result
  seafloor <- rbind(seafloor, result$seafloor)

  # restore original seafloor
  result$seafloor <- temp_result_seafloor
}

seafloor_20m <- seafloor

#### plot model fits for AG:BG 20m ####
# plot model data
plot(bgweight~agweight, data = filter(shayka_sg_data, Distance == 20), col = "black",
     ylab = "BG Biomass (g/m2)",
     xlab = "AG Biomass (g/m2)", cex = 0.5)
# plot empirical data
points(bg_biomass~ag_biomass, data = filter(seafloor, distance == 20),
       col = "red")

# generate and plot model of model data
loglm_mod_20m <- lm(bg_biomass~log(ag_biomass), data = filter(seafloor, distance == 20, ag_biomass < 922))
newd <- data.frame(ag_biomass = seq(min(filter(seafloor, distance == 20, ag_biomass < 922)$ag_biomass), max(filter(seafloor, distance == 20, ag_biomass < 922)$ag_biomass), by = 1))
Xp <- predict(loglm_mod_20m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(loglm_mod_20m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "pink", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="pink", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="pink", lty=2)
title(main = "20 meter AG/BG Regressions")

# generate and plot model of empirical data
loglm_emp_20m <- lm(bgweight~log(agweight), data = filter(shayka_sg_data, Distance == 20))
newd <- data.frame(agweight = seq(min(filter(shayka_sg_data, Distance == 20)$agweight), max(filter(shayka_sg_data, Distance == 20)$agweight), by = 1))
Xp <- predict(loglm_emp_20m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(loglm_emp_20m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "gray", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="gray", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="gray", lty=2)

# print summaries of models
print("EMP LM R2")
print(summary(loglm_emp_20m)$r.squared)
print("MOD LM R2")
print(summary(loglm_mod_20m)$r.squared)


#### plot model fits AG:Production 20m ####
# plot model data
plot(ag_production~ag_biomass, data = filter(seafloor, distance == 20, ag_biomass < 922),
     ylab = "AG Production (g/m2/day)",
     xlab = "AG Biomass (g/m2)",
     col = "red", cex = 0.5)

# generate and plot model of model data
lm_mod_20m <- lm(ag_production~ag_biomass, data = filter(seafloor, distance == 20, ag_biomass < 922))
newd <- data.frame(ag_biomass = seq(min(filter(seafloor, distance == 20, ag_biomass < 922)$ag_biomass),
                                    max(filter(seafloor, distance == 20, ag_biomass < 922)$ag_biomass), by = 1))
Xp <- predict(lm_mod_20m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(lm_mod_20m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "pink", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="pink", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="pink", lty=2)

# plot empirical data
points(coreprodweight~agweight, data = filter(shayka_sg_data, Distance == 20),
       col = "black", cex = 0.5)
title(main = "20 meter AG/Prod Regressions")

# plot model of empirical data
lm_emp_20m <- lm(coreprodweight~(agweight), data = filter(shayka_sg_data, Distance == 20))
newd <- data.frame(agweight = seq(min(filter(shayka_sg_data, Distance == 20)$agweight), max(filter(shayka_sg_data, Distance == 20)$agweight), by = 1))
Xp <- predict(lm_emp_20m, newdata = newd, interval = "confidence", level = 0.95)
betas <- coef(lm_emp_20m)
points(seq(1, nrow(Xp)), Xp[, 1], col = "gray", cex = 0.5)
lines(seq(1, nrow(Xp)), Xp[, 2], col="gray", lty=2)
lines(seq(1, nrow(Xp)), Xp[,3], col="gray", lty=2)

# print summaries of model fits
print("EMP LM R2")
print(summary(loglm_emp_20m)$r.squared)
print("MOD LM R2")
print(summary(loglm_mod_20m)$r.squared)

# write and save the seafloors for 1m and 20m
write.csv(seafloor_1m, "seafloor_1m.csv")
write.csv(seafloor_20m, "seafloor_20m.csv")
