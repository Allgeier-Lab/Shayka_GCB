###########
#This code runs 3 SEMs - near reef fert, near reef no fert, and away from reef - to determine what drives productivity
#run in R version 4.0.0
##########


##Load libraries -------------
library(tidyverse)
library(piecewiseSEM)
library(EnvStats) #for qqPlot()
library(car) #for vif

##Load data -----------------

databydist <- read_csv('processeddata/alldistdata.csv',
                       col_types = cols(Block = "c",))


##Prep data for use in SEMs -----------------

bydist1mfishfert <- databydist %>%
  filter(Distance == 1 & fert == "yes") %>%
  select(Reef, Block, Treatment, Transect, Distance, fish, fert,
         coreprodweight, corebladeweightm2, sheathsm2, rhizomesm2, rootsm2,
         shootsm2, totepweightperarea, corebitesperarea,
         meandC_Blades, meanC_Blades, meanN_Blades, meanP_Blades,
         meandC_Sheaths, meanC_Sheaths, meanN_Sheaths, meanP_Sheaths,
         meandC_Rhizomes, meanC_Rhizomes, meanN_Rhizomes, meanP_Rhizomes,
         meandC_Roots, meanC_Roots, meanN_Roots, meanP_Roots) %>%
  ungroup()

bydist1mfish <- databydist %>%
  filter(Distance == 1 & fert == "no") %>%
  select(Reef, Block, Treatment, Transect, Distance, fish, fert,
         coreprodweight, corebladeweightm2, sheathsm2, rhizomesm2, rootsm2,
         shootsm2, totepweightperarea, corebitesperarea,
         meandC_Blades, meanC_Blades, meanN_Blades, meanP_Blades,
         meandC_Sheaths, meanC_Sheaths, meanN_Sheaths, meanP_Sheaths,
         meandC_Rhizomes, meanC_Rhizomes, meanN_Rhizomes, meanP_Rhizomes,
         meandC_Roots, meanC_Roots, meanN_Roots, meanP_Roots) %>%
  ungroup()

bydist20m <- databydist %>%
  filter(Distance == 20) %>%
  select(Reef, Block, Treatment, Transect, Distance, fish, fert,
         coreprodweight, corebladeweightm2, sheathsm2, rhizomesm2, rootsm2,
         shootsm2, totepweightperarea, corebitesperarea,
         meandC_Blades,  meanC_Blades, meanN_Blades, meanP_Blades,
         meandC_Sheaths, meanC_Sheaths, meanN_Sheaths, meanP_Sheaths,
         meandC_Rhizomes, meanC_Rhizomes, meanN_Rhizomes, meanP_Rhizomes,
         meandC_Roots, meanC_Roots, meanN_Roots, meanP_Roots) %>%
  ungroup()


sem1mfert <- as.data.frame(bydist1mfishfert)
sem1mfish <- as.data.frame(bydist1mfish)
sem20m <- as.data.frame(bydist20m)


#check for normality #1m fish + fert
par(mfrow = c(4,4))
for(i in 8:length(sem1mfert)){
  hist(sem1mfert[,i], main = names(sem1mfert[i]))
  qqPlot(sem1mfert[,i], main = names(sem1mfert[i]), pch=19, cex=1)
  mtext(round(shapiro.test(sem1mfert[,i])$p.value,3), side=3)
} #

par(mfrow = c(1,1))
#transformations
sem1mfert$corebladeweightm2 = log(sem1mfert$corebladeweightm2+10) 
sem1mfert$sheathsm2 = (sem1mfert$sheathsm2)^(1/8) 
sem1mfert$rhizomesm2 = log(sem1mfert$rhizomesm2+10) 
sem1mfert$rootsm2 = log(sem1mfert$rootsm2+10) 
sem1mfert$shootsm2 = log(sem1mfert$shootsm2+1) 
sem1mfert$meanC_Sheaths = log(sem1mfert$meanC_Sheaths) 
sem1mfert$meanP_Sheaths = log(sem1mfert$meanP_Sheaths) 
sem1mfert$meanP_Rhizomes = log(sem1mfert$meanP_Rhizomes+10) 
sem1mfert$meanC_Roots = log(sem1mfert$meanC_Roots+1) 
sem1mfert$meanN_Roots = log(sem1mfert$meanN_Roots+10) 
sem1mfert$meanP_Roots = log(sem1mfert$meanP_Roots) 


#check for normality #1m fish
par(mfrow = c(4,4))
for(i in 8:length(sem1mfish)){
  hist(sem1mfish[,i], main = names(sem1mfish[i]))
  qqPlot(sem1mfish[,i], main = names(sem1mfish[i]), pch=19, cex=1)
  mtext(round(shapiro.test(sem1mfish[,i])$p.value,3), side=3)
} #

par(mfrow = c(1,1))
#transformations
sem1mfish$coreprodweight = log(sem1mfish$coreprodweight) 
sem1mfish$corebladeweightm2 = log(sem1mfish$corebladeweightm2+1) 
sem1mfish$sheathsm2 = log(sem1mfish$sheathsm2+1) 
sem1mfish$rhizomesm2 = log(sem1mfish$rhizomesm2+1) 
sem1mfish$corebitesperarea = sqrt(sem1mfish$corebitesperarea) 
sem1mfish$meandC_Blades = (sem1mfish$meandC_Blades+10)^2 
sem1mfish$meandC_Sheaths = (sem1mfish$meandC_Sheaths+10)^(2) 
sem1mfish$meanC_Sheaths = log(sem1mfish$meanC_Sheaths+10) 
sem1mfish$meanN_Sheaths = sqrt(sem1mfish$meanN_Sheaths+10) 
sem1mfish$meandC_Rhizomes = (sem1mfish$meandC_Rhizomes+10)^(4) 
sem1mfish$meanP_Rhizomes = log(sem1mfish$meanP_Rhizomes) 
sem1mfish$meandC_Roots = (sem1mfish$meandC_Roots+10)^2 
sem1mfish$meanC_Roots = sqrt(sem1mfish$meanC_Roots+10) 
sem1mfish$meanP_Roots = log(sem1mfish$meanP_Roots) 


#check for normality #20m
par(mfrow = c(4,4))
for(i in 8:length(sem20m)){
  hist(sem20m[,i], main = names(sem20m[i]))
  qqPlot(sem20m[,i], main = names(sem20m[i]), pch=19, cex=1)
  mtext(round(shapiro.test(sem20m[,i])$p.value,3), side=3)
} #

par(mfrow = c(1,1))
#sem20m$meanP_Rhizomes[2] <- NA
#transformations
sem20m$coreprodweight = log(sem20m$coreprodweight+1) 
sem20m$corebladeweightm2 = log(sem20m$corebladeweightm2+10) 
sem20m$totepweightperarea = (sem20m$totepweightperarea)^(1/3) 
sem20m$corebitesperarea = (sem20m$corebitesperarea)^(1/3) 
sem20m$meandC_Sheaths = (sem20m$meandC_Sheaths+10)^(2) 
sem20m$meanN_Sheaths = sem20m$meanN_Sheaths^(1/6) 
sem20m$meanP_Sheaths = (sem20m$meanP_Sheaths)^(1/5) 
sem20m$meandC_Rhizomes = (sem20m$meandC_Rhizomes+10)^(3) 
sem20m$meanC_Rhizomes = log(sem20m$meanC_Rhizomes) 
sem20m$meanP_Rhizomes = log(sem20m$meanP_Rhizomes+1) 
sem20m$meandC_Roots = (sem20m$meandC_Roots+10)^3 
sem20m$meanP_Roots = sem20m$meanP_Roots^(1/6) 



##1m fish + fertilizer SEM --------------

prod1lm2 <- lm(coreprodweight ~ corebladeweightm2 + shootsm2 + totepweightperarea + corebitesperarea + meanN_Blades + meanP_Blades, data = sem1mfert) 
blade1lm2 <- lm(corebladeweightm2 ~ totepweightperarea + meandC_Blades + meanC_Blades + meanN_Blades + meanP_Blades, data = sem1mfert) 
sheath1lm2 <- lm(sheathsm2 ~ corebladeweightm2 + meandC_Sheaths + meanC_Sheaths + meanN_Sheaths + meanP_Sheaths, data = sem1mfert) 
rhizome1lm2 <- lm(rhizomesm2 ~ corebladeweightm2 + meandC_Rhizomes + meanC_Rhizomes + meanN_Rhizomes + meanP_Rhizomes, data = sem1mfert) 
root1lm2 <- lm(rootsm2 ~ rhizomesm2 + meandC_Roots + meanC_Roots + meanN_Roots + meanP_Roots, data = sem1mfert) 
eps1lm2 <- lm(totepweightperarea ~ meanN_Blades + meanP_Blades, data = sem1mfert) 
bites1m2 <- lm(corebitesperarea ~ meanN_Blades + meanP_Blades, data = sem1mfert) 


###1m fish + fertilizer SEM V1-----------
sem.fertV1.1 <- psem(
  prod1lm2, blade1lm2, sheath1lm2, rhizome1lm2, root1lm2, eps1lm2, bites1m2,
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfert
)

summary(sem.fertV1.1)
AIC(sem.fertV1.1, aicc = T)
plot(sem.fertV1.1)


sem.fertV1.2 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2 + totepweightperarea + corebitesperarea + meanN_Blades + meanP_Blades, sem1mfert),
  lm(corebladeweightm2 ~ totepweightperarea + meandC_Blades + meanC_Blades + meanN_Blades + meanP_Blades + shootsm2, sem1mfert),
  lm(sheathsm2 ~ corebladeweightm2 + meandC_Sheaths + meanC_Sheaths + meanN_Sheaths + meanP_Sheaths, sem1mfert),
  lm(rhizomesm2 ~ corebladeweightm2 + shootsm2 + sheathsm2 + meandC_Rhizomes + meanC_Rhizomes + meanN_Rhizomes + meanP_Rhizomes, sem1mfert),
  lm(rootsm2 ~ rhizomesm2 + meandC_Roots + meanC_Roots + meanN_Roots + meanP_Roots + meanN_Rhizomes + corebladeweightm2, sem1mfert),
  lm(totepweightperarea ~ meanN_Blades + meanP_Blades, sem1mfert),
  lm(corebitesperarea ~ meanN_Blades + meanP_Blades + meandC_Sheaths + meanN_Sheaths + meandC_Rhizomes + meanN_Roots, sem1mfert),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfert
)

summary(sem.fertV1.2)
AIC(sem.fertV1.2, aicc = T)


sem.fertV1.3 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2, sem1mfert),
  lm(corebladeweightm2 ~ meandC_Blades + shootsm2, sem1mfert),
  lm(sheathsm2 ~ corebladeweightm2, sem1mfert),
  lm(rhizomesm2 ~ sheathsm2 + meandC_Rhizomes, sem1mfert),
  lm(rootsm2 ~ rhizomesm2 + meanN_Roots, sem1mfert),
  lm(totepweightperarea ~ meanP_Blades, sem1mfert),
  lm(corebitesperarea ~ corebladeweightm2 + meanP_Sheaths, sem1mfert),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfert
)

summary(sem.fertV1.3)
AIC(sem.fertV1.3, aicc = T)


sem.fertV1.4 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + meanP_Blades, sem1mfert),
  lm(corebladeweightm2 ~ meandC_Blades + shootsm2, sem1mfert),
  lm(sheathsm2 ~ corebladeweightm2 + coreprodweight + shootsm2, sem1mfert),
  lm(rhizomesm2 ~ sheathsm2 + meandC_Rhizomes + shootsm2, sem1mfert),
  lm(rootsm2 ~ rhizomesm2 + meanN_Roots + corebladeweightm2, sem1mfert),
  lm(totepweightperarea ~ meanP_Blades, sem1mfert),
  lm(corebitesperarea ~ meanN_Roots + meandC_Rhizomes, sem1mfert),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfert
)

summary(sem.fertV1.4)
AIC(sem.fertV1.4, aicc = T)


sem.fertV1.5 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2, sem1mfert),
  lm(corebladeweightm2 ~ meandC_Blades + shootsm2, sem1mfert),
  lm(sheathsm2 ~ coreprodweight, sem1mfert),
  lm(rhizomesm2 ~ sheathsm2 + meandC_Rhizomes + shootsm2, sem1mfert),
  lm(rootsm2 ~ rhizomesm2 + meanN_Roots + corebladeweightm2, sem1mfert),
  lm(corebitesperarea ~ meanN_Roots + corebladeweightm2, sem1mfert),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfert
)

summary(sem.fertV1.5)
AIC(sem.fertV1.5, aicc = T)


sem.fertV1.6 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2 + meanN_Blades, sem1mfert),
  lm(corebladeweightm2 ~ meandC_Blades + shootsm2, sem1mfert), 
  lm(sheathsm2 ~ coreprodweight + meanC_Sheaths, sem1mfert), 
  lm(rhizomesm2 ~ sheathsm2 + meandC_Rhizomes + shootsm2, sem1mfert),
  lm(rootsm2 ~ rhizomesm2 + meanN_Roots + corebladeweightm2 + meandC_Blades + meanC_Blades + meandC_Roots, sem1mfert), #okay assumptions
  lm(corebitesperarea ~ meanN_Roots, sem1mfert),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfert
)

summary(sem.fertV1.6)
AIC(sem.fertV1.6, aicc = T)



#Variance inflation factors check
vif(lm(coreprodweight ~ corebladeweightm2 + shootsm2 + meanN_Blades, sem1mfert)) 
vif(lm(corebladeweightm2 ~ meandC_Blades + shootsm2, sem1mfert)) 
vif(lm(sheathsm2 ~ coreprodweight + meanC_Sheaths, sem1mfert)) 
vif(lm(rhizomesm2 ~ sheathsm2 + meandC_Rhizomes + shootsm2, sem1mfert)) 
vif(lm(rootsm2 ~ rhizomesm2 + meanN_Roots + corebladeweightm2 + meandC_Blades + meanC_Blades + meandC_Roots, sem1mfert)) 







##1m fish SEM --------------

prod1lm <- lm(coreprodweight ~ corebladeweightm2 + shootsm2 + totepweightperarea + corebitesperarea + meanN_Blades + meanP_Blades, data = sem1mfish) 
blade1lm <- lm(corebladeweightm2 ~ totepweightperarea + meandC_Blades + meanC_Blades + meanN_Blades + meanP_Blades, data = sem1mfish) 
sheath1lm <- lm(sheathsm2 ~ corebladeweightm2 + meandC_Sheaths + meanC_Sheaths + meanN_Sheaths + meanP_Sheaths, data = sem1mfish) 
rhizome1lm <- lm(rhizomesm2 ~ corebladeweightm2 + meandC_Rhizomes + meanC_Rhizomes + meanN_Rhizomes + meanP_Rhizomes, data = sem1mfish) 
root1lm <- lm(rootsm2 ~ rhizomesm2 + meandC_Roots + meanC_Roots + meanN_Roots + meanP_Roots, data = sem1mfish) 
eps1lm <- lm(totepweightperarea ~ meanN_Blades + meanP_Blades, data = sem1mfish) 
bites1m <- lm(corebitesperarea ~ meanN_Blades + meanP_Blades, data = sem1mfish) 


###1m fish SEM V1-----------
sem.fishV1.1 <- psem(
  prod1lm, blade1lm, sheath1lm, rhizome1lm, root1lm, eps1lm, bites1m,
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfish
)

summary(sem.fishV1.1)
AIC(sem.fishV1.1, aicc = T)
plot(sem.fishV1.1)


sem.fishV1.2 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2 + totepweightperarea + corebitesperarea + meanN_Blades + meanP_Blades, sem1mfish),
  lm(corebladeweightm2 ~ totepweightperarea + meandC_Blades + meanC_Blades + meanN_Blades + meanP_Blades + meanP_Roots, sem1mfish),
  lm(sheathsm2 ~ corebladeweightm2 + shootsm2 + meandC_Sheaths + meanC_Sheaths + meanN_Sheaths + meanP_Sheaths, sem1mfish),
  lm(rhizomesm2 ~ coreprodweight + corebladeweightm2 + shootsm2 + sheathsm2 + meandC_Rhizomes + meanC_Rhizomes + meanN_Rhizomes + meanP_Rhizomes + meanC_Roots, sem1mfish),
  lm(rootsm2 ~ rhizomesm2 + meandC_Roots + meanC_Roots + meanN_Roots + meanP_Roots + corebitesperarea, sem1mfish),
  lm(totepweightperarea ~ meanN_Blades + meanP_Blades + meandC_Blades, sem1mfish),
  lm(corebitesperarea ~ corebladeweightm2 + meanN_Blades + meanP_Blades + meanC_Rhizomes, data = sem1mfish),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfish
)

summary(sem.fishV1.2)
AIC(sem.fishV1.2, aicc = T)


sem.fishV1.3 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2 + meanN_Blades + meanP_Blades, sem1mfish),
  lm(corebladeweightm2 ~ shootsm2 + meandC_Blades + meanP_Roots + meanN_Roots, sem1mfish),
  lm(sheathsm2 ~ corebladeweightm2 + shootsm2, sem1mfish),
  lm(rhizomesm2 ~ sheathsm2 + meanC_Blades + meanN_Sheaths, sem1mfish),
  lm(rootsm2 ~ rhizomesm2 + meandC_Sheaths, sem1mfish),
  lm(totepweightperarea ~ meanN_Blades + meanP_Blades + meandC_Blades, sem1mfish),
  lm(corebitesperarea ~ meandC_Blades, data = sem1mfish),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfish
)

summary(sem.fishV1.3)
AIC(sem.fishV1.3, aicc = T)


sem.fishV1.4 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + meanN_Blades + meanP_Blades + shootsm2, sem1mfish),
  lm(corebladeweightm2 ~ meandC_Blades + meanP_Roots + meanN_Blades + meanN_Sheaths, sem1mfish),
  lm(sheathsm2 ~ corebladeweightm2 + shootsm2, sem1mfish),
  lm(rhizomesm2 ~ sheathsm2 + meanN_Sheaths + meanC_Blades + corebladeweightm2, sem1mfish),
  lm(rootsm2 ~ rhizomesm2, sem1mfish),
  lm(totepweightperarea ~ meanN_Blades + meanP_Blades + meandC_Blades, sem1mfish),
  lm(corebitesperarea ~ meanC_Rhizomes, data = sem1mfish),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfish
)

summary(sem.fishV1.4)
AIC(sem.fishV1.4, aicc = T)


sem.fishV1.5 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + meanP_Blades + meanN_Blades + shootsm2, sem1mfish),
  lm(corebladeweightm2 ~ shootsm2 + meandC_Blades + meanP_Roots, sem1mfish),
  lm(sheathsm2 ~ corebladeweightm2 + shootsm2, sem1mfish),
  lm(rhizomesm2 ~ sheathsm2 + corebladeweightm2 + meandC_Blades, sem1mfish),
  lm(rootsm2 ~ rhizomesm2, sem1mfish),
  lm(totepweightperarea ~ meanN_Blades + meanP_Blades + meandC_Blades, sem1mfish),
  lm(corebitesperarea ~ meanC_Rhizomes, data = sem1mfish),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfish
)

summary(sem.fishV1.5)
AIC(sem.fishV1.5, aicc = T)


sem.fishV1.6 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2 + meanN_Blades + meanP_Blades, sem1mfish),
  lm(corebladeweightm2 ~ shootsm2 + meandC_Blades + meanP_Roots + meanN_Blades + meanC_Rhizomes, sem1mfish),
  lm(sheathsm2 ~ corebladeweightm2 + shootsm2, sem1mfish),
  lm(rhizomesm2 ~ sheathsm2 + corebladeweightm2, sem1mfish),
  lm(rootsm2 ~ rhizomesm2, sem1mfish),
  lm(totepweightperarea ~ meanN_Blades + meanP_Blades + meandC_Blades, sem1mfish),
  lm(corebitesperarea ~ meanC_Rhizomes, data = sem1mfish),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfish
)

summary(sem.fishV1.6)
AIC(sem.fishV1.6, aicc = T)


sem.fishV1.7 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2 + meanN_Blades + meanP_Blades, sem1mfish),
  lm(corebladeweightm2 ~ shootsm2 + meandC_Blades + meanP_Roots + meanN_Blades + meanC_Rhizomes, sem1mfish),
  lm(sheathsm2 ~ corebladeweightm2 + shootsm2, sem1mfish),
  lm(rhizomesm2 ~ sheathsm2 + meanC_Blades + meanN_Sheaths, sem1mfish),
  lm(rootsm2 ~ rhizomesm2, sem1mfish),
  lm(totepweightperarea ~ meanN_Blades + meanP_Blades + meandC_Blades, sem1mfish),
  lm(corebitesperarea ~ meanC_Rhizomes, data = sem1mfish),
  
  coreprodweight %~~% corebladeweightm2,
  
  data = sem1mfish
)

summary(sem.fishV1.7)
AIC(sem.fishV1.7, aicc = T)


#Variance inflation factors check
vif(lm(coreprodweight ~ corebladeweightm2 + shootsm2 + meanN_Blades + meanP_Blades, sem1mfish)) 
vif(lm(corebladeweightm2 ~ shootsm2 + meandC_Blades + meanP_Roots + meanN_Blades + meanC_Rhizomes, sem1mfish)) 
vif(lm(sheathsm2 ~ corebladeweightm2 + shootsm2, sem1mfish)) 
vif(lm(rhizomesm2 ~ sheathsm2 + meanC_Blades + meanN_Sheaths, sem1mfish)) 
vif(lm(totepweightperarea ~ meanN_Blades + meanP_Blades + meandC_Blades, sem1mfish)) 




##20m ambient SEM -------------

prod1lm3 <- lm(coreprodweight ~ corebladeweightm2 + shootsm2 + totepweightperarea + corebitesperarea + meanN_Blades + meanP_Blades, data = sem20m) 
blade1lm3 <- lm(corebladeweightm2 ~ totepweightperarea + meandC_Blades + meanC_Blades + meanN_Blades + meanP_Blades, data = sem20m) 
sheath1lm3 <- lm(sheathsm2 ~ corebladeweightm2 + meandC_Sheaths + meanC_Sheaths + meanN_Sheaths + meanP_Sheaths, data = sem20m) 
rhizome1lm3 <- lm(rhizomesm2 ~ corebladeweightm2 + meandC_Rhizomes + meanC_Rhizomes + meanN_Rhizomes + meanP_Rhizomes, data = sem20m) 
root1lm3 <- lm(rootsm2 ~ rhizomesm2 + meandC_Roots + meanC_Roots + meanN_Roots + meanP_Roots, data = sem20m) 
eps1lm3 <- lm(totepweightperarea ~ meanN_Blades + meanP_Blades, data = sem20m) 
bites1m3 <- lm(corebitesperarea ~ meanN_Blades + meanP_Blades, data = sem20m) 


###20m ambient SEM V1-----------
sem.ambV1.1 <- psem(
  prod1lm3, blade1lm3, sheath1lm3, rhizome1lm3, root1lm3, eps1lm3, bites1m3,
  
  coreprodweight %~~% corebladeweightm2,
  coreprodweight %~~% shootsm2,
  
  data = sem20m
)

summary(sem.ambV1.1)
AIC(sem.ambV1.1, aicc = T)
plot(sem.ambV1.1)


sem.ambV1.2 <- psem(
  lm(coreprodweight ~ meandC_Blades + corebladeweightm2 + shootsm2 + totepweightperarea + corebitesperarea + meanN_Blades, sem20m),
  lm(corebladeweightm2 ~ meanP_Roots + meanN_Roots + shootsm2 + totepweightperarea + meandC_Blades + meanC_Blades + meanN_Blades, sem20m),
  lm(sheathsm2 ~ coreprodweight + meanP_Roots + meandC_Rhizomes + corebladeweightm2 + meandC_Sheaths + meanN_Sheaths + meanP_Sheaths, sem20m),
  lm(rhizomesm2 ~ corebitesperarea + totepweightperarea + meanN_Roots + meanC_Roots + meandC_Roots + meanN_Sheaths + meanC_Sheaths + meandC_Blades + shootsm2 + meanP_Blades + corebladeweightm2 + meanN_Rhizomes + meanP_Rhizomes, sem20m),
  lm(rootsm2 ~ rhizomesm2 + meandC_Roots, sem20m),
  lm(totepweightperarea ~ meanC_Roots + meanC_Rhizomes + meanC_Sheaths + meanN_Blades + meanP_Blades, sem20m),
  lm(corebitesperarea ~ meanN_Sheaths + meanN_Blades, sem20m),
  
  coreprodweight %~~% corebladeweightm2,
  coreprodweight %~~% shootsm2,
  
  data = sem20m
)

summary(sem.ambV1.2)
AIC(sem.ambV1.2, aicc = T)


sem.ambV1.3 <- psem(
  lm(coreprodweight ~ meandC_Blades + corebladeweightm2 + shootsm2 + meandC_Roots + meandC_Rhizomes, sem20m),
  lm(corebladeweightm2 ~ meanP_Roots + shootsm2, sem20m),
  lm(sheathsm2 ~ coreprodweight + corebladeweightm2 + meanN_Sheaths + meanC_Blades, sem20m),
  lm(rhizomesm2 ~ corebitesperarea + totepweightperarea + shootsm2 + meanP_Blades + meandC_Rhizomes + meanP_Rhizomes, sem20m),
  lm(rootsm2 ~ rhizomesm2 + meanN_Blades + meandC_Roots + sheathsm2, sem20m),
  lm(totepweightperarea ~ meanC_Roots, sem20m),
  lm(corebitesperarea ~ meanN_Sheaths, sem20m),
  
  coreprodweight %~~% corebladeweightm2,
  coreprodweight %~~% shootsm2,
  
  data = sem20m
)

summary(sem.ambV1.3)
AIC(sem.ambV1.3, aicc = T)


sem.ambV1.4 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2 + meanN_Blades, sem20m),
  lm(corebladeweightm2 ~ meanP_Roots + shootsm2, sem20m),
  lm(sheathsm2 ~ coreprodweight + corebladeweightm2 + meanN_Sheaths, sem20m),
  lm(rhizomesm2 ~ corebitesperarea + totepweightperarea + shootsm2 + meanP_Blades + meandC_Rhizomes + meanP_Rhizomes, sem20m),
  lm(rootsm2 ~ rhizomesm2 + meanN_Blades + meandC_Roots + sheathsm2, sem20m),
  lm(totepweightperarea ~ meanC_Roots, sem20m),
  lm(corebitesperarea ~ meanN_Sheaths, sem20m),
  
  coreprodweight %~~% corebladeweightm2,
  coreprodweight %~~% shootsm2,
  
  data = sem20m
)

summary(sem.ambV1.4)
AIC(sem.ambV1.4, aicc = T)


sem.ambV1.5 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2, sem20m),
  lm(corebladeweightm2 ~ meanP_Roots + shootsm2, sem20m),
  lm(sheathsm2 ~ coreprodweight + corebladeweightm2 + meanN_Sheaths, sem20m),
  lm(rhizomesm2 ~ corebitesperarea + totepweightperarea + shootsm2 + meanP_Blades + meandC_Rhizomes + meanP_Rhizomes, sem20m),
  lm(rootsm2 ~ rhizomesm2 + meanN_Blades + meandC_Roots + sheathsm2, sem20m),
  lm(totepweightperarea ~ meanC_Roots, sem20m),
  lm(corebitesperarea ~ meanC_Blades + meanN_Sheaths, sem20m),
  
  coreprodweight %~~% corebladeweightm2,
  coreprodweight %~~% shootsm2,
  
  data = sem20m
)

summary(sem.ambV1.5)
AIC(sem.ambV1.5, aicc = T)


sem.ambV1.6 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2, sem20m),
  lm(corebladeweightm2 ~ meanP_Roots + shootsm2, sem20m),
  lm(sheathsm2 ~ coreprodweight + corebladeweightm2 + meanN_Sheaths, sem20m),
  lm(rhizomesm2 ~ corebitesperarea + totepweightperarea + shootsm2 + meanP_Blades + meandC_Rhizomes + meanP_Rhizomes, sem20m),
  lm(rootsm2 ~ rhizomesm2 + meanN_Blades + meandC_Roots + sheathsm2, sem20m),
  lm(totepweightperarea ~ meanC_Roots, sem20m),
  lm(corebitesperarea ~ meanN_Sheaths, sem20m),
  
  coreprodweight %~~% corebladeweightm2,
  coreprodweight %~~% shootsm2,
  
  data = sem20m
)

summary(sem.ambV1.6)
AIC(sem.ambV1.6, aicc = T)


sem.ambV1.7 <- psem(
  lm(coreprodweight ~ corebladeweightm2 + shootsm2, sem20m),
  lm(corebladeweightm2 ~ meanP_Roots + shootsm2, sem20m),
  lm(sheathsm2 ~ coreprodweight + meanN_Sheaths, sem20m),
  lm(rhizomesm2 ~ corebitesperarea + totepweightperarea + shootsm2 + meanP_Blades + meandC_Rhizomes + meanP_Rhizomes, sem20m),
  lm(rootsm2 ~ rhizomesm2 + meanN_Blades + meandC_Roots + sheathsm2, sem20m),
  lm(totepweightperarea ~ meanC_Roots + meanP_Sheaths, sem20m),
  lm(corebitesperarea ~ meanN_Sheaths, sem20m),
  
  coreprodweight %~~% corebladeweightm2,
  coreprodweight %~~% shootsm2,
  
  data = sem20m
)

summary(sem.ambV1.7)
AIC(sem.ambV1.7, aicc = T)


#Variance inflation factors check
vif(lm(coreprodweight ~ corebladeweightm2 + shootsm2, sem20m)) 
vif(lm(corebladeweightm2 ~ meanP_Roots + shootsm2, sem20m)) 
vif(lm(sheathsm2 ~ coreprodweight + meanN_Sheaths, sem20m)) 
vif(lm(rhizomesm2 ~ corebitesperarea + totepweightperarea + shootsm2 + meanP_Blades + meandC_Rhizomes + meanP_Rhizomes, sem20m)) 
vif(lm(rootsm2 ~ rhizomesm2 + meanN_Blades + meandC_Roots + sheathsm2, sem20m)) 
vif(lm(totepweightperarea ~ meanC_Roots + meanP_Sheaths, sem20m)) 



