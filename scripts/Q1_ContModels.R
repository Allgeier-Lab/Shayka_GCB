###########
#This code analyzes 1 m data with continuous nutrient linear models with Block in the model
##########

##Load libraries -------------
library(tidyverse)
library(nortest) #for checking assumptions of linear models
library(nlme) #for lme function for mixed effects models



##Load functions ------------
stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) }#function


##Load data -----------------

databydist <- read_csv('processeddata/alldistdata.csv',
                       col_types = cols(Block = "c",))

bydist1m <- databydist %>%
  filter(Distance == 1)



##1m models -----------
#everything is standardized (responses and predictors) after transforming the data

###Prod, growth, and biomass models -----------
prodcontmodelN4 <- lme(stand(log(coreprodweight)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(prodcontmodelN4)

#checks model assumptions
quartz()
rm=resid(prodcontmodelN4)
fm=fitted(prodcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


prodcontmodelP4 <- lme(stand(log(coreprodweight)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(prodcontmodelP4)

#checks model assumptions
quartz()
rm=resid(prodcontmodelP4)
fm=fitted(prodcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


prodcontmodelNP4 <- lme(stand(log(coreprodweight+1)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(prodcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(prodcontmodelNP4)
fm=fitted(prodcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(prodcontmodelN4)$AIC #best
summary(prodcontmodelP4)$AIC #w/in 2
summary(prodcontmodelNP4)$AIC



bladeweightcontmodelN4 <- lme(stand(sqrt(corebladeweightm2)^(1/3)) ~ stand(totN), random = ~ 1|Block, data= bydist1m)
summary(bladeweightcontmodelN4) 

#checks model assumptions
quartz()
rm=resid(bladeweightcontmodelN4)
fm=fitted(bladeweightcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeweightcontmodelP4 <- lme(stand(log(corebladeweightm2+1)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(bladeweightcontmodelP4) 

#checks model assumptions
quartz()
rm=resid(bladeweightcontmodelP4)
fm=fitted(bladeweightcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeweightcontmodelNP4 <- lme(stand(sqrt(corebladeweightm2+10)) ~ stand(log(NP+10)), random = ~ 1|Block, data= bydist1m) 
summary(bladeweightcontmodelNP4) 

#checks model assumptions
quartz()
rm=resid(bladeweightcontmodelNP4)
fm=fitted(bladeweightcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(bladeweightcontmodelN4)$AIC #best
summary(bladeweightcontmodelP4)$AIC #w/in 2
summary(bladeweightcontmodelNP4)$AIC 


sheathweightcontmodelN4 <- lme(stand(log(sheathsm2+1)) ~ stand(sqrt(totN+1)), random = ~ 1|Block, data= bydist1m) 
summary(sheathweightcontmodelN4) 

#checks model assumptions
quartz()
rm=resid(sheathweightcontmodelN4)
fm=fitted(sheathweightcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathweightcontmodelP4 <- lme(stand(log(sheathsm2)) ~ stand(totP^(1/3)), random = ~ 1|Block, data= bydist1m) 
summary(sheathweightcontmodelP4) 

#checks model assumptions
quartz()
rm=resid(sheathweightcontmodelP4)
fm=fitted(sheathweightcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathweightcontmodelNP4 <- lme(stand((sheathsm2)^(1/7)) ~ stand((NP)^(1/6)), random = ~ 1|Block, data= bydist1m) 
summary(sheathweightcontmodelNP4) 

#checks model assumptions
quartz()
rm=resid(sheathweightcontmodelNP4)
fm=fitted(sheathweightcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(sheathweightcontmodelN4)$AIC #w/in 2
summary(sheathweightcontmodelP4)$AIC #best



rhzweightcontmodelN4 <- lme(stand(log(rhizomesm2+1)) ~ stand(log(totN+1)), random = ~ 1|Block, data= bydist1m) 
summary(rhzweightcontmodelN4) 

#checks model assumptions
quartz()
rm=resid(rhzweightcontmodelN4)
fm=fitted(rhzweightcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzweightcontmodelP4 <- lme(stand(log(rhizomesm2+1)) ~ stand(log(totP+1)), random = ~ 1|Block, data= bydist1m) 
summary(rhzweightcontmodelP4) 

#checks model assumptions
quartz()
rm=resid(rhzweightcontmodelP4)
fm=fitted(rhzweightcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzweightcontmodelNP4 <- lme(stand(log(rhizomesm2+10)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rhzweightcontmodelNP4) 

#checks model assumptions
quartz()
rm=resid(rhzweightcontmodelNP4)
fm=fitted(rhzweightcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(rhzweightcontmodelP4)$AIC #best


rtweightcontmodelN4 <- lme(stand(sqrt(rootsm2+10)) ~ stand(log(totN+1)), random = ~ 1|Block, data= bydist1m) 
summary(rtweightcontmodelN4) 

#checks model assumptions
quartz()
rm=resid(rtweightcontmodelN4)
fm=fitted(rtweightcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtweightcontmodelP4 <- lme(stand(sqrt(rootsm2+1)) ~ stand(log(totP+1)), random = ~ 1|Block, data= bydist1m) 
summary(rtweightcontmodelP4) 

#checks model assumptions
quartz()
rm=resid(rtweightcontmodelP4)
fm=fitted(rtweightcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtweightcontmodelNP4 <- lme(stand(log(rootsm2+1)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rtweightcontmodelNP4) 

#checks model assumptions
quartz()
rm=resid(rtweightcontmodelNP4)
fm=fitted(rtweightcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(rtweightcontmodelN4)$AIC #same
summary(rtweightcontmodelP4)$AIC #w/in 1





###Shoots, epiphytes, and bites models -----------
shootcontmodelN4 <- lme(stand(sqrt(shootsm2+1)) ~ stand(log(totN+1)), random = ~ 1|Block, data= bydist1m) 
summary(shootcontmodelN4)

#checks model assumptions
quartz()
rm=resid(shootcontmodelN4)
fm=fitted(shootcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


shootcontmodelP4 <- lme(stand(log(shootsm2+10)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(shootcontmodelP4)

#checks model assumptions
quartz()
rm=resid(shootcontmodelP4)
fm=fitted(shootcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


shootcontmodelNP4 <- lme(stand(sqrt(shootsm2+1)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(shootcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(shootcontmodelNP4)
fm=fitted(shootcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(shootcontmodelN4)$AIC #best


epscontmodelN4 <- lme(stand(totepweightperarea^(1/3)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(epscontmodelN4)

#checks model assumptions
quartz()
rm=resid(epscontmodelN4)
fm=fitted(epscontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


epscontmodelP4 <- lme(stand(totepweightperarea^(1/3)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(epscontmodelP4)

#checks model assumptions
quartz()
rm=resid(epscontmodelP4)
fm=fitted(epscontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


epscontmodelNP4 <- lme(stand(totepweightperarea^(1/3)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(epscontmodelNP4)

#checks model assumptions
quartz()
rm=resid(epscontmodelNP4)
fm=fitted(epscontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(epscontmodelN4)$AIC #w/in 2
summary(epscontmodelP4)$AIC #best
summary(epscontmodelNP4)$AIC 




bitescontmodelN4 <- lme(stand((corebitesperarea)^(1/3)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(bitescontmodelN4)

#checks model assumptions
quartz()
rm=resid(bitescontmodelN4)
fm=fitted(bitescontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bitescontmodelP4 <- lme(stand(sqrt(corebitesperarea)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(bitescontmodelP4)

#checks model assumptions
quartz()
rm=resid(bitescontmodelP4)
fm=fitted(bitescontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bitescontmodelNP4 <- lme(stand(sqrt(corebitesperarea)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(bitescontmodelNP4)

#checks model assumptions
quartz()
rm=resid(bitescontmodelNP4)
fm=fitted(bitescontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(bitescontmodelN4)$AIC 
summary(bitescontmodelP4)$AIC
summary(bitescontmodelNP4)$AIC #best




###Blade nutrient models -----------
bladeCcontmodelN4 <- lme(stand(meanC_Blades) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(bladeCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(bladeCcontmodelN4)
fm=fitted(bladeCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeCcontmodelP4 <- lme(stand(meanC_Blades) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(bladeCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(bladeCcontmodelP4)
fm=fitted(bladeCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeCcontmodelNP4 <- lme(stand(meanC_Blades) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(bladeCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(bladeCcontmodelNP4)
fm=fitted(bladeCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeNcontmodelN4 <- lme(stand(meanN_Blades) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(bladeNcontmodelN4)

#checks model assumptions
quartz()
rm=resid(bladeNcontmodelN4)
fm=fitted(bladeNcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeNcontmodelP4 <- lme(stand(log(meanN_Blades+1)) ~ stand(log(totP+1)), random = ~ 1|Block, data= bydist1m) 
summary(bladeNcontmodelP4)

#checks model assumptions
quartz()
rm=resid(bladeNcontmodelP4)
fm=fitted(bladeNcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeNcontmodelNP4 <- lme(stand(meanN_Blades) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(bladeNcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(bladeNcontmodelNP4)
fm=fitted(bladeNcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(bladeNcontmodelN4)$AIC
summary(bladeNcontmodelP4)$AIC #best
summary(bladeNcontmodelNP4)$AIC 



bladePcontmodelN4 <- lme(stand(log(meanP_Blades+1)) ~ stand(log(totN)), random = ~ 1|Block, data= bydist1m) 
summary(bladePcontmodelN4)

#checks model assumptions
quartz()
rm=resid(bladePcontmodelN4)
fm=fitted(bladePcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladePcontmodelP4 <- lme(stand(log(meanP_Blades)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(bladePcontmodelP4)

#checks model assumptions
quartz()
rm=resid(bladePcontmodelP4)
fm=fitted(bladePcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladePcontmodelNP4 <- lme(stand(meanP_Blades) ~ stand(log(NP+1)), random = ~ 1|Block, data= bydist1m) 
summary(bladePcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(bladePcontmodelNP4)
fm=fitted(bladePcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(bladePcontmodelN4)$AIC
summary(bladePcontmodelP4)$AIC #best
summary(bladePcontmodelNP4)$AIC 




bladedCcontmodelN4 <- lme(stand(log(meandC_Blades+100)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(bladedCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(bladedCcontmodelN4)
fm=fitted(bladedCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

bladedCcontmodelP4 <- lme(stand(log(meandC_Blades+100)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
summary(bladedCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(bladedCcontmodelP4)
fm=fitted(bladedCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

bladedCcontmodelNP4 <- lme(stand(log(meandC_Blades+100)) ~ stand(log(NP)), random = ~ 1|Block, data= bydist1m) 
summary(bladedCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(bladedCcontmodelNP4)
fm=fitted(bladedCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(bladedCcontmodelN4)$AIC 
summary(bladedCcontmodelP4)$AIC #best
summary(bladedCcontmodelNP4)$AIC 




###Sheath nutrient models ----------- 
sheathCcontmodelN4 <- lme(stand(meanC_Sheaths) ~ stand(totN), random = ~ 1|Block, data= bydist1m) #bad assumptions
summary(sheathCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(sheathCcontmodelN4)
fm=fitted(sheathCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathCcontmodelP4 <- lme(stand((meanC_Sheaths+1)^(1/3)) ~ stand(totP^(1/3)), random = ~ 1|Block, data= bydist1m) #bad assumptions
summary(sheathCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(sheathCcontmodelP4)
fm=fitted(sheathCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathCcontmodelNP4 <- lme(stand((meanC_Sheaths+1)^(1/3)) ~ stand(NP^(1/3)), random = ~ 1|Block, data= bydist1m) #bad assumptions
summary(sheathCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(sheathCcontmodelNP4)
fm=fitted(sheathCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathNcontmodelN4 <- lme(stand((meanN_Sheaths)^(1/3)) ~ stand((totN)^(1/5)), random = ~ 1|Block, data= bydist1m) 
summary(sheathNcontmodelN4)

#checks model assumptions
quartz()
rm=resid(sheathNcontmodelN4)
fm=fitted(sheathNcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathNcontmodelP4 <- lme(stand(log(meanN_Sheaths)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
summary(sheathNcontmodelP4)

#checks model assumptions
quartz()
rm=resid(sheathNcontmodelP4)
fm=fitted(sheathNcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathNcontmodelNP4 <- lme(stand(meanN_Sheaths) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(sheathNcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(sheathNcontmodelNP4)
fm=fitted(sheathNcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(sheathNcontmodelN4)$AIC #best
summary(sheathNcontmodelP4)$AIC #w/in 1
summary(sheathNcontmodelNP4)$AIC




sheathPcontmodelN4 <- lme(stand(log(meanP_Sheaths)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(sheathPcontmodelN4)

#checks model assumptions
quartz()
rm=resid(sheathPcontmodelN4)
fm=fitted(sheathPcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathPcontmodelP4 <- lme(stand(sqrt(meanP_Sheaths)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(sheathPcontmodelP4)

#checks model assumptions
quartz()
rm=resid(sheathPcontmodelP4)
fm=fitted(sheathPcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathPcontmodelNP4 <- lme(stand((meanP_Sheaths)^(1/4)) ~ stand((NP)^(1/3)), random = ~ 1|Block, data= bydist1m) 
summary(sheathPcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(sheathPcontmodelNP4)
fm=fitted(sheathPcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(sheathPcontmodelN4)$AIC #best
summary(sheathPcontmodelP4)$AIC #w/in 2
summary(sheathPcontmodelNP4)$AIC 




sheathdCcontmodelN4 <- lme(stand(meandC_Sheaths) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(sheathdCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(sheathdCcontmodelN4)
fm=fitted(sheathdCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathdCcontmodelP4 <- lme(stand(meandC_Sheaths) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
summary(sheathdCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(sheathdCcontmodelP4)
fm=fitted(sheathdCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathdCcontmodelNP4 <- lme(stand(meandC_Sheaths) ~ stand(sqrt(NP)), random = ~ 1|Block, data= bydist1m) 
summary(sheathdCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(sheathdCcontmodelNP4)
fm=fitted(sheathdCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(sheathdCcontmodelN4)$AIC #w/in 2
summary(sheathdCcontmodelP4)$AIC #best
summary(sheathdCcontmodelNP4)$AIC




###Rhizome nutrient models ----------- 
rhzCcontmodelN4 <- lme(stand(meanC_Rhizomes) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(rhzCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rhzCcontmodelN4)
fm=fitted(rhzCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzCcontmodelP4 <- lme(stand(log(meanC_Rhizomes+10)) ~ stand(sqrt(totP+1)), random = ~ 1|Block, data= bydist1m) 
summary(rhzCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rhzCcontmodelP4)
fm=fitted(rhzCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzCcontmodelNP4 <- lme(stand(log(meanC_Rhizomes+10)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rhzCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rhzCcontmodelNP4)
fm=fitted(rhzCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(rhzCcontmodelN4)$AIC #best




rhzNcontmodelN4 <- lme(stand(log(meanN_Rhizomes)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(rhzNcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rhzNcontmodelN4)
fm=fitted(rhzNcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzNcontmodelP4 <- lme(stand(sqrt(meanN_Rhizomes)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(rhzNcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rhzNcontmodelP4)
fm=fitted(rhzNcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzNcontmodelNP4 <- lme(stand(log(meanN_Rhizomes)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rhzNcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rhzNcontmodelNP4)
fm=fitted(rhzNcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



rhzPcontmodelN4 <- lme(stand(log(meanP_Rhizomes)) ~ stand(log(totN)), random = ~ 1|Block, data= bydist1m) 
summary(rhzPcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rhzPcontmodelN4)
fm=fitted(rhzPcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

rhzPcontmodelP4 <- lme(stand(log(meanP_Rhizomes)) ~ stand(sqrt(totP+1)), random = ~ 1|Block, data= bydist1m) 
summary(rhzPcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rhzPcontmodelP4)
fm=fitted(rhzPcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzPcontmodelNP4 <- lme(stand(meanP_Rhizomes) ~ stand(log(NP+1)), random = ~ 1|Block, data= bydist1m) 
summary(rhzPcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rhzPcontmodelNP4)
fm=fitted(rhzPcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(rhzPcontmodelN4)$AIC 
summary(rhzPcontmodelP4)$AIC #best
summary(rhzPcontmodelNP4)$AIC 



rhzdCcontmodelN4 <- lme(stand(meandC_Rhizomes) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(rhzdCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rhzdCcontmodelN4)
fm=fitted(rhzdCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzdCcontmodelP4 <- lme(stand(meandC_Rhizomes) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
summary(rhzdCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rhzdCcontmodelP4)
fm=fitted(rhzdCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhzdCcontmodelNP4 <- lme(stand(meandC_Rhizomes) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rhzdCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rhzdCcontmodelNP4)
fm=fitted(rhzdCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(rhzdCcontmodelN4)$AIC #best
summary(rhzdCcontmodelP4)$AIC #w/in 1
summary(rhzdCcontmodelNP4)$AIC




###Root nutrient models ----------- 
rtCcontmodelN4 <- lme(stand(meanC_Roots^(1/5)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) #bad assumptions
summary(rtCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rtCcontmodelN4)
fm=fitted(rtCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtCcontmodelP4 <- lme(stand((meanC_Roots)^(1/6)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) #bad assumptions
summary(rtCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rtCcontmodelP4)
fm=fitted(rtCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtCcontmodelNP4 <- lme(stand(log(meanC_Roots)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) #bad assumptions
summary(rtCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rtCcontmodelNP4)
fm=fitted(rtCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtNcontmodelN4 <- lme(stand(log(meanN_Roots+10)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(rtNcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rtNcontmodelN4)
fm=fitted(rtNcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtNcontmodelP4 <- lme(stand(meanN_Roots) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
summary(rtNcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rtNcontmodelP4)
fm=fitted(rtNcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtNcontmodelNP4 <- lme(stand(meanN_Roots) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rtNcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rtNcontmodelNP4)
fm=fitted(rtNcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



rtPcontmodelN4 <- lme(stand(meanP_Roots^(1/5)) ~ stand(totN^(1/4)), random = ~ 1|Block, data= bydist1m) 
summary(rtPcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rtPcontmodelN4)
fm=fitted(rtPcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtPcontmodelP4 <- lme(stand(log(meanP_Roots)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
summary(rtPcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rtPcontmodelP4)
fm=fitted(rtPcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtPcontmodelNP4 <- lme(stand(log(meanP_Roots+1)) ~ stand(log(NP)), random = ~ 1|Block, data= bydist1m) 
summary(rtPcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rtPcontmodelNP4)
fm=fitted(rtPcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(rtPcontmodelN4)$AIC 
summary(rtPcontmodelP4)$AIC #best
summary(rtPcontmodelNP4)$AIC




rtdCcontmodelN4 <- lme(stand(meandC_Roots) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
summary(rtdCcontmodelN4)

#checks model assumptions
quartz()
rm=resid(rtdCcontmodelN4)
fm=fitted(rtdCcontmodelN4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtdCcontmodelP4 <- lme(stand(sqrt(meandC_Roots+100)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m)
summary(rtdCcontmodelP4)

#checks model assumptions
quartz()
rm=resid(rtdCcontmodelP4)
fm=fitted(rtdCcontmodelP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rtdCcontmodelNP4 <- lme(stand(log(meandC_Roots+10)) ~ stand(log(NP)), random = ~ 1|Block, data= bydist1m) 
summary(rtdCcontmodelNP4)

#checks model assumptions
quartz()
rm=resid(rtdCcontmodelNP4)
fm=fitted(rtdCcontmodelNP4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)

summary(rtdCcontmodelN4)$AIC #best
summary(rtdCcontmodelP4)$AIC #w/in 2
summary(rtdCcontmodelNP4)$AIC


