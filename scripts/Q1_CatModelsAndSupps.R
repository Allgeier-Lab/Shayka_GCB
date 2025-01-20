###########
#This code analyzes 1 m data with categorical nutrient linear models with Block in the model
##########

##Load libraries -------------
library(tidyverse)
library(ggpubr)
library(nortest) #for checking assumptions of linear models
library(nlme) #for lme function for mixed effects models


##Load functions ------------
stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) }#function


##Load data -----------------

databydist <- read_csv('processeddata/alldistdata.csv',
                       col_types = cols(Block = "c",))

bydist1mf <- databydist %>%
  filter(Distance == 1) %>%
  mutate(fish = case_when(fish == "high" ~ "lots",
                          fish == "low" ~ "few")) #flips the alphabetical order of the fish treatments to make interpreting results more intuitive


##1m models -----------
#responses all standardized after transforming the data

###Prod and biomass models -----------

prodmodel1f3 <- lme(stand(log(coreprodweight+1)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) #good assumptions
summary(prodmodel1f3)

#checks model assumptions
quartz()
rm=resid(prodmodel1f3)
fm=fitted(prodmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeweightmodel1f3 <- lme(stand(log(corebladeweightm2)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(bladeweightmodel1f3) 

#checks model assumptions
quartz()
rm=resid(bladeweightmodel1f3)
fm=fitted(bladeweightmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathweightmodel1f3 <- lme(stand(log(sheathsm2)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(sheathweightmodel1f3) 

#checks model assumptions
quartz()
rm=resid(sheathweightmodel1f3)
fm=fitted(sheathweightmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhizomeweightmodel1f3 <- lme(stand(log(rhizomesm2)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(rhizomeweightmodel1f3) 

#checks model assumptions
quartz()
rm=resid(rhizomeweightmodel1f3)
fm=fitted(rhizomeweightmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rootweightmodel1f3 <- lme(stand(rootsm2) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(rootweightmodel1f3) 

#checks model assumptions
quartz()
rm=resid(rootweightmodel1f3)
fm=fitted(rootweightmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###Shoots, epiphytes, and bites models -----------
shootsmodel1f3 <- lme(stand(log(shootsm2+1)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(shootsmodel1f3)

#checks model assumptions
quartz()
rm=resid(shootsmodel1f3)
fm=fitted(shootsmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


epsmodel1af3 <- lme(stand(totepweightperarea) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(epsmodel1af3) 

#checks model assumptions
quartz()
rm=resid(epsmodel1af3)
fm=fitted(epsmodel1af3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bitesmodel1f3 <- lme(stand(sqrt(corebitesperarea)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(bitesmodel1f3) 

#checks model assumptions
quartz()
rm=resid(bitesmodel1f3)
fm=fitted(bitesmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###Blade nutrient models -----------
bladeCmodel1f3 <- lme(stand(meanC_Blades) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(bladeCmodel1f3) 

#checks model assumptions
quartz()
rm=resid(bladeCmodel1f3)
fm=fitted(bladeCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladeNmodel1f3 <- lme(stand(meanN_Blades) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(bladeNmodel1f3) 

#checks model assumptions
quartz()
rm=resid(bladeNmodel1f3)
fm=fitted(bladeNmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladePmodel1f3 <- lme(stand(sqrt(meanP_Blades)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(bladePmodel1f3) 

#checks model assumptions
quartz()
rm=resid(bladePmodel1f3)
fm=fitted(bladePmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


bladedCmodel1f3 <- lme(stand(meandC_Blades) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(bladedCmodel1f3)

#checks model assumptions
quartz()
rm=resid(bladedCmodel1f3)
fm=fitted(bladedCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###Sheath nutrient models ----------- 
sheathCmodel1f3 <- lme(stand(log(meanC_Sheaths+1)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) #bad assumptions
summary(sheathCmodel1f3)

#checks model assumptions
quartz()
rm=resid(sheathCmodel1f3)
fm=fitted(sheathCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathNmodel1f3 <- lme(stand((meanN_Sheaths)^(1/3)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(sheathNmodel1f3) 

#checks model assumptions
quartz()
rm=resid(sheathNmodel1f3)
fm=fitted(sheathNmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathPmodel1f3 <- lme(stand(log(meanP_Sheaths)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(sheathPmodel1f3) 

#checks model assumptions
quartz()
rm=resid(sheathPmodel1f3)
fm=fitted(sheathPmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


sheathdCmodel1f3 <- lme(stand(meandC_Sheaths) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(sheathdCmodel1f3)

#checks model assumptions
quartz()
rm=resid(sheathdCmodel1f3)
fm=fitted(sheathdCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###Rhizome nutrient models ----------- 
rhizomeCmodel1f3 <- lme(stand(log(meanC_Rhizomes+1)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf)
summary(rhizomeCmodel1f3) 

#checks model assumptions
quartz()
rm=resid(rhizomeCmodel1f3)
fm=fitted(rhizomeCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhizomeNmodel1f3 <- lme(stand(log(meanN_Rhizomes)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(rhizomeNmodel1f3) 

#checks model assumptions
quartz()
rm=resid(rhizomeNmodel1f3)
fm=fitted(rhizomeNmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhizomePmodel1f3 <- lme(stand(sqrt(meanP_Rhizomes)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(rhizomePmodel1f3)

#checks model assumptions
quartz()
rm=resid(rhizomePmodel1f3)
fm=fitted(rhizomePmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rhizomedCmodel1f3 <- lme(stand(meandC_Rhizomes) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(rhizomedCmodel1f3) 

#checks model assumptions
quartz()
rm=resid(rhizomedCmodel1f3)
fm=fitted(rhizomedCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


###Root nutrient models ----------- 
rootCmodel1f3 <- lme(stand(log(meanC_Roots)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) #bad assumptions
summary(rootCmodel1f3)

#checks model assumptions
quartz()
rm=resid(rootCmodel1f3)
fm=fitted(rootCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rootNmodel1f3 <- lme(stand(meanN_Roots) ~ fish * fert, random = ~ 1|Block, data= bydist1mf)
summary(rootNmodel1f3) 

#checks model assumptions
quartz()
rm=resid(rootNmodel1f3)
fm=fitted(rootNmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rootPmodel1f3 <- lme(stand(sqrt(meanP_Roots)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(rootPmodel1f3)

#checks model assumptions
quartz()
rm=resid(rootPmodel1f3)
fm=fitted(rootPmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


rootdCmodel1f3 <- lme(stand(meandC_Roots) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
summary(rootdCmodel1f3)

#checks model assumptions
quartz()
rm=resid(rootdCmodel1f3)
fm=fitted(rootdCmodel1f3)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


##FigS2 -----------
prod1m <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=coreprodweight)) + 
  geom_boxplot(outlier.shape = NA) + #this hides the outliers on the boxplot (it doesn't remove them, just hides them so you don't get double points)
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('AG Production ('*~g^-1~m^-2~d^-1*')')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

blades1m <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=corebladeweightm2)) + 
  geom_boxplot(outlier.shape = NA) + #this hides the outliers on the boxplot (it doesn't remove them, just hides them so you don't get double points)
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Blade Biomass ('*g~DW~m^-2*')')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

sheaths1m <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=sheathsm2)) + 
  geom_boxplot(outlier.shape = NA) + #this hides the outliers on the boxplot (it doesn't remove them, just hides them so you don't get double points)
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Sheath Biomass ('*g~DW~m^-2*')')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rhizomes1m <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=rhizomesm2)) + 
  geom_boxplot(outlier.shape = NA) + #this hides the outliers on the boxplot (it doesn't remove them, just hides them so you don't get double points)
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Rhizome Biomass ('*g~DW~m^-2*')')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

roots1m <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=rootsm2)) + 
  geom_boxplot(outlier.shape = NA) + #this hides the outliers on the boxplot (it doesn't remove them, just hides them so you don't get double points)
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Root Biomass ('*g~DW~m^-2*')')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

shoots1m  <- databydist %>% filter(Distance == 1) %>% 
  ggplot(aes(x=labels, y=shootsm2)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Thalassia Shoots'~m^-2)) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

eps1m  <- databydist %>% filter(Distance == 1) %>% 
  ggplot(aes(x=labels, y=totepweightperarea)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Epiphyte Biomass ('*g~mm^-2*')')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

bites1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=corebitesperarea)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Bites'~mm^-2)) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

bladeC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanC_Blades)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Blade %C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

bladeN1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanN_Blades)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Blade %N')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

bladeP1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanP_Blades)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Blade %P')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

bladedC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meandC_Blades)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Blade d13C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

sheathC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanC_Sheaths)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Sheath %C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

sheathN1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanN_Sheaths)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Sheath %N')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

sheathP1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanP_Sheaths)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Sheath %P')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

sheathdC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meandC_Sheaths)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Sheath d13C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rhizomeC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanC_Rhizomes)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Rhizome %C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rhizomeN1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanN_Rhizomes)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Rhizome %N')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rhizomeP1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanP_Rhizomes)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Rhizome %P')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rhizomedC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meandC_Rhizomes)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Rhizome d13C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rootC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanC_Roots)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Root %C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rootN1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanN_Roots)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Root %N')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rootP1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meanP_Roots)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Root %P')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

rootdC1m  <- databydist %>% filter(Distance == 1) %>%
  ggplot(aes(x=labels, y=meandC_Roots)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.2, height=0) +
  stat_summary(fun="mean", geom="point", shape="x", size=6, color="black") + #this puts the mean on as an X
  labs(y = expression('Root d13C')) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"),
        axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

plots1m <- ggpubr::ggarrange(ggarrange(blades1m, sheaths1m, rhizomes1m, roots1m, labels = c("A", "B", "C", "D"), ncol = 4),
                              ggarrange(prod1m, shoots1m, eps1m, bites1m, labels = c("E", "F", "G", "H"), ncol = 4),
                              ggarrange(bladeC1m, bladeN1m, bladeP1m, bladedC1m, labels = c("I", "J", "K", "L"), ncol = 4),
                              ggarrange(sheathC1m, sheathN1m, sheathP1m, sheathdC1m, labels = c("M", "N", "O", "P"), ncol = 4),
                              ggarrange(rhizomeC1m, rhizomeN1m, rhizomeP1m, rhizomedC1m, labels = c("Q", "R", "S", "T"), ncol = 4),
                              ggarrange(rootC1m, rootN1m, rootP1m, rootdC1m, labels = c("U", "V", "W", "X"), ncol = 4),
                              nrow = 6)

ggsave(filename = "FigS2.pdf", path="outputs", plot=plots1m, device = "pdf", width = 18, height = 18, units="in", dpi=300)


