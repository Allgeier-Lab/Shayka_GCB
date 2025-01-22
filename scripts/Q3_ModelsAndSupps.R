###########
#This code runs linear models to determine what drives C productivity
##########


##Load libraries -------------
library(tidyverse)
library(nortest) #for checking assumptions of linear models
library(nlme) #for lme function for mixed effects models
library(MuMIn) #for dredge
library(emmeans) #for post-hoc tests

##Load functions ------------
stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) }

##Load data -----------------

databydist <- read_csv('processeddata/alldistdata.csv',
                       col_types = cols(Block = "c",))

data5nuts <- databydist %>%
  ungroup() %>%
  select(Reef, Block, Treatment, Transect, Distance,
         agthalC, bgthalC, totthalC, agCprod, bgCprod, totCprod, agCturn, bgCturn, totCturn) %>%
  mutate(nuttrt = case_when(Distance == 1 ~ case_when(Treatment == "A" ~ "A",
                                                      Treatment == "B" ~ "B",
                                                      Treatment == "C" ~ "C",
                                                      Treatment == "D" ~ "D"),
                            Distance == 20 ~ "E"))


##5 nutrient treatments ----------
###C stock models -----------------

#AG
agCstockmod <- lme(log(agthalC+1) ~ nuttrt, random = ~ 1|Block, data= data5nuts) #good assumptions
summary(agCstockmod)
emmeans(agCstockmod, pairwise ~ nuttrt) #A=a,c, B=a,b, C=b, D=b, E=c

#checks model assumptions
quartz()
rm=resid(agCstockmod)
fm=fitted(agCstockmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


#BG
bgCstockmod <- lme(bgthalC ~ nuttrt, random = ~ 1|Block, data= data5nuts) #good assumptions
summary(bgCstockmod)
emmeans(bgCstockmod, pairwise ~ nuttrt) #A=a,b, B=a,b, C=a, D=a, E=b

#checks model assumptions
quartz()
rm=resid(bgCstockmod)
fm=fitted(bgCstockmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



#Total
totCstockmod <- lme(totthalC ~ nuttrt, random = ~ 1|Block, data= data5nuts) #good assumptions
summary(totCstockmod)
emmeans(totCstockmod, pairwise ~ nuttrt) #A=a, B=a, C=a, D=a, E=a

#checks model assumptions
quartz()
rm=resid(totCstockmod)
fm=fitted(totCstockmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###C production models -----------------
#AG
agCprodmod <- lme(log(agCprod) ~ nuttrt, random = ~ 1|Block, data= data5nuts) #pretty good assumptions
summary(agCprodmod)
emmeans(agCprodmod, pairwise ~ nuttrt) #A=a, B=b, C=c, D=c, E=a

#checks model assumptions
quartz()
rm=resid(agCprodmod)
fm=fitted(agCprodmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


#BG
bgCprodmod <- lme((bgCprod+1)^2 ~ nuttrt, random = ~ 1|Block, data= data5nuts) #pretty good assumptions
summary(bgCprodmod)
emmeans(bgCprodmod, pairwise ~ nuttrt) #A=a, B=a, C=a, D=a, E=b

#checks model assumptions
quartz()
rm=resid(bgCprodmod)
fm=fitted(bgCprodmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



#Total
totCprodmod <- lme(log(totCprod) ~ nuttrt, random = ~ 1|Block, data= data5nuts) #good assumptions
summary(totCprodmod)
emmeans(totCprodmod, pairwise ~ nuttrt) #A=a, B=b, C=c, D=c, E=d

#checks model assumptions
quartz()
rm=resid(totCprodmod)
fm=fitted(totCprodmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###C turnover models -----------------
#AG
agCturnmod <- lme(log(agCturn) ~ nuttrt, random = ~ 1|Block, data= data5nuts) #pretty good assumptions
summary(agCturnmod)
emmeans(agCturnmod, pairwise ~ nuttrt) #A=a,c, B=a, C=b, D=b, E=c

#checks model assumptions
quartz()
rm=resid(agCturnmod)
fm=fitted(agCturnmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


#BG
bgCturnmod <- lme(log(bgCturn) ~ nuttrt, random = ~ 1|Block, data= data5nuts) #good assumptions
summary(bgCturnmod)
emmeans(bgCturnmod, pairwise ~ nuttrt) #A=a,b, B=a, C=a, D=a, E=b

#checks model assumptions
quartz()
rm=resid(bgCturnmod)
fm=fitted(bgCturnmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



#Total
totCturnmod <- lme(log(totCturn) ~ nuttrt, random = ~ 1|Block, data= data5nuts) #good assumptions
summary(totCturnmod)
emmeans(totCturnmod, pairwise ~ nuttrt) #A=a, B=a,b, C=a,b, D=b, E=c

#checks model assumptions
quartz()
rm=resid(totCturnmod)
fm=fitted(totCturnmod)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###Plots - FigS5 -----------------

#first make a new dataframe with all the posthoc letters; but UPDATE the tukey letters to have the background treatment first
tukeylets <- tibble(nuttrt = c("E","A","B","C","D"),
                    agCstocktuk = c("a","a,b","b,c","c","c"), #E=c, A=a,c, B=a,b, C=b, D=b
                    bgCstocktuk = c("a","a,b","a,b","b","b"), #E=b, A=a,b, B=a,b, C=a, D=a
                    totCstocktuk = c("a","a","a","a","a"), #E=a, A=a, B=a, C=a, D=a
                    agCprodtuk = c("a","a","b","c","c"), #E=a, A=a, B=b, C=c, D=c
                    bgCprodtuk = c("a","b","b","b","b"), #E=b, A=a, B=a, C=a, D=a
                    totCprodtuk = c("a","b","c","d","d"), #E=d, A=a, B=b, C=c, D=c
                    agCturntuk = c("a","a,b","b","c","c"), #E=c, A=a,c, B=a, C=b, D=b
                    bgCturntuk = c("a","a,b","b","b","b"), #E=b, A=a,b, B=a, C=a, D=a
                    totCturntuk = c("a","b","b,c","b,c","c")) #E=c, A=a, B=a,b, C=a,b, D=b


agCstockplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=agthalC)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 220, label=agCstocktuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('AG C Stock ('*gC~m^-2*')'),
       tag = "(a)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))

bgCstockplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=bgthalC)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 450, label=bgCstocktuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('BG C Stock ('*gC~m^-2*')'),
       tag = "(b)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))

totCstockplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=totthalC)) + 
  geom_boxplot(outlier.shape = NA) +
#  geom_text(aes(x=nuttrt, y= 510, label=totCstocktuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('Total C Stock ('*gC~m^-2*')'),
       tag = "(c)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))


agCprodplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=agCprod)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 2.4, label=agCprodtuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('AG C Production ('*gC~m^-2~d^-1*')'),
       tag = "(d)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))

bgCprodplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=bgCprod)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 1.35, label=bgCprodtuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('BG C Production ('*gC~m^-2~d^-1*')'),
       tag = "(e)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))

totCprodplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=totCprod)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 3.2, label=totCprodtuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('Total C Production ('*gC~m^-2~d^-1*')'),
       tag = "(f)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))


agCturnplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=agCturn)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 0.017, label=agCturntuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('AG C Turnover ('*d^-1*')'),
       tag = "(g)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))

bgCturnplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=bgCturn)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 0.021, label=bgCturntuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('BG C Turnover ('*d^-1*')'),
       tag = "(h)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))

totCturnplot <- data5nuts %>%
  mutate(nuttrt = factor(nuttrt, levels=c("E", "A", "B", "C", "D"))) %>%
  ggplot(aes(x=nuttrt, y=totCturn)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(x=nuttrt, y= 0.016, label=totCturntuk), data= tukeylets) +
  theme_classic() +
  geom_point(position = position_jitter(width = 0.15), alpha=0.3) + 
  stat_summary(fun="mean", geom="point", shape="x", size=4, color="black", fill="black") + #this puts the mean on as an X 
  labs(x = "Nutrient treatment",
       y = expression('Total C Turnover ('*d^-1*')'),
       tag = "(i)") +
  scale_x_discrete(labels = c('A' = "-F-A",
                              'B' = "+F-A",
                              'C' = "-F+A",
                              'D' = "+F+A",
                              'E' = "B")) +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13, color = "black"))


allCplots <- ggpubr::ggarrange(agCstockplot, bgCstockplot, totCstockplot, agCprodplot, bgCprodplot, totCprodplot, agCturnplot, bgCturnplot, totCturnplot,
                              ncol = 3, nrow = 3, align = "v") #, common.legend = TRUE, legend="bottom"

ggsave(filename = "figS5.pdf", path="outputs", plot=allCplots, device = "pdf", width = 13, height = 10, units="in", dpi=300)


##Drivers of turnover ----------
###Models -----------
####AG -----------
agCturnpred4 <- lm(log(agCturn) ~ avgheight + shootsm2, data= databydist) #good assumptions
summary(agCturnpred4)

#checks model assumptions
quartz()
rm=resid(agCturnpred4)
fm=fitted(agCturnpred4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


####BG -----------
bgCturnpred4 <- lm(sqrt(bgCturn) ~ avgheight + shootsm2, data= databydist) #good assumptions
summary(bgCturnpred4)

#checks model assumptions
quartz()
rm=resid(bgCturnpred4)
fm=fitted(bgCturnpred4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)


####Total -----------
totCturnpred4 <- lm(sqrt(totCturn) ~ log(avgheight) + log(shootsm2), data= databydist) #good assumptions
summary(totCturnpred4)

#checks model assumptions
quartz()
rm=resid(totCturnpred4)
fm=fitted(totCturnpred4)
model2=lm(rm~fm)
par(mfrow=c(3,2))
plot(model2)
hist(rm)
plot(c(0,5),c(0,5))
text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)



###Supp Plots -----------

####AG C turnover - FigS6 ---------
agCturn_height <- databydist %>%
  ggplot(aes(y=agCturn, x=avgheight)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Blade Height (mm)',
       tag = "(a)")
agCturn_shoots <- databydist %>%
  ggplot(aes(y=agCturn, x=shootsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Thalassia Shoots'~m^-2),
       tag="(b)")
agCturn_agprod <- databydist %>%
  ggplot(aes(y=agCturn, x=coreprodweight)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('AG Production ('*~g^-1~m^-2~d^-1*')'),
       tag="(c)")
agCturn_bgprod <- databydist %>%
  ggplot(aes(y=agCturn, x=bgprod)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('BG Production ('*~g^-1~m^-2~d^-1*')'),
       tag="(d)")
agCturn_blades <- databydist %>%
  ggplot(aes(y=agCturn, x=corebladeweightm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Blade Biomass ('*g~DW~m^-2*')'),
       tag="(e)")
agCturn_sheaths <- databydist %>%
  ggplot(aes(y=agCturn, x=sheathsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Sheath Biomass ('*g~DW~m^-2*')'),
       tag="(f)")
agCturn_rhizomes <- databydist %>%
  ggplot(aes(y=agCturn, x=rhizomesm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Rhizome Biomass ('*g~DW~m^-2*')'),
       tag="(g)")
agCturn_roots <- databydist %>%
  ggplot(aes(y=agCturn, x=rootsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Root Biomass ('*g~DW~m^-2*')'),
       tag="(h)")
agCturn_bladeC <- databydist %>%
  ggplot(aes(y=agCturn, x=meanC_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Blade %C',
       tag="(i)")
agCturn_bladeN <- databydist %>%
  ggplot(aes(y=agCturn, x=meanN_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Blade %N',
       tag="(j)")
agCturn_bladeP <- databydist %>%
  ggplot(aes(y=agCturn, x=meanP_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Blade %P',
       tag="(k)")
agCturn_bladedC <- databydist %>%
  ggplot(aes(y=agCturn, x=meandC_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Blade'~delta^13*'C'),
       tag="(l)")
agCturn_sheathC <- databydist %>%
  ggplot(aes(y=agCturn, x=meanC_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Sheath %C',
       tag="(m)")
agCturn_sheathN <- databydist %>%
  ggplot(aes(y=agCturn, x=meanN_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Sheath %N',
       tag="(n)")
agCturn_sheathP <- databydist %>%
  ggplot(aes(y=agCturn, x=meanP_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Sheath %P',
       tag="(o)")
agCturn_sheathdC <- databydist %>%
  ggplot(aes(y=agCturn, x=meandC_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Sheath'~delta^13*'C'),
       tag="(p)")
agCturn_rhzC <- databydist %>%
  ggplot(aes(y=agCturn, x=meanC_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Rhizome %C',
       tag="(q)")
agCturn_rhzN <- databydist %>%
  ggplot(aes(y=agCturn, x=meanN_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG N Turnover',
       x='Rhizome %N',
       tag="(r)")
agCturn_rhzP <- databydist %>%
  ggplot(aes(y=agCturn, x=meanP_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG P Turnover',
       x='Rhizome %P',
       tag="(s)")
agCturn_rhzdC <- databydist %>%
  ggplot(aes(y=agCturn, x=meandC_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Rhizome'~delta^13*'C'),
       tag="(t)")
agCturn_rtC <- databydist %>%
  ggplot(aes(y=agCturn, x=meanC_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Root %C',
       tag="(u)")
agCturn_rtN <- databydist %>%
  ggplot(aes(y=agCturn, x=meanN_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Root %N',
       tag="(v)")
agCturn_rtP <- databydist %>%
  ggplot(aes(y=agCturn, x=meanP_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x='Root %P',
       tag="(w)")
agCturn_rtdC <- databydist %>%
  ggplot(aes(y=agCturn, x=meandC_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Root'~delta^13*'C'),
       tag="(x)")
agCturn_eps <- databydist %>%
  ggplot(aes(y=agCturn, x=totepweightperarea)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Epiphyte Biomass ('*g~mm^-2*')'),
       tag="(y)")
agCturn_bites <- databydist %>%
  ggplot(aes(y=agCturn, x=corebitesperarea)) +
  geom_point() + 
  theme_classic() +
  labs(y='AG C Turnover',
       x=expression('Bites'~mm^-2),
       tag="(z)")

agCturncorrs <- ggpubr::ggarrange(agCturn_height, agCturn_shoots, agCturn_agprod, agCturn_bgprod,
                                  agCturn_blades, agCturn_sheaths, agCturn_rhizomes, agCturn_roots,
                                  agCturn_bladeC, agCturn_bladeN, agCturn_bladeP, agCturn_bladedC,
                                  agCturn_sheathC, agCturn_sheathN, agCturn_sheathP, agCturn_sheathdC,
                                  agCturn_rhzC, agCturn_rhzN, agCturn_rhzP, agCturn_rhzdC,
                                  agCturn_rtC, agCturn_rtN, agCturn_rtP, agCturn_rtdC,
                                  agCturn_eps, agCturn_bites,
                                  ncol = 4, nrow = 7, align = "v")

ggsave(filename = "figS6.pdf", path="outputs", plot=agCturncorrs, device = "pdf", width = 18, height = 23, units="in", dpi=300)




####BG C turnover - Fig S7 ---------
bgCturn_height <- databydist %>%
  ggplot(aes(y=bgCturn, x=avgheight)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Blade Height (mm)',
       tag = "(a)")
bgCturn_shoots <- databydist %>%
  ggplot(aes(y=bgCturn, x=shootsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Thalassia Shoots'~m^-2),
       tag="(b)")
bgCturn_agprod <- databydist %>%
  ggplot(aes(y=bgCturn, x=coreprodweight)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('AG Production ('*~g^-1~m^-2~d^-1*')'),
       tag="(c)")
bgCturn_bgprod <- databydist %>%
  ggplot(aes(y=bgCturn, x=bgprod)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('BG Production ('*~g^-1~m^-2~d^-1*')'),
       tag="(d)")
bgCturn_blades <- databydist %>%
  ggplot(aes(y=bgCturn, x=corebladeweightm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Blade Biomass ('*g~DW~m^-2*')'),
       tag="(e)")
bgCturn_sheaths <- databydist %>%
  ggplot(aes(y=bgCturn, x=sheathsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Sheath Biomass ('*g~DW~m^-2*')'),
       tag="(f)")
bgCturn_rhizomes <- databydist %>%
  ggplot(aes(y=bgCturn, x=rhizomesm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Rhizome Biomass ('*g~DW~m^-2*')'),
       tag="(g)")
bgCturn_roots <- databydist %>%
  ggplot(aes(y=bgCturn, x=rootsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Root Biomass ('*g~DW~m^-2*')'),
       tag="(h)")
bgCturn_bladeC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanC_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Blade %C',
       tag="(i)")
bgCturn_bladeN <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanN_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Blade %N',
       tag="(j)")
bgCturn_bladeP <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanP_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Blade %P',
       tag="(k)")
bgCturn_bladedC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meandC_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Blade'~delta^13*'C'),
       tag="(l)")
bgCturn_sheathC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanC_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Sheath %C',
       tag="(m)")
bgCturn_sheathN <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanN_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Sheath %N',
       tag="(n)")
bgCturn_sheathP <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanP_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Sheath %P',
       tag="(o)")
bgCturn_sheathdC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meandC_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Sheath'~delta^13*'C'),
       tag="(p)")
bgCturn_rhzC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanC_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Rhizome %C',
       tag="(q)")
bgCturn_rhzN <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanN_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG N Turnover',
       x='Rhizome %N',
       tag="(r)")
bgCturn_rhzP <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanP_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG P Turnover',
       x='Rhizome %P',
       tag="(s)")
bgCturn_rhzdC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meandC_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Rhizome'~delta^13*'C'),
       tag="(t)")
bgCturn_rtC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanC_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Root %C',
       tag="(u)")
bgCturn_rtN <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanN_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Root %N',
       tag="(v)")
bgCturn_rtP <- databydist %>%
  ggplot(aes(y=bgCturn, x=meanP_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x='Root %P',
       tag="(w)")
bgCturn_rtdC <- databydist %>%
  ggplot(aes(y=bgCturn, x=meandC_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Root'~delta^13*'C'),
       tag="(x)")
bgCturn_eps <- databydist %>%
  ggplot(aes(y=bgCturn, x=totepweightperarea)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Epiphyte Biomass ('*g~mm^-2*')'),
       tag="(y)")
bgCturn_bites <- databydist %>%
  ggplot(aes(y=bgCturn, x=corebitesperarea)) +
  geom_point() + 
  theme_classic() +
  labs(y='BG C Turnover',
       x=expression('Bites'~mm^-2),
       tag="(z)")

bgCturncorrs <- ggpubr::ggarrange(bgCturn_height, bgCturn_shoots, bgCturn_agprod, bgCturn_bgprod,
                                  bgCturn_blades, bgCturn_sheaths, bgCturn_rhizomes, bgCturn_roots,
                                  bgCturn_bladeC, bgCturn_bladeN, bgCturn_bladeP, bgCturn_bladedC,
                                  bgCturn_sheathC, bgCturn_sheathN, bgCturn_sheathP, bgCturn_sheathdC,
                                  bgCturn_rhzC, bgCturn_rhzN, bgCturn_rhzP, bgCturn_rhzdC,
                                  bgCturn_rtC, bgCturn_rtN, bgCturn_rtP, bgCturn_rtdC,
                                  bgCturn_eps, bgCturn_bites,
                                  ncol = 4, nrow = 7, align = "v")

ggsave(filename = "figS7.pdf", path="outputs", plot=bgCturncorrs, device = "pdf", width = 18, height = 23, units="in", dpi=300)


####Total C turnover - Fig S8 ---------
totCturn_height <- databydist %>%
  ggplot(aes(y=totCturn, x=avgheight)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Blade Height (mm)',
       tag = "(a)")
totCturn_shoots <- databydist %>%
  ggplot(aes(y=totCturn, x=shootsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Thalassia Shoots'~m^-2),
       tag="(b)")
totCturn_agprod <- databydist %>%
  ggplot(aes(y=totCturn, x=coreprodweight)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('AG Production ('*~g^-1~m^-2~d^-1*')'),
       tag="(c)")
totCturn_bgprod <- databydist %>%
  ggplot(aes(y=totCturn, x=bgprod)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('BG Production ('*~g^-1~m^-2~d^-1*')'),
       tag="(d)")
totCturn_blades <- databydist %>%
  ggplot(aes(y=totCturn, x=corebladeweightm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Blade Biomass ('*g~DW~m^-2*')'),
       tag="(e)")
totCturn_sheaths <- databydist %>%
  ggplot(aes(y=totCturn, x=sheathsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Sheath Biomass ('*g~DW~m^-2*')'),
       tag="(f)")
totCturn_rhizomes <- databydist %>%
  ggplot(aes(y=totCturn, x=rhizomesm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Rhizome Biomass ('*g~DW~m^-2*')'),
       tag="(g)")
totCturn_roots <- databydist %>%
  ggplot(aes(y=totCturn, x=rootsm2)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Root Biomass ('*g~DW~m^-2*')'),
       tag="(h)")
totCturn_bladeC <- databydist %>%
  ggplot(aes(y=totCturn, x=meanC_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Blade %C',
       tag="(i)")
totCturn_bladeN <- databydist %>%
  ggplot(aes(y=totCturn, x=meanN_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Blade %N',
       tag="(j)")
totCturn_bladeP <- databydist %>%
  ggplot(aes(y=totCturn, x=meanP_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Blade %P',
       tag="(k)")
totCturn_bladedC <- databydist %>%
  ggplot(aes(y=totCturn, x=meandC_Blades)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Blade'~delta^13*'C'),
       tag="(l)")
totCturn_sheathC <- databydist %>%
  ggplot(aes(y=totCturn, x=meanC_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Sheath %C',
       tag="(m)")
totCturn_sheathN <- databydist %>%
  ggplot(aes(y=totCturn, x=meanN_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Sheath %N',
       tag="(n)")
totCturn_sheathP <- databydist %>%
  ggplot(aes(y=totCturn, x=meanP_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Sheath %P',
       tag="(o)")
totCturn_sheathdC <- databydist %>%
  ggplot(aes(y=totCturn, x=meandC_Sheaths)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Sheath'~delta^13*'C'),
       tag="(p)")
totCturn_rhzC <- databydist %>%
  ggplot(aes(y=totCturn, x=meanC_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Rhizome %C',
       tag="(q)")
totCturn_rhzN <- databydist %>%
  ggplot(aes(y=totCturn, x=meanN_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total N Turnover',
       x='Rhizome %N',
       tag="(r)")
totCturn_rhzP <- databydist %>%
  ggplot(aes(y=totCturn, x=meanP_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total P Turnover',
       x='Rhizome %P',
       tag="(s)")
totCturn_rhzdC <- databydist %>%
  ggplot(aes(y=totCturn, x=meandC_Rhizomes)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Rhizome'~delta^13*'C'),
       tag="(t)")
totCturn_rtC <- databydist %>%
  ggplot(aes(y=totCturn, x=meanC_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Root %C',
       tag="(u)")
totCturn_rtN <- databydist %>%
  ggplot(aes(y=totCturn, x=meanN_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Root %N',
       tag="(v)")
totCturn_rtP <- databydist %>%
  ggplot(aes(y=totCturn, x=meanP_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x='Root %P',
       tag="(w)")
totCturn_rtdC <- databydist %>%
  ggplot(aes(y=totCturn, x=meandC_Roots)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Root'~delta^13*'C'),
       tag="(x)")
totCturn_eps <- databydist %>%
  ggplot(aes(y=totCturn, x=totepweightperarea)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Epiphyte Biomass ('*g~mm^-2*')'),
       tag="(y)")
totCturn_bites <- databydist %>%
  ggplot(aes(y=totCturn, x=corebitesperarea)) +
  geom_point() + 
  theme_classic() +
  labs(y='Total C Turnover',
       x=expression('Bites'~mm^-2),
       tag="(z)")

totCturncorrs <- ggpubr::ggarrange(totCturn_height, totCturn_shoots, totCturn_agprod, totCturn_bgprod,
                                   totCturn_blades, totCturn_sheaths, totCturn_rhizomes, totCturn_roots,
                                   totCturn_bladeC, totCturn_bladeN, totCturn_bladeP, totCturn_bladedC,
                                   totCturn_sheathC, totCturn_sheathN, totCturn_sheathP, totCturn_sheathdC,
                                   totCturn_rhzC, totCturn_rhzN, totCturn_rhzP, totCturn_rhzdC,
                                   totCturn_rtC, totCturn_rtN, totCturn_rtP, totCturn_rtdC,
                                   totCturn_eps, totCturn_bites,
                                  ncol = 4, nrow = 7, align = "v")

ggsave(filename = "figS8.pdf", path="outputs", plot=totCturncorrs, device = "pdf", width = 18, height = 23, units="in", dpi=300)



