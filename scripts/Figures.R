###########
#This code makes the figures
##########


##Load libraries -------------
library(tidyverse)
library(ggpubr)
library(nortest) #for checking assumptions of linear models
library(nlme) #for lme function for mixed effects models
library(ggpattern) #for nutrient plot


##Load functions ------------
stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) }


##Load data -----------------

databydist <- read_csv('processeddata/alldistdata.csv',
                       col_types = cols(Block = "c",))

bydist1mf <- databydist %>%
  filter(Distance == 1) %>%
  mutate(fish = case_when(fish == "high" ~ "lots",
                          fish == "low" ~ "few")) #flips the alphabetical order of the fish treatments to make interpreting results more intuitive

bydist1m <- databydist %>%
  filter(Distance == 1)

data5nutsLong <- databydist %>%
  ungroup() %>%
  select(Reef, Block, Treatment, Transect, Distance,
         agthalC, bgthalC, agCprod, bgCprod, agCturn, bgCturn) %>%
  mutate(nuttrt = case_when(Distance == 1 ~ case_when(Treatment == "A" ~ "2",
                                                      Treatment == "B" ~ "3",
                                                      Treatment == "C" ~ "4",
                                                      Treatment == "D" ~ "5"),
                            Distance == 20 ~ "1")) %>%
  pivot_longer(cols = c(agthalC:bgthalC, agCprod:bgCprod, agCturn:bgCturn),
               names_to = c("agbg",".value"),
               names_pattern = "(..)(.+)$")
  
data5nutsSummary <- data5nutsLong %>%
  group_by(nuttrt, agbg) %>%
  summarise(Cstockmean = mean(thalC),
            Cstocksd = sd(thalC),
            Cprodmean = mean(Cprod),
            Cprodsd = sd(Cprod),
            Cturnmean = mean(Cturn),
            Cturnsd = sd(Cturn))

reefnutrients <- read_csv('data/reef.fish.nutrients.csv') %>%
  tidyr::separate_wider_delim(cols = reef.date, delim = "-", names = c("reef", "date")) %>%
  tidyr::separate_wider_position(cols = reef, widths = c(Block = 1, Treatment = 1)) %>%
  mutate(Treatment = case_when(Treatment == "a" ~ "A",
                               Treatment == "b" ~ "B",
                               Treatment == "c" ~ "C",
                               Treatment == "d" ~ "D")) %>%
  group_by(Block, Treatment) %>%
  summarise(avgfishN = mean(n), #g/day/reef
            sdfishN = sd(n),
            avgfishP = mean(p), #g/day/reef
            sdfishP = sd(p)) %>% 
  mutate(avgfertN = case_when(Treatment == "A" ~ 0, #g/day/reef
                           Treatment == "B" ~ 0,
                           Treatment == "C" ~ 2.7,
                           Treatment == "D" ~ 2.7),
         sdfertN = 0,
         avgfertP = case_when(Treatment == "A" ~ 0, #g/day/reef
                           Treatment == "B" ~ 0,
                           Treatment == "C" ~ 0.39,
                           Treatment == "D" ~ 0.39),
         sdfertP = 0) %>%
  unite(col = Reef, Block, Treatment, sep = "") %>%
  mutate(Reef = str_c("YM",Reef)) %>%
  mutate(avgfishP = avgfishP*10, avgfertP = avgfertP*10, sdfishP = sdfishP*10, sdfertP = sdfertP*10) %>% #this is for dual axis plotting
  pivot_longer(cols = c(-Reef),
               names_to = c("stat","source",".value"),
               names_pattern = "(.+)(....)(.)$") %>%
  pivot_longer(cols = c(N,P),
               names_to = "nutrient") %>%
  filter(stat=="sd") %>%
  mutate(source = case_when(source == "fish" ~ "avgfish",
                            source == "fert" ~ "fert")) %>%
  rename(stdev = value) %>%
  select(-stat)

##Figure 1 (just the nutrient plot part) -----------------

reefnutsplot <- databydist %>%
  ungroup() %>%
  select(Reef, avgfishN, avgfishP, fertN, fertP) %>%
  mutate(avgfishP = avgfishP * 10, fertP = fertP * 10) %>% #this is for dual axis plotting
  pivot_longer(cols = c(avgfishN,avgfishP,fertN,fertP), 
               names_to = c("source", ".value"), 
               names_pattern = "(.+)(.)$") %>% 
  pivot_longer(cols = c(N,P),
               names_to = "nutrient") %>%
  group_by(Reef, source, nutrient) %>%
  summarise(value = mean(value)) %>%
  left_join(.,reefnutrients) %>% #adds the stdev values to the dataframe
  group_by(Reef,nutrient) %>%
  mutate(total = cumsum(value),
         ordering = max(total)) %>% #this makes a column with max N and P values repeated for each reef to help with sorting x axis by total N
  ggplot(aes(x=forcats::fct_reorder(Reef, ordering), y=total, fill=nutrient)) +
  geom_col(data = . %>% filter(source=="fert"), position = position_dodge(width = 0.9), alpha = 0.5) +  #alpha = 0.5
  geom_col(data = . %>% filter(source=="avgfish"), position = position_dodge(width = 0.9)) + #alpha = 1
  geom_errorbar(data = . %>% filter(source=="avgfish"), aes(ymin=value-stdev, ymax=value+stdev, x=Reef), position = position_dodge(width = 0.9), width=.4, size =.7) + 
  geom_tile(aes(y=NA_integer_, alpha = factor(source))) + 
  scale_alpha_manual(values = c(1,0.4), labels = c("avgfish" = "Fish-derived", "fert" = "Anthropogenic")) +
  guides(alpha = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expression('N'~g~day^-1), sec.axis = sec_axis(~ . / 10, name = expression('P'~g~day^-1)), expand=expansion(mult=0,add=c(0,0.2))) + #removes space at bottom but not top
  theme_classic() + 
  scale_fill_manual(values = c(N = "#EE7624", P = "#00205C")) +
  labs(alpha= "Nutrient Source") +
  guides(fill = "none") +
  theme(axis.title.y = element_text(color = "#EE7624", size=18),
        axis.text.y = element_text(color = "#EE7624", size=18),
        axis.title.y.right = element_text(color = "#00205C", size=18),
        axis.text.y.right = element_text(color = "#00205C", size=18),
        axis.text.x = element_text(size=18, colour = "black", angle = 90, vjust = .5),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "inside",
        legend.position.inside = c(.2,.8),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) + 
  scale_x_discrete(labels = c('YM1A' = '-F-A', 'YM2A'= '-F-A', 'YM3A'= '-F-A', 'YM4A' = '-F-A',
                              'YM1B' = '+F-A', 'YM2B' = '+F-A', 'YM3B' = '+F-A', 'YM4B' = '+F-A',
                              'YM1C' = '-F+A', 'YM2C' = '-F+A', 'YM3C' = '-F+A', 'YM4C' = '-F+A',
                              'YM1D' = '+F+A', 'YM2D' = '+F+A', 'YM3D' = '+F+A', 'YM4D' = '+F+A'))

ggsave(filename = "fig1.pdf", path="outputs", plot=reefnutsplot, device = "pdf", width = 9, height = 4, units="in", dpi=300)


##Figure 2 -----------------
###a - categorical matrix -------------
####load models ---------
prodmodel1f3 <- lme(stand(log(coreprodweight+1)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
bladeweightmodel1f3 <- lme(stand(log(corebladeweightm2)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
sheathweightmodel1f3 <- lme(stand(log(sheathsm2)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rhizomeweightmodel1f3 <- lme(stand(log(rhizomesm2)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rootweightmodel1f3 <- lme(stand(rootsm2) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
shootsmodel1f3 <- lme(stand(log(shootsm2+1)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
bladeCmodel1f3 <- lme(stand(meanC_Blades) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
bladeNmodel1f3 <- lme(stand(meanN_Blades) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
bladePmodel1f3 <- lme(stand(sqrt(meanP_Blades)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
bladedCmodel1f3 <- lme(stand(meandC_Blades) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
sheathNmodel1f3 <- lme(stand((meanN_Sheaths)^(1/3)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
sheathPmodel1f3 <- lme(stand(log(meanP_Sheaths)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
sheathdCmodel1f3 <- lme(stand(meandC_Sheaths) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rhizomeCmodel1f3 <- lme(stand(log(meanC_Rhizomes+1)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rhizomeNmodel1f3 <- lme(stand(log(meanN_Rhizomes)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rhizomePmodel1f3 <- lme(stand(sqrt(meanP_Rhizomes)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rhizomedCmodel1f3 <- lme(stand(meandC_Rhizomes) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rootNmodel1f3 <- lme(stand(meanN_Roots) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rootPmodel1f3 <- lme(stand(sqrt(meanP_Roots)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
rootdCmodel1f3 <- lme(stand(meandC_Roots) ~ fish * fert, random = ~ 1|Block, data= bydist1mf)
epsmodel1af3 <- lme(stand(totepweightperarea) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 
bitesmodel1f3 <- lme(stand(sqrt(corebitesperarea)) ~ fish * fert, random = ~ 1|Block, data= bydist1mf) 

####results matrix compilation -------
mods <- list(prodmodel1f3, bladeweightmodel1f3, sheathweightmodel1f3, rhizomeweightmodel1f3, rootweightmodel1f3,
             shootsmodel1f3, epsmodel1af3, bitesmodel1f3,
             bladeCmodel1f3, bladeNmodel1f3, bladePmodel1f3, bladedCmodel1f3,
             sheathNmodel1f3, sheathPmodel1f3, sheathdCmodel1f3,
             rhizomeCmodel1f3, rhizomeNmodel1f3, rhizomePmodel1f3, rhizomedCmodel1f3,
             rootNmodel1f3, rootPmodel1f3, rootdCmodel1f3)
modresultsmat <- matrix(ncol = 4, nrow = length(mods)) 
modestmat <- matrix(ncol = 4, nrow = length(mods))

for(i in 1:length(mods)){
  modresultsmat[i,] <- summary(mods[[i]])$tTable[,5] 
}

for(i in 1:length(mods)){
  modestmat[i,] <- summary(mods[[i]])$coefficients$fixed 
}


colnames(modresultsmat) <- c('intercept pval', 'fish pval', 'fert pval', 'fish x fert pval')
rownames(modresultsmat) <- c('prod', 'bladeweight', 'sheathweight', 'rhizomeweight', 'rootweight',
                             'shoots', 'epiphytes', 'bites',
                             'bladeC', 'bladeN', 'bladeP', 'bladedC',
                             'sheathN', 'sheathP', 'sheathdC',
                             'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                             'rootN', 'rootP', 'rootdC')

colnames(modestmat) <- c('intercept est', 'fish est', 'fert est', 'fish x fert est')
rownames(modestmat) <- c('prod', 'bladeweight', 'sheathweight','rhizomeweight', 'rootweight',
                         'shoots', 'epiphytes', 'bites',
                         'bladeC', 'bladeN', 'bladeP', 'bladedC',
                         'sheathN', 'sheathP', 'sheathdC',
                         'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                         'rootN', 'rootP', 'rootdC')


####make it in ggplot ---------

mat1f <- data.frame(modresultsmat[,"fish pval"]) %>%
  cbind(modestmat[,"fish est"]) %>%
  cbind("treatment" = "fish") %>%
  cbind("response" = c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 'Epiphytes', 'Bites',
                       'Blade %C', 'Blade %N', 'Blade %P', 'Blade d13C', 'Sheath %N', 'Sheath %P', 'Sheath d13C',
                       'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Rhizome d13C', 'Root %N', 'Root %P', 'Root d13C'))
colnames(mat1f) <- c("pval", "est", "treatment", "response")
rownames(mat1f) <- NULL

mat2f <- data.frame(modresultsmat[,"fert pval"]) %>%
  cbind(modestmat[,"fert est"]) %>%
  cbind("treatment" = "fert") %>%
  cbind("response" = c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 'Epiphytes', 'Bites',
                       'Blade %C', 'Blade %N', 'Blade %P', 'Blade d13C', 'Sheath %N', 'Sheath %P', 'Sheath d13C',
                       'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Rhizome d13C', 'Root %N', 'Root %P', 'Root d13C'))
colnames(mat2f) <- c("pval", "est", "treatment", "response")
rownames(mat2f) <- NULL

mat3f <- data.frame(modresultsmat[,"fish x fert pval"]) %>%
  cbind(modestmat[,"fish x fert est"]) %>%
  cbind("treatment" = "fish x fert") %>%
  cbind("response" = c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 'Epiphytes', 'Bites',
                       'Blade %C', 'Blade %N', 'Blade %P', 'Blade d13C', 'Sheath %N', 'Sheath %P', 'Sheath d13C',
                       'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Rhizome d13C', 'Root %N', 'Root %P', 'Root d13C'))
colnames(mat3f) <- c("pval", "est", "treatment", "response")
rownames(mat3f) <- NULL

long1mf <- rbind(mat1f, mat2f, mat3f) %>%
  mutate("signif" = ifelse(pval<0.05,"*","")) %>%
  mutate(group = case_when(grepl('AG Production|Blade Biomass|Sheath Biomass|Rhizome Biomass|Root Biomass|Shoot Density|Blade d13C|Sheath d13C|Rhizome d13C|Root d13C', response) ~ 1,
                           grepl('Blade %C|Blade %N|Blade %P|Sheath %N|Sheath %P|Rhizome %C|Rhizome %N|Rhizome %P|Root %N|Root %P|Epiphytes', response) ~ 2,
                           grepl('Bites', response) ~ 3,
                           TRUE ~ NA))


long1mf$treatment <- factor(long1mf$treatment, levels = c("fish", "fert", "fish x fert"))
long1mf$response <- factor(long1mf$response, levels=c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 
                                                      'Blade d13C', 'Sheath d13C', 'Rhizome d13C', 'Root d13C',
                                                      'Blade %C', 'Blade %N', 'Blade %P', 'Sheath %N', 'Sheath %P',
                                                      'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Root %N', 'Root %P',
                                                      'Epiphytes', 'Bites'))


catmodelplot1m <- ggplot(long1mf, aes(x=treatment, y=fct_rev(response), fill=est, label=signif)) +
  geom_tile(color = "black",
            lwd = .5,
            linetype = 1,
            size = 2) +
  geom_text(size=8, vjust=0.8) +
  theme_classic() +
  facet_grid(group~., scales = 'free_y', space = 'free_y') +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(5,1,5,5),
        legend.position = "none") + #this hides the legend for the combo figure
  scale_fill_gradient2(low = "#1151FF", mid = "white", high = "#D50E0E", midpoint = 0) +
  scale_y_discrete(labels = c('Blade d13C' = expression('Blade'~delta^13*'C'),
                              'Sheath d13C' = expression('Sheath'~delta^13*'C'),
                              'Rhizome d13C' = expression('Rhizome'~delta^13*'C'),
                              'Root d13C' = expression('Root'~delta^13*'C'))) +
  labs(tag = "(a)") + theme(plot.tag = element_text(face = "bold")) +
  scale_x_discrete(labels = c("fish" = "Fish-derived", "fert" = "Anthropogenic", "fish x fert" = "Fish x Anthro")) 

ggsave(filename = "catmatrix.pdf", path="outputs", plot=catmodelplot1m, device = "pdf", width = 3, height = 6, units="in", dpi=300)


###b - continuous matrix -------------
####load models ---------
prodcontmodelN4 <- lme(stand(log(coreprodweight)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
prodcontmodelP4 <- lme(stand(log(coreprodweight)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
prodcontmodelNP4 <- lme(stand(log(coreprodweight+1)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(prodcontmodelN4)$AIC #best
summary(prodcontmodelP4)$AIC #w/in 2
summary(prodcontmodelNP4)$AIC

bladeweightcontmodelN4 <- lme(stand(sqrt(corebladeweightm2)^(1/3)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
bladeweightcontmodelP4 <- lme(stand(log(corebladeweightm2+1)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
bladeweightcontmodelNP4 <- lme(stand(sqrt(corebladeweightm2+10)) ~ stand(log(NP+10)), random = ~ 1|Block, data= bydist1m) 
summary(bladeweightcontmodelN4)$AIC #best
summary(bladeweightcontmodelP4)$AIC #w/in 2
summary(bladeweightcontmodelNP4)$AIC 

sheathweightcontmodelN4 <- lme(stand(log(sheathsm2+1)) ~ stand(sqrt(totN+1)), random = ~ 1|Block, data= bydist1m) 
sheathweightcontmodelP4 <- lme(stand(log(sheathsm2)) ~ stand(totP^(1/3)), random = ~ 1|Block, data= bydist1m) 
sheathweightcontmodelNP4 <- lme(stand((sheathsm2)^(1/7)) ~ stand((NP)^(1/6)), random = ~ 1|Block, data= bydist1m) 
summary(sheathweightcontmodelN4)$AIC #w/in 2
summary(sheathweightcontmodelP4)$AIC #best

rhzweightcontmodelN4 <- lme(stand(log(rhizomesm2+1)) ~ stand(log(totN+1)), random = ~ 1|Block, data= bydist1m)
rhzweightcontmodelP4 <- lme(stand(log(rhizomesm2+1)) ~ stand(log(totP+1)), random = ~ 1|Block, data= bydist1m) 
rhzweightcontmodelNP4 <- lme(stand(log(rhizomesm2+10)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rhzweightcontmodelP4)$AIC #best

rtweightcontmodelN4 <- lme(stand(sqrt(rootsm2+10)) ~ stand(log(totN+1)), random = ~ 1|Block, data= bydist1m)
rtweightcontmodelP4 <- lme(stand(sqrt(rootsm2+1)) ~ stand(log(totP+1)), random = ~ 1|Block, data= bydist1m)
rtweightcontmodelNP4 <- lme(stand(log(rootsm2+1)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rtweightcontmodelN4)$AIC #best
summary(rtweightcontmodelP4)$AIC #w/in 1

shootcontmodelN4 <- lme(stand(sqrt(shootsm2+1)) ~ stand(log(totN+1)), random = ~ 1|Block, data= bydist1m) 
shootcontmodelP4 <- lme(stand(log(shootsm2+10)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
shootcontmodelNP4 <- lme(stand(sqrt(shootsm2+1)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(shootcontmodelN4)$AIC #best

epscontmodelN4 <- lme(stand(totepweightperarea^(1/3)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
epscontmodelP4 <- lme(stand(totepweightperarea^(1/3)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
epscontmodelNP4 <- lme(stand(totepweightperarea^(1/3)) ~ stand(NP), random = ~ 1|Block, data= bydist1m)
summary(epscontmodelN4)$AIC #w/in 2
summary(epscontmodelP4)$AIC #best
summary(epscontmodelNP4)$AIC 

bitescontmodelN4 <- lme(stand((corebitesperarea)^(1/3)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
bitescontmodelP4 <- lme(stand(sqrt(corebitesperarea)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
bitescontmodelNP4 <- lme(stand(sqrt(corebitesperarea)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(bitescontmodelN4)$AIC 
summary(bitescontmodelP4)$AIC
summary(bitescontmodelNP4)$AIC #best

bladeCcontmodelN4 <- lme(stand(meanC_Blades) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
bladeCcontmodelP4 <- lme(stand(meanC_Blades) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
bladeCcontmodelNP4 <- lme(stand(meanC_Blades) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 

bladeNcontmodelN4 <- lme(stand(meanN_Blades) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
bladeNcontmodelP4 <- lme(stand(log(meanN_Blades+1)) ~ stand(log(totP+1)), random = ~ 1|Block, data= bydist1m) 
bladeNcontmodelNP4 <- lme(stand(meanN_Blades) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(bladeNcontmodelN4)$AIC
summary(bladeNcontmodelP4)$AIC #best
summary(bladeNcontmodelNP4)$AIC 

bladePcontmodelN4 <- lme(stand(log(meanP_Blades+1)) ~ stand(log(totN)), random = ~ 1|Block, data= bydist1m) 
bladePcontmodelP4 <- lme(stand(log(meanP_Blades)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
bladePcontmodelNP4 <- lme(stand(meanP_Blades) ~ stand(log(NP+1)), random = ~ 1|Block, data= bydist1m) 
summary(bladePcontmodelN4)$AIC
summary(bladePcontmodelP4)$AIC #best
summary(bladePcontmodelNP4)$AIC 

bladedCcontmodelN4 <- lme(stand(log(meandC_Blades+100)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
bladedCcontmodelP4 <- lme(stand(log(meandC_Blades+100)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
bladedCcontmodelNP4 <- lme(stand(log(meandC_Blades+100)) ~ stand(log(NP)), random = ~ 1|Block, data= bydist1m) 
summary(bladedCcontmodelN4)$AIC 
summary(bladedCcontmodelP4)$AIC #best
summary(bladedCcontmodelNP4)$AIC 

sheathNcontmodelN4 <- lme(stand((meanN_Sheaths)^(1/3)) ~ stand((totN)^(1/5)), random = ~ 1|Block, data= bydist1m) 
sheathNcontmodelP4 <- lme(stand(log(meanN_Sheaths)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
sheathNcontmodelNP4 <- lme(stand(meanN_Sheaths) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(sheathNcontmodelN4)$AIC #best
summary(sheathNcontmodelP4)$AIC #w/in 1
summary(sheathNcontmodelNP4)$AIC

sheathPcontmodelN4 <- lme(stand(log(meanP_Sheaths)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
sheathPcontmodelP4 <- lme(stand(sqrt(meanP_Sheaths)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
sheathPcontmodelNP4 <- lme(stand((meanP_Sheaths)^(1/4)) ~ stand((NP)^(1/3)), random = ~ 1|Block, data= bydist1m) 
summary(sheathPcontmodelN4)$AIC #best
summary(sheathPcontmodelP4)$AIC #w/in 2
summary(sheathPcontmodelNP4)$AIC 

sheathdCcontmodelN4 <- lme(stand(meandC_Sheaths) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
sheathdCcontmodelP4 <- lme(stand(meandC_Sheaths) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
sheathdCcontmodelNP4 <- lme(stand(meandC_Sheaths) ~ stand(sqrt(NP)), random = ~ 1|Block, data= bydist1m) 
summary(sheathdCcontmodelN4)$AIC #w/in 2
summary(sheathdCcontmodelP4)$AIC #best
summary(sheathdCcontmodelNP4)$AIC

rhzCcontmodelN4 <- lme(stand(meanC_Rhizomes) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
rhzCcontmodelP4 <- lme(stand(log(meanC_Rhizomes+10)) ~ stand(sqrt(totP+1)), random = ~ 1|Block, data= bydist1m) 
rhzCcontmodelNP4 <- lme(stand(log(meanC_Rhizomes+10)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rhzCcontmodelN4)$AIC #best

rhzNcontmodelN4 <- lme(stand(log(meanN_Rhizomes)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
rhzNcontmodelP4 <- lme(stand(sqrt(meanN_Rhizomes)) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
rhzNcontmodelNP4 <- lme(stand(log(meanN_Rhizomes)) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 

rhzPcontmodelN4 <- lme(stand(log(meanP_Rhizomes)) ~ stand(log(totN)), random = ~ 1|Block, data= bydist1m) 
rhzPcontmodelP4 <- lme(stand(log(meanP_Rhizomes)) ~ stand(sqrt(totP+1)), random = ~ 1|Block, data= bydist1m) 
rhzPcontmodelNP4 <- lme(stand(meanP_Rhizomes) ~ stand(log(NP+1)), random = ~ 1|Block, data= bydist1m) 
summary(rhzPcontmodelN4)$AIC 
summary(rhzPcontmodelP4)$AIC #best
summary(rhzPcontmodelNP4)$AIC 

rhzdCcontmodelN4 <- lme(stand(meandC_Rhizomes) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
rhzdCcontmodelP4 <- lme(stand(meandC_Rhizomes) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
rhzdCcontmodelNP4 <- lme(stand(meandC_Rhizomes) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 
summary(rhzdCcontmodelN4)$AIC #best
summary(rhzdCcontmodelP4)$AIC #w/in 1
summary(rhzdCcontmodelNP4)$AIC

rtNcontmodelN4 <- lme(stand(log(meanN_Roots+10)) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
rtNcontmodelP4 <- lme(stand(meanN_Roots) ~ stand(totP), random = ~ 1|Block, data= bydist1m) 
rtNcontmodelNP4 <- lme(stand(meanN_Roots) ~ stand(NP), random = ~ 1|Block, data= bydist1m) 

rtPcontmodelN4 <- lme(stand(meanP_Roots^(1/5)) ~ stand(totN^(1/4)), random = ~ 1|Block, data= bydist1m) 
rtPcontmodelP4 <- lme(stand(log(meanP_Roots)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
rtPcontmodelNP4 <- lme(stand(log(meanP_Roots+1)) ~ stand(log(NP)), random = ~ 1|Block, data= bydist1m) 
summary(rtPcontmodelN4)$AIC 
summary(rtPcontmodelP4)$AIC #best
summary(rtPcontmodelNP4)$AIC

rtdCcontmodelN4 <- lme(stand(meandC_Roots) ~ stand(totN), random = ~ 1|Block, data= bydist1m) 
rtdCcontmodelP4 <- lme(stand(sqrt(meandC_Roots+100)) ~ stand(log(totP)), random = ~ 1|Block, data= bydist1m) 
rtdCcontmodelNP4 <- lme(stand(log(meandC_Roots+10)) ~ stand(log(NP)), random = ~ 1|Block, data= bydist1m) 
summary(rtdCcontmodelN4)$AIC #best
summary(rtdCcontmodelP4)$AIC #w/in 2
summary(rtdCcontmodelNP4)$AIC


####results matrix compilation -------
mods1 <- list(prodcontmodelN4, bladeweightcontmodelN4, sheathweightcontmodelN4, rhzweightcontmodelN4, rtweightcontmodelN4,
              shootcontmodelN4, epscontmodelN4, bitescontmodelN4,
              bladeCcontmodelN4, bladeNcontmodelN4, bladePcontmodelN4, bladedCcontmodelN4,
              sheathNcontmodelN4, sheathPcontmodelN4, sheathdCcontmodelN4,
              rhzCcontmodelN4, rhzNcontmodelN4, rhzPcontmodelN4, rhzdCcontmodelN4,
              rtNcontmodelN4, rtPcontmodelN4, rtdCcontmodelN4)
modresultsmatN <- matrix(ncol = 2, nrow = length(mods1))
modestmatN <- matrix(ncol = 2, nrow = length(mods1))

for(i in 1:length(mods1)){
  modresultsmatN[i,] <- summary(mods1[[i]])$tTable[,5] 
}

for(i in 1:length(mods1)){
  modestmatN[i,] <- summary(mods1[[i]])$coefficients$fixed 
}


colnames(modresultsmatN) <- c('intercept', 'Npval')
rownames(modresultsmatN) <- c('prod', 'bladeweight', 'sheathweight', 'rhizomeweight', 'rootweight',
                              'shoots', 'epiphytes', 'bites',
                              'bladeC', 'bladeN', 'bladeP', 'bladedC',
                              'sheathN', 'sheathP', 'sheathdC',
                              'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                              'rootN', 'rootP', 'rootdC')
colnames(modestmatN) <- c('intercept', 'N')
rownames(modestmatN) <- c('prod', 'bladeweight', 'sheathweight', 'rhizomeweight', 'rootweight',
                          'shoots', 'epiphytes', 'bites',
                          'bladeC', 'bladeN', 'bladeP', 'bladedC',
                          'sheathN', 'sheathP', 'sheathdC',
                          'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                          'rootN', 'rootP', 'rootdC')

modresultsmatN2 <- cbind(modresultsmatN, N=modestmatN[,"N"])



mods2 <- list(prodcontmodelP4, bladeweightcontmodelP4, sheathweightcontmodelP4, rhzweightcontmodelP4, rtweightcontmodelP4,
              shootcontmodelP4, epscontmodelP4, bitescontmodelP4,
              bladeCcontmodelP4, bladeNcontmodelP4, bladePcontmodelP4, bladedCcontmodelP4,
              sheathNcontmodelP4, sheathPcontmodelP4, sheathdCcontmodelP4,
              rhzCcontmodelP4, rhzNcontmodelP4, rhzPcontmodelP4, rhzdCcontmodelP4,
              rtNcontmodelP4, rtPcontmodelP4, rtdCcontmodelP4)
modresultsmatP <- matrix(ncol = 2, nrow = length(mods2)) 
modestmatP <- matrix(ncol = 2, nrow = length(mods2))

for(i in 1:length(mods2)){
  modresultsmatP[i,] <- summary(mods2[[i]])$tTable[,5] 
}

for(i in 1:length(mods2)){
  modestmatP[i,] <- summary(mods2[[i]])$coefficients$fixed 
}


colnames(modresultsmatP) <- c('intercept', 'Ppval')
rownames(modresultsmatP) <- c('prod', 'bladeweight', 'sheathweight', 'rhizomeweight', 'rootweight',
                              'shoots', 'epiphytes', 'bites',
                              'bladeC', 'bladeN', 'bladeP', 'bladedC',
                              'sheathN', 'sheathP', 'sheathdC',
                              'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                              'rootN', 'rootP', 'rootdC')
colnames(modestmatP) <- c('intercept', 'P')
rownames(modestmatP) <- c('prod', 'bladeweight', 'sheathweight', 'rhizomeweight', 'rootweight',
                          'shoots', 'epiphytes', 'bites',
                          'bladeC', 'bladeN', 'bladeP', 'bladedC',
                          'sheathN', 'sheathP', 'sheathdC',
                          'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                          'rootN', 'rootP', 'rootdC')

modresultsmatP2 <- cbind(modresultsmatP, P=modestmatP[,"P"])



mods3 <- list(prodcontmodelNP4, bladeweightcontmodelNP4, sheathweightcontmodelNP4, rhzweightcontmodelNP4, rtweightcontmodelNP4,
              shootcontmodelNP4, epscontmodelNP4, bitescontmodelNP4,
              bladeCcontmodelNP4, bladeNcontmodelNP4, bladePcontmodelNP4, bladedCcontmodelNP4,
              sheathNcontmodelNP4, sheathPcontmodelNP4, sheathdCcontmodelNP4,
              rhzCcontmodelNP4, rhzNcontmodelNP4, rhzPcontmodelNP4, rhzdCcontmodelNP4,
              rtNcontmodelNP4, rtPcontmodelNP4, rtdCcontmodelNP4)
modresultsmatNP <- matrix(ncol = 2, nrow = length(mods3)) 
modestmatNP <- matrix(ncol = 2, nrow = length(mods3))

for(i in 1:length(mods3)){
  modresultsmatNP[i,] <- summary(mods3[[i]])$tTable[,5] 
}

for(i in 1:length(mods3)){
  modestmatNP[i,] <- summary(mods3[[i]])$coefficients$fixed 
}


colnames(modresultsmatNP) <- c('intercept', 'NPpval')
rownames(modresultsmatNP) <- c('prod', 'bladeweight', 'sheathweight', 'rhizomeweight', 'rootweight',
                               'shoots', 'epiphytes', 'bites',
                               'bladeC', 'bladeN', 'bladeP', 'bladedC',
                               'sheathN', 'sheathP', 'sheathdC',
                               'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                               'rootN', 'rootP', 'rootdC')
colnames(modestmatNP) <- c('intercept', 'NP')
rownames(modestmatNP) <- c('prod', 'bladeweight', 'sheathweight', 'rhizomeweight', 'rootweight',
                           'shoots', 'epiphytes', 'bites',
                           'bladeC', 'bladeN', 'bladeP', 'bladedC',
                           'sheathN', 'sheathP', 'sheathdC',
                           'rhizomeC', 'rhizomeN', 'rhizomeP', 'rhizomedC',
                           'rootN', 'rootP', 'rootdC')

modresultsmatNP2 <- cbind(modresultsmatNP, NP=modestmatNP[,"NP"])


####make it in ggplot -----------------

matN2 <- data.frame(modresultsmatN2[,2:3])
matN2 <- cbind(matN2, "treament" = "N")
matN2 <- cbind(matN2, "response" = c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 'Epiphytes', 'Bites',
                                     'Blade %C', 'Blade %N', 'Blade %P', 'Blade d13C', 'Sheath %N', 'Sheath %P', 'Sheath d13C',
                                     'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Rhizome d13C', 'Root %N', 'Root %P', 'Root d13C'))
colnames(matN2) <- c("pval", "est", "treatment", "response")
rownames(matN2) <- NULL
matP2 <- data.frame(modresultsmatP2[,2:3])
matP2 <- cbind(matP2, "treament" = "P")
matP2 <- cbind(matP2, "response" = c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 'Epiphytes', 'Bites',
                                     'Blade %C', 'Blade %N', 'Blade %P', 'Blade d13C', 'Sheath %N', 'Sheath %P', 'Sheath d13C',
                                     'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Rhizome d13C', 'Root %N', 'Root %P', 'Root d13C'))
colnames(matP2) <- c("pval", "est", "treatment", "response")
rownames(matP2) <- NULL
matNP2 <- data.frame(modresultsmatNP2[,2:3])
matNP2 <- cbind(matNP2, "treament" = "NP")
matNP2 <- cbind(matNP2, "response" = c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 'Epiphytes', 'Bites',
                                       'Blade %C', 'Blade %N', 'Blade %P', 'Blade d13C', 'Sheath %N', 'Sheath %P', 'Sheath d13C',
                                       'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Rhizome d13C', 'Root %N', 'Root %P', 'Root d13C'))
colnames(matNP2) <- c("pval", "est", "treatment", "response")
rownames(matNP2) <- NULL

long1m <- rbind(matN2, matP2, matNP2)
long1m <- cbind(long1m, "signif" = ifelse(long1m$pval<0.05,"*",""))


long1m$treatment <- factor(long1m$treatment, levels = c("N", "P", "NP"))
long1m$response <- factor(long1m$response, levels=c('AG Production', 'Blade Biomass', 'Sheath Biomass', 'Rhizome Biomass', 'Root Biomass', 'Shoot Density', 
                                                    'Blade d13C', 'Sheath d13C', 'Rhizome d13C', 'Root d13C',
                                                    'Blade %C', 'Blade %N', 'Blade %P', 'Sheath %N', 'Sheath %P',
                                                    'Rhizome %C', 'Rhizome %N', 'Rhizome %P', 'Root %N', 'Root %P',
                                                    'Epiphytes', 'Bites'))

long1m <- cbind(long1m, "bestmodel" = c(1,1,1,0,1,1,1,0,0,0,0,0,1,1,1,1,0,0,1,0,0,1,
                                        1,1,1,1,1,0,1,0,0,1,1,1,1,1,1,0,0,1,1,0,1,1,
                                        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) #these need to be in order of data, not order of levels
long1m <- cbind(long1m, "group" = case_when(grepl('AG Production|Blade Biomass|Sheath Biomass|Rhizome Biomass|Root Biomass|Shoot Density|Blade d13C|Sheath d13C|Rhizome d13C|Root d13C', long1m$response) ~ 1,
                                            grepl('Blade %C|Blade %N|Blade %P|Sheath %N|Sheath %P|Rhizome %C|Rhizome %N|Rhizome %P|Root %N|Root %P|Epiphytes', long1m$response) ~ 2,
                                            grepl('Bites', long1m$response) ~ 3,
                                            TRUE ~ NA))

contmodelplot1m <- ggplot(long1m, aes(x=treatment, y=fct_rev(response), fill=est, label=signif)) +
  geom_tile(color = "black",
            lwd = .5,
            linetype = 1,
            size = 2) +
  geom_text(size=8, vjust=0.8, aes(color=factor(bestmodel))) +
  scale_color_manual(values = c("black", "white")) +
  theme_classic() +
  facet_grid(group~., scales = 'free_y', space = 'free_y') +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text.y = element_blank(), #hides y axis labels for combo figure
        plot.margin = margin(5,20,5,5)) +
  scale_fill_gradient2(low = "#1151FF", mid = "white", high = "#D50E0E", midpoint = 0) +
  guides(fill = guide_colourbar(title = "Effect Size", size = 12),
         color = "none") + 
  labs(tag = "(b)") + theme(plot.tag = element_text(face = "bold")) +
  scale_y_discrete(labels = c('Blade d13C' = expression('Blade'~delta^13*'C'),
                              'Sheath d13C' = expression('Sheath'~delta^13*'C'),
                              'Rhizome d13C' = expression('Rhizome'~delta^13*'C'),
                              'Root d13C' = expression('Root'~delta^13*'C')))

ggsave(filename = "contmatrix.pdf", path="outputs", plot=contmodelplot1m, device = "pdf", width = 3, height = 5.5, units="in", dpi=300)


###c - summary matrix -------------
combosummary <- tibble(treatment = c("fish","fish","fish", "fert","fert","fert", "fish x fert","fish x fert","fish x fert",
                                     "N","N","N","P","P","P","NP","NP","NP"),
                       group = c("Quantity","Quality","Herbivory","Quantity","Quality","Herbivory","Quantity","Quality","Herbivory",
                                 "Quantity","Quality","Herbivory","Quantity","Quality","Herbivory","Quantity","Quality","Herbivory"),
                       total = c(2,4,0,8,7,0,1,0,1,
                                 8,4,0,9,7,0,0,0,1),
                       xgroup = case_when(grepl('fish|fert|fish x fert', treatment) ~ 1,
                                          grepl('N|P|NP', treatment) ~ 2,
                                          TRUE ~ NA)) %>%
  mutate(treatment = factor(treatment, levels = c("fish", "fert", "fish x fert", "N", "P", "NP")),
         group = factor(group, levels = c("Quantity", "Quality", "Herbivory"))) %>%
  ggplot(aes(x=treatment, y=fct_rev(group), label=total)) +
  geom_tile(color = "black",
            lwd = .5,
            linetype = 1,
            fill = "white") +
  geom_text(size=4.5) + 
  theme_classic() +
  facet_grid(.~xgroup, scales = 'free_x', space = 'free_x') +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1, size = 12), 
        axis.text.y = element_text(size = 12),
        plot.margin = margin(5,98,89,50.5)) +
  theme(panel.spacing.x = unit(1.6, "lines")) +
  labs(tag = "(c)") + theme(plot.tag = element_text(face = "bold")) +
  scale_x_discrete(labels = c("fish" = "Fish-derived", "fert" = "Anthropogenic", "fish x fert" = "Fish x Anthro", "NP" = "N:P")) 

ggsave(filename = "matrixsummary.pdf", path="outputs", plot=combosummary, device = "pdf", width = 3.7, height = 2.5, units="in", dpi=300)


###Fig 2 all together -------------

fig2long <- ggpubr::ggarrange(ggarrange(catmodelplot1m,contmodelplot1m, ncol=2, nrow=1),
                              ggarrange(combosummary,ncol=1,nrow = 1), nrow=2, heights = c(2,1))

ggsave(filename = "fig2.pdf", path="outputs", plot=fig2long, device = "pdf", width = 5, height = 9.5, units="in", dpi=300)

##Figure 4 (addition to the scaled seagrass images) -----------------

Cstockcomboplot <- ggplot() + 
  geom_point(aes(x=nuttrt, y=thalC, color=agbg), data= data5nutsLong,
             alpha=0.3, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), size = 3) + 
  geom_pointrange(aes(x=nuttrt, y=Cstockmean, 
                      ymin=Cstockmean-Cstocksd, ymax=Cstockmean+Cstocksd,
                      color=agbg), data= data5nutsSummary,
                  position = position_dodge(width = 0.75), fatten = 9, linewidth = 1) +
  theme_classic() +
  labs(x = "Nutrient treatment",
       y = expression('C Stock ('*gC~m^-2*')'),
       tag = "(b)") + theme(plot.tag = element_text(face = "bold")) +
  scale_x_discrete(labels = c('1' = "B",
                              '2' = "-F-A",
                              '3' = "+F-A",
                              '4' = "-F+A",
                              '5' = "+F+A")) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17, color = "black"),
        legend.title= element_blank(),
        legend.text = element_text(size=17),
        plot.tag = element_text(size=17)) +
  scale_color_manual(values = c("#4C6444", "#8A6240"), labels = c("Aboveground", "Belowground")) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

Cprodcomboplot <- ggplot() + 
  geom_point(aes(x=nuttrt, y=Cprod, color=agbg), data= data5nutsLong,
             alpha=0.3, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), size = 3) + 
  geom_pointrange(aes(x=nuttrt, y=Cprodmean, 
                      ymin=Cprodmean-Cprodsd, ymax=Cprodmean+Cprodsd,
                      color=agbg), data= data5nutsSummary,
                  position = position_dodge(width = 0.75), fatten = 9, linewidth = 1) +
  theme_classic() +
  labs(x = "Nutrient treatment",
       y = expression('C Production ('*gC~m^-2~d^-1*')'),
       tag = "(c)") + theme(plot.tag = element_text(face = "bold")) +
  scale_x_discrete(labels = c('1' = "B",
                              '2' = "-F-A",
                              '3' = "+F-A",
                              '4' = "-F+A",
                              '5' = "+F+A")) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17, color = "black"),
        legend.title= element_blank(),
        legend.text = element_text(size=17),
        plot.tag = element_text(size=17)) +
  scale_color_manual(values = c("#4C6444", "#8A6240"), labels = c("Aboveground", "Belowground")) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

Cturncomboplot <- ggplot() + 
  geom_point(aes(x=nuttrt, y=Cturn, color=agbg), data= data5nutsLong,
             alpha=0.3, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), size = 3) + 
  geom_pointrange(aes(x=nuttrt, y=Cturnmean, 
                      ymin=Cturnmean-Cturnsd, ymax=Cturnmean+Cturnsd,
                      color=agbg), data= data5nutsSummary,
                  position = position_dodge(width = 0.75), fatten = 9, linewidth = 1) +
  theme_classic() +
  labs(x = "Nutrient treatment",
       y = expression('C Turnover ('*d^-1*')'),
       tag = "(d)") + theme(plot.tag = element_text(face = "bold")) +
  scale_x_discrete(labels = c('1' = "B",
                              '2' = "-F-A",
                              '3' = "+F-A",
                              '4' = "-F+A",
                              '5' = "+F+A")) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=17),
        plot.tag = element_text(size=17)) +
  scale_color_manual(values = c("#4C6444", "#8A6240"), labels = c("Aboveground", "Belowground")) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

plotsforfig4 <- ggpubr::ggarrange(Cstockcomboplot, Cprodcomboplot, Cturncomboplot,
                  ncol = 3, nrow = 1,
                  common.legend = T, legend = "bottom", widths = c(1,1,1.05))

ggsave(filename = "fig4bcd.pdf", path="outputs", plot=plotsforfig4, device = "pdf", width = 13.5, height = 5, units="in", dpi=300)


##Figure 5 -----------------
agCturn_comboplot <- databydist %>%
  ggplot(aes(x=avgheight, y=agCturn, fill=shootsm2)) +
  geom_point(pch=21, colour="white", stroke=.2, size=5) + 
  theme_classic() +
  labs(y=expression('AG C Turnover ('*day^-1*')'),
       x='Blade Height (mm)',
       fill=expression('Shoots'~m^-2),
       tag = "(a)") + theme(plot.tag = element_text(face = "bold")) +
  theme(legend.position="bottom",
        axis.text = element_text(size = 13.5), 
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.tag = element_text(size=15)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA))

bgCturn_comboplot <- databydist %>%
  ggplot(aes(x=avgheight, y=bgCturn, fill=shootsm2)) +
  geom_point(pch=21, colour="white", stroke=.2, size=5) + 
  theme_classic() +
  labs(y=expression('BG C Turnover ('*day^-1*')'),
       x='Blade Height (mm)',
       fill=expression('Shoots'~m^-2),
       tag = "(b)") + theme(plot.tag = element_text(face = "bold")) +
  theme(legend.position="bottom",
        axis.text = element_text(size = 13.5), 
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.tag = element_text(size=15)) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

totCturn_comboplot <- databydist %>%
  ggplot(aes(x=avgheight, y=totCturn, fill=shootsm2)) +
  geom_point(pch=21, colour="white", stroke=.2, size=5) + #stroke changes size of border
  theme_classic() +
  labs(y=expression('Total C Turnover ('*day^-1*')'),
       x='Blade Height (mm)',
       fill=expression('Shoots'~m^-2),
       tag = "(c)") + theme(plot.tag = element_text(face = "bold")) +
  theme(legend.position="bottom",
        axis.text = element_text(size = 13.5), 
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.tag = element_text(size=15)) +
  theme(panel.border = element_rect(colour = "black", fill=NA))

corrsforfig5 <- ggpubr::ggarrange(agCturn_comboplot, bgCturn_comboplot, totCturn_comboplot,
                                   ncol = 3, nrow = 1, common.legend = T, legend="bottom")

ggsave(filename = "fig5.pdf", path="outputs", plot=corrsforfig5, device = "pdf", width = 14, height = 4.5, units="in", dpi=300)





