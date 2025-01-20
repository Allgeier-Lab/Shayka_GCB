###########
#This code takes initial data and calculates biomass, growth, and productivity
##########

##Load libraries -------------
library(tidyverse)


##Load data -----------------

agdata <- read_csv('data/AGdata.csv') %>%
  select(Days, CoreID, Reef, Transect, Distance, Core, Shoot, Blade, Growth, Length, Width, Blade_dry_weight, Epiphyte_dry_weight, Bitten, Bites) %>%
  filter(!is.na(Reef)) %>%
  dplyr::mutate(Epiphyte_dry_weight = replace_na(Epiphyte_dry_weight, 0)) #this replaces NAs in epiphyte weight column with zeros so that it accurately reflects 0 grams of epiphytes

bgdata <- read_csv('data/BGdata.csv') %>%
  filter(!is.na(Reef))

agnutdata <- read_csv('data/AGnutrientdata.csv') %>%
  tidyr::separate(CoreID, c('Reef', 'Transect', 'Distance', 'Core'))

bgnutdata <- read_csv('data/BGnutrientdata.csv') %>%
  tidyr::separate(SampleID, c('Reef', 'Transect', 'Distance', 'Core', 'ContentAbrev')) #some samples have an additional number after the "ContentAbrev"; R will discard that; that's ok!

agpdata <- read_csv('data/AGPdata.csv')

bgpdata <- read_csv('data/BGPdata.csv') %>%
  tidyr::separate(Contents, c('ContentAbrev')) #this discards the additional number after the "ContentAbrev", just like in the bgnutdata code above

bgproddata <- read_csv('data/bladeareaweight_BG.csv') %>%
  select(Reef, Transect, Distance, Core, bgprod)

syrdata <- read_csv('data/Syrdata.csv') %>%
  select(CoreID, Contents, Dry_weight) %>%
  tidyr::separate(CoreID, c('Reef', 'Transect', 'Distance', 'Core'))

##Data prep -------------

###nutrient data -------------
agcn <- agnutdata %>%
  mutate(Transect = as.numeric(Transect),
         Distance = as.numeric(Distance),
         Core = as.numeric(Core)) %>%
  group_by(Reef, Transect, Distance, Core) %>%
  summarise(meanC_Blades = mean(C,na.rm = T),
            meandC_Blades = mean(d13C,na.rm = T),
            meanN_Blades = mean(N,na.rm = T),
            meandN_Blades = mean(d15N,na.rm = T))

agp <- agpdata %>%
  group_by(Reef, Transect, Distance, Core) %>%
  summarise(meanP_Blades = mean(P,na.rm = T))
  

agnuts <- left_join(agcn, agp, by = c("Reef", "Transect", "Distance", "Core"))

bgcn <- bgnutdata %>%
  mutate(Transect = as.numeric(Transect),
         Distance = as.numeric(Distance),
         Core = as.numeric(Core)) %>%
  mutate(Contents = case_when(ContentAbrev == "S" ~ "Sheaths",
                              ContentAbrev == "Rz" ~ "Rhizomes",
                              ContentAbrev == "Rt" ~ "Roots")) %>%
  group_by(Reef, Transect, Distance, Core, Contents) %>%
  summarise(meanC = mean(C,na.rm = T),
            meandC = mean(d13C,na.rm = T),
            meanN = mean(N,na.rm = T),
            meandN = mean(d15N,na.rm = T))

bgp <- bgpdata %>%
  mutate(Contents = case_when(ContentAbrev == "S" ~ "Sheaths",
                              ContentAbrev == "Rz" ~ "Rhizomes",
                              ContentAbrev == "Rt" ~ "Roots")) %>%
  group_by(Reef, Transect, Distance, Core, Contents) %>%
  summarise(meanP = mean(P,na.rm = T))

bgnuts <- left_join(bgp, bgcn, by = c("Reef", "Transect", "Distance", "Core", "Contents"))

###growth and morphology data -------------
#for first two cores with epiphytes by core
growthdata1 <- agdata %>%
  filter(CoreID %in% c("YM4B_1_1_3", "YM4A_1_20_3")) %>% #the %in% function instead of == here lets you find rows that match either of the cores
  mutate(garea = Growth * Width, #in mm2 #area of new growth during Days time period
         totarea = Length * Width) %>% #in mm2
  group_by(Reef, Transect, Distance, Core, Shoot) %>%
  summarise(gshootarea = sum(garea), #in mm2 #don't remove NAs
            totshootarea = sum(totarea), #in mm2 #don't remove NAs; search for them in L and W columns above to make sure there aren't any
            totshootweight = mean(Blade_dry_weight, na.rm=T), #in g #this is mean b/c all blades were weighed together from a shoot and same shoot weight value is recorded for every blade
            days = mean(Days, na.rm=T), #removes NA values where there are some numbers and some NAs; otherwise just get NAs
            epweight = mean(Epiphyte_dry_weight), #in g #total epiphyte weight by shoot; use mean b/c epiphyte weights are copied across all blades of a core for these two cores (DO NOT SUM)
            shootheight = max(Length), #in mm #height of tallest blade in each shoot (even if bitten)
            notop = sum(Bitten == "y", na.rm=T), #number of blades missing top (by shoot)
            blades = n(), #number of blades (by shoot)
            bites = mean(Bites) #number of bites (in each shoot); use mean b/c bites are copied across all blades of a shoot (DO NOT SUM)
            ) %>% 
  mutate(growthratearea = gshootarea/days, #in mm2/day for each shoot
         gshootweight = gshootarea* (totshootweight/totshootarea), #in g #weight of new growth during Days time period for each shoot
         growthrateweight = gshootweight/days,  #in g/day for each shoot
         notopratio = notop/blades,
         bitesperarea = bites/totshootarea) #in number/mm2

databycore1 <- growthdata1 %>%
  group_by(Reef, Transect, Distance, Core) %>%
  summarise(totgrowtharea = sum(growthratearea, na.rm=T), # mm2/day by core #THIS ONLY GIVES THE MEASURED GROWTH, NOT THE TOTAL CORE'S GROWTH; TOTAL CORE GROWTH CALCULATED BELOW
            avggrowtharea = mean(growthratearea, na.rm=T), #avg mm2/day/shoot within a core
            totgrowthweight = sum(growthrateweight, na.rm=T), # g/day by core #SAME AS COMMENT ABOVE
            avggrowthweight = mean(growthrateweight, na.rm=T), #avg g/day/shoot within a core
            totepweightperarea = mean(epweight, na.rm=T)/sum(totshootarea, na.rm=T), # g/mm2 by core #using mean b/c eps were grouped by core for these two cores
            totepweightperweight = mean(epweight, na.rm=T)/sum(totshootweight, na.rm=T), #g/g by core #using mean b/c eps were grouped by core for these two cores
            avgbites = mean(bitesperarea), #avg bites/mm2 by core
            totbites = sum(bites), #total bites by core
            avgnotop = mean(notopratio), #average percent of blades missing top per shoot (aka ratio of top missing to full blades per shoot) per core
            avgbladearea = mean(totshootarea), #average blade area per shoot, by core
            totbladearea = sum(totshootarea), #total blade area by core
            avgbladeweight = mean(totshootweight), #average shoot weight (blades only), by core
            totbladeweight = sum(totshootweight), #in g #total blade weight by core aka blade biomass for each core
            avgheight = mean(shootheight), #average height of tallest blade in each shoot, by core
            shoots = sum(Shoot != "F"), #of shoots (by core) not including floaters (F)
            .groups = "drop") #explicitly tells summarize to drop groups at the end


#for all other cores with epiphytes by shoot
growthdata2 <- agdata %>%
  filter(CoreID != "YM4B_1_1_3", CoreID != "YM4A_1_20_3") %>% 
  mutate(garea = Growth * Width, #in mm2 #area of new growth during Days time period
         totarea = Length * Width) %>% #in mm2
  group_by(Reef, Transect, Distance, Core, Shoot) %>%
  summarise(gshootarea = sum(garea), #in mm2 #don't remove NAs
            totshootarea = sum(totarea), #in mm2 #don't remove NAs; search for them in L and W columns above to make sure there aren't any
            totshootweight = mean(Blade_dry_weight, na.rm=T), #in g #this is mean b/c all blades were weighed together from a shoot and same shoot weight value is recorded for every blade
            days = mean(Days, na.rm=T), #removes NA values where there are some numbers and some NAs; otherwise just get NAs
            epweight = mean(Epiphyte_dry_weight), #in g #total epiphyte weight by shoot; use mean b/c epiphyte weights are copied across all blades of a shoot (DO NOT SUM)
            shootheight = max(Length), #in mm #height of tallest blade in each shoot (even if bitten)
            notop = sum(Bitten == "y", na.rm=T), #number of blades missing top (by shoot)
            blades = n(), #number of blades (by shoot)
            bites = mean(Bites) #number of bites (in each shoot); use mean b/c bites are copied across all blades of a shoot (DO NOT SUM)
            ) %>% 
  mutate(growthratearea = gshootarea/days, #in mm2/day for each shoot
         gshootweight = gshootarea* (totshootweight/totshootarea), #in g #weight of new growth during Days time period for each shoot
         growthrateweight = gshootweight/days,  #in g/day for each shoot
         notopratio = notop/blades,
         epweightperarea = epweight/totshootarea,  #in g/mm2
         epweightperweight = epweight/totshootweight, #in g/g
         bitesperarea = bites/totshootarea) #in number/mm2



databycore <- growthdata2 %>%
  group_by(Reef, Transect, Distance, Core) %>%
  summarise(totgrowtharea = sum(growthratearea, na.rm=T), # mm2/day by core #THIS ONLY GIVES THE MEASURED GROWTH, NOT THE TOTAL CORE'S GROWTH; TOTAL CORE GROWTH CALCULATED BELOW
            avggrowtharea = mean(growthratearea, na.rm=T), #avg mm2/day/shoot within a core
            totgrowthweight = sum(growthrateweight, na.rm=T), # g/day by core #SAME AS COMMENT ABOVE
            avggrowthweight = mean(growthrateweight, na.rm=T), #avg g/day/shoot within a core
            totepweightperarea = sum(epweight, na.rm=T)/sum(totshootarea, na.rm=T), # g/mm2 by core
            totepweightperweight = sum(epweight, na.rm=T)/sum(totshootweight, na.rm=T), #g/g by core 
            avgbites = mean(bitesperarea), #avg bites/mm2 by core
            totbites = sum(bites), #total bites by core
            avgnotop = mean(notopratio), #average percent of blades missing top per shoot (aka ratio of top missing to full blades per shoot) per core
            avgbladearea = mean(totshootarea), #average blade area per shoot, by core
            totbladearea = sum(totshootarea), #total blade area by core
            avgbladeweight = mean(totshootweight), #average shoot weight (blades only), by core
            totbladeweight = sum(totshootweight), #in g #total blade weight by core aka blade biomass for each core
            avgheight = mean(shootheight), #average height of tallest blade (mm) in each shoot, by core
            shoots = sum(Shoot != "F"), #of shoots (by core) not including floaters (F)
            .groups = "drop") %>% #explicitly tells summarize to drop groups at the end
  rbind(databycore1) %>% #adds in the rows for the two cores that were pulled out at the beginning
  tidyr::separate(col=Reef, into=c(NA, NA, NA, "Block","Treatment"), sep = "", remove=F) %>% #this separates a character column into multiple character columns, sep="" separates every character, remove=F keeps input column
  mutate(shoots = replace_na(shoots, 0), #this replaces NAs in shoots column with zeros so that it accurately reflects 0 shoots
         coregrowtharea = avggrowtharea*shoots, #(mm2/day/shoot)*shoot = mm2/day for each core
         coregrowthweight = avggrowthweight*shoots,#(g/day/shoot)*shoot = g/day for each core
         coreprodarea = coregrowtharea*100, #100cm2=0.01m2 #mm2/m2/day
         coreprodweight = coregrowthweight*100, #g/m2/day
         corebladeaream2 = totbladearea*100, #blade area per m2 by core #mm2/m2
         corebladeweightm2 = totbladeweight*100, #blade weight per m2 by core #g/m2
         shootsm2 = shoots*100, #shoots/m2 by core
         corebitesperarea = totbites/totbladearea) %>% #bites/mm2 - can compare to avgbites to see how diff
  mutate(prodtoblade = coreprodweight/corebladeweightm2) %>%
  left_join(agnuts, by = c("Reef", "Transect", "Distance", "Core"))
         


bgweightdata <- bgdata %>%
  select(CoreID, Reef, Transect, Distance, Core, Contents, Dry_weight) %>%
  group_by(Reef, Transect, Distance, Core, Contents) %>% 
  summarise(totweight = sum(Dry_weight, na.rm=T)) %>% #sum of weight of each component by core #in g
  left_join(bgnuts, by = c("Reef", "Transect", "Distance", "Core", "Contents")) %>%
  pivot_wider(names_from = Contents, values_from = c(totweight,meanC,meandC,meanN,meandN,meanP)) %>%
  dplyr::mutate(totweight_Flowers = replace_na(totweight_Flowers, 0),
                totweight_Detritus = replace_na(totweight_Detritus, 0),
                detritusm2 = totweight_Detritus*100, #detritus weight per m2 by core #g/m2
                rhizomesm2 = totweight_Rhizomes*100, #rhizome weight per m2 by core #g/m2
                rootsm2 = totweight_Roots*100, #root weight per m2 by core #g/m2
                sheathsm2 = totweight_Sheaths*100) #sheath weight per m2 by core #g/m2

syrweightdata <- syrdata %>%
  mutate(Transect = as.numeric(Transect),
         Distance = as.numeric(Distance),
         Core = as.numeric(Core)) %>%
  group_by(Reef, Transect, Distance, Core, Contents) %>%
  summarise(syrtotweight = sum(Dry_weight, na.rm=T)) %>% #sum of weight of each component by core #in g
  pivot_wider(names_from = Contents, values_from = c(syrtotweight)) %>%
  dplyr::mutate(syrtotweight_Blades = replace_na(Blades, 0),
                syrtotweight_Rhizomes = replace_na(Rhizomes, 0),
                syrtotweight_Roots = replace_na(Roots, 0),
                syrbladesm2 = syrtotweight_Blades*100, #blade weight per m2 by core #g/m2
                syrrhizomesm2 = syrtotweight_Rhizomes*100, #rhizome weight per m2 by core #g/m2
                syrrootsm2 = syrtotweight_Roots*100) %>% #root weight per m2 by core #g/m2
  select(Reef, Transect, Distance, Core, syrtotweight_Blades, syrtotweight_Rhizomes, syrtotweight_Roots, syrbladesm2, syrrhizomesm2, syrrootsm2)
  



allcoredata <- left_join(databycore, bgweightdata, by = c("Reef", "Transect", "Distance", "Core")) %>%
  left_join(syrweightdata, by = c("Reef", "Transect", "Distance", "Core")) %>%
  left_join(bgproddata, by = c("Reef", "Transect", "Distance", "Core")) %>%
  mutate(fish = case_when(Treatment == "A" ~ "low",
                          Treatment == "B" ~ "high",
                          Treatment == "C" ~ "low",
                          Treatment == "D" ~ "high"),
         fert = case_when(Treatment == "A" ~ "no",
                          Treatment == "B" ~ "no",
                          Treatment == "C" ~ "yes",
                          Treatment == "D" ~ "yes")) %>%
  mutate(bgweight = rhizomesm2 + rootsm2, #g/m2 #total biomass of roots and rhizomes
         agweight = sheathsm2 + corebladeweightm2, #g/m2 #total biomass of blades and sheaths
         bgshthweight = sheathsm2 + rhizomesm2 + rootsm2, #g/m2 #total biomass of roots, rhizomes, and sheaths
         agtobg = agweight/bgweight, #aboveground to belowground biomass ratio
         bladetobg = corebladeweightm2/bgweight, #blade to belowground biomass ratio
         bladetoelse = corebladeweightm2/(bgweight + sheathsm2), #ratio of blades to everything else biomass
         agtorz = agweight/rhizomesm2, #ratio of aboveground to rhizome biomass
         bladetorz = corebladeweightm2/rhizomesm2, #ratio of blade to rhizome biomass
         agtoroot = agweight/rootsm2, #ratio of aboveground to root biomass
         bladetoroot = corebladeweightm2/rootsm2, #ratio of blade to root biomass
         roottorz = rootsm2/rhizomesm2, #ratio of root to rhizome biomass
         bladetosheath = corebladeweightm2/sheathsm2) %>% #ratio of blade to sheath biomass
  mutate(totbio = corebladeweightm2+sheathsm2+rhizomesm2+rootsm2, #g/m2 #total biomass
         bladeC = corebladeweightm2 * meanC_Blades/100, #gC/m2 in blades
         sheathC = sheathsm2 * meanC_Sheaths/100, #gC/m2 in sheaths
         rzC = rhizomesm2 * meanC_Rhizomes/100, #gC/m2 in rhizomes
         rtC = rootsm2 * meanC_Roots/100, #gC/m2 in roots
         bladeN = corebladeweightm2 * meanN_Blades/100, #gN/m2 in blades
         sheathN = sheathsm2 * meanN_Sheaths/100, #gN/m2 in sheaths
         rzN = rhizomesm2 * meanN_Rhizomes/100, #gN/m2 in rhizomes
         rtN = rootsm2 * meanN_Roots/100, #gN/m2 in roots
         bladeP = corebladeweightm2 * meanP_Blades/100, #gP/m2 in blades
         sheathP = sheathsm2 * meanP_Sheaths/100, #gP/m2 in sheaths
         rzP = rhizomesm2 * meanP_Rhizomes/100, #gP/m2 in rhizomes
         rtP = rootsm2 * meanP_Roots/100) %>% #gP/m2 in Thalassia
  mutate(agthalC = bladeC + sheathC, #gC/m2 in AG Thalassia 
         agthalN = bladeN + sheathN, #gN/m2 in AG Thalassia 
         agthalP = bladeP + sheathP, #gP/m2 in AG Thalassia
         bgthalC = rzC + rtC, #gC/m2 in BG Thalassia
         bgthalN = rzN + rtN, #gN/m2 in BG Thalassia
         bgthalP = rzP + rtP, #gP/m2 in BG Thalassia
         totthalC = bladeC + sheathC + rzC + rtC, #gC/m2 in Thalassia
         totthalN = bladeN + sheathN + rzN + rtN, #gN/m2 in Thalassia
         totthalP = bladeP + sheathP + rzP + rtP) %>% #gP/m2 in Thalassia
  mutate(agCprod = coreprodweight * (meanC_Blades/100), #gC m-2 d-1 #AG C production
         bgCprod = bgprod * ((rootsm2/bgweight)*(meanC_Roots/100) + (rhizomesm2/bgweight)*(meanC_Rhizomes/100))) %>% #gC m-2 d-1 #BG C production
  mutate(totCprod = agCprod + bgCprod, #gC m-2 d-1 #total C production
         agCturn = agCprod/agthalC, #d-1 #C turnover rate for AG
         bgCturn = bgCprod/bgthalC) %>% #d-1 #C turnover rate for BG
  mutate(totCturn = totCprod/totthalC) #d-1 #C turnover rate for Thalassia




###reef nutrient treatments -------------
reefnuts <- read_csv('data/reef.fish.nutrients.csv') %>%
  tidyr::separate_wider_delim(cols = reef.date, delim = "-", names = c("reef", "date")) %>%
  tidyr::separate_wider_position(cols = reef, widths = c(Block = 1, Treatment = 1)) %>%
  mutate(Treatment = case_when(Treatment == "a" ~ "A",
                               Treatment == "b" ~ "B",
                               Treatment == "c" ~ "C",
                               Treatment == "d" ~ "D")) %>%
  group_by(Block, Treatment) %>%
  summarise(avgfishN = mean(n), #units
            avgfishP = mean(p)) %>% #units
  mutate(fertN = case_when(Treatment == "A" ~ 0, #g/day/reef
                           Treatment == "B" ~ 0,
                           Treatment == "C" ~ 2.7,
                           Treatment == "D" ~ 2.7),
         fertP = case_when(Treatment == "A" ~ 0, #g/day/reef
                           Treatment == "B" ~ 0,
                           Treatment == "C" ~ 0.39,
                           Treatment == "D" ~ 0.39)) %>%
  mutate(totN = avgfishN + fertN, #g/day/reef
         totP = avgfishP + fertP, #g/day/reef
         NP = (totN/14)/(totP/30.97)) #molar ratio

allcoredata2 <- left_join(allcoredata,reefnuts,by = c("Block", "Treatment"))


###all data summarized/combined by distance -------------
databydist <- allcoredata %>%
  group_by(Reef, Block, Treatment, Transect, Distance, fish, fert) %>% #below are all the averages of core values by distance
  summarise(coregrowtharea = mean(coregrowtharea, na.rm=T), #avg of whole core mm2/day at dist
            coregrowthweight = mean(coregrowthweight, na.rm=T), #avg of whole core g/day at dist
            coreprodarea = mean(coreprodarea, na.rm=T), #avg of whole core mm2/m2/day at dist
            coreprodweight = mean(coreprodweight, na.rm=T), #avg of whole core g/m2/day at dist
            corebladeaream2 = mean(corebladeaream2, na.rm=T), #avg mm2/m2
            corebladeweightm2 = mean(corebladeweightm2, na.rm=T), #avg g/m2
            bgprod = mean(bgprod, na.rm=T), #avg of whole core g/m2/day at dist
            avgheight = mean(avgheight, na.rm=T), #avg btw cores of avg height of tallest blade (mm) w/in cores
            shootsm2 = mean(shootsm2, na.rm=T), #avg shoots/m2
            totepweightperarea = mean(totepweightperarea, na.rm=T), #avg g/g
            totepweightperweight = mean(totepweightperweight, na.rm=T), #avg g/mm2
            corebitesperarea = mean(corebitesperarea, na.rm=T), #avg of core level bites/mm2 by dist
            sheathsm2 = mean(sheathsm2, na.rm=T), #avg g/m2 by dist
            rhizomesm2 = mean(rhizomesm2, na.rm=T), #g/m2
            rootsm2 = mean(rootsm2, na.rm=T), #g/m2
            detritusm2 = mean(detritusm2, na.rm=T), #g/m2
            syrbladesm2 = mean(syrbladesm2, na.rm=T), #g/m2
            syrrhizomesm2 = mean(syrrhizomesm2, na.rm=T), #g/m2
            syrrootsm2 = mean(syrrootsm2, na.rm=T), #g/m2
            bgweight = mean(bgweight, na.rm=T), #g/m2 #avg biomass of roots and rhizomes
            agweight = mean(agweight, na.rm=T), #g/m2 #avg biomass of blades and sheaths
            bgshthweight = mean(bgshthweight, na.rm=T), #g/m2 #avg biomass of roots, rhizomes, and sheaths
            totweight = mean(totbio, na.rm=T), #g/m2 #avg biomass of Thalassia
            meanC_Blades = mean(meanC_Blades, na.rm=T),
            meandC_Blades = mean(meandC_Blades, na.rm=T),
            meanN_Blades = mean(meanN_Blades, na.rm=T),
            meandN_Blades = mean(meandN_Blades, na.rm=T),
            meanP_Blades = mean(meanP_Blades, na.rm=T),
            meanC_Sheaths = mean(meanC_Sheaths, na.rm=T),
            meandC_Sheaths = mean(meandC_Sheaths, na.rm=T),
            meanN_Sheaths = mean(meanN_Sheaths, na.rm=T),
            meandN_Sheaths = mean(meandN_Sheaths, na.rm=T),
            meanP_Sheaths = mean(meanP_Sheaths, na.rm=T),
            meanC_Rhizomes = mean(meanC_Rhizomes, na.rm=T),
            meandC_Rhizomes = mean(meandC_Rhizomes, na.rm=T),
            meanN_Rhizomes = mean(meanN_Rhizomes, na.rm=T),
            meandN_Rhizomes = mean(meandN_Rhizomes, na.rm=T),
            meanP_Rhizomes = mean(meanP_Rhizomes, na.rm=T),
            meanC_Roots = mean(meanC_Roots, na.rm=T),
            meandC_Roots = mean(meandC_Roots, na.rm=T),
            meanN_Roots = mean(meanN_Roots, na.rm=T),
            meandN_Roots = mean(meandN_Roots, na.rm=T),
            meanP_Roots = mean(meanP_Roots, na.rm=T),
            bladeC = mean(bladeC, na.rm=T), #gC/m2 #avg of total C in blades
            sheathC = mean(sheathC, na.rm=T), #gC/m2 #avg of total C in sheaths
            rzC = mean(rzC, na.rm=T), #gC/m2 #avg of total C in rhizomes
            rtC = mean(rtC, na.rm=T), #gC/m2 #avg of total C in roots
            bladeN = mean(bladeN, na.rm=T), #gN/m2 #avg of total N in blades
            sheathN = mean(sheathN, na.rm=T), #gN/m2 #avg of total N in sheaths
            rzN = mean(rzN, na.rm=T), #gN/m2 #avg of total N in rhizomes
            rtN = mean(rtN, na.rm=T), #gN/m2 #avg of total N in roots
            bladeP = mean(bladeP, na.rm=T), #gP/m2 #avg of total P in blades
            sheathP = mean(sheathP, na.rm=T), #gP/m2 #avg of total P in sheaths
            rzP = mean(rzP, na.rm=T), #gP/m2 #avg of total P in rhizomes
            rtP = mean(rtP, na.rm=T), #gP/m2 #avg of total P in roots
            agthalC = mean(agthalC, na.rm=T), #gC/m2 #avg of total C in AG Thalassia
            agthalN = mean(agthalN, na.rm=T), #gN/m2 #avg of total C in AG Thalassia
            agthalP = mean(agthalP, na.rm=T), #gP/m2 #avg of total C in AG Thalassia
            bgthalC = mean(bgthalC, na.rm=T), #gC/m2 #avg of total C in BG Thalassia
            bgthalN = mean(bgthalN, na.rm=T), #gN/m2 #avg of total C in BG Thalassia
            bgthalP = mean(bgthalP, na.rm=T), #gP/m2 #avg of total C in BG Thalassia
            totthalC = mean(totthalC, na.rm=T), #gC/m2 #avg of total C in Thalassia
            totthalN = mean(totthalN, na.rm=T), #gN/m2 #avg of total N in Thalassia
            totthalP = mean(totthalP, na.rm=T), #gP/m2 #avg of total P in Thalassia
            agCprod = mean(agCprod, na.rm=T), #gC m-2 d-1 #avg of AG C production
            bgCprod = mean(bgCprod, na.rm=T), #gC m-2 d-1 #avg of BG C production
            totCprod = mean(totCprod, na.rm=T), #gC m-2 d-1 #avg of total C production
            agCturn = mean(agCturn, na.rm=T), #d-1 #avg of C turnover rate for AG (blades only)
            bgCturn = mean(bgCturn, na.rm=T), #d-1 #avg of C turnover rate for BG
            totCturn = mean(totCturn, na.rm=T) #d-1 #avg of C turnover rate for Thalassia
  ) %>%
  left_join(reefnuts, by = c("Block", "Treatment")) %>%
  mutate(labels = case_when(Treatment == "A" ~ "-F-A",
                            Treatment == "B" ~ "+F-A",
                            Treatment == "C" ~ "-F+A",
                            Treatment == "D" ~ "+F+A"))



## Exporting data -------

write.csv(databydist, file = "processeddata/alldistdata.csv", row.names = FALSE)

