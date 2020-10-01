



require(tidyr)
require(dplyr)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

## Set parameters
res <- 0.4 #leaf resorption 0.4 for normal, 1 for upper bound, res is amount lost here
resr <- 0.73 #root resporption 0.73 for normal, 1 for upper bound, res is amount lost here

res <- 1
resr <- 1
source("scripts/ndfa_bootstrapping.R") #go to this file, regular for normal, change to 1s for upper bound

spdata <- read.csv("ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")
spdata <- dplyr::select(spdata, SPCD, GENUS, JENKINS_SPGRPCD)


# using RGR from same species ------------------------------------

# Load datafrom ndfa_GRnearby.R
agb_1_2 <- readRDS("output/agb_1_2.RDS") #biomass measured and predicted for t2, fixers only

agb_1_2 <- agb_1_2[agb_1_2$dbh2pred>=2,] #jenkins equations only work for trees dbh>2

agb_1_2 <- merge(agb_1_2, spdata, by.x="spcd", by.y="SPCD", all.x=T, all.y=F)

# assign evergreen or deciduous
agb_1_2$EV0DEC1 <- NA
agb_1_2[agb_1_2$JENKINS_SPGRPCD <= 5, ]$EV0DEC1 <- 0 #no evergreen fixers (here only have fixers)
agb_1_2[agb_1_2$JENKINS_SPGRPCD > 5, ]$EV0DEC1 <- 1

# assign foliar resorption for evergreen/deciduous fixers/nonfixers
agb_1_2$resf <- NA
agb_1_2[agb_1_2$EV0DEC1==0 & agb_1_2$FIX==0, ]$resf <- 0.545 #from Vergutz 2012
agb_1_2[agb_1_2$EV0DEC1==0 & agb_1_2$FIX==1, ]$resf <- 0.418 #from Vergutz 2012
agb_1_2[agb_1_2$EV0DEC1==1 & agb_1_2$FIX==0, ]$resf <- 0.597 #from Vergutz 2012
agb_1_2[agb_1_2$EV0DEC1==1 & agb_1_2$FIX==1, ]$resf <- 0.495 #from Vergutz 2012


## Get biomass allocations
# hardwood parameters (Jenkins, 2003)
B0.foliage <- -4.0813
B1.foliage <- 5.8816
B0.roots <- -1.6911
B1.roots <- 0.8160
B0.bark <- -2.0129
B1.bark <- -1.6805
B0.wood <- -0.3065
B1.wood <- -5.4240

# agb fraction from Jenssen
b.foliage <- exp(B0.foliage+(B1.foliage/agb_1_2$dbh2pred))
b.bark <- exp(B0.bark+(B1.bark/agb_1_2$dbh2pred))
b.wood <- exp(B0.wood+(B1.wood/agb_1_2$dbh2pred))
b.root <- exp(B0.roots+(B1.roots/agb_1_2$dbh2pred))
b.other <- 1 - (b.foliage+b.bark+b.wood)

agb_1_2[agb_1_2$agb2pred < 0.5,]$agb2pred <- 0.05 #to prevent 0 in denominator
agb_1_2$agbd <- agb_1_2$agb2pred - agb_1_2$agb #get difference in agb
agb_1_2[agb_1_2$agbd < 0,]$agbd <- 0.05 #if no change assign small increment

#get change in biomass (kgC in foliage and wood)
plot_agb <- agb_1_2 %>% 
  mutate(agb.foliage = exp(B0.foliage+(B1.foliage/dbh2pred))*agb2pred,
         agb.wood = (agbd),
         bgb.root = exp(B0.roots+(B1.roots/dbh2pred))*agb2pred)
# should agb.wood = ((agbd)/t)*(b.wood+b.bark)
plot_agb <- plot_agb[plot_agb$tpha > 0.5,]
plot_agb <- plot_agb[plot_agb$dbh < 200,] #weird 2 outliers driving calcs
plot_agb <- plot_agb[!is.na(plot_agb$dbh2pred),]

#use C:N ratios since biomass is in KgC
#cn.foliage <- 15  #Ferlian et al, 2017
#cn.wood <- 465    #Johnson et al, 2014 (from sugar maple)
cn.fr <- 43       #Gordon & Jackson, 2000
cn.foliage <- 35.1 #McGrody 2004
cn.wood <- 350     #Du and Vries

# foliage (full crown of foliage)
fol <- plot_agb %>% 
  mutate(foliagen = ifelse(Genus %in% c("Acacia"), agb.foliage*ndfas$aca*(1/cn.foliage)*tpha*resf, 
                           ifelse(Genus %in% c("Albizia"), agb.foliage*ndfas$alb*(1/cn.foliage)*tpha*resf,
                                  ifelse(Genus %in% c("Alnus"), agb.foliage*ndfas$ald*(1/cn.foliage)*tpha*resf,
                                         ifelse(Genus %in% c("Cercocarpus"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha*resf,
                                                ifelse(Genus %in% c("Prosopis"), agb.foliage*ndfas$pro*(1/cn.foliage)*tpha*resf, 
                                                       ifelse(Genus %in% c("Casuarina"), agb.foliage*ndfas$cas*(1/cn.foliage)*tpha*resf,
                                                              ifelse(Genus %in% c("Robinia"), agb.foliage*ndfas$rob*(1/cn.foliage)*tpha*resf,
                                                                     ifelse(Genus %in% c("Olneya"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha*resf,
                                                                            ifelse(Genus %in% c("Elaeagnus"), agb.foliage*ndfas$ela*(1/cn.foliage)*tpha*resf,
                                                                                   0))))))))))



woodc <- fol %>% 
  mutate(woodn = ifelse(Genus %in% c("Acacia"), agb.wood*ndfas$aca*(1/cn.wood)*tpha, 
                        ifelse(Genus %in% c("Albizia"), agb.wood*ndfas$alb*(1/cn.wood)*tpha,
                               ifelse(Genus %in% c("Alnus"), agb.wood*ndfas$ald*(1/cn.wood)*tpha,
                                      ifelse(Genus %in% c("Cercocarpus"), agb.wood*ndfas$other*(1/cn.wood)*tpha,
                                             ifelse(Genus %in% c("Prosopis"), agb.wood*ndfas$pro*(1/cn.wood)*tpha, 
                                                    ifelse(Genus %in% c("Casuarina"), agb.wood*ndfas$cas*(1/cn.wood)*tpha,
                                                           ifelse(Genus %in% c("Robinia"), agb.wood*ndfas$rob*(1/cn.wood)*tpha,
                                                                  ifelse(Genus %in% c("Olneya"), agb.wood*ndfas$other*(1/cn.wood)*tpha,
                                                                         ifelse(Genus %in% c("Elaeagnus"), agb.wood*ndfas$ela*(1/cn.wood)*tpha,
                                                                                0))))))))))



#belowground biomass (Li, 2003)
agb_bgb <- woodc %>%
  mutate(rb = (1.576*((agb2pred+agb)/2))^0.615) %>% #calc root biomass
  mutate(pfr = 0.072 + 0.354*exp(-0.06*rb)) %>%     #calc frac root biomass as FR
  mutate(fr = rb*pfr) %>%                           #calc FR biomass
  mutate(rootn = (1/cn.fr)*fr) %>% #does this rootn not include turnover (rootn = (1/cn.fr)*frt*fr)
  mutate(totaln = rootn*resr + foliagen + woodn)

agb_bgb <- agb_bgb[agb_bgb$totaln < 100,]

saveRDS(agb_bgb, "output/agb_bgb_perakis.RDS")
