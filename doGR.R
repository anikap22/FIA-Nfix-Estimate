
# doGR.R
# processes FIA data (growth rates)
# by Anika Petach 
# 6/16/17

require(tidyr)
require(dplyr)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

## Set parameters
res <- 0.4 #leaf resorption 0.4 for normal, 1 for upper bound, res is amount lost here
resr <- 0.73 #root resporption 0.73 for normal, 1 for upper bound, res is amount lost here
res <- 1
resr <- 1
source("scripts/ndfa_bootstrapping.R") #go to this file, regular for normal, change to 1s for upper bound



# using RGR from same species ------------------------------------

# Load datafrom ndfa_GRnearby.R
agb_1_2 <- readRDS("output/agb_1_2.RDS") #biomass measured and predicted for t2, fixers only

agb_1_2 <- agb_1_2[agb_1_2$dbh2pred>=2,] #jenkins equations only work for trees dbh>2
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
cn.foliage <- 15  #Ferlian et al, 2017
cn.wood <- 465    #Johnson et al, 2014 (from sugar maple)
cn.fr <- 43       #Gordon & Jackson, 2000

# foliage (full crown of foliage)
fol <- plot_agb %>% 
  mutate(foliagen = ifelse(Genus %in% c("Acacia"), agb.foliage*ndfas$aca*(1/cn.foliage)*tpha, 
                           ifelse(Genus %in% c("Albizia"), agb.foliage*ndfas$alb*(1/cn.foliage)*tpha,
                                  ifelse(Genus %in% c("Alnus"), agb.foliage*ndfas$ald*(1/cn.foliage)*tpha,
                                         ifelse(Genus %in% c("Cercocarpus"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha,
                                                ifelse(Genus %in% c("Prosopis"), agb.foliage*ndfas$pro*(1/cn.foliage)*tpha, 
                                                       ifelse(Genus %in% c("Casuarina"), agb.foliage*ndfas$cas*(1/cn.foliage)*tpha,
                                                              ifelse(Genus %in% c("Robinia"), agb.foliage*ndfas$rob*(1/cn.foliage)*tpha,
                                                                     ifelse(Genus %in% c("Olneya"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha,
                                                                            ifelse(Genus %in% c("Elaeagnus"), agb.foliage*ndfas$ela*(1/cn.foliage)*tpha,
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
  mutate(rb = (1.576*((agb2pred+agb)/2))^0.615) %>%
  mutate(pfr = 0.072 + 0.354*exp(-0.06*rb)) %>%
  mutate(fr = rb*pfr) %>%
  mutate(rootn = (1/cn.fr)*fr) %>% #does this rootn not include turnover (rootn = (1/cn.fr)*frt*fr)
  mutate(totaln = rootn*resr + foliagen*res + woodn)

agb_bgb <- agb_bgb[agb_bgb$totaln < 100,]

#saveRDS(agb_bgb, "output/agb_bgb.RDS")
saveRDS(agb_bgb, "output/agb_bgb_upper.RDS")




# Ndfa state scaling (no neighbors) OLD --------------------------------
############ Step1: Load cleaned data witout neighbors
gdata <- readRDS("output/gdata.RDS") #from clean.R

grdata <- readRDS("output/gdataGR1.RDS") #grdata has fixers and non-fixers from plots with a fixer present

#recode -999.00 as na in tpha
grdata[grdata$tpha1 == -999,]$tpha1 <- NA
grdata[grdata$tpha2 == -999,]$tpha2 <- NA

#remove islands (I think this isn't necessary, PR is still in the factor list but not present in data)
grdata <- grdata[grdata$state != "PR",]

#if dbh2=1.00 then set it to dbh1 value
grdata[grdata$dbh2==1.00,]$dbh2 <- grdata[grdata$dbh2==1.00,]$dbh1

#state scaling
gplot <- gdata[!duplicated(gdata$pcn),]       #get unique plots in gdata
grplot <- grdata[!duplicated(grdata$pcn),]    #get unique plots in grdata
grstate <- count(grplot, state)
gstate <- count(gplot, state)
scale <- merge(gstate, grstate, by = "state")
colnames(scale) <- c("state","gstate","grstate")
scale <- mutate(scale, scaling = scale$gstate/scale$grstate)

# Get biomass allocations
#hardwood parameters (Jenkins, 2003)
B0.foliage <- -4.0813
B1.foliage <- 5.8816
B0.roots <- -1.6911
B1.roots <- 0.8160
B0.bark <- -2.0129
B1.bark <- -1.6805
B0.wood <- -0.3065
B1.wood <- -5.4240

# agb fraction from Jenssen
b.foliage <- exp(B0.foliage+(B1.foliage/grdata$dbh2))
b.bark <- exp(B0.bark+(B1.bark/grdata$dbh2))
b.wood <- exp(B0.wood+(B1.wood/grdata$dbh2))
b.root <- exp(B0.roots+(B1.roots/grdata$dbh2))
b.other <- 1 - (b.foliage+b.bark+b.wood)

grdata[grdata$agb2 < 0.1,]$agb2 <- grdata[grdata$agb2 < 0.1,]$agb1
grdata$agbd <- grdata$agb2 - grdata$agb1
grdata[grdata$agbd < 0,]$agbd <- 0.05

#get change in biomass (kgC in foliage and wood)
plot_agb <- grdata %>% 
  mutate(agb.foliage = exp(B0.foliage+(B1.foliage/dbh2))*agb2,
         agb.wood = (agbd)/t,
         bgb.root = exp(B0.roots+(B1.roots/grdata$dbh2))*agb2)
# should agb.wood = ((agbd)/t)*(b.wood+b.bark)

########### Step: Get % N by tissue and species
# Acacia koa (from Pearson & Vitousek, 2001), %N as frac
# n.branch <- 0.0085
# n.twig <- 0.017
# n.dead <- 0.0104
# n.stem <- 0.0053
# n.foliage <- 0.0289
# n.bark <- n.dead
# n.other <- mean(n.branch, n.twig)

#use C:N ratios since biomass is in KgC
cn.foliage <- 15  #Ferlian et al, 2017
cn.wood <- 465    #Johnson et al, 2014 (from sugar maple)
cn.fr <- 43       #Gordon & Jackson, 2000

########### Step: Calculate N used
# foliage (full crown of foliage)
fol <- plot_agb %>% 
  mutate(foliage = ifelse(Genus %in% c("Acacia"), agb.foliage*ndfas$aca*(1/cn.foliage)*tpha1, 
                          ifelse(Genus %in% c("Albizia"), agb.foliage*ndfas$alb*(1/cn.foliage)*tpha1,
                                 ifelse(Genus %in% c("Alnus"), agb.foliage*ndfas$ald*(1/cn.foliage)*tpha1,
                                        ifelse(Genus %in% c("Cercocarpus"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha1,
                                               ifelse(Genus %in% c("Prosopis"), agb.foliage*ndfas$pro*(1/cn.foliage)*tpha1, 
                                                      ifelse(Genus %in% c("Casuarina"), agb.foliage*ndfas$cas*(1/cn.foliage)*tpha1,
                                                             ifelse(Genus %in% c("Ebenopsis"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha1,
                                                                    ifelse(Genus %in% c("Sophora"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha1,
                                                                           ifelse(Genus %in% c("Robinia"), agb.foliage*ndfas$rob*(1/cn.foliage)*tpha1,
                                                                                  ifelse(Genus %in% c("Olneya"), agb.foliage*ndfas$other*(1/cn.foliage)*tpha1,
                                                                                         ifelse(Genus %in% c("Elaeagnus"), agb.foliage*ndfas$ela*(1/cn.foliage)*tpha1,
                                                                                                0))))))))))))



woodc <- fol %>% 
  mutate(wood = ifelse(Genus %in% c("Acacia"), agb.wood/t*ndfas$aca*(1/cn.wood)*tpha1, 
                       ifelse(Genus %in% c("Albizia"), agb.wood/t*ndfas$alb*(1/cn.wood)*tpha1,
                              ifelse(Genus %in% c("Alnus"), agb.wood/t*ndfas$ald*(1/cn.wood)*tpha1,
                                     ifelse(Genus %in% c("Cercocarpus"), agb.wood/t*ndfas$other*(1/cn.wood)*tpha1,
                                            ifelse(Genus %in% c("Prosopis"), agb.wood/t*ndfas$pro*(1/cn.wood)*tpha1, 
                                                   ifelse(Genus %in% c("Casuarina"), agb.wood/t*ndfas$cas*(1/cn.wood)*tpha1,
                                                          ifelse(Genus %in% c("Ebenopsis"), agb.wood/t*ndfas$other*(1/cn.wood)*tpha1,
                                                                 ifelse(Genus %in% c("Sophora"), agb.wood/t*ndfas$other*(1/cn.wood)*tpha1,
                                                                        ifelse(Genus %in% c("Robinia"), agb.wood/t*ndfas$rob*(1/cn.wood)*tpha1,
                                                                               ifelse(Genus %in% c("Olneya"), agb.wood/t*ndfas$other*(1/cn.wood)*tpha1,
                                                                                      ifelse(Genus %in% c("Elaeagnus"), agb.wood/t*ndfas$ela*(1/cn.wood)*tpha1,
                                                                                             0))))))))))))



biom <- woodc %>%
  group_by(state) %>%
  summarize(foli = sum(foliage)*res, woo = sum(wood)) %>%
  mutate(total = foli+woo) #total gives kgN/year

biomass <- merge(biom,scale,by="state")
kgN <- sum((biomass$total*biomass$scaling),na.rm=T)*2428.114
kgN #total N fixed with %ndfa method without roots

#belowground biomass (Li, 2003)
biot <- woodc %>%
  mutate(rb = (1.576*((agb2+agb1)/2))^0.615) %>%
  mutate(pfr = 0.072 + 0.354*exp(-0.06*rb)) %>%
  mutate(fr = rb*pfr) %>%
  mutate(frt = 0.641*fr)

# biofr <- biot %>%
#   mutate(rootn = (1/cn.fr)*frt*fr)

biofr <- biot %>%
  mutate(rootn = ifelse(FIX==1, (1/cn.fr)*fr,
                        0))

# biofr <- biot %>%
#   mutate(rootn = (1/cn.fr)*fr) #does not include turnover rate

#incorporate resorption and remove nonfixers
biofr$root <- biofr$rootn*resr
biofr$total <- biofr$wood+biofr$foliage+biofr$root
biofr$foliage <- biofr$foliage*res

biofr[!is.na(biofr$pcn),]

#saveRDS(biofr, "output/GR_ndfa.RDS")

#remove nonfixers
biomr <- biofr %>%
  group_by(state) %>%
  summarize(foli = sum(foliage,na.rm=T), woo = sum(wood,na.rm=T), roo = sum(rootn,na.rm=T)) %>%
  mutate(total = foli+woo+roo) #total gives kgN/year

biomp <- biofr %>%
  group_by(pcn) %>%
  summarize(tot=sum(total), tottpha=sum(total*tpha2))

(sum((biomp$tot),na.rm=T))/length(unique(grdata$pcn)) #kgN/ha/yr with %ndfa method

biomassr <- merge(biomr,scale,by="state")
kgN <- sum((biomassr$total*biomassr$scaling),na.rm=T)*2428.114
kgN

test <- biofr %>%
  group_by(pcn) %>%
  summarize(pfix=sum(fix.r,na.rm=T))

# P deltoides & p. trichocarpa ------------------------------------------
# Get all poplars to make growth curve
# No poplars in time series data
grdata <- readRDS("output/gdataGR1.RDS") # from time series data
allstems <- readRDS("output/allstems.RDS")
spdata <- read.csv("ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")

allstems_pop <- dplyr::filter(allstems, spcd %in% c(742, 745))
rm(allstems)

allstems.allo <- merge(allstems_pop, 
                       spdata[,c("SPCD","JENKINS_TOTAL_B1","JENKINS_TOTAL_B2")], 
                       by.x="spcd", by.y="SPCD", 
                       all.x=T)

populus <- grdata[grdata$spcd %in% c(752,749,747,741,745,742,748,743,744,753,740,746), ]

spdata <- dplyr::select(spdata, SPCD, GENUS, STOCKING_SPGRPCD)

hardwood <- merge(grdata, spdata, by.x="spcd1", by.y="SPCD", all.x=T, all.y=F)
hardwood <- dplyr::filter(hardwood, STOCKING_SPGRPCD==37)

model.hardwood <- loess(hardwood$LGRp ~ hardwood$dbh1,
                        control = loess.control(surface="direct"))

s0 <- data.table(allstems.allo) #make data table
s0[, predfit := predict(model.hardwood, dbh, se=T)$fit] #models[[i]]
s0[, predse := predict(model.hardwood, dbh, se=T)$se.fit] #models[[i]]
for(k in 1:nrow(s0)){
    s0$RGRpred[k] <- rnorm(1, mean=s0$predfit[k], sd=s0$predse[k])
}
s0[, dbh2pred:=exp(RGRpred+log(dbh))]
s0[, agb2pred:=exp(JENKINS_TOTAL_B1+JENKINS_TOTAL_B2*log(dbh2pred))]


agb_1_2 <- s0 

agb_1_2 <- agb_1_2[agb_1_2$dbh2pred>=2,] #jenkins equations only work for trees dbh>2
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
cn.foliage <- 15  #Ferlian et al, 2017
cn.wood <- 465    #Johnson et al, 2014 (from sugar maple)
cn.fr <- 43       #Gordon & Jackson, 2000

# foliage (full crown of foliage)
fol <- plot_agb %>% 
  mutate(foliagen = agb.foliage*.65*(1/cn.foliage)*tpha)

woodc <- fol %>% 
  mutate(woodn = agb.wood*0.65*(1/cn.wood)*tpha)

#belowground biomass (Li, 2003)
agb_bgb <- woodc %>%
  mutate(rb = (1.576*((agb2pred+agb)/2))^0.615) %>%
  mutate(pfr = 0.072 + 0.354*exp(-0.06*rb)) %>%
  mutate(fr = rb*pfr) %>%
  mutate(rootn = (1/cn.fr)*fr) %>% #does this rootn not include turnover (rootn = (1/cn.fr)*frt*fr)
  mutate(totaln = rootn*resr + foliagen*res + woodn)

agb_bgb <- agb_bgb[agb_bgb$totaln < 100,]

# get agb_bgb
# then line 27 in map.R gfdata <- agb_bgb
# then line 1052 in map.R

sum(agb_bgb$totaln, na.rm=T)

