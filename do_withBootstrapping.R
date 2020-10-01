
# do_withBootstrapping.R
# processes FIA data (accretion data) with bootstrapped F rates
# by Anika Petach 
# 5/24/18

require(tidyr)
require(dplyr)
library(gridExtra)
require(ggplot2)
library(lattice)
library(ggmap)
library(maps)
library(mapdata)
require(RColorBrewer)
library(mapproj)
require(MASS)
require(raster)
require(bit64)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

############ Step1: Load cleaned data

gdata <- readRDS("output/gdata.RDS") #trees only
sdata_pl <- readRDS("output/saplingdata.RDS") #saplings only
speciesdata <- read.csv("ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")

############ Step2: Bootstrapping fixation rates
source("scripts/genus_regressions_new.R")

#plot distribution
plot(density(ac.samp),main=expression(italic(Acacia)),xlab=expression('fixation rate (kg N ha'^-1*' yr'^-1*') per stand basal area'))
plot(density(al.samp),main=expression(italic(Alnus)),xlab=expression('fixation rate (kg N ha'^-1*' yr'^-1*') per stand basal area'))
plot(density(pr.samp),main=expression(italic(Prosopis)),xlab=expression('fixation rate (kg N ha'^-1*' yr'^-1*') per stand basal area'))
plot(density(ro.samp),main=expression(italic(Robinia)),xlab=expression('fixation rate (kg N ha'^-1*' yr'^-1*') per stand basal area'))
plot(density(ot.samp),main=expression(italic(Others)),xlab=expression('fixation rate (kg N ha'^-1*' yr'^-1*') per stand basal area'))

table(gdata$Genus)

nfixers <- gdata %>% 
  filter(FIX==1) %>% 
  summarise(n())
gdata %>% 
  filter(FIX == 1) %>% 
  group_by(Genus) %>% 
  summarise(cnt=n()) %>% 
  mutate(pctfixerstems = cnt/nfixers[[1]]*100)

# combine trees and saplings -----------------------------------
#select cols from gdata for merging with saplings
mature <- gdata %>% dplyr::select(pcn,fixer_present,meastime,lat,lon,elev,stdage,
                                  state,xmh,spcd,tcn,dbh,cclcd,agb,tpha,Genus,
                                  COMMON_NAME,FIX,ACT1vsRIZ0,BAm2)
#add BAm2ha, TREECOUNT to mature
mature$BAm2ha <- mature$BAm2*mature$tpha
mature$TREECOUNT <- 1

# select cols from sdata for merging with mature trees
sap <- sdata_pl %>% dplyr::select(pcn,fixer_present,meastime,lat,lon,elev,stdage,
                                  state,xmh,SPCD,CN,tpha,Genus,COMMON_NAME,TREECOUNT,
                                  FIX,ACT1vsRIZ0,BAm2ha)
#add cclcd, agb, BAm2 to sap
sap$cclcd <- -1
dbhcm <- sqrt(2.54*12.7)
sap$dbh <- dbhcm
speciesdata.small <- speciesdata %>% dplyr::select(SPCD,JENKINS_TOTAL_B1,JENKINS_TOTAL_B2)
sap0 <- merge(sap, speciesdata.small, by="SPCD", all.x=TRUE, all.y=FALSE)
sap <- sap0 %>%
  mutate(agb = exp(JENKINS_TOTAL_B1 + JENKINS_TOTAL_B2*log(dbhcm)))
sap <- sap %>% dplyr::select(-JENKINS_TOTAL_B1, -JENKINS_TOTAL_B2)
sap$BAm2 <- sap$BAm2ha/sap$tpha

#rename CN to tcn, SPCD to spcd in sap
names(sap)[names(sap) == 'CN'] <- 'tcn'
names(sap)[names(sap) == 'SPCD'] <- 'spcd'

# combine mature and saplings
allstems <- rbind(mature, sap)

#saveRDS(allstems, "output/allstems.RDS")

# using mvrnorm -----------------------
# Group by genus and pcn
source("scripts/genus_regressions_new.R")
allstems <- readRDS("output/allstems.RDS") #read in if don't have it from prev step

stems0 <- allstems %>%
  group_by(pcn, Genus) %>%
  summarise(pBA = sum(BAm2ha),
            FIX = first(FIX),
            ACT1vsRIZ0 = first(ACT1vsRIZ0))

#while loop to get CI around estimate
totalt <- vector()
i <- 1
while(i<2){
  n <- 10000
  plot_sens <- stems0 %>% 
    group_by(pcn, Genus) %>%
    mutate(fix = ifelse(Genus %in% c("Acacia"), sum(ac.samp[sample(1:n, 1),]*c(1,pBA)), 
                        ifelse(Genus %in% c("Albizia"), sum(ot.samp[sample(1:n, 1),]*c(1,pBA)),
                               ifelse(Genus %in% c("Alnus"), al.samp[sample(1:n, 1)],
                                      ifelse(Genus %in% c("Cercocarpus"), sum(pr.samp[sample(1:n, 1),]*c(1,pBA)),
                                             ifelse(Genus %in% c("Prosopis"), sum(pr.samp[sample(1:n, 1),]*c(1,pBA)), 
                                                    ifelse(Genus %in% c("Ebenopsis"), sum(ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                           ifelse(Genus %in% c("Robinia"),ro.samp[sample(1:n, 1)],
                                                                  ifelse(Genus %in% c("Olneya"), sum(ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                                         ifelse(Genus %in% c("Elaeagnus"), sum(ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                                                0)))))))))) 

  #saveRDS(plot_sens, "output/stemsfix.RDS")
  
  
  pfix <- plot_sens %>%
    group_by(pcn) %>%
    summarize(pfix = sum(fix, na.rm=T))
  
  totalt[i] <- mean(pfix$pfix)
  i<-i+1
  print(i)
}

t.test(totalt)

numplots <- length(unique(stems0$pcn))


mean(pfix$pfix) #mean F rate in plots with fixers
sum(pfix$pfix)/length(unique(stems0$pcn)) #mean F rate in all plots

BA_summarystats <- plot_sens %>%
  group_by(Genus) %>%
  summarize(BAmean = mean(pBA, na.rm=T),
            BAmin = min(pBA, na.rm=T),
            BAmax = max(pBA, na.rm=T),
            BAmed = median(pBA, na.rm=T))

fix_summarystats <- plot_sens %>%
  group_by(Genus) %>%
  summarize(FIXmean = mean(fix, na.rm=T),
            FIXmin = min(fix, na.rm=T),
            FIXmax = max(fix, na.rm=T),
            FIXmed = median(fix, na.rm=T))

mean(plot_sens[plot_sens$Genus=="Alnus",]$pBA) #avg plot BA of Alnus

#Break down BA by genus
plot_sens %>%
  group_by(Genus) %>%
  summarize(meanBA=mean(pBA,na.rm=T), maxBA=max(pBA,na.rm=T), minBA=min(pBA,na.rm=T))

# Bootstrap mvrnorm -----------------------------------
#while loop to get CI around estimate
allstems_bs <- stems0 %>% 
  group_by(pcn, Genus) %>%
  mutate(fix = ifelse(Genus %in% c("Acacia"), sum(ac.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                      ifelse(Genus %in% c("Albizia"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                             ifelse(Genus %in% c("Alnus"), sum(al.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                    ifelse(Genus %in% c("Cercocarpus"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                           ifelse(Genus %in% c("Prosopis"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                                                  ifelse(Genus %in% c("Ebenopsis"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                         ifelse(Genus %in% c("Robinia"),rnorm(1,mean(robinia$rate),var(robinia$rate)^0.5),
                                                                ifelse(Genus %in% c("Olneya"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                       ifelse(Genus %in% c("Elaeagnus"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                              0)))))))))) 
i <- 1
while(i<999){
  plot_sens <- stems0 %>% 
    group_by(pcn, Genus) %>%
    mutate(fix = ifelse(Genus %in% c("Acacia"), sum(ac.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                        ifelse(Genus %in% c("Albizia"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                               ifelse(Genus %in% c("Alnus"), sum(al.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                      ifelse(Genus %in% c("Cercocarpus"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                             ifelse(Genus %in% c("Prosopis"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                                                    ifelse(Genus %in% c("Ebenopsis"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                           ifelse(Genus %in% c("Robinia"),rnorm(1,mean(robinia$rate),var(robinia$rate)^0.5),
                                                                  ifelse(Genus %in% c("Olneya"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                         ifelse(Genus %in% c("Elaeagnus"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                                0)))))))))) 
  


  allstems_bs <- bind_cols(allstems_bs,plot_sens[,"fix"])
  i<-i+1
  print(i)
}


# Bootstrap with mvrnorm using function from functions.R -----------------------
fixcut <- 0.9 # 1 for confirmed fixers
latcut <- 35
dlat <- 1 # grid cell size
dlon <- 1 # grid cell size

# load data
FOR_COV_FOR_ALIGN <- raster("output/FOR_COV_FOR_ALIGN.grd") #forest area in km2 per grid cell, used compareRaster to check
FOR_COV_ALL_ALIGN <- raster("output/FOR_COV_ALL_ALIGN.grd") #ground area in km2 per grid cell, used compareRaster to check
FOR_COV_ALIGN <- raster("output/FOR_COV_ALIGN.grd")
p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv", header=T)
allstems <- readRDS("output/allstems.RDS") #for bootstrapping

# prep lat long list
p <- p[p$state != 'PR' & p$state != 'VI' & p$state != 'HI', ]
p[p$lat == -999, ]$lat <- NA
p[p$lon == -999, ]$lon <- NA
latmin <- min(p$lat, na.rm = TRUE)
latmax <- max(p$lat, na.rm = TRUE)
lonmin <- min(p$lon, na.rm = TRUE)
lonmax <- max(p$lon, na.rm = TRUE)
latminint <- floor(latmin)
latmaxint <- ceiling(latmax)
lonminint <- floor(lonmin)
lonmaxint <- ceiling(lonmax)
lat.list <- seq(latminint, latmaxint, dlat)
lon.list <- seq(lonminint, lonmaxint, dlon)
nlat <- length(lat.list)
nlon <- length(lon.list)

source("scripts/functions.R")

bootlength <- 10
bootstrapped <- matrix(nrow=bootlength, ncol=10)
bootstrapped <- as.data.frame(bootstrapped)
colnames(bootstrapped) <- c("tgn", "acapct", "albpct", "alnpct", "cerpct",
                            "elapct", "olnpct", "propct", "robpct")
for(i in 1:bootlength){
  bootstrapped[i,] <- plotfix_new(allstems, 6, normalrates) #canpos=6 for regular
  print(i)
}




# Light downregulation, no boot ----------------------------------

stems0l <- allstems %>%
  dplyr::filter(cclcd < 4) %>%
  group_by(pcn, Genus) %>%
  summarise(pBA = sum(BAm2ha),
            FIX = first(FIX),
            ACT1vsRIZ0 = first(ACT1vsRIZ0))

allstems_lightlim <- stems0l %>%
  group_by(pcn, Genus) %>%
  mutate(fix = ifelse(Genus %in% c("Acacia"), sum(ac.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                      ifelse(Genus %in% c("Albizia"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                             ifelse(Genus %in% c("Alnus"), sum(al.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                    ifelse(Genus %in% c("Cercocarpus"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                           ifelse(Genus %in% c("Prosopis"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                                                  ifelse(Genus %in% c("Ebenopsis"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                         ifelse(Genus %in% c("Robinia"),rnorm(1,mean(robinia$rate),var(robinia$rate)^0.5),
                                                                ifelse(Genus %in% c("Olneya"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                       ifelse(Genus %in% c("Elaeagnus"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                              0)))))))))) 
saveRDS(allstems_lightlim, "output/stemsfix_lightlim.RDS")

# Light limitation (bootstrap) -----------------------------------
fixcut <- 0.9 # 1 for confirmed fixers
latcut <- 35
dlat <- 1 # grid cell size
dlon <- 1 # grid cell size

# load data
FOR_COV_FOR_ALIGN <- raster("output/FOR_COV_FOR_ALIGN.grd") #forest area in km2 per grid cell, used compareRaster to check
FOR_COV_ALL_ALIGN <- raster("output/FOR_COV_ALL_ALIGN.grd") #ground area in km2 per grid cell, used compareRaster to check
FOR_COV_ALIGN <- raster("output/FOR_COV_ALIGN.grd")
p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv", header=T)
allstems <- readRDS("output/allstems.RDS") #for bootstrapping

# prep lat long list
p <- p[p$state != 'PR' & p$state != 'VI' & p$state != 'HI', ]
p[p$lat == -999, ]$lat <- NA
p[p$lon == -999, ]$lon <- NA
latmin <- min(p$lat, na.rm = TRUE)
latmax <- max(p$lat, na.rm = TRUE)
lonmin <- min(p$lon, na.rm = TRUE)
lonmax <- max(p$lon, na.rm = TRUE)
latminint <- floor(latmin)
latmaxint <- ceiling(latmax)
lonminint <- floor(lonmin)
lonmaxint <- ceiling(lonmax)
lat.list <- seq(latminint, latmaxint, dlat)
lon.list <- seq(lonminint, lonmaxint, dlon)
nlat <- length(lat.list)
nlon <- length(lon.list)

# # data for fixation rate regression
source("scripts/genus_regressions_new.R")

# regression, get mean and SE
# ac <- lm(rate ~ BA, data = acacia)
# al <- lm(rate ~ BA, data = alnus)
# pr <- lm(rate ~ BA, data = prosopis)
# ot <- lm(rate ~ BA, data = other)
n = 10000

# get forest area
foresta.grid <- ground.grid <- fracfor.grid <- array(dim = c(nlon, nlat))
for(i in 1:nlon){
  for(j in 1:nlat){
    coords <- data.frame(lon = lon.list[i], 
                         lat = lat.list[j])
    foresta.grid[i,j] <- raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
    ground.grid[i,j] <- raster::extract(FOR_COV_ALL_ALIGN, coords)*100 #
    fracfor.grid[i,j] <- raster::extract(FOR_COV_ALIGN, coords)
  }
}

#for getrates function, plotfix_new function
source("scripts/functions.R")

#genus should be in lowercase inside quotes
#type can be max or min, inside quotes
#for all genera at normal values use genus="standard" and leave type blank
# getrates <- function(genus, type){
#   n = 10000
#   #standard rates
#   ac.samp <- mvrnorm(n, coef(ac)[1:2], vcov(ac))
#   pr.samp <- mvrnorm(n, coef(pr)[1:2], vcov(pr))
#   al.samp <- mvrnorm(n, coef(al)[1:2], vcov(al))
#   ot.samp <- mvrnorm(n, coef(ot)[1:2], vcov(ot))
#   ro.samp <- rnorm(n, mean(robinia$rate), var(robinia$rate)^0.5) #gets random distribution of constant robinia rate
#   
#   #modified rates
#   if(genus=="acacia"){
#     if(type=="max"){
#       ac.samp <- mvrnorm(n, coef(ac)[1:2]*c(1,1.5), vcov(ac)) #slope 150%, int same
#     }else if(type=="min"){
#       ac.samp <- mvrnorm(n, coef(ac)[1:2]*c(1,0.5), vcov(ac)) #slope 50%, int same
#     }
#   }else if(genus=="albizia" | genus=="ebenopsis" | genus=="olneya" | genus=="elaeagnus" | genus=="other"){
#     if(type=="max"){
#       ot.samp <- mvrnorm(n, coef(ot)[1:2]*c(1,1.5), vcov(ot))
#     }else if(type=="min"){
#       ot.samp <- mvrnorm(n, coef(ot)[1:2]*c(1,0.5), vcov(ot))
#     }
#   }else if(genus=="alnus"){
#     if(type=="max"){
#       al.samp <- mvrnorm(n, coef(al)[1:2]*c(1,1.5), vcov(al))
#     }else if(type=="min"){
#       al.samp <- mvrnorm(n, coef(al)[1:2]*c(1,0.5), vcov(al))
#     }
#   }else if(genus=="cercocarpus" | genus=="prosopis"){
#     if(type=="max"){
#       pr.samp <- mvrnorm(n, coef(pr)[1:2]*c(1,1.5), vcov(pr))
#     }else if(type=="min"){
#       pr.samp <- mvrnorm(n, coef(pr)[1:2]*c(1,0.5), vcov(pr))
#     }
#   }else if(genus=="robinia"){
#     if(type=="max"){
#       ro.samp <- rnorm(n, mean(robinia$rate)*1.5, var(robinia$rate)^0.5)
#     }else if(type=="min"){
#       ro.samp <- rnorm(n, mean(robinia$rate)*0.5, var(robinia$rate)^0.5)
#     }
#   }else if(genus=="standard"){
#     ac.samp <- mvrnorm(n, coef(ac)[1:2], vcov(ac))
#     pr.samp <- mvrnorm(n, coef(pr)[1:2], vcov(pr))
#     al.samp <- mvrnorm(n, coef(al)[1:2], vcov(al))
#     ot.samp <- mvrnorm(n, coef(ot)[1:2], vcov(ot))
#     ro.samp <- rnorm(n, mean(robinia$rate), var(robinia$rate)^0.5)
#   }
#   r <- list(ac.samp = ac.samp,
#             pr.samp = pr.samp,
#             al.samp = al.samp,
#             ot.samp = ot.samp,
#             ro.samp = ro.samp)
#   return(r)
# }

# plotdata should be fed allstems, 
# r should be fed output from getrates function, 
# canpos = 6 for all trees or 4 for canopy only
# plotfix <- function(plotdata, canpos, r) {
#   
#   allstemsf <- plotdata[plotdata$cclcd < canpos, ] 
# 
#   stems00 <- allstemsf %>%
#     group_by(pcn, Genus) %>%
#     summarise(pBA = sum(BAm2ha),
#               FIX = first(FIX),
#               ACT1vsRIZ0 = first(ACT1vsRIZ0))
# 
#   n <- 10000
#   # get fixation rates for plots
#   plot_sens <- stems00 %>% 
#     group_by(pcn, Genus) %>%
#     mutate(fix = ifelse(Genus %in% c("Acacia"), sum(r$ac.samp[sample(1:n, 1),]*c(1,pBA)), 
#                         ifelse(Genus %in% c("Albizia"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
#                                ifelse(Genus %in% c("Alnus"), sum(r$al.samp[sample(1:n, 1),]*c(1,pBA)),
#                                       ifelse(Genus %in% c("Cercocarpus"), sum(r$pr.samp[sample(1:n, 1),]*c(1,pBA)),
#                                              ifelse(Genus %in% c("Prosopis"), sum(r$pr.samp[sample(1:n, 1),]*c(1,pBA)), 
#                                                     ifelse(Genus %in% c("Ebenopsis"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
#                                                            ifelse(Genus %in% c("Robinia"),r$ro.samp[sample(1:n, 1)],
#                                                                   ifelse(Genus %in% c("Olneya"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
#                                                                          ifelse(Genus %in% c("Elaeagnus"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
#                                                                                 0)))))))))) 
# 
#   
#   
#   an <- as.data.frame.table(by(plot_sens$fix, plot_sens$pcn, sum)) 
#   an <- data.frame(pcn = as.numeric(as.character(an[ ,1])), 
#                    fixr = as.numeric(an[ ,2]))
# 
#   p <- merge(p, an, 
#              by = "pcn", 
#              all.x = TRUE, all.y = FALSE)
#   
#   bygenus <- plot_sens %>%
#     group_by(Genus) %>%
#     #summarise_all(funs(sum(., na.rm=T)))
#     summarise(fix = sum(fix, na.rm = T))
#   bygenus <- bygenus[c(1:5,7,8,10:11),] #remove Ebenopsis (too few stems), Piscidia (1 stem), Sophora (1 stem), NAs
#   
#   
#   fix.grid <- fixfor.grid <- array(dim = c(nlon, nlat))
#   
#   print("starting iteration over grid")
#   for(i in 1:nlon){
#     lon.a <- lon.list[i] - dlon/2
#     lon.b <- lon.list[i] + dlon/2
#     
#     for(j in 1:nlat){
#       lat.a <- lat.list[j] - dlat/2
#       lat.b <- lat.list[j] + dlat/2
#       inds <- which(p$lon>=lon.a & p$lon<lon.b & p$lat>=lat.a & p$lat<lat.b)
#       
#       if(sum(inds) > 0){
#         fixfor.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)
#         fix.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)*foresta.grid[i,j] #fix rate times forest
#       }else{
#         fixfor.grid[i,j] <- NA
#         fix.grid[i,j] <- NA
#       }
#       
#     }
#   }
#   
#   fixedn <- sum(fix.grid, na.rm = T) #total N fixed across US from accretion method
#   totalk <- sum(bygenus[ ,"fix"])
#   acan <- (bygenus[[1,"fix"]]/totalk)*100 #*fixedn to get kgN instead of %
#   albn <- (bygenus[[2,"fix"]]/totalk)*100
#   alnn <- (bygenus[[3,"fix"]]/totalk)*100
#   casn <- (bygenus[[4,"fix"]]/totalk)*100
#   cern <- (bygenus[[5,"fix"]]/totalk)*100
#   elan <- (bygenus[[6,"fix"]]/totalk)*100
#   olnn <- (bygenus[[7,"fix"]]/totalk)*100
#   pron <- (bygenus[[8,"fix"]]/totalk)*100
#   robn <- (bygenus[[9,"fix"]]/totalk)*100
#   
#   totalacc <- sum(fix.grid, na.rm = T)/1e9 #total N fixed across US from accretion method 
# 
#   output <- list(tgn=totalacc, acapct=acan, albpct=albn, alnpct=alnn, caspct=casn, 
#                  cerpct=cern, ealpct=elan, olnpct=olnn, propct=pron, robpct=robn)
#   return(output)
# }

normalrates <- getrates("standard")
orig_fixedn <- plotfix_new(allstems, 6, normalrates)
mean(orig_fixedn$tgn, na.rm=T)*1e9/sum(foresta.grid, na.rm=T) #kgN/ha forest/yr for accretion

bootlength <- 10
bootstrapped <- matrix(nrow=bootlength, ncol=11)
bootstrapped <- as.data.frame(bootstrapped)
colnames(bootstrapped) <- c("tgn", "acapct", "albpct", "alnpct", "cerpct",
                            "elapct", "olnpct", "propct", "robpct", "fixperfor", 
                            "fixpergr")
for(i in 1:bootlength){
  bootstrapped[i,] <- plotfix_new(allstems, 6, normalrates) #canopos=4 for light lim
  print(i)
}

saveRDS(bootstrapped, "output/bootstrapped_acc.RDS")
mean(bootstrapped$tgn, na.rm=T)*1e9/sum(foresta.grid, na.rm=T) #kgN/ha forest/yr for light lim
mean(bootstrapped$tgn, na.rm=T)*1e9/sum(ground.grid, na.rm=T) #kgN/ha ground/yr for light lim


bootlength <- 10
bootstrapped <- matrix(nrow=bootlength, ncol=11)
bootstrapped <- as.data.frame(bootstrapped)
colnames(bootstrapped) <- c("tgn", "acapct", "albpct", "alnpct", "cerpct",
                            "elapct", "olnpct", "propct", "robpct", "fixperfor", 
                            "fixpergr")
for(i in 1:bootlength){
  bootstrapped[i,] <- plotfix_new(allstems, 4, normalrates) #canopos=4 for light lim
  print(i)
}

saveRDS(bootstrapped, "output/bootstrapped_lightlim.RDS")
mean(bootstrapped$tgn, na.rm=T)*1e9/sum(foresta.grid, na.rm=T) #kgN/ha forest/yr for light lim
mean(bootstrapped$tgn, na.rm=T)*1e9/sum(ground.grid, na.rm=T) #kgN/ha ground/yr for light lim


# %Ndfa bootstrap by genus ------------------------------------------------------
agb_1_2 <- readRDS("output/agb_1_2.RDS") #biomass measured and predicted for t2, fixers only
FOR_COV_FOR_ALIGN <- raster("output/FOR_COV_FOR_ALIGN.grd") #forest area in km2 per grid cell, used compareRaster to check
FOR_COV_ALL_ALIGN <- raster("output/FOR_COV_ALL_ALIGN.grd") #ground area in km2 per grid cell, used compareRaster to check
p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv", header=T)

fixcut <- 0.9 # 1 for confirmed fixers
latcut <- 35
dlat <- 1 # grid cell size
dlon <- 1 # grid cell size

p[p$lat == -999, ]$lat <- NA
p[p$lon == -999, ]$lon <- NA
latmin <- min(p$lat, na.rm = TRUE)
latmax <- max(p$lat, na.rm = TRUE)
lonmin <- min(p$lon, na.rm = TRUE)
lonmax <- max(p$lon, na.rm = TRUE)
latminint <- floor(latmin)
latmaxint <- ceiling(latmax)
lonminint <- floor(lonmin)
lonmaxint <- ceiling(lonmax)
lat.list <- seq(latminint, latmaxint, dlat)
lon.list <- seq(lonminint, lonmaxint, dlon)
nlat <- length(lat.list)
nlon <- length(lon.list)

agb_1_2 <- agb_1_2[agb_1_2$dbh2pred>=2,] #jenkins equations only work for trees dbh>2
agb_1_2[agb_1_2$agb2pred < 0.5,]$agb2pred <- 0.05 #to prevent 0 in denominator
agb_1_2$agbd <- agb_1_2$agb2pred - agb_1_2$agb #get difference in agb
agb_1_2[agb_1_2$agbd < 0,]$agbd <- 0.05 #if no change assign small increment

# get forest area
foresta.grid <- ground.grid <- array(dim = c(nlon, nlat))
for(i in 1:nlon){
  for(j in 1:nlat){
    coords <- data.frame(lon = lon.list[i], 
                         lat = lat.list[j])
    foresta.grid[i,j] <- raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
    ground.grid[i,j] <- raster::extract(FOR_COV_ALL_ALIGN, coords)*100 
  }
}

# get Ndfa data
ndfa <- read.csv("data/ndfa_litreview.csv", header=TRUE)

# separate genera
aca <- ndfa[ndfa$Genus=="Acacia" & ndfa$Location != "pot",]
alb <- ndfa[ndfa$Genus=="Albizia" & ndfa$Location != "pot",]
ald <- ndfa[ndfa$Genus=="Alnus" & ndfa$Location != "pot",]
pro <- ndfa[ndfa$Genus=="Prosopis" & ndfa$Location != "pot",]
ela <- ndfa[ndfa$Genus=="Elaeagnus" & ndfa$Location != "pot",]
cas <- ndfa[ndfa$Genus=="Casuarina" & ndfa$Location != "pot",]
rob <- ndfa[ndfa$Genus=="Robinia" & ndfa$Location != "pot",]

# get mean and standard error of genera
aca.mean <- mean(aca$Value, na.rm = T)
aca.se <- sd(aca$Value, na.rm = T)/nrow(aca)
alb.mean <- mean(alb$Value, na.rm = T)
alb.se <- sd(alb$Value, na.rm = T)/nrow(alb)
ald.mean <- mean(ald$Value, na.rm = T)
ald.se <- sd(ald$Value, na.rm = T)/nrow(ald)
pro.mean <- mean(pro$Value, na.rm = T)
pro.se <- sd(pro$Value, na.rm = T)/nrow(pro)
ela.mean <- mean(ela$Value, na.rm = T)
ela.se <- sd(ela$Value, na.rm = T)/nrow(ela)
cas.mean <- mean(cas$Value, na.rm = T)
cas.se <- sd(cas$Value, na.rm = T)/nrow(cas)
rob.mean <- mean(rob$Value, na.rm = T)
rob.se <- sd(rob$Value, na.rm = T)/nrow(rob)
other.mean <- mean(ndfa$Value, na.rm = T)
other.se <- sd(ndfa$Value, na.rm = T)/nrow(ndfa)


#plotdata should be agb_1_2
plotfixndfa <- function(plotdata, 
                        res = 0.4, 
                        resr = 0.73, 
                        cn.foliage = 35.1, 
                        cn.wood = 350, 
                        cn.fr = 43, 
                        ndfatype="realistic"){
  res <- res #leaf resorption 0.4 for normal, 1 for upper bound, res is amount lost here
  resr <- resr #root resporption 0.73 for normal, 1 for upper bound, res is amount lost here
  
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
  
  #use C:N ratios since biomass is in KgC
  cn.foliage <- cn.foliage  #Ferlian et al, 2017
  cn.wood <- cn.wood    #Johnson et al, 2014 (from sugar maple)
  cn.fr <- cn.fr       #Gordon & Jackson, 2000

  #ac.samp[sample(1:n, 1, replace=T),]
  
  if(ndfatype=="realistic"){
    ndfas <- list(aca = rnorm(1, aca.mean, aca.se)/100, 
                  alb = rnorm(1, alb.mean, alb.se)/100, 
                  ald = rnorm(1, ald.mean, ald.se)/100, 
                  pro = rnorm(1, pro.mean, pro.se)/100,
                  ela = rnorm(1, ela.mean, ela.se)/100, 
                  cas = rnorm(1, cas.mean, cas.se)/100, 
                  rob = rnorm(1, rob.mean, rob.se)/100, 
                  other = rnorm(1, other.mean, other.se)/100)
  }else if(ndfatype=="upperbound"){
    ndfas <- list(aca = 1, 
                  alb = 1, 
                  ald = 1, 
                  pro = 1,
                  ela = 1, 
                  cas = 1, 
                  rob = 1, 
                  other = 1)
  }
  
  
  #ndfas <- list(aca=1, alb=1, ald=1, pro=1, ela=1, cas=1, rob=1, other=1)

  # agb fraction from Jenssen
  # b.foliage <- exp(B0.foliage+(B1.foliage/agb_1_2$dbh2pred))
  # b.bark <- exp(B0.bark+(B1.bark/agb_1_2$dbh2pred))
  # b.wood <- exp(B0.wood+(B1.wood/agb_1_2$dbh2pred))
  # b.root <- exp(B0.roots+(B1.roots/agb_1_2$dbh2pred))
  # b.other <- 1 - (b.foliage+b.bark+b.wood)
  # print("calc agb fractions")
  
  #get change in biomass (kgC in foliage and wood)
  plot_agb <- plotdata %>% 
    mutate(agb.foliage = exp(B0.foliage+(B1.foliage/dbh2pred))*agb2pred,
           agb.wood = (agbd),
           bgb.root = exp(B0.roots+(B1.roots/dbh2pred))*agb2pred)
  plot_agb <- plot_agb[plot_agb$tpha > 0.5,]
  plot_agb <- plot_agb[plot_agb$dbh < 200,] #weird 2 outliers driving calcs
  plot_agb <- plot_agb[!is.na(plot_agb$dbh2pred),]

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
    mutate(rb = (1.576*((agb2pred+agb)/2))^0.615) %>% #root biomass
    mutate(pfr = 0.072 + 0.354*exp(-0.06*rb)) %>%     #percent root biomass that's fine root
    mutate(fr = rb*pfr) %>%                           #fine root biomass
    mutate(rootn = (1/cn.fr)*fr) %>% #does this rootn not include turnover (rootn = (1/cn.fr)*frt*fr)
    mutate(totaln = rootn*resr + foliagen*res + woodn)

  agb_bgb <- agb_bgb[agb_bgb$totaln < 100,]
  
  gf <- as.data.frame.table(by(agb_bgb$totaln, agb_bgb$pcn, sum))
  gf <- data.frame(pcn = as.numeric(as.character(gf[ ,1])), 
                   ndem = as.numeric(gf[ ,2]))
  
  pw <- merge(p, gf, 
              by.x = "pcn", by.y = "pcn", 
              all.x = TRUE, all.y = TRUE)
  pw[is.na(pw$ndem), ]$ndem <- 0
  pw$llu <- pw$lat/100 + pw$lon*100
  wu <- pw[!duplicated(pw$llu), ]

  bygenus <- agb_bgb %>%
    group_by(Genus) %>%
    summarise(fix = sum(totaln, na.rm = T))
  
  fixgr.grid <- array(dim = c(nlon, nlat))
  
  print("Starting iteration over grid")
  for(i in 1:nlon){
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(pw$lon>=lon.a & pw$lon<lon.b & pw$lat>=lat.a & pw$lat<lat.b)
      
      if(sum(inds) > 0){
        fixgr.grid[i,j] <- mean(pw[inds, ]$ndem, na.rm=TRUE)*foresta.grid[i,j]
      }else{
        fixgr.grid[i,j] <- NA
      }
      
    }
  }

  fixedn <- sum(fixgr.grid, na.rm = T)/1e9 #total N fixed across US from accretion method
  totalk <- sum(bygenus[ ,"fix"])
  acan <- (bygenus[[1,"fix"]]/totalk)*100 #*fixedn to get kgN instead of %
  albn <- (bygenus[[2,"fix"]]/totalk)*100
  alnn <- (bygenus[[3,"fix"]]/totalk)*100
  casn <- (bygenus[[4,"fix"]]/totalk)*100
  cern <- (bygenus[[5,"fix"]]/totalk)*100
  elan <- (bygenus[[6,"fix"]]/totalk)*100
  olnn <- (bygenus[[7,"fix"]]/totalk)*100
  pron <- (bygenus[[8,"fix"]]/totalk)*100
  robn <- (bygenus[[9,"fix"]]/totalk)*100
  
  output <- list(tgn=fixedn, acapct=acan, albpct=albn, alnpct=alnn, caspct=casn, 
                 cerpct=cern, ealpct=elan, olnpct=olnn, propct=pron, robpct=robn)
  return(output)
}

orig_fixedndfa <- plotfixndfa(agb_1_2, ndfatype="realistic")

# for realistic %Ndfa
bootlength <- 1000
bootstrapped2 <- matrix(nrow=1, ncol=10)
bootstrapped2 <- as.data.frame(bootstrapped2)
colnames(bootstrapped2) <- c("tgn", "acapct", "albpct", "alnpct", "caspct", "cerpct",
                            "elapct", "olnpct", "propct", "robpct")
for(i in 1:bootlength){
  bootstrapped2[i,] <- plotfixndfa(agb_1_2, ndfatype="realistic")
  print(i)
}

#saveRDS(bootstrapped2, "output/bootstrapped_ndfa.RDS")
mean(bootstrapped2$tgn, na.rm=T)*1e9/sum(foresta.grid, na.rm=T) #kgN/ha forest/yr for %ndfa
mean(bootstrapped2$tgn, na.rm=T)*1e9/sum(ground.grid, na.rm=T) #kgN/ha ground/yr for %ndfa

# for upper bound
bootstrapped3 <- plotfixndfa(agb_1_2, ndfatype="upperbound")

#saveRDS(bootstrapped3, "output/bootstrapped_upper.RDS")
bootstrapped3$tgn*1e9/sum(foresta.grid, na.rm=T) #kgN/ha forest/yr with upper bound
bootstrapped3$tgn*1e9/sum(ground.grid, na.rm=T) #kgN/ha ground/yr with upper bound

xtabs(~ Genus, agb_1_2)


# P contorta ---------------------------------
################
# P. contorta
# endophytic N fixer? (see Paul papers)
# Is it important to the total N fixed?
###############

# filter out P. contorta
pcontorta <- gdata %>% 
  filter(spcd == 1082 | spcd == 1081 | spcd == 108)
nrow(pcontorta) #192154
#saveRDS(pcontorta, "output/pcontorta.RDS")

# At Alnus Rate
pcon_fix <- pcontorta %>% 
  group_by(pcn) %>%
  summarize(plotBA=sum(BAm2*tpha,na.rm=T), cnt=n()) %>%
  mutate(fix = sum(al.samp[sample(1:n, 1, replace=T),]*c(1,plotBA)))

sum(pcon_fix$fix) #1266878 kgN (if fixing at same rate as Alnus)

totalacc <- 383778385

(totalacc+sum(pcon_fix$fix))/totalacc #percent change in value for fixed N

# 17 % change in the amount of N fixed
# Assumes that P. contorta fixes at same rate as Alnus
# probably too high so let's check with a lower rate

# At 10% Alnus rate
ch <- vector()
i <- 1
while(i<100){
  pcon_fix <- pcontorta %>% 
    group_by(pcn) %>%
    summarize(plotBA=sum(BAm2*tpha,na.rm=T), cnt=n()) %>%
    mutate(fix = sum(al.samp[sample(1:n, 1, replace=T),]*0.1*c(1,plotBA)))
  
  sum(pcon_fix$fix) #1346508398 kgN (if fixing at same rate as Alnus)
  
  ch[i] <- (totalacc+sum(pcon_fix$fix))/totalacc #percent change in value for fixed N
  i <- i+1
}
t.test(ch)

#for mapping at 10% alnus rate
source("scripts/genus_regressions_new.R")

pconto <- pcontorta %>% 
  group_by(pcn) %>%
  summarize(plotBA = sum(BAm2*tpha,na.rm=T), 
            cnt = n()) 
pcon_fix <- pconto %>%
  group_by(pcn) %>%
  mutate(fix = sum(pcon.samp[sample(1:n, 1, replace=T),]*c(1, plotBA)))
#saveRDS(pcon_fix, "output/pcon_fix.RDS")
sum(pcon_fix$fix) #123212.9

# 5.05 % change in the amount of N fixed
# Assumes that P. contorta fixes at 10% Alnus 
# This is what the Paul paper suggests

# At 1% Alnus rate
ch <- vector()
i <- 1
while(i<100){
  pcon_fix <- pcontorta %>% 
    group_by(pcn) %>%
    summarize(plotBA=sum(BAm2*tpha,na.rm=T), cnt=n()) %>%
    mutate(fix = sum(al.samp[sample(1:n, 1, replace=T),]*0.01*c(1,plotBA)))
  
  sum(pcon_fix$fix) #114322.4 kgN (if fixing at same rate as Alnus)
  
  ch[i] <- (totalacc+sum(pcon_fix$fix))/totalacc #percent change in value for fixed N
  i <- i+1
}
t.test(ch)


# 1.404 % change +/- 0.052 in the amount of N fixed
# Assumes that P. contorta fixes at 10% Alnus 
# This is what the Paul paper suggests

# At 0.5% Alnus rate
# while loop to bootstrap percent change for CI
ch <- vector()
i <- 1
while(i<100){
  pcon_fix <- pcontorta %>% 
    group_by(pcn) %>%
    summarize(plotBA=sum(BAm2*tpha,na.rm=T), cnt=n()) %>%
    mutate(fix = sum(al.samp[sample(1:n, 1, replace=T),]*0.005*c(1,plotBA)))
  
  sum(pcon_fix$fix) #63321393 kgN (if fixing at same rate as Alnus)
  
  ch[i] <- (totalacc+sum(pcon_fix$fix))/totalacc #percent change in value for fixed N
  i <- i+1
}
t.test(ch)

# 1.18 % +/- 0.033 change in the amount of N fixed
# Assumes that P. contorta fixes at 0.5% Alnus 



