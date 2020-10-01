
# sensitivity.R
# sensitivity analysis
# by Anika Petach 
# 10/8/18

# fators for %Ndfa method: 
#    fraction N in foliage, C:N of each tissue/genus, N fixer abundance, GR of each genus

# factors for accretion method: 
#    fixation rate (rate/BA slope), controls on fixation, N fixer abundance

require(dplyr)
library(gridExtra)
require(ggplot2)
require(MASS)
require(raster)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

set.seed(123) #to get consistent answers for checking code, for bootstrapping remove
# set.seed(.Random.seed) #if you want to remove set seed

# %ndfa method -----------------------------------------------------------------------------
#orig_fixedndfa <- 0.30 #TgN/yr

agb_1_2 <- readRDS("output/agb_1_2.RDS") #biomass measured and predicted for t2, fixers only
#FOR_COV_FOR_ALIGN <- raster("output/FOR_COV_FOR_ALIGN.grd") #forest area in km2 per grid cell, used compareRaster to check
#p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv", header=T)

agb_1_2 <- agb_1_2[agb_1_2$dbh2pred >= 2,] #jenkins equations only work for trees dbh>2
agb_1_2[agb_1_2$agb2pred < 0.5,]$agb2pred <- 0.05 #to prevent 0 in denominator
agb_1_2$agbd <- agb_1_2$agb2pred - agb_1_2$agb #get difference in agb
agb_1_2[agb_1_2$agbd < 0,]$agbd <- 0.05 #if no change assign small increment


# get Ndfa data
ndfa <- read.csv("data/ndfa_litreview.csv", header=TRUE)

# separate genera
aca <- ndfa[ndfa$Genus == "Acacia" & ndfa$Location != "pot",]
alb <- ndfa[ndfa$Genus == "Albizia" & ndfa$Location != "pot",]
ald <- ndfa[ndfa$Genus == "Alnus" & ndfa$Location != "pot",]
pro <- ndfa[ndfa$Genus == "Prosopis" & ndfa$Location != "pot",]
ela <- ndfa[ndfa$Genus == "Elaeagnus" & ndfa$Location != "pot",]
cas <- ndfa[ndfa$Genus == "Casuarina" & ndfa$Location != "pot",]
rob <- ndfa[ndfa$Genus == "Robinia" & ndfa$Location != "pot",]

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

# make gamma distributions
n = 10000
# aca.samp <- rnorm(n, aca.mean, aca.se)
# alb.samp <- rnorm(n, alb.mean, alb.se)
# ald.samp <- rnorm(n, ald.mean, ald.se)
# pro.samp <- rnorm(n, pro.mean, pro.se)
# ela.samp <- rnorm(n, ela.mean, ela.se)
# cas.samp <- rnorm(n, cas.mean, cas.se)
# rob.samp <- rnorm(n, rob.mean, rob.se)
# other.samp <- rnorm(n, other.mean, other.se)

source("scripts/functions.R") #to get newest version of plotfixndfa function

#orig_fixedndfa <- plotfixndfa(agb_1_2, ndfatype="realistic")

#orig_fixedndfa <- readRDS("output/bootstrapped_ndfa_new.RDS")
#orig_fixedndfa <- mean(orig_fixedndfa$tgn)

orig_fixedndfa <- 0.304

bootlength <- 10

## 1. fraction N in foliage

## 2. C:N of each tissue
bootstrap <- vector()
for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, cn.foliage = 35.1*0.5, ndfatype="realistic")$tgn
}
cnfoliage.max <- mean(bootstrap)
cnfoliage.max.sd <- sd(bootstrap)
cnfolup <- (cnfoliage.max - orig_fixedndfa) / orig_fixedndfa * 100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, cn.foliage = 35.1*1.5, ndfatype="realistic")$tgn
}
cnfoliage.min <- mean(bootstrap)
cnfoliage.min.sd <- sd(bootstrap)
cnfoldown <- (cnfoliage.min - orig_fixedndfa) / orig_fixedndfa * 100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, cn.wood = 350*0.5, ndfatype="realistic")$tgn
}
cnwood.max <- mean(bootstrap)
cnwood.max.sd <- sd(bootstrap)
cnwoodup <- (cnwood.max - orig_fixedndfa) / orig_fixedndfa * 100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, cn.wood = 350*1.5, ndfatype="realistic")$tgn
}
cnwood.min <- mean(bootstrap)
cnwood.min.sd <- sd(bootstrap)
cnwooddown <- (cnwood.min - orig_fixedndfa) / orig_fixedndfa * 100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, cn.fr = 43*0.5, ndfatype="realistic")$tgn
}
cnfr.max <- mean(bootstrap)
cnfr.max.sd <- sd(bootstrap)
cnfrup <- (cnfr.max - orig_fixedndfa) / orig_fixedndfa * 100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, cn.fr = 43*1.5, ndfatype="realistic")$tgn
}
cnfr.min <- mean(bootstrap)
cnfr.min.sd <- sd(bootstrap)
cnfrdown <- (cnfr.min - orig_fixedndfa) / orig_fixedndfa * 100

## 3. C:N of each genus

## 4. N fixer abundance
agb_upabundance <- agb_1_2
agb_upabundance$agb2pred <- agb_upabundance$agb2pred*1.5
for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_upabundance, ndfatype="realistic")$tgn
}
abundance.max <- mean(bootstrap)
abundance.max.sd <- sd(bootstrap)
abundupn <- (abundance.max - orig_fixedndfa) / orig_fixedndfa * 100

agb_downabundance <- agb_1_2
agb_downabundance$agb2pred <- agb_downabundance$agb2pred*0.5
for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_downabundance, ndfatype="realistic")$tgn
}
abundance.min <- mean(bootstrap)
abundance.min.sd <- sd(bootstrap)
abunddownn <- (abundance.min - orig_fixedndfa) / orig_fixedndfa * 100

## 5. GR splines

## 6. Resorption (of leaves and fine roots)
for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, res = 0.4*1.5, ndfatype="realistic")$tgn
}
res.max <- mean(bootstrap)
res.max.sd <- sd(bootstrap)
folresup <- (res.max - orig_fixedndfa) / orig_fixedndfa * 100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, res = 0.4*0.5, ndfatype="realistic")$tgn
}
res.min <- mean(bootstrap)
res.min.sd <- sd(bootstrap)
folresdown <- (res.min-orig_fixedndfa)/orig_fixedndfa*100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, resr = 0.73*1.5, ndfatype="realistic")$tgn
}
resfr.max <- mean(bootstrap)
resfr.max.sd <- sd(bootstrap)
frresup <- (resfr.max-orig_fixedndfa)/orig_fixedndfa*100

for(i in 1:bootlength){
  bootstrap[i] <- plotfixndfa(agb_1_2, resr = 0.73*0.5, ndfatype="realistic")$tgn
}
resfr.min <- mean(bootstrap)
resfr.min.sd <- sd(bootstrap)
frresdown <- (resfr.min-orig_fixedndfa)/orig_fixedndfa*100

# upper bound ------------------------------------------------------------------------
## 1. fraction N in foliage

## 2. C:N of each tissue
cnfoliage.max <- plotfixndfa(agb_1_2, cn.foliage=35.1*0.5, ndfatype="upperbound")$tgn
(cnfoliage.max-orig_fixedndfa)/orig_fixedndfa*100

cnfoliage.min <- plotfixndfa(agb_1_2, cn.foliage=35.1*1.5, ndfatype="upperbound")$tgn
(cnfoliage.min-orig_fixedndfa)/orig_fixedndfa*100

cnwood.max <- plotfixndfa(agb_1_2, cn.wood=350*0.5, ndfatype="upperbound")$tgn
(cnwood.max-orig_fixedndfa)/orig_fixedndfa*100

cnwood.min <- plotfixndfa(agb_1_2, cn.wood=350*1.5, ndfatype="upperbound")$tgn
(cnwood.min-orig_fixedndfa)/orig_fixedndfa*100

cnfr.max <- plotfixndfa(agb_1_2, cn.fr=43*0.5, ndfatype="upperbound")$tgn
(cnfr.max-orig_fixedndfa)/orig_fixedndfa*100

cnfr.min <- plotfixndfa(agb_1_2, cn.fr=43*1.5, ndfatype="upperbound")$tgn
(cnfr.min-orig_fixedndfa)/orig_fixedndfa*100

## 3. C:N of each genus

## 4. N fixer abundance
agb_upabundance <- agb_1_2
agb_upabundance$agb2pred <- agb_upabundance$agb2pred*1.5
abundance.max <- plotfixndfa(agb_upabundance, ndfatype="upperbound")$tgn
abundupn <- (abundance.max-orig_fixedndfa)/orig_fixedndfa*100

agb_downabundance <- agb_1_2
agb_downabundance$agb2pred <- agb_downabundance$agb2pred*0.5
abundance.min <- plotfixndfa(agb_downabundance, ndfatype="upperbound")$tgn
abunddownn <- (abundance.min-orig_fixedndfa)/orig_fixedndfa*100

## 5. GR splines

## 6. Resorption (of leaves and fine roots)
res.max <- plotfixndfa(agb_1_2, res=0.4*1.5, ndfatype="upperbound")$tgn
(res.max-orig_fixedndfa)/orig_fixedndfa*100

res.min <- plotfixndfa(agb_1_2, res=0.4*0.5, ndfatype="upperbound")$tgn
(res.min-orig_fixedndfa)/orig_fixedndfa*100



# accretion method -------------------------------------------------------------------------
# run once
orig_fixedn <- 0.878 #TgN/yr

# load data
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


source("scripts/genus_regressions_new.R")

n = 10000

#getrates function from functions.R

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

## plotfix function from functions.R

# plotdata should be fed allstems, 
# r should be fed output from getrates function, 
# canpos = 6 for all trees or 4 for canopy only
# plotfix <- function(plotdata, canpos, r) {
#   
#   allstemsf <- plotdata[plotdata$cclcd < canpos, ] 
#   print("filtered out understory")
#   
#   stems00 <- allstemsf %>%
#     group_by(pcn, Genus) %>%
#     summarise(pBA = sum(BAm2ha),
#               FIX = first(FIX),
#               ACT1vsRIZ0 = first(ACT1vsRIZ0))
#   print("agg to plot level")
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
#   print("assigned fixation rates")
#   
#   
#   
#   plot_sens$pcn <- as.integer64(plot_sens$pcn)
#   ps <- merge(plot_sens, psnap,
#               by = "pcn",
#               all.x = T, all.y = F)
#   
#   ps$fixedn <- ps$fix*ps$EXPCURR*haperacre #multiply N fixed in plot by area plot represents and convert to ha
#   ps2 <- ps[!is.na(ps$Genus),]
#   ps2[is.na(ps2$fixedn),]$fixedn <- 0
#   
#   anb <- as.data.frame.table(by(ps2$fixedn, ps2$Genus, sum)) #if not working do e$fix
#   anb <- data.frame(Genus = anb[ ,1], 
#                     fixboot = as.numeric(anb[ ,2]))
#   
#   print("aggregated to genus")
#   
#   
#   totalacc <- sum(anb$fixboot, na.rm = T)/1e9 #total N fixed across US from accretion method 
#   return(totalacc)
# }

source("scripts/functions.R")

normalrates <- getrates("standard")
orig_fixedn <- plotfix_new(allstems, 6, normalrates)
orig_fixedn <- readRDS("output/bootstrapped_acc.RDS")
orig_fixedn <- mean(orig_fixedn$tgn)

orig_fixedn <- 0.886

bootlength <- 3

output <- NULL
for(i in 1:100){
  output[i] <- plotfix_new(allstems, 6, normalrates)[[1]]
}

## 1. fixation regression slope of each genus
aca.max <- getrates("acacia", "max")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 6, aca.max)
  output[i] <- temp[[1]]
}
aca.max.fix <- mean(output)
aca.max.fix.sd <- sd(output)
acup <- (aca.max.fix - orig_fixedn)/orig_fixedn*100

for(i in 1:bootlength){
  aca.min <- getrates("acacia", "min")
  temp <- plotfix_new(allstems, 6, aca.min)
  output[i] <- temp[[1]]
}
aca.min.fix <- mean(output)
aca.min.fix.sd <- sd(output)
acdown <- (aca.min.fix - orig_fixedn)/orig_fixedn*100

pro.max <- getrates("prosopis", "max")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 6, pro.max)
  output[i] <- temp[[1]]
}
pro.max.fix <- mean(output)
pro.max.fix.sd <- sd(output)
proup <- (pro.max.fix - orig_fixedn)/orig_fixedn

for(i in 1:bootlength){
  pro.min <- getrates("prosopis", "min")
  temp <- plotfix_new(allstems, 6, pro.min)
  output[i] <- temp[[1]]
}
pro.min.fix <- mean(output)
pro.min.fix.sd <- sd(output)
prodown <- (pro.min.fix - orig_fixedn)/orig_fixedn*100

aln.max <- getrates("alnus", "max")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 6, aln.max)
  output[i] <- temp[[1]]
}
aln.max.fix <- mean(output)
aln.max.fix.sd <- sd(output)
alnup <- (aln.max.fix - orig_fixedn)/orig_fixedn*100

aln.min <- getrates("alnus", "min")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 6, aln.min)
  output[i] <- temp[[1]]
}
aln.min.fix <- mean(output)
aln.min.fix.sd <- sd(output)
alndown <-( aln.min.fix - orig_fixedn)/orig_fixedn*100

ot.max <- getrates("other", "max")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 6, ot.max)
  output[i] <- temp[[1]]
}
ot.max.fix <- mean(output)
ot.max.fix.sd <- sd(output)
otup <- (ot.max.fix - orig_fixedn)/orig_fixedn*100

for(i in 1:bootlength){
  ot.min <- getrates("other", "min")
  temp <- plotfix_new(allstems, 6, ot.min)
  output[i] <- temp[[1]]
}
ot.min.fix <- mean(output)
ot.min.fix.sd <- sd(output)
otdown <- (ot.min.fix - orig_fixedn)/orig_fixedn*100

rob.max <- getrates("robinia", "max")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 6, rob.max)
  output[i] <- temp[[1]]
}
rob.max.fix <- mean(output)
rob.max.fix.sd <- sd(output)
roup <- (rob.max.fix - orig_fixedn)/orig_fixedn*100

rob.min <- getrates("robinia", "min")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 6, rob.min)
  output[i] <- temp[[1]]
}
rob.min.fix <- mean(output)
rob.min.fix.sd <- sd(output)
rodown <- (rob.min.fix - orig_fixedn)/orig_fixedn*100


## 2. controls on fixation (canopy position)
normalrates <- getrates("standard")
for(i in 1:bootlength){
  temp <- plotfix_new(allstems, 4, normalrates)
  output[i] <- temp[[1]]
}
lightlim.fix <- mean(output)
lightlim.fix.sd <- sd(output)
lightdown <- (lightlim.fix - orig_fixedn)/orig_fixedn*100

## 3. N fixer abundance
normalrates <- getrates("standard")
allstems_upabundance <- allstems
allstems_upabundance$BAm2ha <- allstems_upabundance$BAm2ha*1.5
for(i in 1:bootlength){
  temp <- plotfix_new(allstems_upabundance, 6, normalrates)
  output[i] <- temp[[1]]
}
abundance.max.fix <- mean(output)
abundance.max.fix.sd <- sd(output)
abundup <- (abundance.max.fix - orig_fixedn)/orig_fixedn*100

allstems_downabundance <- allstems
allstems_downabundance$BAm2ha <- allstems_downabundance$BAm2ha*0.5
for(i in 1:bootlength){
  temp <- plotfix_new(allstems_downabundance, 6, normalrates)
  output[i] <- temp[[1]]
}
abundance.min.fix <- mean(output)
abundance.min.fix.sd <- sd(output)
abunddown <- (abundance.min.fix - orig_fixedn)/orig_fixedn*100

# Figure prep ---------------------------------------------------------
#Fig 5a prep
accsen <- data.frame(ac=c(acup, acdown), aln=c(alnup, alndown), pro=c(proup, prodown),
                     ot=c(otup, otdown), ro=c(roup, rodown), abund=c(abundup, abunddown),
                     light=c(NA, lightdown)) #no light
# accsen <- data.frame(ac=c(0.91, -0.66), aln=c(18.27, -17.63), pro=c(0.01, -0.79),
#                      ot=c(0.26, -0.72), ro=c(25.87, -25.91), abund=c(0.72, -0.33))
rownames(accsen) <- c("Upper", "Lower")
accsen <- as.matrix(accsen)
saveRDS(accsen, "output/forfigs/accsen_new.RDS")

#Fig 5b prep
ndfasen <- data.frame(cnfol=c(cnfolup, cnfoldown), cnwood=c(cnwoodup, cnwooddown),
                      cnfr=c(cnfrup, cnfrdown), abund=c(abundupn, abunddownn),
                      folres=c(folresup, folresdown), frres=c(frresup, frresdown))
# ndfasen <- data.frame(cnfol=c(34.33, -11.48), cnwood=c(63.58, -21.36),
#                       cnfr=c(1.42, -0.47), abund=c(17.25, -17.36),
#                       folres=c(17.14, -17.21), frres=c(0.71, -0.71))
rownames(ndfasen) <- c("Upper", "Lower")
ndfasen <- as.matrix(ndfasen)
saveRDS(ndfasen, "output/forfigs/ndfasen_new.RDS")

