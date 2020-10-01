
# map.R
# present processed FIA data
# maps of BNF by accretion method, %Ndfa method, BNF by genus, basal area
# by Anika Petach 
# 7/13/17

#rm(list = ls())

library(fields)
require(rgdal)
require(raster)
require(sp)
require(dplyr)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

# load cleaned data ----------------
# data cleaned in clean.R, sapling.R, cleanGR.R
# data processed in do_withBootstrapping.R, doGR.R
p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv", header=T)
p <- p[p$state!='PR'&p$state!='VI'&p$state!='HI', ]
p[p$lat == -999, ]$lat <- NA
p[p$lon == -999, ]$lon <- NA

#for sapling and tree data (growth rate %Ndfa)
gfdata <- readRDS("output/agb_bgb_perakis.RDS") #from doGR.R
#gfdata <- readRDS("output/agb_bgb_upper.RDS") #from doGR.R

#for saplings and tree data (accretion)
sf <- readRDS("output/stemsfix.RDS")
sf0 <- merge(sf, p, by="pcn", all.x=TRUE, all.y=FALSE)

#for light limitation
#sf <- readRDS("output/stemsfix_lightlim.RDS")

#saplings and tree data before fixation calculated
allstems <- readRDS("output/allstems.RDS") #for bootstrapping

#for assigning FIA region
regions <- read.csv("data/FIA_regions.csv", header=TRUE)

#forest cover data
FOR_COV_FOR_ALIGN <- raster("output/FOR_COV_FOR_ALIGN2.grd") #forest area in km2 per grid cell, used compareRaster to check
FOR_COV_ALL_ALIGN <- raster("output/FOR_COV_ALL_ALIGN2.grd") #ground area in km2 per grid cell, used compareRaster to check
FOR_COV_ALIGN <- raster("output/FOR_COV_ALIGN2.grd") #frac forest (forest/ground)

fixcut <- 0.9 # 1 for confirmed fixers
latcut <- 35

p$llu <- p$lat/100 + p$lon*100
pu <- p[!duplicated(p$llu), ]  # Number of plot record and unique plots

nrow(sf[sf$FIX == 1, ])  #Number of fixers
nrow(sf[sf$FIX == 0, ])  #Number of non-fixers

# map of plot density ------------------
# Get latitude and longitude bounds of AK, coterminous US, Mexico, PR, VI
latmin <- min(p$lat, na.rm=TRUE)
latmax <- max(p$lat, na.rm=TRUE)
lonmin <- min(p$lon, na.rm=TRUE)
lonmax <- max(p$lon, na.rm=TRUE)
latminint <- floor(latmin)
latmaxint <- ceiling(latmax)
lonminint <- floor(lonmin)
lonmaxint <- ceiling(lonmax)
dlat <- 1 # grid cell size
dlon <- 1 # grid cell size
lat.list <- seq(latminint, latmaxint, dlat)
lon.list <- seq(lonminint, lonmaxint, dlon)
nlat <- length(lat.list)
nlon <- length(lon.list)
#fiascaling <- 2428.114 #constant forest each plot represents
fia.avg.plotarea <- (7.3^2)*pi*4/10000

# Make matrix of plot numbers
nplot.grid <- plotarea.grid <- array(dim=c(nlon, nlat)) # plot area in ha
nuniqueplot.grid <- uniqueplotarea.grid <- array(dim=c(nlon, nlat)) # unique plot number and area in ha

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    nplot.grid[i,j] <- sum(p$lon>=lon.a & p$lon<lon.b & 
                             p$lat>=lat.a & p$lat<lat.b,na.rm=TRUE)
    nuniqueplot.grid[i,j] <- sum(pu$lon>=lon.a & pu$lon<lon.b & 
                                   pu$lat>=lat.a & pu$lat<lat.b,na.rm=TRUE)
    plotarea.grid[i,j] <- sum(p$lon>=lon.a & p$lon<lon.b & 
                                p$lat>=lat.a & p$lat<lat.b,na.rm=TRUE)*fia.avg.plotarea
    uniqueplotarea.grid[i,j] <- sum(pu$lon>=lon.a & pu$lon<lon.b & 
                                      pu$lat>=lat.a & pu$lat<lat.b,na.rm=TRUE)*fia.avg.plotarea 
    if(nplot.grid[i,j]==0){
      nplot.grid[i,j] <- NA
      nuniqueplot.grid[i,j] <- NA
      plotarea.grid[i,j] <- NA
      uniqueplotarea.grid[i,j] <- NA
    }
  }
}

nlev <- 64
ticks <- 2^(0:13)
image.plot(lon.list, lat.list, log(nplot.grid), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list( at=log(ticks),labels=ticks),
           xlab="longitude",ylab="latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main="Plot records")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

nlev <- 64
ticks <- 2^(0:13)
image.plot(lon.list, lat.list, log(nuniqueplot.grid), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list( at=log(ticks),labels=ticks),
           xlab="longitude",ylab="latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main="Number of plots")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)


# Plot area as well to make sure it worked ok
ticks <- 2^(0:13)
image.plot(lon.list, lat.list, log(uniqueplotarea.grid), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list( at=log(ticks),labels=ticks),
           xlab="longitude",ylab="latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main="Sampling area (ha)")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# Relative basal area, fixed N by accretion  -------------- 
# Need plot-level, grid-cell-level, and latitude-level values
# Need pba of all fixers, pba of rhizobials, pba of actinorhizals

# plot-level basal area
#d$pBA <- d$BAm2*d$tpha #don't need if using tree/sapling combined data (already calc)
d <- sf

fixcut <- 0.9
a <- as.data.frame.table(by(d$pBA, d$pcn, sum))
a <- data.frame(pcn = as.numeric(as.character(a[ ,1])), 
                tba = as.numeric(a[ ,2]))
p <- merge(p, a, 
           by = "pcn", 
           all.x = TRUE, all.y = FALSE) #if tba is 0 that means no fixers in plot

af <- as.data.frame.table(by(d[d$FIX >= fixcut, ]$pBA,
                             d[d$FIX >= fixcut, ]$pcn,
                             sum))
af <- data.frame(pcn=as.numeric(as.character(af[ ,1])),
                 fba=as.numeric(af[ ,2]))
ar <- as.data.frame.table(by(d[d$FIX >= fixcut & d$ACT1vsRIZ0 == 0, ]$pBA,
                             d[d$FIX >= fixcut & d$ACT1vsRIZ0 == 0, ]$pcn,
                             sum))
ar <- data.frame(pcn=as.numeric(as.character(ar[ ,1])),
                 rba=as.numeric(ar[ ,2]))
aa <- as.data.frame.table(by(d[d$FIX >= fixcut & d$ACT1vsRIZ0 == 1, ]$pBA,
                             d[d$FIX >= fixcut & d$ACT1vsRIZ0 == 1, ]$pcn,
                             sum))
aa <- data.frame(pcn=as.numeric(as.character(aa[ ,1])),
                 aba=as.numeric(aa[ ,2]))
p <- merge(p, af, by="pcn", all.x=TRUE, all.y=FALSE)
p <- merge(p, ar, by="pcn", all.x=TRUE, all.y=FALSE)
p <- merge(p, aa, by="pcn", all.x=TRUE, all.y=FALSE)
p[is.na(p$fba), ]$fba <- 0
p[is.na(p$rba), ]$rba <- 0
p[is.na(p$aba), ]$aba <- 0
p$pfba <- p$fba/p$tba*100 ##some tba are NA no fixers in plot
p$prba <- p$rba/p$tba*100
p$paba <- p$aba/p$tba*100

an <- as.data.frame.table(by(sf$fix, sf$pcn, sum)) #was mean before ... wrong
an <- data.frame(pcn = as.numeric(as.character(an[ ,1])), 
                 fixr = as.numeric(an[ ,2]))
p <- merge(p, an, 
           by = "pcn", 
           all.x = TRUE, all.y = FALSE)

#saveRDS(p, "output/p_from_map.RDS")

# pfba = plot basal area all fixers
# prba = plot basal area of rhizobials
# paba = plot basal area of actinorhizals

# grid-cell-level across continent
pfba.grid <- prba.grid <- paba.grid <- pfs.grid <- prs.grid <- pas.grid <- fba.grid <- fbagr.grid <- fn.grid <- fngr.grid <- forest.grid <- ground.grid <- forestfrac.grid <- array(dim=c(nlon,nlat))

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    lat.l <- 90.767  #length of a lat in km
    lon.l <- cos((lat.a+dlat/2)*3.14/180)*110.567 #length of a lon in km
    inds <- which(p$lon>=lon.a & p$lon<lon.b & p$lat>=lat.a & p$lat<lat.b)
    coords <- data.frame(lon=lon.list[i], lat=lat.list[j])
    
    if(sum(inds) > 0){
      pfba.grid[i,j] <- weighted.mean(p[inds,]$pfba, p[inds,]$tba, na.rm=TRUE)
      prba.grid[i,j] <- weighted.mean(p[inds,]$prba, p[inds,]$tba, na.rm=TRUE)
      paba.grid[i,j] <- weighted.mean(p[inds,]$paba, p[inds,]$tba, na.rm=TRUE)
      fba.grid[i,j] <- mean(p[inds, ]$fba, na.rm=TRUE)
      fbagr.grid[i,j] <- mean(p[inds, ]$fba, na.rm=TRUE)*(raster::extract(FOR_COV_FOR_ALIGN, coords)/raster::extract(FOR_COV_ALL_ALIGN, coords)) #Fixer BA per ground area
      fn.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)
      fngr.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)*raster::extract(FOR_COV_ALIGN, coords)
      
      forest.grid[i,j] <- raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
      ground.grid[i,j] <- raster::extract(FOR_COV_ALL_ALIGN, coords)*100 #convert from km2 to ha
      #forestfrac.grid[i,j] <- raster::extract(FOR_COV_ALIGN, coords)
    }else{
      pfba.grid[i,j] <- NA
      prba.grid[i,j] <- NA
      paba.grid[i,j] <- NA
      fba.grid[i,j] <- NA
      fbagr.grid[i,j] <- NA
      fn.grid[i,j] <- NA
     # fngr.grid[i,j] <- NA
      forest.grid[i,j] <- NA
      ground.grid[i,j] <- NA
     # forestfrac.grid[i,j] <- NA
    }
  }
}

mean(p$fixr, na.rm = TRUE) #average fixed N rate (plot scale)
mean(fn.grid, na.rm = TRUE) #averaged fixed N rate (grid scale)

max(fba.grid, na.rm=T) #max rate of N-fixer plot level BA
min(fba.grid, na.rm=T) #min rate of N-fixer plot level BA
mean(fba.grid, na.rm=T) #avg rate of N-fixer plot level BA

saveRDS(fn.grid,"output/fngrid_new.RDS") #accretion fixation rate per forest
#saveRDS(fba.grid, "output/fbagrid.RDS") #accretion fixation rate per ground
#saveRDS(fbagr.grid, "output/fbagrgrid.RDS")
saveRDS(fngr.grid, "output/fngrgrid_new.RDS")

##map of % BA here!
nlev <- 64
ticks <- 2^(0:13)
image.plot(lon.list, lat.list, pfba.grid, nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=(ticks),labels=ticks),
           xlab="longitude",ylab="latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main="Percent Fixing BA", cex.main=1.5)
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

######map of stand level BA here!
nlev <- 64
ticks <- 2^(-2:13)
image.plot(lon.list, lat.list, log(fba.grid+0.1), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=log(ticks),labels=ticks),
           xlab="longitude",ylab="latitude", 
           legend.lab="Basal area",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim=c(latminint-1,60+1),xlim=c(-140-2,lonmaxint+1),zlim=c(log(0.1),log(8)))
title(main=expression(paste("Fixer Basal Area (m"^"2"," ha"^"-1",")")), cex.main=1.5)
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

##map of fixed nitrogen rate average in cell
# fixed N per forest area
nlev <- 64
ticks <- 2^(-2:5)
image.plot(lon.list, lat.list, log(fn.grid+0.1), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=log(ticks),labels=ticks),
           xlab="longitude",ylab="latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1),
           zlim=c(log(0.1),log(35)))
#N Fixation Rate (kg N ha forest area-1 yr-1)
title(main=expression(paste("N Fixation Rate (kg N ha forest"^"-1"," yr"^"-1",")")), cex.main=1.5)
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# fixed N per ground area
nlev <- 64
ticks <- 2^(-2:5)
image.plot(lon.list, lat.list, log(fngr.grid+0.1), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=log(ticks),labels=ticks),
           xlab="longitude",ylab="latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1),zlim=c(log(0.1),log(35)))
title(main=expression(paste("N Fixation Rate (kg N ha ground"^"-1"," yr"^"-1",")")), cex.main=1.5)
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# forest fraction (0-1)
nlev <- 64
ticks <- 2^(-5:2)
image.plot(lon.list, lat.list, log(forestfrac.grid+0.01), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=log(ticks),labels=ticks),
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1),zlim=c(log(0.01),log(1)))

# forest fraction (not log transformed)
nlev <- 64
ticks <- c(0,0.25,0.50,0.75,1)
image.plot(lon.list, lat.list, forestfrac.grid, nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=(ticks),labels=ticks),
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1),zlim=c(0,1))
title(main=expression(paste("Fraction of grid cell forested")), cex.main=1.5)
map('world',regions='usa',add=TRUE)

nlev <- 64
ticks <- 2^(0:13)
image.plot(lon.list, lat.list, log(forest.grid+0.01), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=log(ticks),labels=ticks),
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1))

ticks <- c(0,100000,200000,300000,400000,500000,600000,700000)
image.plot(lon.list, lat.list, (forest.grid), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=(ticks),labels=ticks),
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1))
title(main=expression(paste("Forest area (ha)")), cex.main=1.5)
map('world',regions='usa',add=TRUE)

# Total understory, asymbiotic -------------------------------------------
sum(forest.grid,na.rm=T) #total forest area (209 million ha, FIA says should be 302 million but that includes all of AK)
sum(ground.grid,na.rm=T) #total land area ()
1.4*sum(forest.grid,na.rm=T) #total asymbiotic fixation (kgN/ha forest/yr)*(ha forest), 292923222
7.32*sum(forest.grid,na.rm=T) #total understory fixation (10% cover), (kgN/ha forest/yr)*(ha forest), 1531569988
2.20*sum(forest.grid,na.rm=T) #total understory fixation (3% cover), (kgN/ha forest/yr)*(ha forest), 460307920
12.44*sum(forest.grid,na.rm=T) #total understory fixation (17% cover),

hist(forestfrac.grid, na.rm=T, breaks=20,
     main="Histogram of fraction grid cell forested",
     xlab="Fraction cell forested")

# Sum total from accretion & regional breakdown ----------------
pr <- merge(p, regions, by.x="state", by.y="State.Abb", all.x=T)

ordering <- c("Interior West","Northern","Pacific Northwest","Southern")
pr$regionsFactor <- factor(pr$FIA.Region, levels=ordering)

fix.grid <- array(dim=c(nlon,nlat))
fixsp.acc <- data.frame()

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    inds <- which(pr$lon>=lon.a & pr$lon<lon.b & pr$lat>=lat.a & pr$lat<lat.b)
    coords <- data.frame(lon=lon.list[i], lat=lat.list[j])
    
    if(sum(inds) > 0){
      fix.grid[i,j] <- mean(pr[inds, ]$fixr, na.rm=TRUE)*raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
      reg <- pr[inds[1],]$regionsFactor
    }else{
      fix.grid[i,j] <- NA
      reg <- NA
    }
    
    out <- c(fix.grid[i,j], lon.list[i], lat.list[j], reg)
    fixsp.acc <- rbind(fixsp.acc, out)
  }
}

totalacc <- sum(fix.grid, na.rm=T) #total N fixed across US from accretion method (with bootstrap)

colnames(fixsp.acc) <- c("fixedn","lon","lat","fiaregion") #for making figure 4

saveRDS(fixsp.acc, "output/fixsp.acc.RDS")


# Accretion bootstrapping (at plot scale, save by genus) ----------------------------------------------------------
#allstems_bs <- readRDS("output/allstems_bs_396.RDS")
pb <- p

# Total sum from accretion bootstrapped
#source("scripts/genus_regressions.R") #need to actually open file and run it, source not working
source("scripts/genus_regressions_new.R") #need to actually open file and run it

stems0 <- allstems %>%
  filter(dbh > 0 & BAm2ha >0) %>%
  group_by(pcn, Genus) %>%
  summarise(pBA = sum(BAm2ha),
            FIX = first(FIX),
            ACT1vsRIZ0 = first(ACT1vsRIZ0))

fixedn <- acan <- albn <- alnn <- casn <- cern <- elan <- olnn <- pron <- robn <- NULL
fix.grid <- array(dim=c(nlon,nlat))

k <- 1
bootlength <- 100
while(k < bootlength){
  
  plot_sens <- stems0 %>%
    #filter(FIX == 1) %>%
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
  

  anb <- as.data.frame.table(by(plot_sens$fix, plot_sens$pcn, sum)) #if not working do e$fix
  anb <- data.frame(pcn = as.numeric(as.character(anb[ ,1])), 
                    fixboot = as.numeric(anb[ ,2]))
  pb <- merge(pb, anb, 
              by = "pcn", 
              all.x = TRUE, all.y = FALSE)
  
  bygenus <- plot_sens %>%
    group_by(Genus) %>%
    #summarise_all(funs(sum(., na.rm=T)))
    summarise(fix = sum(fix, na.rm = T))
  bygenus <- bygenus[c(1:5,7,8,10:11),] #remove Ebenopsis (too few stems), Piscidia (1 stem), Sophora (1 stem), NAs
  
  for(i in 1:nlon){
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(pb$lon>=lon.a & pb$lon<lon.b & pb$lat>=lat.a & pb$lat<lat.b)
      coords <- data.frame(lon = lon.list[i], 
                           lat = lat.list[j])
      
      if(sum(inds) > 0){
        fix.grid[i,j] <- mean(pb[inds, ]$fixboot, na.rm=TRUE)*raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
      }else{
        fix.grid[i,j] <- NA
      }
    }
  }
  
  fixedn[k] <- sum(fix.grid, na.rm = T) #total N fixed across US from accretion method
  totalk <- sum(bygenus[ ,"fix"])
  acan[k] <- (bygenus[[1,"fix"]]/totalk)*fixedn[k]
  albn[k] <- (bygenus[[2,"fix"]]/totalk)*fixedn[k]
  alnn[k] <- (bygenus[[3,"fix"]]/totalk)*fixedn[k]
  casn[k] <- (bygenus[[4,"fix"]]/totalk)*fixedn[k]
  cern[k] <- (bygenus[[5,"fix"]]/totalk)*fixedn[k]
  elan[k] <- (bygenus[[6,"fix"]]/totalk)*fixedn[k]
  olnn[k] <- (bygenus[[7,"fix"]]/totalk)*fixedn[k]
  pron[k] <- (bygenus[[8,"fix"]]/totalk)*fixedn[k]
  robn[k] <- (bygenus[[9,"fix"]]/totalk)*fixedn[k]
  
  pb <- subset(pb, select = -c(fixboot))
  
  k <- k+1
  print(k/bootlength*100)
}

bootstrap <- cbind(fixedn, acan, albn, alnn, casn, cern, elan, olnn, pron, robn)
bootstrap <- as.data.frame(bootstrap)
#saveRDS(bootstrap, "output/accretion_boot_700.RDS")
saveRDS(bootstrap, "output/accretion_boot_100.RDS")

#convert all values from kg to tg
bootstraptg <- bootstrap/1e9

#get mean and sd for each column
# do we want sd, se, or mean-95% CI? mean(bootstrap$fixedn/1e9)-quantile(bootstrap$fixedn/1e9,0.025)
gs <- bootstraptg %>% summarise_all(funs(sd(., na.rm = TRUE),mean(., na.rm=TRUE)))



# Accretion bootstrapping (at plot scale, save by grid cell) ----------------------------------------------------------
#allstems_bs <- readRDS("output/allstems_bs_396.RDS")
pb <- p

# Total sum from accretion bootstrapped
#source("scripts/genus_regressions.R") #need to actually open file and run it, source not working
source("scripts/genus_regressions_new.R")

stems0 <- allstems %>%
  group_by(pcn, Genus) %>%
  summarise(pBA = sum(BAm2ha),
            FIX = first(FIX),
            ACT1vsRIZ0 = first(ACT1vsRIZ0))

bootlength <- 10
fixedn <- NULL
fix.grid <- array(dim=c(nlon,nlat,bootlength))

for(k in 1:bootlength){
  
  plot_sens <- stems0 %>%
    #filter(FIX == 1) %>%
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
  
  
  anb <- as.data.frame.table(by(plot_sens$fix, plot_sens$pcn, sum)) #if not working do e$fix
  anb <- data.frame(pcn = as.numeric(as.character(anb[ ,1])), 
                    fixboot = as.numeric(anb[ ,2]))
  pb <- merge(pb, anb, 
              by = "pcn", 
              all.x = TRUE, all.y = FALSE)
  
  for(i in 1:nlon){
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(pb$lon>=lon.a & pb$lon<lon.b & pb$lat>=lat.a & pb$lat<lat.b)
      coords <- data.frame(lon = lon.list[i], 
                           lat = lat.list[j])
      
      if(sum(inds) > 0){
        fix.grid[i,j,k] <- mean(pb[inds, ]$fixboot, na.rm = TRUE)
        
        #fix.grid[i,j,k] <- mean(pb[inds, ]$fixboot, na.rm=TRUE)*raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
      }else{
        fix.grid[i,j,k] <- NA
      }
    }
  }
  
  fixedn[k] <- sum(fix.grid, na.rm = T) #total N fixed across US from accretion method
  
  pb <- subset(pb, select = -c(fixboot))
  print(k/bootlength*100)
}

means <- apply(fix.grid, c(1,2), mean)
sds <- apply(fix.grid, c(1,2), sd)

bootstrap <- cbind(fixedn, acan, albn, alnn, casn, cern, elan, olnn, pron, robn)
bootstrap <- as.data.frame(bootstrap)
#saveRDS(bootstrap, "output/accretion_boot_700.RDS")
#convert all values from kg to tg
bootstraptg <- bootstrap/1e9

image.plot(lon.list, lat.list, log(sds+0.001))

nlev <- 64
ticks <- 2^(-2:5)
image.plot(lon.list, lat.list, log(sds+0.1), nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=log(ticks),labels=ticks),
           xlab="",ylab="",
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1),zlim=c(log(0.1),log(35)))
title(main=expression('Accretion Uncertainty'), cex.main=1.5)
mtext("(d)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# Find grid cells for testing -----------------------------------------------------------------------
save.grid <- array(dim = c(nlon, nlat))
for(i in 1:nlon){
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    inds <- which(pb$lon>=lon.a & pb$lon<lon.b & pb$lat>=lat.a & pb$lat<lat.b)
    
    save.grid[i,j] <- length(inds)
  }
}




# Accretion bootstrapping (at grid cell scale) ------------------------------------------------------
pb <- p

# Make forest area grid

foresta.grid <- array(dim = c(nlon, nlat))

for(i in 1:nlon){
  for(j in 1:nlat){
    coords <- data.frame(lon = lon.list[i], 
                         lat = lat.list[j])
    foresta.grid[i,j] <- raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
  }
}

# Total sum from accretion bootstrapped
#source("scripts/genus_regressions.R") #need to actually open file and run it, source not working
source("scripts/genus_regressions.R")

#get total plot BA for each genus in each plot
stems0 <- allstems %>%
  filter(FIX == 1) %>%
  filter(!Genus %in% c("Casuarina","Sophora","Piscidia")) %>%
  group_by(pcn, Genus) %>%
  summarise(pBA = sum(BAm2ha),
            FIX = first(FIX),
            ACT1vsRIZ0 = first(ACT1vsRIZ0))

fixedn <- acan <- albn <- alnn <- casn <- cern <- elan <- olnn <- pron <- robn <- NULL
fix.grid <- acaf.grid <- albf.grid <- alnf.grid <- casf.grid <- cerf.grid <- elaf.grid <- olnf.grid <- prof.grid <- robf.grid <- array(dim = c(nlon, nlat))

k <- 1
bootlength <- 5


while(k < bootlength){

  for(i in 1:nlon){
    print(i/nlon*100)
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(pb$lon>=lon.a & pb$lon<lon.b & pb$lat>=lat.a & pb$lat<lat.b)
      pcns <- pb[inds,]$pcn
      
      #get mean plot BA for each genus
      s0 <- stems0 %>% 
        filter(pcn %in% pcns) %>%
        group_by(Genus) %>%
        summarize(gridBA = mean(pBA, na.rm=T)) 
      
      if(nrow(s0) > 0){
        #get BNF rate for each genus in the grid cell
        grid_sens <- s0 %>%
          mutate(fix = ifelse(Genus %in% c("Acacia"), sum(ac.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)), 
                              ifelse(Genus %in% c("Albizia"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)),
                                     ifelse(Genus %in% c("Alnus"), sum(al.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)),
                                            ifelse(Genus %in% c("Cercocarpus"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)),
                                                   ifelse(Genus %in% c("Prosopis"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)), 
                                                          ifelse(Genus %in% c("Ebenopsis"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)),
                                                                 ifelse(Genus %in% c("Robinia"),rnorm(1,mean(robinia$rate),var(robinia$rate)^0.5),
                                                                        ifelse(Genus %in% c("Olneya"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)),
                                                                               ifelse(Genus %in% c("Elaeagnus"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,gridBA)),
                                                                                      0)))))))))) 
        
        #make all genera present so don't need if statements
        Genus <- c("Acacia","Albizia","Alnus","Cercocarpus","Prosopis","Ebenopsis","Robinia","Olneya","Elaeagnus")
        gridBA <- rep(0,9)
        fix <- rep(0,9)
        grid_fix <- data.frame(Genus, gridBA, fix)
        
        grid_fix$gridBA[match(grid_sens$Genus,grid_fix$Genus)] <- grid_sens$gridBA
        grid_fix$fix[match(grid_sens$Genus,grid_fix$Genus)] <- grid_sens$fix
        
        fix.grid[i,j] <- sum(grid_fix$fix, na.rm=TRUE)*forest.grid[i,j]
        acaf.grid[i,j] <- grid_fix[Genus=="Acacia",]$fix*forest.grid[i,j]
        albf.grid[i,j] <- grid_fix[Genus=="Albizia",]$fix*forest.grid[i,j]
        alnf.grid[i,j] <- grid_fix[Genus=="Alnus",]$fix*forest.grid[i,j]
        cerf.grid[i,j] <- grid_fix[Genus=="Cercocarpus",]$fix*forest.grid[i,j]
        elaf.grid[i,j] <- grid_fix[Genus=="Elaeagnus",]$fix*forest.grid[i,j]
        olnf.grid[i,j] <- grid_fix[Genus=="Olneya",]$fix*forest.grid[i,j]
        prof.grid[i,j] <- grid_fix[Genus=="Prosopis",]$fix*forest.grid[i,j]
        robf.grid[i,j] <- grid_fix[Genus=="Robinia",]$fix*forest.grid[i,j]
        
      }else{
        fix.grid[i,j] <- NA
        acaf.grid[i,j] <- NA
        albf.grid[i,j] <- NA
        alnf.grid[i,j] <- NA
        cerf.grid[i,j] <- NA
        elaf.grid[i,j] <- NA
        olnf.grid[i,j] <- NA
        prof.grid[i,j] <- NA
        robf.grid[i,j] <- NA
      }
    }
  }
    
  fixedn[k] <- sum(fix.grid, na.rm=T) #total N fixed across US from accretion method
  acan[k] <- sum(acaf.grid, na.rm=T)/fixedn[k]
  albn[k] <- sum(albf.grid, na.rm=T)/fixedn[k]
  alnn[k] <- sum(alnf.grid, na.rm=T)/fixedn[k]
  casn[k] <- sum(casf.grid, na.rm=T)/fixedn[k]
  cern[k] <- sum(cerf.grid, na.rm=T)/fixedn[k]
  elan[k] <- sum(elaf.grid, na.rm=T)/fixedn[k]
  olnn[k] <- sum(olnf.grid, na.rm=T)/fixedn[k]
  pron[k] <- sum(prof.grid, na.rm=T)/fixedn[k]
  robn[k] <- sum(robf.grid, na.rm=T)/fixedn[k]
  
  k <- k+1
  print(k/bootlength*100)
}


bootstrap <- cbind(fixedn, acan, albn, alnn, casn, cern, elan, olnn, pron, robn)
bootstrap <- as.data.frame(bootstrap)
#saveRDS(bootstrap, "output/accretion_boot_700.RDS")
#convert all values from kg to tg
bootstraptg <- bootstrap/1e9

#get mean and sd for each column
bootstraptg_plot <- readRDS("output/accretion_boot_700.RDS")
# do we want sd, se, or mean-95% CI? mean(bootstrap$fixedn/1e9)-quantile(bootstrap$fixedn/1e9,0.025)
gs <- bootstraptg_plot %>% summarise_all(funs(sd(., na.rm = TRUE),mean(., na.rm=TRUE)))




# Accretion bootstrapping (at continent scale) -------------------------------------------------
# change line 715 depending if light limited or not
pb <- p

# Total sum from accretion bootstrapped
#source("scripts/genus_regressions.R") #need to actually open file and run it, source not working
source("scripts/genus_regressions_new.R")

stems0 <- allstems %>%
  filter(cclcd < 4) %>% #keep when doing light limited
  group_by(pcn, Genus) %>%
  summarise(pBA = sum(BAm2ha),
            FIX = first(FIX),
            ACT1vsRIZ0 = first(ACT1vsRIZ0))

fixedn <- acan <- albn <- alnn <- casn <- cern <- elan <- olnn <- pron <- robn <- NULL
fix.grid <- array(dim = c(nlon, nlat))

k <- 1
bootlength <- 100

while(k < bootlength){
  
  #define BNF slope and int for whole continent
  ac <- ac.samp[sample(1:n, 1, replace=T),]
  ot <- ot.samp[sample(1:n, 1, replace=T),]
  al <- al.samp[sample(1:n, 1, replace=T),]
  pr <- pr.samp[sample(1:n, 1, replace=T),]
  ro <- rnorm(1, mean(robinia$rate), var(robinia$rate)^0.5)
  
  plot_sens <- stems0 %>%
    #filter(FIX == 1) %>%
    group_by(pcn, Genus) %>%
    mutate(fix = ifelse(Genus %in% c("Acacia"), sum(ac*c(1,pBA)), 
                        ifelse(Genus %in% c("Albizia"), sum(ot*c(1,pBA)),
                               ifelse(Genus %in% c("Alnus"), sum(al*c(1,pBA)),
                                      ifelse(Genus %in% c("Cercocarpus"), sum(pr*c(1,pBA)),
                                             ifelse(Genus %in% c("Prosopis"), sum(pr*c(1,pBA)), 
                                                    ifelse(Genus %in% c("Ebenopsis"), sum(ot*c(1,pBA)),
                                                           ifelse(Genus %in% c("Robinia"),ro,
                                                                  ifelse(Genus %in% c("Olneya"), sum(ot*c(1,pBA)),
                                                                         ifelse(Genus %in% c("Elaeagnus"), sum(ot*c(1,pBA)),
                                                                                0)))))))))) 
  
  
  anb <- as.data.frame.table(by(plot_sens$fix, plot_sens$pcn, sum)) #if not working do e$fix
  anb <- data.frame(pcn = as.numeric(as.character(anb[ ,1])), 
                    fixboot = as.numeric(anb[ ,2]))
  pb <- merge(pb, anb, 
              by = "pcn", 
              all.x = TRUE, all.y = FALSE)
  
  bygenus <- plot_sens %>%
    group_by(Genus) %>%
    #summarise_all(funs(sum(., na.rm=T)))
    summarise(fix = sum(fix, na.rm = T))
  bygenus <- bygenus[c(1:5,7,8,10:11),] #remove Ebenopsis (too few stems), Piscidia (1 stem), Sophora (1 stem), NAs
  
  for(i in 1:nlon){
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(pb$lon>=lon.a & pb$lon<lon.b & pb$lat>=lat.a & pb$lat<lat.b)
      coords <- data.frame(lon = lon.list[i], 
                           lat = lat.list[j])
      
      if(sum(inds) > 0){
        fix.grid[i,j] <- mean(pb[inds, ]$fixboot, na.rm=TRUE)*raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
      }else{
        fix.grid[i,j] <- NA
      }
    }
  }
  
  fixedn[k] <- sum(fix.grid, na.rm = T) #total N fixed across US from accretion method
  totalk <- sum(bygenus[ ,"fix"])
  acan[k] <- (bygenus[[1,"fix"]]/totalk)*fixedn[k]
  albn[k] <- (bygenus[[2,"fix"]]/totalk)*fixedn[k]
  alnn[k] <- (bygenus[[3,"fix"]]/totalk)*fixedn[k]
  casn[k] <- (bygenus[[4,"fix"]]/totalk)*fixedn[k]
  cern[k] <- (bygenus[[5,"fix"]]/totalk)*fixedn[k]
  elan[k] <- (bygenus[[6,"fix"]]/totalk)*fixedn[k]
  olnn[k] <- (bygenus[[7,"fix"]]/totalk)*fixedn[k]
  pron[k] <- (bygenus[[8,"fix"]]/totalk)*fixedn[k]
  robn[k] <- (bygenus[[9,"fix"]]/totalk)*fixedn[k]
  
  pb <- subset(pb, select = -c(fixboot))
  
  k <- k+1
  print(k/bootlength*100)
}

bootstrap <- cbind(fixedn, acan, albn, alnn, casn, cern, elan, olnn, pron, robn)
bootstrap <- as.data.frame(bootstrap)
#saveRDS(bootstrap, "output/accretion_boot_continent.RDS")
#convert all values from kg to tg
bootstraptg <- bootstrap/1e9
saveRDS(bootstrap, "output/accretion_boot_lightlim_continent.RDS")

#get mean and sd for each column
# do we want sd, se, or mean-95% CI? mean(bootstrap$fixedn/1e9)-quantile(bootstrap$fixedn/1e9,0.025)
gs <- bootstraptg %>% summarise_all(funs(sd(., na.rm = TRUE),mean(., na.rm=TRUE)))




# Fig S1: Genus specific fixation --------------
# Not necessary because get the same results from bootstrapped
# acacia, albizia, olneya, prosopis, robinia, alnus, casuarina, cercocarpus, elaeagnus
# aca.grid = acacia, acba
# alb.grid = albizia, alba
# oln.grid = olneya, olba
# pro.grid = prosopis, poba
# rob.grid = robinia, roba
# aln.grid = alnus, anba
# cas.grid = casuarina, caba
# cer.grid = cercocarpus, ceba
# ela.grid = elaeagnus, elba

# acacia, albizia, olneya, prosopis, robinia, alnus, casuarina, cercocarpus, elaeagnus
aca <- sf0 %>% filter(Genus %in% c("Acacia"))
alb <- sf0 %>% filter(Genus %in% c("Albizia"))
oln <- sf0 %>% filter(Genus %in% c("Olneya"))
pro <- sf0 %>% filter(Genus %in% c("Prosopis"))
rob <- sf0 %>% filter(Genus %in% c("Robinia"))
aln <- sf0 %>% filter(Genus %in% c("Alnus"))
cas <- sf0 %>% filter(Genus %in% c("Casuarina"))
cer <- sf0 %>% filter(Genus %in% c("Cercocarpus"))
ela <- sf0 %>% filter(Genus %in% c("Elaeagnus"))
pco <- readRDS("output/pcontorta.RDS")

#get mean fixation rate by genus
se <- function(x) sd(x)/sqrt(length(x))
mean(aca$fix.r)
se(aca$fix.r) #standard error
mean(alb$fix.r)
se(alb$fix.r)
mean(oln$fix.r)
se(oln$fix.r)
mean(pro$fix.r)
se(pro$fix.r)
mean(rob$fix.r)
se(rob$fix.r)
mean(aln$fix.r)
se(aln$fix.r)
mean(cas$fix.r)
se(cas$fix.r)
mean(cer$fix.r)
se(cer$fix.r)
mean(ela$fix.r)
se(ela$fix.r)

# grid-cell-level across continent
aca.grid <- alb.grid <- oln.grid <- pro.grid <- rob.grid <- aln.grid <- cas.grid <- cer.grid <- ela.grid <- array(dim=c(nlon,nlat))
dfList <- list(aca=aca, alb=alb, aln=aln, pro=pro, rob=rob, aln=aln, cas=cas, cer=cer, ela=ela)

for(k in 1:length(dfList)){
  m <- get(names(dfList)[k])
  q <- paste(names(dfList)[k], ".grid", sep="")
  out <- assign(q, array(dim=c(nlon,nlat)))
  
  for(i in 1:nlon){
    print(i/nlon*100)
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(m$lon>=lon.a & m$lon<lon.b & m$lat>=lat.a & m$lat<lat.b)
      
      if(k == 1 & sum(inds) > 0){
        aca.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 1 & sum(inds) <= 0){
        aca.grid[i,j] <- NA
      }else if(k == 2 & sum(inds) > 0){
        alb.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 2 & sum(inds) <= 0){
        alb.grid[i,j] <- NA
      }else if(k == 3 & sum(inds) > 0){
        oln.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 3 & sum(inds) <= 0){
        oln.grid[i,j] <- NA
      }else if(k == 4 & sum(inds) > 0){
        pro.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 4 & sum(inds) <= 0){
        pro.grid[i,j] <- NA
      }else if(k == 5 & sum(inds) > 0){
        rob.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 5 & sum(inds) <= 0){
        rob.grid[i,j] <- NA
      }else if(k == 6 & sum(inds) > 0){
        aln.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 6 & sum(inds) <= 0){
        aln.grid[i,j] <- NA
      }else if(k == 7 & sum(inds) > 0){
        cas.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 7 & sum(inds) <= 0){
        cas.grid[i,j] <- NA
      }else if(k == 8 & sum(inds) > 0){
        cer.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 8 & sum(inds) <= 0){
        cer.grid[i,j] <- NA
      }else if(k == 9 & sum(inds) > 0){
        ela.grid[i,j] <- mean(m[inds, ]$fix.r, na.rm=TRUE)
      }else if(k == 9 & sum(inds) <= 0){
        ela.grid[i,j] <- NA
      }
    }
  }
}

######map of genus specific fixation
png("output/fixationmap_riz.png", width = 12, height = 6, units = 'in', res = 300)
par(mfrow = c(2, 3))
par(mar = c(3.1, 5.1, 4.1, 2.1)+0.1,
    oma = c(1, 1, 1, 4),
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
# Acacia
nlev <- 64
#ticks <- 2^(0:13)
ticks <- c(2, 32, 64, 128)
image.plot(lon.list, lat.list, aca.grid, nlevel=nlev, col=tim.colors(nlev),
           axis.args = list(at=(ticks), labels=ticks, las=3),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Acacia))), cex.main = 1.5)
mtext("(a)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Albizia
image.plot(lon.list, lat.list, alb.grid, nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=(ticks),labels=ticks),
           xlab="",ylab="",
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1),zlim=c(0,155))
title(main=substitute(paste("", italic(Albizia))), cex.main=1.5)
mtext("(b)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# Olneya
image.plot(lon.list, lat.list, oln.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           legend.args = list(text = expression(paste("BNF Rate (kg N ha"^"-1"," yr"^"-1",")")),
                              cex = 1,
                              side = 4,
                              line = 4),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Olneya))), cex.main = 1.5)
mtext("(c)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Prosopis
image.plot(lon.list, lat.list, pro.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Prosopis))), cex.main = 1.5)
mtext("(d)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Robinia
image.plot(lon.list, lat.list, rob.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           legend.args = list(text = expression(paste("BNF Rate (kg N ha"^"-1"," yr"^"-1",")")),
                              cex = 1,
                              side = 4,
                              line = 4),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Robinia))), cex.main = 1.5)
mtext("(e)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)
mtext("BNF Rate by Genus", outer = TRUE, cex = 1.3, line = -0.5)
dev.off()

# Actinorhizal
png("output/fixationmap_act.png", width = 12, height = 6, units = 'in', res = 300)
par(mfrow = c(2, 3))
par(mar = c(3.1, 5.1, 4.1, 2.1)+0.1,
    oma = c(1, 1, 1, 4),
    cex.main = 1.5,
    cex.lab = 1.5, cex.axis = 1.3,
    bty = "n")
# Alnus
image.plot(lon.list, lat.list, aln.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks),labels = ticks),
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Alnus))), cex.main = 1.5)
mtext("(f)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Casuarina
image.plot(lon.list, lat.list, cas.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Casuarina))), cex.main = 1.5)
mtext("(g)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Cercocarpus
image.plot(lon.list, lat.list, cer.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           legend.args = list(text=expression(paste("BNF Rate (kg N ha"^"-1"," yr"^"-1",")")), 
                            cex = 1, 
                            side = 4, 
                            line = 4),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Cercocarpus))), cex.main = 1.5)
mtext("(h)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Elaeagnus
image.plot(lon.list, lat.list, ela.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           legend.args = list(text = expression(paste("BNF Rate (kg N ha"^"-1"," yr"^"-1",")")), cex=1, side=4, line=4),           
           xlab="",ylab="",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
title(main = substitute(paste("", italic(Elaeagnus))), cex.main = 1.5)
mtext("(i)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)
dev.off()

# %ndfa -------------
# Remove extreme outliers
#gfdata <- gfdata[gfdata$wood < 2,]
#gfdata <- gfdata[gfdata$total < 10 & gfdata$total > 0,]
# Calculate demand per area
#gfdata$ndem <- gfdata$total*gfdata$tpha1
gfdata <- gfdata[gfdata$totaln < 100, ] #totaln is kgN/ha

#gf <- as.data.frame.table(by(gfdata[gfdata$FIX>=fixcut,]$totaln,gfdata[gfdata$FIX>=fixcut,]$pcn,sum))
gf <- as.data.frame.table(by(gfdata$totaln, gfdata$pcn, sum))
gf <- data.frame(pcn = as.numeric(as.character(gf[ ,1])), 
                 ndem = as.numeric(gf[ ,2]))

p$pcn <- as.integer64(p$pcn)
gf$pcn <- as.integer64(gf$pcn)

pw <- merge(p, gf, 
            by.x = "pcn", by.y = "pcn", 
            all.x = TRUE, all.y = TRUE)
pw[is.na(pw$ndem),]$ndem <- 0
pw$llu <- pw$lat/100 + pw$lon*100
wu <- pw[!duplicated(pw$llu), ]
# Number of plot record and unique plots
print(nrow(pw))
print(nrow(wu))

# need to sum ndfa demand to plot level
# then in loop take the mean value
# divide mean value by size of fia plot in ha

grn.grid <- grngr.grid <- array(dim = c(nlon, nlat))
for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    
    inds <- which(pw$lon >= lon.a & pw$lon < lon.b & pw$lat >= lat.a & pw$lat < lat.b) #normal data
    coords <- data.frame(lon = lon.list[i], 
                         lat = lat.list[j])
    
    if(sum(inds) > 0){
      grn.grid[i,j] <- mean(pw[inds, ]$ndem, na.rm=T)
      grngr.grid[i,j] <- mean(pw[inds, ]$ndem, na.rm=T)*raster::extract(FOR_COV_ALIGN, coords) #was FOR_COV_ALL_ALIGN but I think that's wrong
    }else{
      grn.grid[i,j] <- NA
      grngr.grid[i,j] <- NA
    }
  }
}

mean(grn.grid, na.rm = T) #mean fixation rate per forest area
sum(grn.grid, na.rm=T)

# save output 
#saveRDS(grn.grid, "output/grn_grid_new.RDS")
#saveRDS(grngr.grid, "output/grngr_grid_new.RDS")

##map of N demand based on ndfa
nlev <- 64
ticks <- 2^(-2:13)
image.plot(lon.list, lat.list, log(grn.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), labels = round(ticks, digits = 2)),
           xlab = "longitude", ylab = "latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1))
title(main = expression(paste("N demand by fixers (kg N ha forest"^"-1"," yr"^"-1",")")),
      cex.main = 1.5)
map('world', regions = 'usa', add = TRUE)

##map of N demand based on ndfa
nlev <- 64
ticks <- 2^(-2:13)
image.plot(lon.list, lat.list, log(grngr.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), labels = round(ticks, digits = 2)),
           xlab = "longitude", ylab = "latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1))
title(main = expression(paste("N demand by fixers (kg N ha ground"^"-1"," yr"^"-1",")")),
      cex.main = 1.5)
map('world', regions = 'usa', add = TRUE)

##Fig S2: map of N demand for Populus deltoides ------------
# run section above with gfdata <- agb_bgb from the populus section in doGR.R
grn.grid.1 <- grn.grid
grn.grid.1[grn.grid.1==0] <- NA

png("output/figureS2b.png", width = 12, height = 8, units = 'in', res = 300)
par(bty = "n")
nlev <- 64
ticks <- 2^(-2:13)
image.plot(lon.list, lat.list, log(grn.grid.1+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), labels = round(ticks, digits = 2)),
           legend.args = list(text = expression(paste("BNF Rate (kg N ha"^"-1"," yr"^"-1",")")), 
                              cex=1, side=4, line=3),
           xlab = "longitude", ylab = "latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1))
title(main = substitute(paste("", italic(Populus))), cex.main = 1.5)
map('world', regions = 'usa', add = TRUE)
    
dev.off()

## Fig S2a: P contorta map --------------------------------
# grid-cell-level across continent
pcon_fix <- readRDS("output/pcon_fix.RDS")
p.fn.grid <- p.fngr.grid <- forest.grid <- ground.grid <- forestfrac.grid <- array(dim=c(nlon,nlat))

p <- merge(p, pcon_fix, by="pcn", all=T)

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    lat.l <- 90.767  #length of a lat in km
    lon.l <- cos((lat.a+dlat/2)*3.14/180)*110.567 #length of a lon in km
    inds <- which(p$lon>=lon.a & p$lon<lon.b & p$lat>=lat.a & p$lat<lat.b)
    coords <- data.frame(lon=lon.list[i], lat=lat.list[j])
    
    if(sum(inds) > 0){
      p.fn.grid[i,j] <- mean(p[inds, ]$fix, na.rm=TRUE)
      p.fngr.grid[i,j] <- mean(p[inds, ]$fix, na.rm=TRUE)*raster::extract(FOR_COV_ALIGN, coords)
      
      #forest.grid[i,j] <- raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
      #ground.grid[i,j] <- raster::extract(FOR_COV_ALL_ALIGN, coords)*100 #convert from km2 to ha
      #forestfrac.grid[i,j] <- raster::extract(FOR_COV_ALIGN, coords)
    }else{
      p.fn.grid[i,j] <- NA
      p.fngr.grid[i,j] <- NA
      #forest.grid[i,j] <- NA
      #ground.grid[i,j] <- NA
      # forestfrac.grid[i,j] <- NA
    }
  }
}

png("output/figureS2a.png", width = 12, height = 8, units = 'in', res = 300)
par(bty = "n")
nlev <- 64
ticks <- 2^(-2:13)
image.plot(lon.list, lat.list, log(p.fngr.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), labels = round(ticks, digits = 2)),
           legend.args = list(text = expression(paste("BNF Rate (kg N ha"^"-1"," yr"^"-1",")")), 
                              cex=1, side=4, line=3),
           xlab = "longitude", ylab = "latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1))
title(main = substitute(paste("", italic(Pinus~contorta))), cex.main = 1.5)
map('world', regions = 'usa', add = TRUE)

dev.off()

# Sum total from %ndfa ----------------
pi <- merge(pw, regions, 
            by.x = "state", by.y = "State.Abb", 
            all.x = T)

ordering <- c("Interior West","Northern","Pacific Northwest","Southern")
pi$regionsFactor <- factor(pi$FIA.Region, levels = ordering)

fixgr.grid <- array(dim = c(nlon, nlat))
fixsp.ndfa <- data.frame()

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    inds <- which(pi$lon>=lon.a & pi$lon<lon.b & pi$lat>=lat.a & pi$lat<lat.b)
    coords <- data.frame(lon = lon.list[i],
                         lat = lat.list[j])
    
    if(sum(inds) > 0){
      fixgr.grid[i,j] <- mean(pi[inds, ]$ndem, na.rm=TRUE)*raster::extract(FOR_COV_FOR_ALIGN, coords)*100 #convert from km2 to ha
      reg <- pi[inds[1],]$regionsFactor
    }else{
      fixgr.grid[i,j] <- NA
      reg <- NA
    }
    
    out <- c(fixgr.grid[i,j],lon.list[i],lat.list[j],reg)
    fixsp.ndfa <- rbind(fixsp.ndfa, out)
  }
}

totalndfa <- sum(fixgr.grid,na.rm=T)
totalndfa/1e9

colnames(fixsp.ndfa) <- c("fixedn","lon","lat","fiaregion")

saveRDS(fixsp.ndfa, "output/fixsp.ndfa.RDS")

# Figure 1 N-fixer BA ------------------

# Figure 1a: N fixer presence map 
png("output/figure1a.png", width = 12, height = 8, units = 'in', res = 300)
par(bty = "n")
usa <- map_data("usa")
world <- map_data("world", xlim = c(-140, -65), ylim = c(23, 63))

ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
               alpha = 0.01, 
               colour = 'black', 
               fill = NA) + 
  coord_fixed(1.3) +
  geom_point(aes(x = gfdata$lon, y = gfdata$lat, color = as.factor(gfdata$Genus))) +
  geom_jitter() +
  labs(x = " ", y = " ", title = "Nitrogen Fixing Trees in FIA") +
  labs(color = 'Genus') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.text = element_text(face = "italic")) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  scale_colour_brewer(palette = "Paired")  

dev.off()

# Figure 1b: N-fixer basal area
fba.grid <- readRDS("output/fbagrid.RDS")
png("output/figure1b.png", width = 12, height = 8, units = 'in', res = 300)
par(bty = "n",
    cex.lab = 1.5,
    cex.axis = 1.5,
    oma = c(0, 0, 0, 3.5), #for right side of plot in view
    mar = c(5.1, 5, 4.1, 2.1)) #for left side of plot in view
nlev <- 64
ticks <- 2^(-2:13)
image.plot(lon.list, lat.list, log(fba.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), labels = ticks, cex.axis = 2),
           legend.args = list(text = expression(paste("Basal Area (m"^"2"," ha"^"-1",")")), 
                            cex = 2, 
                            side = 4, 
                            line = 5, #moves legend text away from numbers
                            mar = 5),
           xlab = "Longitude", ylab = "Latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140-2, lonmaxint+1),
           zlim = c(log(0.1), log(12)))
title(main = expression(paste("N-Fixing tree basal area per forest area")), 
      cex.main = 2.5)
map('world', regions = 'usa', add = TRUE, col = 'black', lwd = 2)

dev.off()

# Figure 1c: N-fixer basal area per ground area
fbagr.grid <- readRDS("output/fbagrgrid.RDS")
png("output/figure1c.png", width = 12, height = 8, units = 'in', res = 300)
par(bty = "n",
    cex.lab = 1.5,
    cex.axis = 1.5,
    oma = c(0, 0, 0, 3.5), #for right side of plot in view
    mar = c(5.1, 5, 4.1, 2.1)) #for left side of plot in view
nlev <- 64
ticks <- 2^(-2:13)
image.plot(lon.list, lat.list, log(fbagr.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), labels = ticks, cex.axis = 2),
           legend.args = list(text = expression(paste("Basal Area (m"^"2"," ha"^"-1",")")), 
                            cex = 2, 
                            side = 4, 
                            line = 5,
                            mar = 5),
           xlab = "Longitude", ylab = "Latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140-2, lonmaxint+1),
           zlim = c(log(0.1), log(12)))
title(main = expression(paste("N-Fixing tree basal area per ground area")), 
      cex.main = 2.5)
map('world', regions = 'usa', add = TRUE, col = 'black')

dev.off()

# Figure 2, Fixation Rate Maps -----------------------------------------
fn.grid <- readRDS("output/fngrid.RDS") #check correct version was saved
fngr.grid <- readRDS("output/fngrgrid.RDS")
grn.grid <- readRDS("output/grn_grid.RDS")
grngr.grid <- readRDS("output/grngr_grid.RDS")

png("output/figure2_new.png", width = 12, height = 8, units = 'in', res = 300)

par(bty = "n",
    cex.lab = 1.5,
    cex.axis = 1.5,
    oma = c(0, 0, 0, 4),
    mai = c(0.8, 1, 0.8, 1))
set.panel(2, 2)

  #  mar = c(5.1, 5, 4.1, 2.1)) #for left side of plot in view

mtext("Title for Two Plots", outer = TRUE, cex = 1.5)

##map of fixed nitrogen (per forest area) rate average in cell
nlev <- 64
ticks <- 2^(-2:5)
image.plot(lon.list, lat.list, log(fn.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), 
                            labels = ticks,
                            cex.axis = 1.5),
           legend.args = list(text = expression('kg N ha forest'^-1*' yr'^-1), 
                              side = 4, 
                              font = 2, 
                              line = 3.5, 
                              cex = 1,
                              mar = 5),
           xlab = "", ylab = "latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1),
           zlim = c(log(0.1), log(53)))
title(main = "N Fixation Rate \n Per forest area \n Accretion", 
      cex.main = 1.5)
mtext("(a)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE, col = 'black')
#expression('N Fixation Rate (kg N ha forest'^-1*' yr'^-1*'), Accretion')
#zlim was max a log(35)

# fixed N per ground area rate accretion
nlev <- 64
ticks <- 2^(-2:5)
image.plot(lon.list, lat.list, log(fngr.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), 
                            labels = ticks,
                            cex.axis = 1.5),
           legend.args = list(text = expression('kg N ha ground'^-1*' yr'^-1), 
                              side = 4, 
                              font = 2, 
                              line = 3.5, 
                              cex = 1,
                              mar = 5),
           xlab = "", ylab = "",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1),
           zlim = c(log(0.1), log(53)))
title(main = "N Fixation Rate \n Per ground area \n Accretion",
      cex.main = 1.5)
mtext("(b)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions  = 'usa', add = TRUE)
#expression('N Fixation Rate (kg N ha ground'^-1*' yr'^-1*'), Accretion')

##map of N demand based on ndfa (per forest area)
nlev <- 64
ticks <- 2^(-2:5)
image.plot(lon.list, lat.list, log(grn.grid+0.1),
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), labels = round(ticks, digits = 2),
                            cex.axis = 1.5),
           legend.args = list(text = expression('kg N ha forest'^-1*' yr'^-1), 
                              side = 4, 
                              font = 2, 
                              line = 3.5, 
                              cex = 1,
                              mar = 5),
           xlab = "longitude", ylab = "latitude",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1),
           zlim = c(log(0.1), log(53)))
title(main = "N Fixation Rate \n per forest area \n %Ndfa", 
      cex.main = 1.5)
mtext("(c)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)
#expression('N Fixation Rate (kg N ha forest'^-1*' yr'^-1*'), %Ndfa')

##map of N demand based on ndfa
nlev <- 64
ticks <- 2^(-2:5)
image.plot(lon.list, lat.list, log(grngr.grid+0.1), 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = log(ticks), 
                            labels = round(ticks, digits = 2),
                            cex.axis = 1.5),
           legend.args = list(text = expression('kg N ha ground'^-1*' yr'^-1), 
                              side = 4, 
                              font = 2, 
                              line = 3.5, 
                              cex = 1,
                              mar = 5),
           xlab = "longitude", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1),
           zlim = c(log(0.1), log(53)))
title(main = 'N Fixation Rate \n per ground area \n %Ndfa', 
      cex.main = 1.5)
mtext("(d)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)
#title(main = expression('N Fixation Rate (kg N ha ground'^-1*' yr'^-1*'), %Ndfa')

dev.off()

# figure 3 - %BNF from each genus --------------------- 
# Figure 3

## ACCRETION from bootstrap (continent scale)
bootstrap <- readRDS("output/accretion_boot_continent.RDS")
rows <- c("total","Acacia","Albizia","Alnus","Casuarina","Cercocarpus","Elaeagnus",
          "Olneya","Prosopis","Robinia")
means <- c(mean(bootstrap$fixedn), 
           mean(bootstrap$acan/bootstrap$fixedn*100), 
           mean(bootstrap$albn/bootstrap$fixedn*100),
           mean(bootstrap$alnn/bootstrap$fixedn*100),
           mean(bootstrap$casn/bootstrap$fixedn*100),
           mean(bootstrap$cern/bootstrap$fixedn*100),
           mean(bootstrap$elan/bootstrap$fixedn*100),
           mean(bootstrap$olnn/bootstrap$fixedn*100),
           mean(bootstrap$pron/bootstrap$fixedn*100),
           mean(bootstrap$robn/bootstrap$fixedn*100))
sds <- c(sd(bootstrap$fixedn),
         sd(bootstrap$acan/bootstrap$fixedn*100),
         sd(bootstrap$albn/bootstrap$fixedn*100),
         sd(bootstrap$alnn/bootstrap$fixedn*100),
         sd(bootstrap$casn/bootstrap$fixedn*100),
         sd(bootstrap$cern/bootstrap$fixedn*100),
         sd(bootstrap$elan/bootstrap$fixedn*100),
         sd(bootstrap$olnn/bootstrap$fixedn*100),
         sd(bootstrap$pron/bootstrap$fixedn*100),
         sd(bootstrap$robn/bootstrap$fixedn*100))
genus.acc.boot <- dplyr::bind_cols(data.frame(rows),
                                   data.frame(means),
                                   data.frame(sds))
genus.acc.boot <- genus.acc.boot[-1,]
# genus.acc.boot tells you how much was fixed by each genus too
colMeans(bootstrap)/1e9
apply(bootstrap, 2, sd, na.rm = TRUE)/1e9


## LIGHT
#bootstrapped_light <- readRDS("output/bootstrapped_lightlim.RDS")
bootstrapped_light <- readRDS("output/accretion_boot_lightlim_continent.RDS")
totalfix <- mean(bootstrapped_light$fixedn, na.rm = T)
bootlight_pct <- bootstrapped_light[ ,2:ncol(bootstrapped_light)]/bootstrapped_light[ ,1]*100
means <- colMeans(bootlight_pct, na.rm = T)
sds <- apply(bootlight_pct, 2, sd)
genus.light.boot <- dplyr::bind_cols(data.frame(rows[2:10]), 
                                    data.frame(means), 
                                    data.frame(sds))
colnames(genus.light.boot) <- c("rows", "means", "sds")

## NDFA
bootstrapped_ndfa <- readRDS("output/bootstrapped_ndfa.RDS")
totalfix <- mean(bootstrapped_ndfa$tgn, na.rm = T)
means <- colMeans(bootstrapped_ndfa, na.rm = T)
sds <- apply(bootstrapped_ndfa, 2, sd)
genus.ndfa.boot <- dplyr::bind_cols(data.frame(rows), 
                                    data.frame(means), 
                                    data.frame(sds))
genus.ndfa.boot <- genus.ndfa.boot[-1,]
#get fixed n by genus with sd
Ggn_bygen <- (bootstrapped_ndfa[ ,2:ncol(bootstrapped_ndfa)]/100*bootstrapped_ndfa[ ,1])*1e3
colMeans(Ggn_bygen)
apply(Ggn_bygen, 2, sd)

## NDFA upper bound
bootstrapped_upper <- readRDS("output/bootstrapped_upper.RDS")
totalfix <- bootstrapped_upper$tgn
genus.upper.boot <- bootstrapped_upper[-1]
pctupper <- as.vector(unlist(genus.upper.boot))
genus.upper.boot <- dplyr::bind_cols(data.frame(rows[-1]),
                                     data.frame(pctupper))
colnames(genus.upper.boot) <- c("rows","means")


png("output/figure3.png", width = 18, height = 8, units = 'in', res = 300)

plot.acc <- ggplot(genus.acc.boot, aes(x = rows, y = means)) +
  geom_bar(stat="identity", position = position_dodge(0.9), alpha=0.7) +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black") +
  geom_errorbar(aes(ymax = pmin(means+sds, 100), ymin = pmax(means-sds, 0)),
  #geom_errorbar(aes(ymin = means-sds, ymax = means+sds), 
                colour = "black", 
                width = 0.02,
                position = position_dodge(0.9),
                size = 1) +
  geom_hline(yintercept = 0, cex = 0.9) +
  #geom_vline(xintercept=which(genus$genus == 0)) + 
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(subtitle = "Accretion Method", 
      x = "Genus", 
      y = "% of tree-based BNF") +
  ylim(0,100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face =  "italic",
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color =  "#666666", 
                                     face =   "bold", 
                                     size =   32, 
                                     hjust =  0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color =  "#666666", 
                                  face =   "bold", 
                                  size =   22)) 

plot.light <- ggplot(genus.light.boot, aes(x = rows, y = means)) +
  geom_bar(stat="identity", position = position_dodge(0.9), alpha=0.7) +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black") +
  geom_errorbar(aes(ymin = pmax(means-sds,0), ymax = pmin(means+sds, 100)), 
                colour = "black", 
                width = 0.02,
                position = position_dodge(0.9),
                size = 1) +
  geom_hline(yintercept = 0, cex = 0.9) +
  #geom_vline(xintercept=which(genus$genus == 0)) + 
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(subtitle = "Light Limitation", 
       x = "Genus", 
       y = "") +
  ylim(0,100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face =  "italic",
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color =  "#666666", 
                                     face =   "bold", 
                                     size =   32, 
                                     hjust =  0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color =  "#666666", 
                                  face =   "bold", 
                                  size =   22)) 


plot.ndfa <- ggplot(genus.ndfa.boot, aes(x = rows, y = means)) +
  geom_bar(stat="identity", position = position_dodge(0.9), alpha=0.7) +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black") +
  geom_errorbar(aes(ymin = means-sds, ymax = means+sds), 
                colour = "black", 
                width = 0.02,
                position = position_dodge(0.9),
                size = 1) +
  geom_hline(yintercept = 0, cex = 0.9) +
  #geom_vline(xintercept=which(genus$genus == 0)) + 
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(subtitle = "%Ndfa Method", 
       x = "Genus", 
       y = "") +
  ylim(0,100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face =  "italic",
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color =  "#666666", 
                                     face =   "bold", 
                                     size =   32, 
                                     hjust =  0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color =  "#666666", 
                                  face =   "bold", 
                                  size =   22)) 


plot.ndfaup <- ggplot(genus.upper.boot, aes(x = rows, y = means)) +
  # draw the bar plot
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black") +
  geom_hline(yintercept = 0, cex = 0.9) +
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(subtitle = "Upper Bound", x = "Genus", y = "") +
  ylim(0, 100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face="italic",
                                   angle=90, 
                                   hjust=1, 
                                   vjust=0.5),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color="#666666", 
                                     face="bold", 
                                     size=32, 
                                     hjust=0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color="#666666", 
                                  face="bold", 
                                  size=22)) 

grid.arrange(plot.acc, plot.light, plot.ndfa, plot.ndfaup, ncol=4, nrow = 1)
dev.off()

# Total N fixed by each genus ----------------------
#by accretion method
genus.acc %>% mutate(totalfix = (totalacc*genfixpct/100)/1e9)

#by %ndfa method
genus.ndfa %>% mutate(totalfix = (totalndfa*genfixpct/100)/1e9)


# figure 4: BNF by region ------------------------------------
# N by region
regionNames <- c("Interior West","Northern","Pacific Northwest","Southern")

#totalndfa is the total N fixed by the %Ndfa method
#fixsp.ndfa is the N fixed in each grid cell
fixsp.ndfa$fiaregion <- factor(fixsp.ndfa$fiaregion)
fixsp.ndfa <- fixsp.ndfa[!is.na(fixsp.ndfa$fiaregion), ]
regions.ndfa <- fixsp.ndfa %>%
  group_by(fiaregion) %>%
  summarise(fix.region = sum(fixedn, na.rm = TRUE)) %>%
  mutate(fixpct = fix.region/totalndfa*100)
fixsp.ndfa <- cbind(regions.ndfa, regionNames)

fixsp.acc$fiaregion <- factor(fixsp.acc$fiaregion)
fixsp.acc <- fixsp.acc[!is.na(fixsp.acc$fiaregion), ]
regions.acc <- fixsp.acc %>%
  group_by(fiaregion) %>%
  summarise(fix.region = sum(fixedn, na.rm = TRUE)) %>%
  mutate(fixpct = fix.region/totalacc*100)
fixsp.acc <- cbind(regions.acc, regionNames)

plot.reg.ndfa <- ggplot(regions.ndfa, aes(x = regionNames, y = fixpct)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(fixpct, digits = 1)), 
            vjust = -0.3, 
            colour = "grey22") +
  geom_hline(yintercept = 0, cex = 0.9) +
  geom_vline(xintercept = 0.4, cex = 1.5) +
  ylim(0, 70) +
  labs(title = "%Ndfa Method", 
       x = "FIA Region", 
       y = "") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot.reg.acc <- ggplot(regions.acc, aes(x = regionNames, y = fixpct)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(fixpct, digits = 1)), vjust = -0.3, colour = "grey22") +
  geom_hline(yintercept = 0, cex = 0.9) +
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(title = "Accretion Method", 
       x = "FIA Region", 
       y = "% BNF") +
  ylim(0, 70) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 20))

png("output/figure4.png", width = 14, height = 8, units = 'in', res = 300)
grid.arrange(plot.reg.acc, plot.reg.ndfa, 
             ncol=2, nrow = 1, 
             top = textGrob("Percent of total BNF in FIA Region",
                            gp = gpar(fontsize = 25)))
dev.off()
