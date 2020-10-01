

###############################
#
# %ndfa from nearby plots GR
#
# Anika Petach
# 9/12/18
#
################################

library(fields)
library(dplyr)
library(ggplot2)
require(data.table)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

#----------- Load processed data
# plot data
p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv",header=T)
p <- p[p$state!='PR'&p$state!='VI'&p$state!='HI', ]
p[p$lat==-999,]$lat <- NA
p[p$lon==-999,]$lon <- NA

# species data
speciesdata <- read.csv("ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")

# census data
grdata <- readRDS("output/gdataGR1.RDS") # from time series data
allstems <- readRDS("output/allstems.RDS") # from all plots

# merge to get lat/lon
grdata.loc <- merge(grdata, p, by="pcn", all.x=T)

# merge to get species data
allstems.allo <- merge(allstems, 
                       speciesdata[,c("SPCD","JENKINS_TOTAL_B1","JENKINS_TOTAL_B2")], 
                       by.x="spcd", by.y="SPCD", 
                       all.x=T)
allstems.fix <- allstems.allo %>% dplyr::filter(FIX==1, !Genus %in% c("Ebenopsis","Piscidia","Sophora"))

# change factor levels for Genus
grdata.loc$Genus <- factor(grdata.loc$Genus)

#------------ Set up lat lon lists
# Get latitude and longitude bounds of AK, coterminous US, Mexico, PR, VI
latmin <- min(p$lat,na.rm=TRUE)
latmax <- max(p$lat,na.rm=TRUE)
lonmin <- min(p$lon,na.rm=TRUE)
lonmax <- max(p$lon,na.rm=TRUE)
latminint <- floor(latmin)
latmaxint <- ceiling(latmax)
lonminint <- floor(lonmin)
lonmaxint <- ceiling(lonmax)
dlat <- 10 # grid cell size
dlon <- 10 # grid cell size
lat.list <- seq(latminint,latmaxint,dlat)
lon.list <- seq(lonminint,lonmaxint,dlon)
nlat <- length(lat.list)
nlon <- length(lon.list)
fiascaling <- 2428.114

#------------Build RGR vs dbh1 distributions (genus only)
# graphs of all plots data pooled by fixer genus
grdata %>% filter(FIX==1) %>% ggplot(aes(x=dbh1, y=LGRp))+geom_smooth(method="loess")+facet_wrap(~Genus,ncol=3)


#### Break down regression by genus only
grdata.loc$Genus <- factor(grdata.loc$Genus) #get rid of unused factor levels
temp <- split(grdata.loc,grdata.loc$Genus)

model.prosopis <- model.acacia <- loess(temp$Prosopis$LGRp~temp$Prosopis$dbh1,control=loess.control(surface="direct"))
model.albizia <- loess(temp$Albizia$LGRp~temp$Albizia$dbh1,control=loess.control(surface="direct"))
model.robinia <- model.olneya <- loess(temp$Robinia$LGRp~temp$Robinia$dbh1,control=loess.control(surface="direct"))
model.alnus <- loess(temp$Alnus$LGRp~temp$Alnus$dbh1,control=loess.control(surface="direct"))
model.cercocarpus <- model.casuarina <- loess(temp$Cercocarpus$LGRp~temp$Cercocarpus$dbh1,control=loess.control(surface="direct"))
model.elaeagnus <- loess(temp$Elaeagnus$LGRp~temp$Elaeagnus$dbh1,control=loess.control(surface="direct"))

gens <- as.vector(unique(allstems.fix$Genus))
gens <- c("Acacia","Albizia","Alnus","Cercocarpus","Prosopis","Casuarina",
          "Olneya","Elaeagnus") #removed Robinia because too large to run in loop
models <- list(model.acacia,model.albizia,model.alnus,model.cercocarpus,
               model.prosopis,model.casuarina,
               model.olneya,model.elaeagnus)
#models order: prosopis (acacia), albizia, alnus, cercocarpus, prosopis, cercocarpus, robinia, elaeagnus

# dplyr implementation
# requires too much memory -- go to data.table implementation
out <- NULL
for(i in 1:length(gens)){
  outn <- allstems.allo %>% 
    dplyr::filter(Genus %in% c(gens[i])) %>%
    mutate(predfit=predict(models[[i]],dbh,se=T)$fit,
           predse=predict(models[[i]],dbh,se=T)$se.fit) %>%
    mutate(RGRpred=rnorm(1,mean=predfit,sd=predse)) %>%
    mutate(dbh2pred=exp(RGRpred+log(dbh))) %>%
    mutate(agb2pred=JENKINS_TOTAL_B1+JENKINS_TOTAL_B2*log(dbh2pred))
  out <- rbind(out,outn)
}

#data.table implementation
s1 <- data.table(allstems.allo) #make data table
agb_1_2 <- NULL
for(i in 1:length(gens)){
  #s2 <- s1[Genus %in% c(gens[i])] #subset rows (gens[i])
  s0 <- s1[Genus %in% c("Robinia")]
  s2 <- s0[45001:nrow(s0),]
  s2[, predfit := predict(model.robinia,dbh,se=T)$fit] #models[[i]]
  s2[, predse := predict(model.robinia,dbh,se=T)$se.fit] #models[[i]]
  for(k in 1:nrow(s2)){
    s2$RGRpred[k] <- rnorm(1,mean=s2$predfit[k],sd=s2$predse[k])
  }
  #s2[, RGRpred := rnorm(1,mean=predfit,sd=predse)]
  s2[, dbh2pred:=exp(RGRpred+log(dbh))]
  s2[, agb2pred:=exp(JENKINS_TOTAL_B1+JENKINS_TOTAL_B2*log(dbh2pred))]
  agb_1_2 <- rbind(agb_1_2, s2)
}
## Robinia is too big to do predict in one step
## Manually replace gens[i] with "Robinia"
## Manually replace models[[i]] with model.robinia in 2 lines
## Subset s2 <- s0[1:10000,] then s0[10001:20000,] then s0[20001:30000,] 
## then s0[30001:40000,] then s0[40001:45000,] then s0[45001:nrow(s0),] on line 108

#saveRDS(agb_1_2,"output/agb_1_2.RDS")

# Check whether predictions are in model range ----------------
par(mfrow=c(3,3))
#Acacia
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Acacia"))
plot(model.acacia, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Acacia (from Prosopis")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.acacia)[[4]]),font=2,adj=0.01,padj=-0.5)
#Albizia
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Albizia"))
plot(model.albizia, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Albizia")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.albizia)[[4]]),font=2,adj=0.01,padj=-0.5)
#Alnus
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Alnus"))
plot(model.alnus, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Alnus")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.alnus)[[4]]),font=2,adj=0.01,padj=-0.5)
#Cercocarpus
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Cercocarpus"))
plot(model.cercocarpus, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Cercocarpus")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.cercocarpus)[[4]]),font=2,adj=0.01,padj=-0.5)
#Prosopis
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Prosopis"))
plot(model.prosopis, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Prosopis")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.prosopis)[[4]]),font=2,adj=0.01,padj=-0.5)
#Robinia
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Robinia"))
plot(model.robinia, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Robinia")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.robinia)[[4]]),font=2,adj=0.01,padj=-0.5)
#Casuarina
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Casuarina"))
plot(model.cercocarpus, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Casuarina (from Cercocarpus)")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.casuarina)[[4]]),font=2,adj=0.01,padj=-0.5)
#Olneya
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Olneya"))
plot(model.robinia, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Olneya (from Robinia)")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.olneya)[[4]]),font=2,adj=0.01,padj=-0.5)
#Elaeagnus
temp1 <- allstems.allo %>% dplyr::filter(Genus %in% c("Elaeagnus"))
plot(model.elaeagnus, xlim=c(min(temp1$dbh),max(temp1$dbh)), main="Elaeagnus")
lines(density(temp1$dbh), col='red')
mtext(paste("Equivalent parameters ", summary(model.elaeagnus)[[4]]),font=2,adj=0.01,padj=-0.5)

#DBH before and after
ggplot() + 
  geom_smooth(data=filter(grdata.loc,FIX==1),aes(x=dbh1, y=dbh2), color='black') + #color=Genus
  geom_smooth(data=agb_bgb,aes(x=dbh, y=dbh2pred), shape=21, color='red') + #fill=Genus
  scale_color_discrete("Train model") +
  scale_fill_discrete("Predicted") +
  geom_abline(intercept=0, slope=1, col='blue') +
  ggtitle("DBH1 and DBH2 for training data (black) and predictions (red)")

#--------------------- Build RGR distributions (genus and locaion), didn't use
#### Break down regressions by region AND genus
grdata.loc$RGRpred <- NA # to store predicted output
for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    
    # develop loess model
    temp <- grdata.loc %>% dplyr::filter(FIX==1 & lon>lon.a & lon<lon.b & lat>lat.a & lat<lat.a)
    templ <- nrow(temp)
    temp$Genus <- factor(temp$Genus) #get rid of unused factor levels
    temp <- split(temp, temp$Genus) #make sep df for each genus
    
    if(is.null(temp$Prosopis)==FALSE){
      model.prosopis <- model.acacia <- loess(temp$Prosopis$LGRp~temp$Prosopis$dbh1,control=loess.control(surface="direct"))
      
      # predict values for RGR based on dbh1
      inds.pro <- which(grdata.loc$Genus %in% c("Prosopis","Acacia") & grdata.loc$FIX==1 & grdata.loc$lon>=lon.a & grdata.loc$lon<lon.b & grdata.loc$lat>=lat.a & grdata.loc$lat<lat.b) # get rows with Nfixer inside lat lon bounds
      if(length(inds.pro)>0){
        for(k in 1:length(inds.pro)){
          pred <- predict(model.prosopis, grdata.loc[inds.pro[k],]$dbh1,se=T)
          grdata.loc[inds.pro[k],]$RGRpred <- rnorm(1, mean=pred$fit, sd=pred$se.fit)
        }
      }
    }
    if(is.null(temp$Albizia)==FALSE){
      model.albizia <- loess(temp$Albizia$LGRp~temp$Albizia$dbh1,control=loess.control(surface="direct"))
      inds.alb <- which(grdata.loc$Genus %in% c("Albizia") & grdata.loc$FIX==1 & grdata.loc$lon>=lon.a & grdata.loc$lon<lon.b & grdata.loc$lat>=lat.a & grdata.loc$lat<lat.b) # get rows with Nfixer inside lat lon bounds
      if(length(inds.alb)>0){
        for(k in 1:length(inds.alb)){
          pred <- predict(model.albizia, grdata.loc[inds.alb[k],]$dbh1,se=T)
          grdata.loc[inds.alb[k],]$RGRpred <- rnorm(1, mean=pred$fit, sd=pred$se.fit)
        }
      }
    }
    if(is.null(temp$Robinia)==FALSE){
      model.robinia <- model.olneya <- loess(temp$Robinia$LGRp~temp$Robinia$dbh1,control=loess.control(surface="direct"))
      inds.rob <- which(grdata.loc$Genus %in% c("Robinia","Olneya") & grdata.loc$FIX==1 & grdata.loc$lon>=lon.a & grdata.loc$lon<lon.b & grdata.loc$lat>=lat.a & grdata.loc$lat<lat.b) # get rows with Nfixer inside lat lon bounds
      if(length(inds.rob)>0){
        for(k in 1:length(inds.rob)){
          pred <- predict(model.robinia, grdata.loc[inds.rob[k],]$dbh1,se=T)
          grdata.loc[inds.rob[k],]$RGRpred <- rnorm(1, mean=pred$fit, sd=pred$se.fit)
        }
      }
    }
    if(is.null(temp$Alnus)==FALSE){
      model.alnus <- loess(temp$Alnus$LGRp~temp$Alnus$dbh1,control=loess.control(surface="direct"))
      inds.aln <- which(grdata.loc$Genus %in% c("Alnus") & grdata.loc$FIX==1 & grdata.loc$lon>=lon.a & grdata.loc$lon<lon.b & grdata.loc$lat>=lat.a & grdata.loc$lat<lat.b) # get rows with Nfixer inside lat lon bounds
      if(length(inds.aln)>0){
        for(k in 1:length(inds.aln)){
          pred <- predict(model.alnus, grdata.loc[inds.aln[k],]$dbh1,se=T)
          grdata.loc[inds.aln[k],]$RGRpred <- rnorm(1, mean=pred$fit, sd=pred$se.fit)
        }
      }
    }
    if(is.null(temp$Casuarina)==FALSE){
      model.casuarina <- loess(temp$Casuarina$LGRp~temp$Casuarina$dbh1,control=loess.control(surface="direct"))
      inds.cas <- which(grdata.loc$Genus %in% c("Casuarina") & grdata.loc$FIX==1 & grdata.loc$lon>=lon.a & grdata.loc$lon<lon.b & grdata.loc$lat>=lat.a & grdata.loc$lat<lat.b) # get rows with Nfixer inside lat lon bounds
      if(length(inds.cas)>0){
        for(k in 1:length(inds.cas)){
          pred <- predict(model.casuarina, grdata.loc[inds.cas[k],]$dbh1,se=T)
          grdata.loc[inds.cas[k],]$RGRpred <- rnorm(1, mean=pred$fit, sd=pred$se.fit)
        }
      }
    }
    if(is.null(temp$Cercocarpus)==FALSE){
      model.cercocarpus <- loess(temp$Cercocarpus$LGRp~temp$Cercocarpus$dbh1,control=loess.control(surface="direct"))
      inds.cer <- which(grdata.loc$Genus %in% c("Cercocarpus") & grdata.loc$FIX==1 & grdata.loc$lon>=lon.a & grdata.loc$lon<lon.b & grdata.loc$lat>=lat.a & grdata.loc$lat<lat.b) # get rows with Nfixer inside lat lon bounds
      if(length(inds.cer)>0){
        for(k in 1:length(inds.cer)){
          pred <- predict(model.cercocarpus, grdata.loc[inds.cer[k],]$dbh1,se=T)
          grdata.loc[inds.cer[k],]$RGRpred <- rnorm(1, mean=pred$fit, sd=pred$se.fit)
        }
      }
    }
    if(is.null(temp$Elaeagnus)==FALSE){
      model.elaeagnus <- loess(temp$Elaeagnus$LGRp~temp$Elaeagnus$dbh1,control=loess.control(surface="direct"))
      inds.ela <- which(grdata.loc$Genus %in% c("Elaeagnus") & grdata.loc$FIX==1 & grdata.loc$lon>=lon.a & grdata.loc$lon<lon.b & grdata.loc$lat>=lat.a & grdata.loc$lat<lat.b) # get rows with Nfixer inside lat lon bounds
      if(length(inds.ela)>0){
        for(k in 1:length(inds.ela)){
          pred <- predict(model.elaeagnus, grdata.loc[inds.ela[k],]$dbh1,se=T)
          grdata.loc[inds.ela[k],]$RGRpred <- rnorm(1, mean=pred$fit, sd=pred$se.fit)
        }
      }
    }
  }
}

# RGR by genus ----------------------

agb_df <- as.data.frame(agb_1_2)
agb_df %>% group_by(Genus) %>% hist(RGRpred)

# RGR
ggplot(agb_df,aes(x=RGRpred, y=..density..))+
  geom_histogram()+
  facet_wrap(~Genus, ncol=3)+
  theme_bw() +
  xlim(c(-0.1,0.2)) +
  xlab("Relative Growth Rate") +
  geom_vline(xintercept=0)

# 
ggplot(agb_df,aes(x=dbh2pred, y=..density..))+
  geom_histogram()+
  facet_wrap(~Genus, ncol=3)+
  theme_bw() +
  xlab("Predicted DBH 2") +
  geom_vline(xintercept=0)

