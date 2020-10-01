

#functions

#newest version of plotfixndfa with:
# Perakis constants, new accretion data
#plotdata should be agb_1_2
plotfixndfa <- function(plotdata, 
                        res = 0.4, 
                        resr = 0.73, 
                        cn.foliage = 35.1, 
                        cn.wood = 350, 
                        cn.fr = 43, 
                        ndfatype = "realistic"){
  
  res <- res #leaf resorption 0.4 for normal, 1 for upper bound, res is amount lost here
  resr <- resr #root resporption 0.73 for normal, 1 for upper bound, res is amount lost here
  ## Get biomass allocations
  
  if(ndfatype=="realistic"){
    ndfas <- list(aca = rnorm(1, aca.mean, aca.se)/100, 
                  alb = rnorm(1, alb.mean, alb.se)/100, 
                  ald = rnorm(1, ald.mean, ald.se)/100, 
                  pro = rnorm(1, pro.mean, pro.se)/100,
                  ela = rnorm(1, ela.mean, ela.se)/100, 
                  cas = rnorm(1, cas.mean, cas.se)/100, 
                  rob = rnorm(1, rob.mean, rob.se)/100, 
                  other = rnorm(1, other.mean, other.se)/100)
    lightcut <- 6
    
  }else if(ndfatype=="upperbound"){
    ndfas <- list(aca = 1, 
                  alb = 1, 
                  ald = 1, 
                  pro = 1,
                  ela = 1, 
                  cas = 1, 
                  rob = 1, 
                  other = 1)
    lightcut <- 6
    
  }else if(ndfatype=="lightlim"){
    ndfas <- list(aca = rnorm(1, aca.mean, aca.se)/100, 
                  alb = rnorm(1, alb.mean, alb.se)/100, 
                  ald = rnorm(1, ald.mean, ald.se)/100, 
                  pro = rnorm(1, pro.mean, pro.se)/100,
                  ela = rnorm(1, ela.mean, ela.se)/100, 
                  cas = rnorm(1, cas.mean, cas.se)/100, 
                  rob = rnorm(1, rob.mean, rob.se)/100, 
                  other = rnorm(1, other.mean, other.se)/100)
    lightcut <- 4
  }
  
  # hardwood parameters (Jenkins, 2003)
  B0.foliage <- -4.0813
  B1.foliage <- 5.8816
  B0.roots <- -1.6911
  B1.roots <- 0.8160
  B0.bark <- -2.0129
  B1.bark <- -1.6805
  B0.wood <- -0.3065
  B1.wood <- -5.4240
  
  #get change in biomass (kgC in foliage and wood)
  plot_agb <- plotdata %>% 
    mutate(agb.foliage = exp(B0.foliage+(B1.foliage/dbh2pred))*agb2pred,
           agb.wood = (agbd),
           bgb.root = exp(B0.roots+(B1.roots/dbh2pred))*agb2pred)
  plot_agb <- plot_agb[plot_agb$tpha > 0.5,]
  plot_agb <- plot_agb[plot_agb$dbh < 200,] #weird 2 outliers driving calcs
  plot_agb <- plot_agb[!is.na(plot_agb$dbh2pred),]
  
  # filter for light limitation
  plot_agb <- filter(plot_agb, cclcd <= lightcut)
  
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
  
  output <- list(tgn=fixedn)
  return(output)
}




# newest version of plotfixacc
plotfixacc <- function(plotdata){
  fixedn <- acan <- albn <- alnn <- cern <- elan <- olnn <- pron <- robn <- NULL
  
  #define BNF slope and int for whole continent
  ac <- ac.samp[sample(1:n, 1, replace=T),]
  ot <- ot.samp[sample(1:n, 1, replace=T),]
  al <- al.samp[sample(1:n, 1, replace=T)]
  pr <- pr.samp[sample(1:n, 1, replace=T),]
  #ro <- rnorm(1, mean(robinia$rate), var(robinia$rate)^0.5)
  ro <- ro.samp[sample(1:n, 1, replace=T)]
  
  plot_sens <- stems0 %>%
    #filter(FIX == 1) %>%
    filter(!Genus %in% c("Casuarina","Sophora","Piscidia")) %>%
    group_by(pcn, Genus) %>%
    mutate(fix = ifelse(Genus %in% c("Acacia"), sum(ac.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                        ifelse(Genus %in% c("Albizia"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                               ifelse(Genus %in% c("Alnus"), al.samp[sample(1:n, 1, replace=T)],
                                      ifelse(Genus %in% c("Cercocarpus"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                             ifelse(Genus %in% c("Prosopis"), sum(pr.samp[sample(1:n, 1, replace=T),]*c(1,pBA)), 
                                                    ifelse(Genus %in% c("Ebenopsis"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                           ifelse(Genus %in% c("Robinia"),rnorm(1,mean(robinia$rate),var(robinia$rate)^0.5),
                                                                  ifelse(Genus %in% c("Olneya"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                         ifelse(Genus %in% c("Elaeagnus"), sum(ot.samp[sample(1:n, 1, replace=T),]*c(1,pBA)),
                                                                                0)))))))))) 
  
  
  #add pb onto plot_sens to get expcurr
  plot_sens$pcn <- as.integer64(plot_sens$pcn)
  ps <- merge(plot_sens, psnap,
              by = "pcn",
              all.x = T, all.y = F)
  
  ps$fixedn <- ps$fix*ps$EXPCURR*haperacre #multiply N fixed in plot by area plot represents and convert to ha
  ps2 <- ps[!is.na(ps$Genus),]
  ps2[is.na(ps2$fixedn),]$fixedn <- 0
  
  anb <- as.data.frame.table(by(ps2$fixedn, ps2$Genus, sum)) #if not working do e$fix
  anb <- data.frame(Genus = anb[ ,1], 
                    fixboot = as.numeric(anb[ ,2]))
  
  bygenus <- anb[!is.na(anb$fixboot),]
  
  fixedn <- sum(bygenus$fixboot, na.rm = T)/1e9 #total N fixed across US from accretion method (in TgN/yr)
  acan <- bygenus[bygenus$Genus %in% c("Acacia"),"fixboot"]/1e9
  albn <- bygenus[bygenus$Genus %in% c("Albizia"),"fixboot"]/1e9
  alnn <- bygenus[bygenus$Genus %in% c("Alnus"),"fixboot"]/1e9
  cern <- bygenus[bygenus$Genus %in% c("Cercocarpus"),"fixboot"]/1e9
  elan <- bygenus[bygenus$Genus %in% c("Elaeagnus"),"fixboot"]/1e9
  olnn <- bygenus[bygenus$Genus %in% c("Olneya"),"fixboot"]/1e9
  pron <- bygenus[bygenus$Genus %in% c("Prosopis"),"fixboot"]/1e9
  robn <- bygenus[bygenus$Genus %in% c("Robinia"),"fixboot"]/1e9
  
  # output <- list(tgn=fixedn, acan=acan, albn=albn, alnn=alnn, 
  #                cern=cern, elan=elan, olnn=olnn, pron=pron, robn=robn)
  output <- cbind(fixedn, acan, albn, alnn, cern, elan, olnn, pron, robn)
  return(output)
  
}


# plotdata should be fed allstems, 
# r should be fed output from getrates function, 
# canpos = 6 for all trees or 4 for canopy only
plotfix_new <- function(plotdata, canpos, r) {
  
  allstemsf <- plotdata[plotdata$cclcd < canpos, ] 
  print("filtered out understory")
  
  stems00 <- allstemsf %>%
    group_by(pcn, Genus) %>%
    summarise(pBA = sum(BAm2ha),
              FIX = first(FIX),
              ACT1vsRIZ0 = first(ACT1vsRIZ0))
  print("agg to plot level")
  
  n <- 10000
  # get fixation rates for plots
  plot_sens <- stems00 %>% 
    group_by(pcn, Genus) %>%
    mutate(fix = ifelse(Genus %in% c("Acacia"), sum(r$ac.samp[sample(1:n, 1),]*c(1,pBA)), 
                        ifelse(Genus %in% c("Albizia"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                               ifelse(Genus %in% c("Alnus"), r$al.samp[sample(1:n, 1)],
                                      ifelse(Genus %in% c("Cercocarpus"), sum(r$pr.samp[sample(1:n, 1),]*c(1,pBA)),
                                             ifelse(Genus %in% c("Prosopis"), sum(r$pr.samp[sample(1:n, 1),]*c(1,pBA)), 
                                                    ifelse(Genus %in% c("Ebenopsis"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                           ifelse(Genus %in% c("Robinia"),r$ro.samp[sample(1:n, 1)],
                                                                  ifelse(Genus %in% c("Olneya"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                                         ifelse(Genus %in% c("Elaeagnus"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                                                0)))))))))) 
  print("assigned fixation rates")
  
  
  
  an <- as.data.frame.table(by(plot_sens$fix, plot_sens$pcn, sum)) 
  an <- data.frame(pcn = as.numeric(as.character(an[ ,1])), 
                   fixr = as.numeric(an[ ,2]))
  print("aggregated plots")
  
  p <- merge(p, an, 
             by = "pcn", 
             all.x = TRUE, all.y = FALSE)
  
  bygenus <- plot_sens %>%
    group_by(Genus) %>%
    #summarise_all(funs(sum(., na.rm=T)))
    summarise(fix = sum(fix, na.rm = T))
  bygenus <- bygenus[c(1:5,7,8,10:11),] #remove Ebenopsis (too few stems), Piscidia (1 stem), Sophora (1 stem), NAs
  
  fix.grid <- fixfor.grid <- fngr.grid <- array(dim = c(nlon, nlat))
  
  for(i in 1:nlon){
    print(i/nlon*100)
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(p$lon>=lon.a & p$lon<lon.b & p$lat>=lat.a & p$lat<lat.b)
      
      if(sum(inds) > 0){
        fixfor.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)
        fix.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)*foresta.grid[i,j] #fix rate times forest
        fngr.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)*fracfor.grid[i,j] #fix rate times fraction of ground forested
      }else{
        fixfor.grid[i,j] <- NA
        fix.grid[i,j] <- NA
        fngr.grid[i,j] <- NA
      }
      
    }
  }
  
  fixedn <- sum(fix.grid, na.rm = T) #total N fixed across US from accretion method
  totalk <- sum(bygenus[ ,"fix"])
  acan <- (bygenus[[1,"fix"]]/totalk)*100 #*fixedn to get kgN instead of %
  albn <- (bygenus[[2,"fix"]]/totalk)*100
  alnn <- (bygenus[[3,"fix"]]/totalk)*100
  #casn <- (bygenus[[4,"fix"]]/totalk)*100
  cern <- (bygenus[[5,"fix"]]/totalk)*100
  elan <- (bygenus[[6,"fix"]]/totalk)*100
  olnn <- (bygenus[[7,"fix"]]/totalk)*100
  pron <- (bygenus[[8,"fix"]]/totalk)*100
  robn <- (bygenus[[9,"fix"]]/totalk)*100
  
  #totalacc <- sum(fix.grid, na.rm = T)/1e9 #total N fixed across US from accretion method 
  fixfor <- mean(fixfor.grid, na.rm=T)
  fixgr <- mean(fngr.grid, na.rm=T)
  
  output <- list(tgn=fixedn/1e9, acapct=acan, albpct=albn, alnpct=alnn,
                 cerpct=cern, ealpct=elan, olnpct=olnn, propct=pron, robpct=robn,
                 fixperfor=fixfor, fixpergr=fixgr)
  return(output)
  
}



#plot_sens step only
plotsensonly <- function(plotdata, canpos, r, focalgenus) {
  
  allstemsf <- plotdata[plotdata$cclcd < canpos, ] 
  print("filtered out understory")
  
  stems00 <- allstemsf %>%
    group_by(pcn, Genus) %>%
    summarise(pBA = sum(BAm2ha))
  print("agg to plot level")
  
  n <- 10000
  # get fixation rates for plots
  plot_sens <- stems00 %>% 
    group_by(pcn, Genus) %>%
    mutate(fix = ifelse(Genus %in% c("Acacia"), sum(r$ac.samp[sample(1:n, 1),]*c(1,pBA)), 
                        ifelse(Genus %in% c("Albizia"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                               ifelse(Genus %in% c("Alnus"), r$al.samp[sample(1:n, 1)],
                                      ifelse(Genus %in% c("Cercocarpus"), sum(r$pr.samp[sample(1:n, 1),]*c(1,pBA)),
                                             ifelse(Genus %in% c("Prosopis"), sum(r$pr.samp[sample(1:n, 1),]*c(1,pBA)), 
                                                    ifelse(Genus %in% c("Ebenopsis"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                           ifelse(Genus %in% c("Robinia"),r$ro.samp[sample(1:n, 1)],
                                                                  ifelse(Genus %in% c("Olneya"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                                         ifelse(Genus %in% c("Elaeagnus"), sum(r$ot.samp[sample(1:n, 1),]*c(1,pBA)),
                                                                                0)))))))))) 
  print("assigned fixation rates")
  plot_sens2 <- plot_sens %>% 
    mutate(fix = replace(fix, Genus != focalgenus, 0))
  
  an <- as.data.frame.table(by(plot_sens2$fix, plot_sens2$pcn, sum)) 
  an <- data.frame(pcn = as.numeric(as.character(an[ ,1])), 
                   fixr = as.numeric(an[ ,2]))
  an$pcn <- as.integer64(an$pcn)
  print("aggregated plots")
  
  p$pcn <- as.integer64(p$pcn)
  p <- merge(p, an, 
             by = "pcn", 
             all.x = TRUE, all.y = FALSE)
  
  fix.grid <- fixfor.grid <- fngr.grid <- array(dim = c(nlon, nlat))
  
  for(i in 1:nlon){
    print(i/nlon*100)
    lon.a <- lon.list[i] - dlon/2
    lon.b <- lon.list[i] + dlon/2
    
    for(j in 1:nlat){
      lat.a <- lat.list[j] - dlat/2
      lat.b <- lat.list[j] + dlat/2
      inds <- which(p$lon>=lon.a & p$lon<lon.b & p$lat>=lat.a & p$lat<lat.b)
      
      if(sum(inds) > 0){
        fixfor.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)
        fix.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)*foresta.grid[i,j] #fix rate times forest
        fngr.grid[i,j] <- mean(p[inds, ]$fixr, na.rm=TRUE)*fracfor.grid[i,j] #fix rate times fraction of ground forested
      }else{
        fixfor.grid[i,j] <- NA
        fix.grid[i,j] <- NA
        fngr.grid[i,j] <- NA
      }
      
    }
  }
  
  fixedn <- sum(fix.grid, na.rm = T) #total N fixed across US from accretion method
 
  output <- list(fixfor.grid = fixfor.grid)
  return(fixfor.grid)
  
}


#genus should be in lowercase inside quotes
#type can be max or min, inside quotes
#for all genera at normal values use genus="standard" and leave type blank
getrates <- function(genus, type){
  n = 10000
  #standard rates
  # al.samp <- mvrnorm(n, coef(al)[1:2], vcov(al))
  # ot.samp <- mvrnorm(n, coef(ot)[1:2], vcov(ot))
  ac.samp <- mvrnorm(n,coef(ac)[1:2],vcov(ac))
  pr.samp <- mvrnorm(n,coef(pr)[1:2],vcov(pr))
  ro.samp <- rnorm(n, ro.mu, ro.se) #gets random distribution of constant robinia rate
  al.samp <- rnorm(n, al.mu, al.se)
  
  BA <- rnorm(n, other.mu, other.se) #slope
  Intercept <- rnorm(n, other.int, other.intse) #intercept
  ot.samp <- cbind(Intercept, BA)
  
  #modified rates
  if(genus=="acacia"){
    if(type=="max"){
      ac.samp <- mvrnorm(n, coef(ac)[1:2]*c(1.5,1.5), vcov(ac)) #slope 150%, int same
    }else if(type=="min"){
      ac.samp <- mvrnorm(n, coef(ac)[1:2]*c(0.5,0.5), vcov(ac)) #slope 50%, int same
    }
  }else if(genus=="albizia" | genus=="ebenopsis" | genus=="olneya" | genus=="elaeagnus" | genus=="cercocarpus" | genus=="other"){
    if(type=="max"){
      #ot.samp <- mvrnorm(n, coef(ot)[1:2]*c(1,1.5), vcov(ot))
      BA <- rnorm(n, other.mu*1.5, other.se) #slope
      Intercept <- rnorm(n, other.int, other.intse) #intercept
      ot.samp <- cbind(Intercept, BA)
    }else if(type=="min"){
      #ot.samp <- mvrnorm(n, coef(ot)[1:2]*c(1,0.5), vcov(ot))
      BA <- rnorm(n, other.mu*0.5, other.se) #slope
      Intercept <- rnorm(n, other.int, other.intse) #intercept
      ot.samp <- cbind(Intercept, BA)
    }
  }else if(genus=="alnus"){
    if(type=="max"){
      #al.samp <- mvrnorm(n, coef(al)[1:2]*c(1,1.5), vcov(al))
      al.samp <- rnorm(n, al.mu*1.5, al.se)
    }else if(type=="min"){
      #al.samp <- mvrnorm(n, coef(al)[1:2]*c(1,0.5), vcov(al))
      al.samp <- rnorm(n, al.mu*0.5, al.se)
    }
  }else if(genus=="prosopis"){
    if(type=="max"){
      pr.samp <- mvrnorm(n, coef(pr)[1:2]*c(1.5,1.5), vcov(pr))
    }else if(type=="min"){
      pr.samp <- mvrnorm(n, coef(pr)[1:2]*c(0.5,0.5), vcov(pr))
    }
  }else if(genus=="robinia"){
    if(type=="max"){
      ro.samp <- rnorm(n, ro.mu*1.5, ro.se)
    }else if(type=="min"){
      ro.samp <- rnorm(n, ro.mu*0.5, ro.se)
    }
  }else if(genus=="standard"){
    # ac.samp <- mvrnorm(n, coef(ac)[1:2], vcov(ac))
    # pr.samp <- mvrnorm(n, coef(pr)[1:2], vcov(pr))
    # al.samp <- mvrnorm(n, coef(al)[1:2], vcov(al))
    # ot.samp <- mvrnorm(n, coef(ot)[1:2], vcov(ot))
    ac.samp <- mvrnorm(n,coef(ac)[1:2],vcov(ac))
    pr.samp <- mvrnorm(n,coef(pr)[1:2],vcov(pr))
    ro.samp <- rnorm(n, ro.mu, ro.se) #gets random distribution of constant robinia rate
    al.samp <- rnorm(n, al.mu, al.se)
    
    BA <- rnorm(n, other.mu, other.se) #slope
    Intercept <- rnorm(n, other.int, other.intse) #intercept
    ot.samp <- cbind(Intercept, BA)
  }
  r <- list(ac.samp = ac.samp,
            pr.samp = pr.samp,
            al.samp = al.samp,
            ot.samp = ot.samp,
            ro.samp = ro.samp)
  return(r)
}

