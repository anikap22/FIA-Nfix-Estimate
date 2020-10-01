

# 
# bootstrap %NDFA values to get species specific 
# by Anika Petach 
# 4/20/18

require(tidyr)
require(dplyr)

############ Step1: Load cleaned data

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")
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
aca.mean <- mean(aca$Value, na.rm=T)
aca.se <- sd(aca$Value, na.rm=T)/nrow(aca)
alb.mean <- mean(alb$Value, na.rm=T)
alb.se <- sd(alb$Value, na.rm=T)/nrow(alb)
ald.mean <- mean(ald$Value, na.rm=T)
ald.se <- sd(ald$Value, na.rm=T)/nrow(ald)
pro.mean <- mean(pro$Value, na.rm=T)
pro.se <- sd(pro$Value, na.rm=T)/nrow(pro)
ela.mean <- mean(ela$Value, na.rm=T)
ela.se <- sd(ela$Value, na.rm=T)/nrow(ela)
cas.mean <- mean(cas$Value, na.rm=T)
cas.se <- sd(cas$Value, na.rm=T)/nrow(cas)
rob.mean <- mean(rob$Value, na.rm=T)
rob.se <- sd(rob$Value, na.rm=T)/nrow(rob)
other.mean <- mean(ndfa$Value, na.rm=T)
other.se <- sd(ndfa$Value, na.rm=T)/nrow(ndfa)

# make gamma distributions
n = 10000
aca.samp <- rnorm(n,aca.mean, aca.se)
alb.samp <- rnorm(n,alb.mean, alb.se)
ald.samp <- rnorm(n,ald.mean, ald.se)
pro.samp <- rnorm(n,pro.mean, pro.se)
ela.samp <- rnorm(n,ela.mean, ela.se)
cas.samp <- rnorm(n,cas.mean, cas.se)
rob.samp <- rnorm(n,rob.mean, rob.se)
other.samp <- rnorm(n,other.mean, other.se)

# get mean and 95% CI for each species
mean(aca.samp)
mean(alb.samp)
mean(ald.samp)
mean(pro.samp)
mean(ela.samp)
mean(cas.samp)
mean(other.samp)

ndfas <- list(aca=mean(aca.samp)/100, alb=mean(alb.samp)/100, ald=mean(ald.samp)/100, pro=mean(pro.samp)/100,
              ela=mean(ela.samp)/100, cas=mean(cas.samp)/100, rob=mean(rob.samp)/100, other=mean(other.samp)/100)

#ndfas <- list(aca=1, alb=1, ald=1, pro=1, ela=1, cas=1, rob=1, other=1)

