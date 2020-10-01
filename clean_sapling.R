
# 
# FIA Sapling data
# clean data, add fixation rates, calc N fixed by saplings
# by Anika Petach 
# 5/21/18

library(dplyr)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

sap <- read.csv("data/SEEDLING.csv", header=T)
#sap <- sapt[sapt$state!='PR'&sapt$state!='VI'&sapt$state!='HI', ]

#---------- Assign fixing status
fixer <- read.csv("data/FIA_Fixer_spcd.csv", header=T) 
sap_spcd <- merge(sap, fixer, 
                  by.x = "SPCD", by.y = "SPCD", 
                  all.x = TRUE, all.y = FALSE)
sap_spcd[is.na(sap_spcd$FIX), ]$FIX <- 0
sap_spcd <- sap_spcd[, c(1:5, 6, 7, 13, 14, 22, 33:36)]

#---------- Get the plots with fixer present
# load in plot data.
p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv",header=T)
ffp <- p[, c(1, 2, 3, 5, 6, 7, 8, 9)] # pcn (plot number), meastime, lat, lon, elev, stdage(stand age), xmh
sdata <- merge(ffp, sap_spcd, 
               by.x = "pcn", by.y = "PLT_CN", 
               all.x = FALSE, all.y = TRUE) 

sdata_fp_pcn <- unique(sdata[sdata$FIX==1, ]$pcn) #unique plot number with fixer present
sdata_fp <- data.frame(pcn = sdata_fp_pcn, 
                       fixer_present = 1)
sdata_fp1 <- merge(sdata_fp, sdata, 
                   by.x = "pcn", by.y = "pcn", 
                   all.x = TRUE, all.y = TRUE)
sdata_fp2 <- sdata_fp1[!is.na(sdata_fp1$fixer_present), ]

sdata <- sdata_fp1[sdata_fp1$state != 'PR' & 
                     sdata_fp1$state != 'VI' & 
                     sdata_fp1$state != 'HI', ] #why didn't we use gdata_fp2?
sdata[sdata$tpha == -999, ]$TPA_UNADJ <- median(sdata[sdata$tpha != -999, ]$TPA_UNADJ)

#--------- Get total number of fixers (for some reason there are tons of na rows)
notna <- sdata[!is.na(sdata$Genus), ]
nrow(notna)

# get total number of fixers using treecount
notna2 <- notna[!is.na(notna$TREECOUNT), ]
counts <- notna %>%
  group_by(Genus) %>%
  summarize(sum(TREECOUNT, na.rm=T))
colnames(counts) <- c("Genus","numstem")
sum(counts$numstem)

dbhcm <- sqrt(2.54*12.7) #Assume all trees at max seedling size
dbhm <- dbhcm/100

sdatat <- sdata[!is.na(sdata$TREECOUNT), ]
sdatat[sdatat$TREECOUNT==999, ]$TREECOUNT <- NA

sdata_pl <- sdatat %>%
  mutate(tpha=TPA_UNADJ*0.404686, BAm2ha=(dbhm/2)^2*3.14*TREECOUNT*tpha)

################## intercept method
# run first part of do_bootstrapping.R to get data for lm
#regression, get mean and SE
ac <- lm(rate ~ BA, data=acacia)
ac.mu <- summary(ac)$coefficients[2,1]
ac.se <- summary(ac)$coefficients[2,2]
ac.int <- summary(ac)$coefficients[1,1]
ac.ints <- summary(ac)$coefficients[1,2]

al <- lm(rate ~ BA, data=alnus)
al.mu <- summary(al)$coefficients[2,1]
al.se <- summary(al)$coefficients[2,2]
al.int <- summary(al)$coefficients[1,1]
al.ints <- summary(al)$coefficients[1,2]

pr <- lm(rate ~ BA, data=prosopis)
pr.mu <- summary(pr)$coefficients[2,1]
pr.se <- summary(pr)$coefficients[2,2]
pr.int <- summary(pr)$coefficients[1,1]
pr.ints <- summary(pr)$coefficients[1,2]

ro <- lm(rate ~ BA, data=robinia)
#ro.mu <- mean(robinia$rate)
#ro.se <- sd(robinia$rate)/sqrt(nrow(robinia))
ro.mu <- summary(ro)$coefficients[2,1]
ro.se <- summary(ro)$coefficients[2,2]
ro.int <- summary(ro)$coefficients[1,1]
ro.ints <- summary(ro)$coefficients[1,2]

other.mu <- mean(ac.mu, al.mu, pr.mu, ro.mu)
other.se <- sqrt(ro.se^2+pr.se^2+al.se^2+ac.se^2)
other.int <- (ac.int + al.int + pr.int + ro.int)/4
other.ints <- sqrt(ro.ints^2+pr.ints^2+al.ints^2+ac.ints^2)

################## bootstrap F
# Robinia fixed rate: rnorm(1,mean(robinia$rate),var(robinia$rate)^(0.5)) 
# Robinia varying rate: rnorm(1,ro.mu,ro.se)*plotBA+ro.int
# filter(notna2, FIX==1)
plot_fixs <- sdata_pl %>% 
  group_by(pcn, Genus) %>%
  summarize(plotBA=sum(BAm2ha,na.rm=T), cnt=n()) %>%
  mutate(fix.s = ifelse(Genus %in% c("Acacia"), rnorm(1,ac.mu,ac.se)*plotBA+ac.int, 
                        ifelse(Genus %in% c("Albizia"), rnorm(1,other.mu,other.se)*plotBA+other.int,
                               ifelse(Genus %in% c("Alnus"), rnorm(1,al.mu,al.se)*plotBA+al.int,
                                      ifelse(Genus %in% c("Cercocarpus"), rnorm(1,other.mu,other.se)*plotBA+other.int,
                                             ifelse(Genus %in% c("Prosopis"), rnorm(1,pr.mu,pr.se)*plotBA+pr.int, 
                                                    ifelse(Genus %in% c("Ebenopsis"), rnorm(1,other.mu,other.se)*plotBA+other.int,
                                                           ifelse(Genus %in% c("Robinia"), rnorm(1,mean(robinia$rate),var(robinia$rate)^(0.5)),
                                                                  ifelse(Genus %in% c("Olneya"), rnorm(1,other.mu,other.se)*plotBA+other.int,
                                                                         ifelse(Genus %in% c("Elaeagnus"), rnorm(1,other.mu,other.se)*plotBA+other.int,
                                                                                0))))))))))

numplots <- length(unique(sdata_pl$pcn))
test <- plot_fixs %>% 
  group_by(pcn) %>%
  summarize(pfix=sum(fix.s,na.rm=T), sumBA=sum(plotBA,na.rm=T)) #add total basal area to this as well

test1 <- plot_fixs %>% 
  filter(!is.na(Genus)) %>%
  group_by(pcn) %>%
  summarize(fixerBA=sum(plotBA,na.rm=T)) #add total basal area to this as well

totals <- mean(test$pfix)*2428.114*numplots
totals

hist(plot_fixs$plotBA,breaks=60,xlim=c(0,200))

sum(sdata_pl$TREECOUNT*sdata_pl$FIX, na.rm=T) #number of nfixing stems

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/output")
saveRDS(test, "saplingfix.RDS")
saveRDS(test1, "saplingba.RDS")
saveRDS(sdata_pl, "saplingdata.RDS")
saveRDS(plot_fixs, "saplingfixgenus.RDS")

# get stem counts of each genus
fixers <- filter(sdata_pl, FIX==1)
fixers %>%
  group_by(Genus) %>%
  summarize(cnt=sum(TREECOUNT,na.rm=T))
