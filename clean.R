
# clean.R
# loads and cleans FIA data (all data)
# by Anika Petach (parts adapted from Wenying Liao)
# 6/5/17

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

# load and prep tree data -----------------------
ft <- readRDS("output/t0.RDS")

# remove islands, assign regions
regions <- read.csv("data/FIA_regions.csv", header=TRUE)
ft <- ft[ft$state != 'PR' &
           ft$state != 'VI' &
           ft$state != 'HI', ]
ft <- merge(ft, regions[,c(1,2)], 
            by.x = "state", by.y = "State.Abb", 
            all.x = TRUE, all.y = FALSE, 
            sort = FALSE)

# Assign fixing status
fixer <- read.csv("data/FIA_Fixer_spcd.csv", header=T) 
tree_spcd <- merge(ft, fixer, 
                   by.x = "spcd", by.y = "SPCD", 
                   all.x = TRUE, all.y = FALSE)
tree_spcd[is.na(tree_spcd$FIX), ]$FIX <- 0 #species not in fixer list get FIX=0
tree_spcd[tree_spcd$tpha == -999, ]$tpha <- median(tree_spcd[tree_spcd$tpha != -999, ]$tpha) #if missing use median value

# Get the plots with fixer present
# load in plot data.
p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv", header=T)
ffp <- p[ ,c(2,3,5,6,7,8,9)] # pcn (plot number), meastime, lat, lon, elev, stdage(stand age), xmh
gdata <- merge(ffp, tree_spcd, 
               by.x = "pcn", by.y = "pcn", 
               all.x = FALSE, all.y = TRUE) 

gdata_fp_pcn <- unique(gdata[gdata$FIX==1, ]$pcn) #unique plot number with fixer present
gdata_fp <- data.frame(pcn = gdata_fp_pcn, 
                       fixer_present = 1)
rm(t, gdata_fp_pcn, p, ft) #free up some memory
gdata_fp1 <- merge(gdata_fp, gdata, 
                   by.x = "pcn", by.y = "pcn", 
                   all.x = TRUE, all.y = TRUE)
gdata_fp2 <- gdata_fp1[!is.na(gdata_fp1$fixer_present), ]
gdata <- gdata_fp1 # to get all trees (not just plots with a fixer)

gdata <- gdata[gdata$dbh >= 0,] #only take positive dbh (will delete others) 1.96% trees
gdata$BAm2 <- ((gdata$dbh/200)^2*3.14)

rm(gdata_fp1)

#recode NAs
gdata[gdata$meastime < 0, ]$meastime <- NA
gdata[gdata$elev < 0, ]$elev <- NA
gdata[gdata$stdage < 0, ]$stdage <- NA

# Save the cleaned data --------------------
saveRDS(gdata, "output/gdata.RDS")
