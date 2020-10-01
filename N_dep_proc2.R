
##################
#
# Nitrogen inputs from raster data (SNF trees, understory, asymbiotic, deposition)
# Gets N deposition
# Author: Anika Petach
# Date: 4/6/18
#
##################

#rm(list=ls())

require(rgdal)
require(raster)
require(sp)
require(fields)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")

# Nitrogen Deposition (Load and project) -----------------

#dry hno3
N_HNO3 <- raster("Ndep_clim_nitrogen_730/data/NDDN_dry_deposition_hno3conc_hno3vd_0.5x0.5_grid_annual.txt")
proj4string(N_HNO3) <- crs("+init=epsg:4326") #defining coordinate system to the WGS84 lat/long
myCRS <- CRS(proj4string(N_HNO3))
N_HNO3_1 <- projectRaster(N_HNO3, res = 1, crs = myCRS)
ncol <- 16
plot(N_HNO3_1, 
     xlab = "lon", 
     ylab = "lat", 
     main = "N from Nitrate Deposition", 
     xlim = c(-126, -66), 
     ylim = c(24, 51),
     breaks=seq(0, ncol, by = 1), 
     col = terrain.colors(ncol))
extent(N_HNO3_1)
e <- extent(-135.82, -67+1, 24, 59.38)
N_HNO3_1 <- extend(N_HNO3_1, e) #match extent

#dry no3
N_NO3 <- raster("Ndep_clim_nitrogen_730/data/NDDN_dry_deposition_no3conc_no3vd_0.5x0.5_grid_annual.txt")
proj4string(N_NO3) <- crs("+init=epsg:4326") #defining coordinate system to the WGS84 lat/long
N_NO3_1 <- projectRaster(N_NO3, res = 1, crs = myCRS)
e <- extent(-135.82, -67+1, 24, 59.38)
N_NO3_1 <- extend(N_NO3_1, e)
plot(N_NO3_1, xlab = "lon", ylab = "lat", main = "N from Dry Nitrate Deposition")

#wet no3
N_NO3w <- raster("Ndep_clim_nitrogen_730/data/NADP_wet_deposition_no3_0.5x0.5_grid_annual_R1.txt")
proj4string(N_NO3w) <- crs("+init=epsg:4326") #defining coordinate system to the WGS84 lat/long
N_NO3w_1 <- projectRaster(N_NO3w, res = 1, crs = myCRS)
e <- extent(-135.82, -67+1, 24, 59.38)
N_NO3w_1 <- extend(N_NO3w_1, e)
plot(N_NO3w_1, xlab = "lon", ylab = "lat", main = "N from Wet Nitrate Deposition")

#total nh4
N_NO3_T <- N_HNO3_1 + N_NO3_1 + N_NO3w_1
plot(N_NO3_T, xlab = "lon", ylab = "lat", main = "N from Nitrate Deposition")
map('world', regions = 'usa', add = TRUE)

#dry nh4
N_NH4 <- raster("Ndep_clim_nitrogen_730/data/NDDN_dry_deposition_nh4conc_particulatevd_0.5x0.5_grid_annual.txt")
proj4string(N_NH4) <- crs("+init=epsg:4326") #defining coordinate system to the WGS84 lat/long
image(N_NH4)
myCRS <- CRS(proj4string(N_NH4))
N_NH4_1 <- projectRaster(N_NH4, res = 1, crs = myCRS)
e <- extent(-135.82,-67+1,24,59.38)
N_NH4_1 <- extend(N_NH4_1, e)
plot(N_NH4_1, xlab = "lon", ylab = "lat", main = "N from Dry Ammonium Deposition")
map('world', regions = 'usa', add = TRUE)
extent(N_NH4_1)

#wet nh4
N_NH4w <- raster("Ndep_clim_nitrogen_730/data/NADP_wet_deposition_nh4_0.5x0.5_grid_annual_R1.txt")
proj4string(N_NH4w) <- crs("+init=epsg:4326") #defining coordinate system to the WGS84 lat/long
N_NH4w_1 <- projectRaster(N_NH4w, res = 1, crs = myCRS)
e <- extent(-135.82, -67+1, 24, 59.38)
N_NH4w_1 <- extend(N_NH4w_1, e)
plot(N_NH4w_1, xlab = "lon", ylab = "lat", main = "N from Wet Ammonium Deposition")
map('world', regions = 'usa', add = TRUE)

#total nh4
N_NH4_T <- N_NH4w_1 + N_NH4_1
plot(N_NH4_T, xlab = "lon", ylab = "lat", main = "N from Ammonium Deposition")
map('world', regions = 'usa', add = TRUE)

# Forest Cover -----------------------------------------
#Import and project forest cover
FOR_COV_ORIG <- raster("foresti020l/foresti020l.tif")
plot(FOR_COV_ORIG)
proj4string(FOR_COV_ORIG) <- crs("+proj=laea +lat_0=45 +lon_0=-100 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Make masks to isolate forest
# forest = 2-22, 25
# ground = 23
# water = 0, 24
# non-US = 1
FOR_MASK <- FOR_COV_ORIG #make mask for SNF data
FOR_MASK[(FOR_MASK >= 2 & FOR_MASK < 23) | FOR_MASK == 25] <- 10
FOR_MASK_1 <- FOR_MASK
FOR_MASK_1[(FOR_MASK_1 < 10)|(FOR_MASK_1 >= 23)] <- 0 #non-forest
FOR_MASK_2 <- FOR_MASK_1
FOR_MASK_2[FOR_MASK_2 == 10] <- 1 #forest

# Reproject
FOR_MASK_2r <- projectRaster(FOR_MASK_2, crs = crs(N_NH4w_1)) #reproject based on dep crs

# Clip extent
FOR_MASK_2s <- crop(FOR_MASK_2r, extent(N_NH4w_1))
FOR_MASK_2_1s <- FOR_MASK_2s
FOR_MASK_2_1s[FOR_MASK_2_1s > 0] <- 1

# Make masks to isolate ground
GR_MASK <- FOR_COV_ORIG #make mask for SNF data
GR_MASK[(GR_MASK >= 2 & GR_MASK < 24) | GR_MASK == 25] <- 10 #ground
GR_MASK_1 <- GR_MASK
GR_MASK_1[(GR_MASK_1 < 2)|(GR_MASK_1 == 24)] <- 0 #non-ground
GR_MASK_1[GR_MASK_1 == 10] <- 1 #ground
#remove unnecessary rasters
# Reproject
GR_MASK_1r <- projectRaster(GR_MASK_1, crs = crs(N_NH4w_1)) #reproject based on dep crs

# Clip extent
GR_MASK_1s <- crop(GR_MASK_1r, extent(N_NH4w_1))

#aggregate
#aggregate(x, fact=2, fun=mean, expand=TRUE, na.rm=TRUE, filename='', ...) where  two integers (horizontal and vertical aggregation factor)
Nhr <- dim(FOR_MASK_2s) # resolution of high-res raster
Nlr <- dim(N_NH4w_1) # resolution of low-res raster
Nratio <- c(as.integer(Nhr[2]/Nlr[2])+1, as.integer(Nhr[1]/Nlr[1])+1) # ratio of high to low resolutions, to nearest integer value for aggregation
FOR_COV_FOR <- aggregate(FOR_MASK_2s, 
                         fact = Nratio, 
                         fun = sum) #fact=c(62,28) #gives raster of number of km^2 of forest within 1 degree cell
FOR_COV_ALL <- aggregate(GR_MASK_1s, #was FOR_MASK_2_1s before
                         fact = Nratio, 
                         fun = sum) #get total number km2 cells per degree2 cell
FOR_COV <- FOR_COV_FOR/FOR_COV_ALL #frac of cell covered in forest

#projectRaster(from, to, res, crs, method="bilinear", alignOnly=FALSE, over=FALSE, filename="", ...) 
FOR_COV_ALIGN <- projectRaster(from = FOR_COV, 
                               to = N_NH4w_1, 
                               res = res(N_NH4w_1), 
                               crs = crs(N_NH4w_1), 
                               method = "bilinear", 
                               alignOnly = FALSE, 
                               over = FALSE) 

FOR_COV_FOR_ALIGN <- projectRaster(from = FOR_COV_FOR, 
                               to = N_NH4w_1, 
                               res = res(N_NH4w_1), 
                               crs = crs(N_NH4w_1), 
                               method = "bilinear", 
                               alignOnly = FALSE, 
                               over = FALSE) 

FOR_COV_ALL_ALIGN <- projectRaster(from = FOR_COV_ALL, 
                               to = N_NH4w_1, 
                               res = res(N_NH4w_1), 
                               crs = crs(N_NH4w_1), 
                               method = "bilinear", 
                               alignOnly = FALSE, 
                               over = FALSE) 

writeRaster(FOR_COV_FOR_ALIGN, "output/FOR_COV_FOR_ALIGN2.grd")
writeRaster(FOR_COV_ALIGN, "output/FOR_COV_ALIGN2.grd")
writeRaster(FOR_COV_ALL_ALIGN, "output/FOR_COV_ALL_ALIGN2.grd")

res(N_NO3)
res(FOR_COV_ORIG)

# N-Fixing Trees --------------------------------
#Import and project SNF
SNFgr_GRID <- readRDS("output/fngrgrid_new.RDS")
#SNF_GRID <- readRDS("output/SNF_grid_mat.RDS")
SNF_GRID <- readRDS("output/fngrid_new.RDS")

SNF_RAS_0 <- raster(SNF_GRID) #convert to raster
SNFgr_RAS_0 <- raster(SNFgr_GRID) #convert to raster

SNF_RAS <- t(flip(SNF_RAS_0, 1)) #rotate because loads sideways
SNFgr_RAS <- t(flip(SNFgr_RAS_0, 1)) #rotate because loads sideways

SNF_extent <- extent(-155, -64, 17, 62) #lat,lon bounds used to create SNF_grid_mat.RDS
extent(SNF_RAS) <- SNF_extent
extent(SNFgr_RAS) <- SNF_extent

# Plot output
plot(SNF_RAS, 
     xlab = "lon", 
     ylab = "lat", 
     main = "N from SNF Trees", 
     xlim = c(-136, -66), 
     ylim = c(24, 60))
map('world',
    regions = 'usa',
    add = TRUE) #map is ok

# Resolution (want 1x1 so can mix with other data)
res(SNF_RAS) # [1] 0.9891304 0.9782609
res(SNFgr_RAS)
SNF_RES <- resample(SNF_RAS, N_NH4_T, method = 'bilinear') #get resolution from N_NH4_T
SNFgr_RES <- resample(SNFgr_RAS, N_NH4_T, method = 'bilinear')

plot(SNF_RES, 
     xlab = "lon", 
     ylab = "lat", 
     main = "N from SNF Trees", 
     xlim = c(-136, -66), 
     ylim = c(24, 60))
map('world', regions = 'usa', add = TRUE)

# Fill SNF_RAS missing values with 0 
# (make mask of FOR_COV_ALIGN, set all values to 0, add to SNF_RAS)
FOR_COV_SNF <- FOR_COV_ALIGN
FOR_COV_SNF[FOR_COV_SNF >= 0.01] <- 0.01
SNF_RASf <- SNF_RES
plot(SNF_RASf)

# Make masks to check intersection
SNF_MASK <- SNF_RASf #make mask for SNF data
SNF_MASK[SNF_MASK > 0] <- 1
plot(SNF_MASK, main = "SNF Mask")

NO3_MASK <- N_NO3_1 #make mas for nitrate data
NO3_MASK[NO3_MASK > 0] <- 1
plot(NO3_MASK, main = "NO3 Mask")

SNF_NO3_MASK <- SNF_MASK * NO3_MASK #where we have data for both data rasters
plot(SNF_NO3_MASK)

# Calculate maps --------------------------------------------------
# Understory
FOR_COV_MASK <- FOR_COV_ALIGN * SNF_NO3_MASK #where we have data for both rasters
N_UNDER <- 2.20*FOR_COV_MASK
#N_UNDER <- 7.32*FOR_COV_MASK
plot(N_UNDER, xlab = "lon", ylab = "lat", main = "N from SNF Understory", zlim = c(0,24))
map('world', regions = 'usa', add = TRUE)

# Asymbiotic
N_ASYM <- 1.70*FOR_COV_ALIGN*SNF_NO3_MASK
#N_ASYM <- 1.70*FOR_COV_ALIGN
plot(N_ASYM, xlab = "lon", ylab = "lat", main = "N from BNF Asymbiotic", zlim = c(0,24))
map('world', regions = 'usa', add = TRUE)

# Make mask where 1 is when there should be forest
FOR_COV_ALIGN_MASK <- FOR_COV_ALIGN
FOR_COV_ALIGN_MASK[FOR_COV_ALIGN_MASK > 0.00001] <- 1
FOR_COV_ALIGN_MASK[FOR_COV_ALIGN_MASK <= 0.00001] <- NA
plot(FOR_COV_ALIGN_MASK)

# Deposition FOR_COV_ALIGN=frac gric cell forested, FOR_COV_ALIGN_MASK=0 over water 1 over land
N_dep <- N_NO3_T*FOR_COV_ALIGN*FOR_COV_ALIGN_MASK + 
  N_NH4_T*FOR_COV_ALIGN*FOR_COV_ALIGN_MASK
plot(N_dep, 
     xlab = "lon", 
     ylab = "lat", 
     main = "N from N Deposition", 
     xlim = c(-126, -66), 
     ylim = c(24, 51), 
     zlim = c(0, 24))
map('world', regions = 'usa', add = TRUE)

# Calculate total N across grid cells
N_NO3_Ty <- reclassify(N_NO3_T, cbind(NA, 0)) 
N_NH4_Ty <- reclassify(N_NH4_T, cbind(NA, 0))
N_TOT <- SNF_RASf + 
  N_NO3_Ty*FOR_COV_ALIGN_MASK + 
  N_NH4_Ty*FOR_COV_ALIGN_MASK + 
  N_UNDER + 
  N_ASYM
plot(N_TOT, xlab = "lon", ylab = "lat", main = "Total N input", zlim = c(0, 24))
map('world',regions='usa',add=TRUE)

# Percent SNF across grid cells
N_PCT <- SNF_RASf/N_TOT
plot(N_PCT*100, xlab = "lon", ylab = "lat", main = "Percent N Input from SNF Trees")
map('world', regions = 'usa', add = TRUE)

#Percent symbiotic
N_PCT_SNF <- (SNF_RASf*FOR_COV_ALIGN_MASK + N_UNDER*FOR_COV_ALIGN_MASK)/N_TOT*FOR_COV_ALIGN_MASK
plot(N_PCT_SNF*100, xlab = "lon", ylab = "lat", main = "Percent N Input from all SNF")
map('world', regions = 'usa', add = TRUE)

#Percent biological
N_PCT_BIO <- (SNF_RASf + N_UNDER + N_ASYM)/N_TOT
colfunc <- colorRampPalette(c("white", "black"))
plot(N_PCT_SNF*100, 
     xlab = "lon", 
     ylab = "lat", 
     main = "Percent N Input from Biological Fixation", 
     col = colfunc(10))
#col=cm.colors(10)
map('world', regions = 'usa', add = TRUE)

#Nitrate deposition
N_NO3_masked <- N_NO3_T*FOR_COV_ALIGN
plot(N_NO3_masked, 
     xlab = "lon", 
     ylab = "lat", 
     main = "N from Nitrate Deposition", 
     xlim = c(-126, -66), 
     ylim = c(24, 51))
map('world', regions = 'usa', add = TRUE)

N_PCT_tree <- (SNF_RASf*FOR_COV_ALIGN_MASK)/N_TOT*FOR_COV_ALIGN_MASK


# save outputs

writeRaster(SNF_RASf, "output/SNF_RASf.grd")
writeRaster(N_ASYM, "output/N_ASYM.grd")
writeRaster(N_dep, "output/N_dep.grd")
writeRaster(N_TOT, "output/N_TOT.grd")
writeRaster(N_UNDER, "output/N_UNDER.grd")
writeRaster(N_PCT_SNF, "output/N_PCT_SNF.grd")
writeRaster(N_PCT_tree, "output/N_PCT_tree.grd")

# FIGURE 5: map figures -------------------------------------
label <- function(px, py, lab, ..., adj=c(0, 1)) {
  usr <- par("usr")
  text(usr[1] + px*(usr[2] - usr[1]),
       usr[3] + py*(usr[4] - usr[3]),
       lab, adj=adj, ...)
}

SNF_RASf <- raster("output/SNF_RASf.grd")
N_ASYM <- raster("output/N_ASYM.grd")
N_dep <- raster("output/N_dep.grd")
N_TOT <- raster("output/N_TOT.grd")
N_UNDER <- raster("output/N_UNDER.grd")
N_PCT_SNF <- raster("output/N_PCT_SNF.grd")

png('output/figure6.png', width = 8, height = 9, units = 'in', res = 300)
#par(mfrow=c(3, 2), mar=c(4.1, 2.1, 1.5, 0.5), oma=c(0, 2, 0, 0), mgp=c(2.25, 1, 0)) #keep all y-axes
par(mfrow = c(3, 2), 
    mar = c(4.1, 2.2, 2.2, 1), 
    oma = c(0, 3, 3, 3), 
    mgp = c(2.25, 1, 0),
    bty="n") #drop some y-axes, c(b,l,t,r)
#trees 
ticks <- c(0, 5, 20, 35) #*FOR_COV_ALIGN_MASK
plot(sqrt((SNF_RASf)),
     xlim = c(-136, -66), 
     ylim = c(24, 59), 
     zlim = c(sqrt(0), sqrt(35)), 
     cex.main = 1.5, 
     cex.lab = 1.5, 
     las = 1,
     axis.args = list(at = sqrt(ticks),
                      labels = ticks),
     legend.args = list(text = expression(paste("N Fixed (kg N ha"^"-1"," yr"^"-1",")")), 
                        side = 4, 
                        line = 3))
map('world', 
    regions = 'usa', 
    add = TRUE)
mtext("N from Tree BNF", 
      3, 
      outer = FALSE, 
      line = 0.5, 
      cex = 1.5)
label(.08, 0.96, "a)", cex=1.5)
#understory 
plot(sqrt(N_UNDER), 
     xlim = c(-136, -66), 
     ylim = c(24, 59), 
     zlim = c(sqrt(0), sqrt(35)), 
     cex.main = 1.5, 
     las = 1, 
     yaxt = "n",
     axis.args = list(at = sqrt(ticks),
                      labels = ticks),
     legend.args = list(text = expression(paste("N Fixed (kg N ha"^"-1"," yr"^"-1",")")), 
                        side = 4, 
                        line = 3))
map('world',
    regions = 'usa',
    add = TRUE)
mtext("N from Understory BNF", 3, outer=FALSE, line=0.5, cex=1.5)
label(.08, .96, "b)", cex=1.5)
#axis(2, labels=FALSE) #if you still want y tick marks
#asymbiotic 
plot(sqrt(N_ASYM), 
     xlim = c(-136,-66), 
     ylim = c(24,59), 
     zlim = c(sqrt(0),sqrt(35)), 
     cex.main = 1.5, 
     las = 1,
     axis.args = list(at = sqrt(ticks),
                      labels = ticks),
     legend.args = list(text = expression(paste("N Fixed (kg N ha"^"-1"," yr"^"-1",")")), 
                        side = 4, 
                        line = 3))
map('world',
    regions = 'usa',
    add = TRUE)
mtext("N from Asymbiotic BNF", 3, 
      outer = FALSE, 
      line = 0.5, 
      cex = 1.5)
label(.08, .96, "c)", cex = 1.5)
#deposition 
plot(sqrt(N_dep), 
     xlim = c(-136,-66), 
     ylim = c(24,59), 
     zlim = c(sqrt(0),sqrt(35)), 
     cex.main = 1.5, 
     las = 1, 
     yaxt = "n",
     axis.args = list(at = sqrt(ticks),
                      labels = ticks),
     legend.args = list(text = expression(paste("N (kg N ha"^"-1"," yr"^"-1",")")), 
                        side = 4, 
                        line = 3))
map('world',
    regions = 'usa',
    add = TRUE)
mtext("N from N Deposition", 
      3, 
      outer = FALSE, 
      line = 0.5, 
      cex = 1.5)
label(.08, .96, "d)", cex = 1.5)
#total N
plot(sqrt(N_TOT), 
     xlim = c(-136, -66), 
     ylim = c(24, 59), 
     zlim = c(sqrt(0), sqrt(35)), 
     cex.main = 1.5, 
     las = 1,
     axis.args = list(at = sqrt(ticks),
                      labels = ticks),
     legend.args = list(text = expression(paste("N (kg N ha"^"-1"," yr"^"-1",")")), 
                        side = 4, 
                        line = 3))
map('world',
    regions = 'usa',
    add = TRUE)
mtext("Total N Input", 3, outer = FALSE, line = 0.5, cex = 1.5)
label(.08, .96, "e)", cex=1.5)
#percent biological, colfunc(10)
colfunc <- colorRampPalette(c("white", "black"))
plot(N_PCT_SNF*100, 
     xlim = c(-136, -66), 
     ylim = c(24, 59), 
     zlim = c(0, 100),
     col = tim.colors(10), 
     cex.main = 1.5, 
     las = 1, 
     yaxt = "n",
     legend.args = list(text = "% of total", 
                        side = 4, 
                        line = 3))
map('world',
    regions = 'usa',
    add = TRUE)
mtext("% N Input from Biological N-Fixation", 
      3, 
      outer = FALSE, 
      line = 0.5, 
      yaxt = "n", 
      cex = 1.3)
label(.08, .96, "f)", cex=1.5)
dev.off()

# figure 7
# percent from trees
N_PCT_tree <- raster("output/N_PCT_tree.grd")

png('output/figure7.png', 
    width = 10, 
    height = 8, 
    units = 'in', 
    res = 300)
par(bty="n")
colfunc <- colorRampPalette(c("grey98", "black"))
plot(N_PCT_tree*100, 
     xlim = c(-136, -66), 
     ylim = c(24, 59), 
     zlim = c(0, 100),
     xlab = "longitude", 
     ylab = "latitude", 
     col = tim.colors(100), 
     cex.main = 1.5, 
     main = "Percent N Input from Tree BNF",
     legend.args = list(text = "% of total", 
                        side = 4, 
                        line = 2))
map('world',
    regions = 'usa',
    add = TRUE)
dev.off()


############## Misc ##############
#SNF Trees + Deposition
N_Tree_Dep <- SNF_RAS + N_NO3_T + N_NH4_T 
plot(N_Tree_Dep, xlab="lon", ylab="lat", main="N from SNF Trees & Deposition")
map('world',regions='usa',add=TRUE)

#Forest Percent
plot(FOR_COV_MASK*100, xlab="lon", ylab="lat", main="Percent Forest Cover")
map('world',regions='usa',add=TRUE)

#check values
click(FOR_COV_MASK,n=1)
