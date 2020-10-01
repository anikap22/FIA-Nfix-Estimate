



####################
#
#
# figs_new.R
#
# 3/4/19
# Anika Petach
#
# present processed FIA data
# maps and barcharts
#
# uses new accretion data, Perakis constants, FIA expns scaling
# 
#
#####################

library(fields)
require(rgdal)
require(raster)
require(sp)
require(dplyr)
require(gridExtra)

setwd("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/")


# Fig S1: Genus specific fixation --------------

source("scripts/functions.R")
source("scripts/genus_regressions_new.R")

p <- read.csv("ACGCA_FIAv51_from_JWL/ACGCA_FIAv51_plot.csv", header=T)
p <- p[p$state!='PR'&p$state!='VI'&p$state!='HI', ]
p[p$lat == -999, ]$lat <- NA
p[p$lon == -999, ]$lon <- NA
sf <- readRDS("output/stemsfix.RDS")
sf0 <- merge(sf, p, by="pcn", all.x=TRUE, all.y=FALSE)

# acacia, albizia, olneya, prosopis, robinia, alnus, casuarina, cercocarpus, elaeagnus
pco <- readRDS("output/pcontorta.RDS")

normalrates <- getrates("standard")

aca.grid <- plotsensonly(allstems, 6, normalrates, "Acacia")
aca.grid[aca.grid == 0] <- NA
alb.grid <- plotsensonly(allstems, 6, normalrates, "Albizia")
alb.grid[alb.grid == 0] <- NA
oln.grid <- plotsensonly(allstems, 6, normalrates, "Olneya")
oln.grid[oln.grid == 0] <- NA
pro.grid <- plotsensonly(allstems, 6, normalrates, "Prosopis")
pro.grid[pro.grid == 0] <- NA
rob.grid <- plotsensonly(allstems, 6, normalrates, "Robinia")
rob.grid[rob.grid == 0] <- NA
aln.grid <- plotsensonly(allstems, 6, normalrates, "Alnus")
aln.grid[aln.grid == 0] <- NA
cer.grid <- plotsensonly(allstems, 6, normalrates, "Cercocarpus")
cer.grid[cer.grid == 0] <- NA
ela.grid <- plotsensonly(allstems, 6, normalrates, "Elaeagnus")
ela.grid[ela.grid == 0] <- NA

saveRDS(aca.grid, "output/aca.grid.RDS")
saveRDS(alb.grid, "output/alb.grid.RDS")
saveRDS(oln.grid, "output/oln.grid.RDS")
saveRDS(pro.grid, "output/pro.grid.RDS")
saveRDS(rob.grid, "output/rob.grid.RDS")
saveRDS(aln.grid, "output/aln.grid.RDS")
saveRDS(cer.grid, "output/cer.grid.RDS")
saveRDS(ela.grid, "output/ela.grid.RDS")

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
#ticks <- 2^(0:6)
ticks <- c(1, 8, 16, 32, 64)
image.plot(lon.list, lat.list, (aca.grid), nlevel=nlev, col=tim.colors(nlev),
           axis.args = list(at=(ticks), labels=ticks, las=3),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), 
           zlim = c(0, 71))
title(main = substitute(paste("", italic(Acacia))), cex.main = 1.5)
mtext("(a)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Albizia
image.plot(lon.list, lat.list, alb.grid, nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=(ticks),labels=ticks),
           xlab="",ylab="",
           ylim=c(latminint-1,60+1),xlim=c(-140,lonmaxint+1),zlim=c(0, 71))
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
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 71))
title(main = substitute(paste("", italic(Olneya))), cex.main = 1.5)
mtext("(c)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Prosopis
image.plot(lon.list, lat.list, pro.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 71))
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
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 71))
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
           axis.args = list(at = (ticks), labels = ticks),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 71))
title(main = substitute(paste("", italic(Alnus))), cex.main = 1.5)
mtext("(f)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# # Casuarina
# image.plot(lon.list, lat.list, cas.grid, 
#            nlevel = nlev, col = tim.colors(nlev),
#            axis.args = list(at = (ticks), labels = ticks),
#            xlab = "", ylab = "",
#            ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 155))
# title(main = substitute(paste("", italic(Casuarina))), cex.main = 1.5)
# mtext("(g)", font = 2, adj = 0.01, padj = -0.5)
# map('world', regions = 'usa', add = TRUE)

# Cercocarpus
image.plot(lon.list, lat.list, cer.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           xlab = "", ylab = "",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 71))
title(main = substitute(paste("", italic(Cercocarpus))), cex.main = 1.5)
mtext("(g)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)

# Elaeagnus
image.plot(lon.list, lat.list, ela.grid, 
           nlevel = nlev, col = tim.colors(nlev),
           axis.args = list(at = (ticks), labels = ticks),
           legend.args = list(text = expression(paste("BNF Rate (kg N ha"^"-1"," yr"^"-1",")")), cex=1, side=4, line=4),           
           xlab="",ylab="",
           ylim = c(latminint-1, 60+1), xlim = c(-140, lonmaxint+1), zlim = c(0, 71))
title(main = substitute(paste("", italic(Elaeagnus))), cex.main = 1.5)
mtext("(h)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)
dev.off()


##Fig S2: map of N demand for Populus deltoides ------------
# run section above with gfdata <- agb_bgb from the populus section in doGR.R
pdel.grid <- readRDS("output/pdel_fix.RDS")
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


#Fig S3: BNF by region ------------------------------------
# N by region
regionNames <- c("Interior West","Northern","Pacific Northwest","Southern")

fixsp.ndfa <- readRDS("output/fixsp.ndfa.RDS")
fixsp.acc <- readRDS("output/fixsp.acc.RDS")

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
fn.grid <- readRDS("output/fngrid_new.RDS") #check correct version was saved
fngr.grid <- readRDS("output/fngrgrid_new.RDS")
grn.grid <- readRDS("output/grn_grid_new.RDS")
grngr.grid <- readRDS("output/grngr_grid_new.RDS")

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
ticks <- 2^(-2:6)
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
           zlim = c(log(0.1), log(64)))
title(main = "N Fixation Rate \n Per forest area \n Accretion", 
      cex.main = 1.5)
mtext("(a)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE, col = 'black')
#expression('N Fixation Rate (kg N ha forest'^-1*' yr'^-1*'), Accretion')
#zlim was max a log(35)

# fixed N per ground area rate accretion
nlev <- 64
ticks <- 2^(-2:6)
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
           zlim = c(log(0.1), log(64)))
title(main = "N Fixation Rate \n Per ground area \n Accretion",
      cex.main = 1.5)
mtext("(b)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions  = 'usa', add = TRUE)
#expression('N Fixation Rate (kg N ha ground'^-1*' yr'^-1*'), Accretion')

##map of N demand based on ndfa (per forest area)
nlev <- 64
ticks <- 2^(-2:6)
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
           zlim = c(log(0.1), log(64)))
title(main = "N Fixation Rate \n per forest area \n %Ndfa", 
      cex.main = 1.5)
mtext("(c)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)
#expression('N Fixation Rate (kg N ha forest'^-1*' yr'^-1*'), %Ndfa')

##map of N demand based on ndfa
nlev <- 64
ticks <- 2^(-2:6)
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
           zlim = c(log(0.1), log(64)))
title(main = 'N Fixation Rate \n per ground area \n %Ndfa', 
      cex.main = 1.5)
mtext("(d)", font = 2, adj = 0.01, padj = -0.5)
map('world', regions = 'usa', add = TRUE)
#title(main = expression('N Fixation Rate (kg N ha ground'^-1*' yr'^-1*'), %Ndfa')

dev.off()

# figure 3 - %BNF from each genus --------------------- 
# Figure 3

## ACCRETION from bootstrap (continent scale)
#bootstrap <- readRDS("output/accretion_boot_100.RDS")
bootstrap <- readRDS("output/bootstrapped_acc.RDS")
#bootstrap <- bootstrap/1e9
bootstrap <- bootstrap[,1:11] #only because first time I made this extra col at edge
#bootstrap <- readRDS("output/accretion_boot_continent.RDS")
rows <- c("total","Acacia","Albizia","Alnus","Cercocarpus","Elaeagnus",
          "Olneya","Prosopis","Robinia")
means <- c(mean(bootstrap$tgn), 
           mean(bootstrap$acapct), 
           mean(bootstrap$albpct),
           mean(bootstrap$alnpct),
           mean(bootstrap$cerpct),
           mean(bootstrap$elapct),
           mean(bootstrap$olnpct),
           mean(bootstrap$propct),
           mean(bootstrap$robpct))
sds <- c(sd(bootstrap$tgn),
         sd(bootstrap$acapct),
         sd(bootstrap$albpct),
         sd(bootstrap$alnpct),
         sd(bootstrap$cerpct),
         sd(bootstrap$elapct),
         sd(bootstrap$olnpct),
         sd(bootstrap$propct),
         sd(bootstrap$robpct))
genus.acc.boot <- dplyr::bind_cols(data.frame(rows),
                                   data.frame(means),
                                   data.frame(sds))
genus.acc.boot <- genus.acc.boot[-1,]
# genus.acc.boot tells you how much was fixed by each genus too
colMeans(bootstrap)
apply(bootstrap, 2, sd, na.rm = TRUE)/1e9
genus.acc.boot$means_orig <- genus.light.boot$means
genus.acc.boot$means[genus.acc.boot$means<1] <- round(genus.acc.boot$means[genus.acc.boot$means<1], digits=1)
genus.acc.boot$means[genus.acc.boot$means>=1] <- trunc(round(genus.acc.boot$means[genus.acc.boot$means>=1], digits=0))



## LIGHT (from ndfa)
bootstrapped_light <- readRDS("output/bootstrapped_lightlim.RDS")
#bootstrapped_light <- readRDS("output/accretion_boot_lightlim_continent.RDS")
bootstrapped_light <- bootstrapped_light[,1:9] #remove fixperfor and fixpergr
totalfix <- mean(bootstrapped_light$tgn, na.rm = T)
means <- colMeans(bootstrapped_light, na.rm = T)
sds <- apply(bootstrapped_light, 2, sd)
genus.light.boot <- dplyr::bind_cols(data.frame(rows), 
                                    data.frame(means), 
                                    data.frame(sds))
genus.light.boot <- genus.light.boot[-1,]
#get fixed n by genus with sd
#Ggn_bygen <- (bootstrapped_light[ ,2:ncol(bootstrapped_light)]/100*bootstrapped_light[ ,1])*1e3
#colMeans(Ggn_bygen)
#apply(Ggn_bygen, 2, sd)
genus.light.boot$means_orig <- genus.light.boot$means
genus.light.boot$means[genus.light.boot$means<1] <- round(genus.light.boot$means[genus.light.boot$means<1], digits=1)
genus.light.boot$means[genus.light.boot$means>=1] <- trunc(round(genus.light.boot$means[genus.light.boot$means>=1], digits=0))


## NDFA
bootstrapped_ndfa <- readRDS("output/bootstrapped_ndfa.RDS")
totalfix <- mean(bootstrapped_ndfa$tgn, na.rm = T)
bootstrapped_ndfa <- bootstrapped_ndfa[,-5]
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
genus.ndfa.boot$means_orig <- genus.ndfa.boot$means
genus.ndfa.boot$means[genus.ndfa.boot$means<1] <- round(genus.ndfa.boot$means[genus.ndfa.boot$means<1], digits=1)
genus.ndfa.boot$means[genus.ndfa.boot$means>=1] <- trunc(round(genus.ndfa.boot$means[genus.ndfa.boot$means>=1], digits=0))


## NDFA upper bound
bootstrapped_upper <- readRDS("output/bootstrapped_upper_new.RDS")
totalfix <- bootstrapped_upper$tgn
genus.upper.boot <- bootstrapped_upper
pctupper <- as.vector(unlist(genus.upper.boot))
genus.upper.boot <- dplyr::bind_cols(data.frame(rows),
                                     data.frame(pctupper))
colnames(genus.upper.boot) <- c("rows","means")
genus.upper.boot <- genus.upper.boot[-1,]
genus.upper.boot$means_orig <- genus.upper.boot$means
genus.upper.boot$means[genus.upper.boot$means<1] <- round(genus.upper.boot$means[genus.upper.boot$means<1], digits=1)
genus.upper.boot$means[genus.upper.boot$means>=1] <- trunc(round(genus.upper.boot$means[genus.upper.boot$means>=1], digits=0))


png("output/figure3_new.png", width = 18, height = 8, units = 'in', res = 300)

plot.acc <- ggplot(genus.acc.boot, aes(x = rows, y = means)) +
  geom_bar(stat="identity", position = position_dodge(0.9), alpha=0.7, fill="black") +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black", size=8) +
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
       x = "", 
       y = "% of tree-based BNF") +
  ylim(0,100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face =  "italic",
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   color = "black"),
        axis.text.y = element_text(color="black"),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color =  "black", 
                                     face =   "bold", 
                                     size =   32, 
                                     hjust =  0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color =  "black", 
                                  face =   "bold", 
                                  size =   22)) 

plot.light <- ggplot(genus.light.boot, aes(x = rows, y = means)) +
  geom_bar(stat="identity", position = position_dodge(0.9), alpha=0.7, fill="black") +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black", size=8) +
  geom_errorbar(aes(ymin = pmax(means-sds,0), ymax = pmin(means+sds, 100)), 
                colour = "black", 
                width = 0.02,
                position = position_dodge(0.9),
                size = 1) +
  geom_hline(yintercept = 0, cex = 0.9) +
  #geom_vline(xintercept=which(genus$genus == 0)) + 
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(subtitle = "Light Limitation", 
       x = "", 
       y = "") +
  ylim(0,100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face =  "italic",
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   color = "black"),
        axis.text.y = element_text(color="black"),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color =  "black", 
                                     face =   "bold", 
                                     size =   32, 
                                     hjust =  0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color =  "black", 
                                  face =   "bold", 
                                  size =   22)) 


plot.ndfa <- ggplot(genus.ndfa.boot, aes(x = rows, y = means)) +
  geom_bar(stat="identity", position = position_dodge(0.9), alpha=0.7, fill="black") +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black", size=8) +
  geom_errorbar(aes(ymin = means-sds, ymax = means+sds), 
                colour = "black", 
                width = 0.02,
                position = position_dodge(0.9),
                size = 1) +
  geom_hline(yintercept = 0, cex = 0.9) +
  #geom_vline(xintercept=which(genus$genus == 0)) + 
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(subtitle = "%Ndfa Method", 
       x = "", 
       y = "") +
  ylim(0,100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face =  "italic",
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   color = "black"),
        axis.text.y = element_text(color="black"),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color =  "black", 
                                     face =   "bold", 
                                     size =   32, 
                                     hjust =  0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color =  "black", 
                                  face =   "bold", 
                                  size =   22)) 


plot.ndfaup <- ggplot(genus.upper.boot, aes(x = rows, y = means)) +
  # draw the bar plot
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7, fill = "black") +
  geom_text(aes(label = round(means, digits = 1)), vjust = -2, colour = "black", size=8) +
  geom_hline(yintercept = 0, cex = 0.9) +
  geom_vline(xintercept = 0.4, cex = 1.5) +
  labs(subtitle = "Upper Bound", x = "", y = "") +
  ylim(0, 100) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'),
        text = element_text(size=20), 
        axis.text.x = element_text(face="italic",
                                   angle=90, 
                                   hjust=1, 
                                   vjust=0.5,
                                   color="black"),
        axis.text.y = element_text(color="black"),
        plot.subtitle = element_text(family = "Trebuchet MS", 
                                     color="black", 
                                     face="bold", 
                                     size=32, 
                                     hjust=0.5),
        axis.title = element_text(family = "Trebuchet MS", 
                                  color="black", 
                                  face="bold", 
                                  size=22)) 

grid.arrange(plot.acc, plot.light, plot.ndfa, plot.ndfaup, ncol=4, nrow = 1)
dev.off()

# Total N fixed by each genus ----------------------
#by accretion method
genus.acc %>% mutate(totalfix = (totalacc*genfixpct/100)/1e9)

#by %ndfa method
genus.ndfa %>% mutate(totalfix = (totalndfa*genfixpct/100)/1e9)


# Figure 4: Plot sensivity for accretion method ------------------------------------------
png("output/figure4a_new.png", width = 12, height = 6, units = 'in', res = 300)

accsen <- readRDS("output/forfigs/accsen_new.RDS")
accsen <- accsen[,-7]

par(mar = c(5, 5.5, 4, 2)) #(bottom, left, top, right), default=5,4,4,2
senbara <- barplot(accsen,
                   col = c("#999999", "#E69F00"),
                   ylim = c(-45, 45),
                   beside = T,
                   names.arg = c(expression(atop(italic("Acacia"),~"BNF Rate")),
                                 bquote(atop(italic("Alnus"),~ "BNF Rate")),
                                 bquote(atop(italic("Prosopis"),~ "BNF Rate")),
                                 "Other \ngenera \nBNF Rate",
                                 bquote(atop(italic("Robinia"),~ "BNF Rate")),
                                 "N-fixer \nAbundance"),
                   ylab = expression(paste("% ", Delta, " N fixed by trees")),
                   cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5, cex.names=1.3)

abline(h = 0)
pctch1 <- round(accsen, 2)
pctch1 <- pctch1 + rep(c(4, -4), 6)
text(senbara, 
     pctch1, 
     round(accsen, digits = 2), 
     cex = 1)
mtext("(a)", font = 2, adj = 0.01, padj = -0.5)
legend("topright", c("Upper", "Lower"), 
       fill = c("#999999","#E69F00"), bty = "n")
title("Sensitivity Analysis for Accretion Method", adj = 0.5, line = 0, cex.main=1.5)

dev.off()

png("output/figure4b_new.png", width = 12, height = 6, units = 'in', res = 300)

ndfasen <- readRDS("output/forfigs/ndfasen.RDS") #use ndfasen, ndfasen_new is from wrong scaling method

par(mar = c(5, 5.5, 4, 2)) #(bottom, left, top, right), default=5,4,4,2

senbar <- barplot(ndfasen, col=c("#999999","#E69F00"),
                  ylim=c(-50,95),
                  beside=T,
                  names.arg=c("C:N \nFoliage","C:N \nWood", "C:N \nFine Roots", "N-fixer \nAbundance",
                              "Foliar N \nResorption", "Root N \nResorption"),
                  ylab = expression(paste("% ", Delta, " N fixed by trees")),
                  cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.3)
abline(h = 0)
pctch1 <- round(ndfasen, 2)
#pctch1 <- pctch1 + c(5, -5, 5, -5, 5, -5, 5, -5, 5, -5, 5, -5)
pctch1 <- pctch1 + rep(c(5, -5), 6)
text(senbar, 
     pctch1, 
     round(ndfasen, digits = 2), 
     cex = 1) 
mtext("(b)", font = 2, adj = 0.01, padj = -0.5)
legend("topright", c("Upper", "Lower"), 
       fill = c("#999999","#E69F00"), bty = "n")
title("Sensitivity Analysis for N Demand Method", adj = 0.5, line = 0, cex.main = 1.5)

dev.off()
