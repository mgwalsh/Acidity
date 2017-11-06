# AfSIS cation exchange reference data setup
# M. Walsh, October 2017

# Required packages
# install.packages(c("downloader","rgdal","raster")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("CEC_dat", showWarnings=F)
setwd("./CEC_dat")

# download reference data
download("https://www.dropbox.com/s/4sdz7hj8o4p3quf/CEC_ref_data.zip?raw=1", "CEC_ref_data.zip", mode="wb")
unzip("CEC_ref_data.zip", overwrite=T)
prof <- read.table("Profiles.csv", header=T, sep=",")
samp <- read.table("Samples.csv", header=T, sep=",")
cecr <- merge(prof, samp, by="SSN")
cecr <- na.omit(cecr)

# download Gtifs and stack in raster (note this is a big 550+ Mb download)
download("https://www.dropbox.com/s/uyzsr8oqrwfju19/AF_EEV_grids.zip?raw=1", "AF_EEV_grids.zip", mode="wb")
unzip("AF_EEV_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Data setup ---------------------------------------------------------------
# project GeoSurvey coords to grid CRS
cec.proj <- as.data.frame(project(cbind(cecr$Lon, cecr$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(cec.proj) <- c("x","y")
cecr <- cbind(cecr, cec.proj)
coordinates(cecr) <- ~x+y
projection(cecr) <- projection(grids)

# extract gridded variables at sentinel site locations
cecgrid <- extract(grids, cecr)
cecdat <- as.data.frame(cbind(cecr, cecgrid)) 

# Write output file -------------------------------------------------------
dir.create("Results", showWarnings=F)
write.csv(cecdat, "./Results/AF_cec_dat.csv", row.names = FALSE)

# Site-level glmer --------------------------------------------------------
require(arm)
hp.glmer <- glmer(I(Hp>1)~pH+(1|Site), family=binomial, data=cecdat)
summary(hp.glmer)
fixefs <- fixef(hp.glmer)

# Case definition
p <- 0.01 ## specify false negative probability
CD <- (log(p/(1-p))-fixefs[1])/fixefs[2] ## pH at which p/1-p of observations are expected below 1 cmol/kg Hp

# pH vs Hp plot
par(mfrow=c(1,1), mar=c(4.5,4.5,1,1), pty="s")
plot(Hp~pH, cecdat, ylab = expression("Exch. Acidity" ~ (cmol[c] ~ kg^{-1})), xlab="pH (Water)", cex.lab=1.3)
abline(h=1, v=CD, col="red", lwd=1.5)
text(round(CD,2), 4, paste("Case definition: pH <", round(CD,2)), pos=4, col="red", cex=1.2)

# EEV plots ---------------------------------------------------------------
par(mfrow=c(2,2), mar=c(4.5,4.5,1,1), pty="s")
plot(pH~MAP, cecdat, xlab = expression("Mean annual precipitation" ~ (mm ~ yr^{-1})), ylab="Soil pH (Water)", cex.lab=1.3)
lines(lowess(cecdat$MAP,cecdat$pH), lwd=1.5, col="blue")
abline(h=round(CD,2), lwd=1.5, col="red")
plot(pH~LSTD, cecdat, xlab = expression("Mean land surface temp." ~ (C^{o})), ylab="Soil pH (Water)", cex.lab=1.3)
lines(lowess(cecdat$LSTD,cecdat$pH), lwd=1.5, col="blue")
abline(h=round(CD,2), lwd=1.5, col="red")
plot(pH~PARA, cecdat, xlab = "Mean fAPAR (%)", ylab="Soil pH (Water)", xlim=c(0,100), cex.lab=1.3)
lines(lowess(cecdat$PARA,cecdat$pH), lwd=1.5, col="blue")
abline(h=round(CD,2), lwd=1.5, col="red")
plot(pH~NPPA, cecdat, xlab = expression("Mean NPP" ~ (kg ~ ha^{-1} ~ yr^{-1})), ylab="Soil pH (Water)", xlim=c(0,2), cex.lab=1.3)
lines(lowess(cecdat$NPPA,cecdat$pH), lwd=1.5, col="blue")
abline(h=round(CD,2), lwd=1.5, col="red")
