# AfSIS cation exchange reference data setup
# M. Walsh, October 2017

# Required packages
# install.packages(c("downloader","rgdal","raster")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(lattice)
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
download("https://www.dropbox.com/s/a3zjx6m72ya9vc6/AF_prod_grids.zip?raw=1", "AF_prod_grids.zip", mode="wb")
unzip("AF_prod_grids.zip", overwrite=T)
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
hp.glmer <- glmer(I(Hp<1)~pH+(1|Site), family=binomial, data=cecdat)
summary(hp.glmer)
fixefs <- fixef(hp.glmer)
ED <- (log(0.99/0.01)-fixefs[1])/fixefs[2] ## pH at which 99% of observations are expected below 1 cmol/kg Hp

# pH vs Hp plot

plot(Hp~pH, cecdat, ylab = expression("Exch. acidity" ~ (cmol[c] ~ kg^{-1})), xlab="pH (Water)")
abline(h=1, v=ED, col="red")

