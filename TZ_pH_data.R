# Tanzania pH data setup
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
dir.create("TZ_data", showWarnings=F)
setwd("./TZ_data")

# download reference data
download("https://www.dropbox.com/s/7vzb6fhrvga3lh2/TZ_pH.csv.zip?raw=1", "TZ_pH.csv.zip", mode="wb")
unzip("TZ_pH.csv.zip", overwrite=T)
phdat <- read.table("TZ_pH.csv", header=T, sep=",")

# download Gtifs and stack in raster (note this is a big 800+ Mb download)
download("https://www.dropbox.com/s/pshrtvjf7navegu/TZ_250m_2017.zip?raw=1", "TZ_250m_2017.zip", mode="wb")
unzip("TZ_250m_2017.zip", overwrite=T)

# download GeoSurvey prediction grids
download("https://www.dropbox.com/s/3px2xh9l4a6b38g/TZ_GS_preds.zip?raw=1", "TZ_GS_preds.zip", mode="wb")
unzip("TZ_GS_preds.zip", overwrite=T)

# stack grids
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Data setup ---------------------------------------------------------------
# project GeoSurvey coords to grid CRS
ph.proj <- as.data.frame(project(cbind(phdat$lon, phdat$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(ph.proj) <- c("x","y")
phdat <- cbind(phdat, ph.proj)
coordinates(phdat) <- ~x+y
projection(phdat) <- projection(grids)

# extract gridded variables at sentinel site locations
phgrid <- extract(grids, phdat)
phdat <- as.data.frame(cbind(phdat, phgrid))

# Write output file --------------------------------------------------------
dir.create("Results", showWarnings=F)
write.csv(phdat, "./Results/TZ_phdat.csv", row.names = FALSE)
