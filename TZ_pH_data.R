# Tanzania topsoil (0-20 cm) pH & OC data setup
# M. Walsh, April 2019

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
sdat <- read.table("TZ_pH.csv", header=T, sep=",")

# download Gtifs and stack in raster (note this is a big 800MB+ download)
download("https://www.dropbox.com/s/ejl3h62hojnhh3a/TZ_250m_2019.zip?raw=1", "TZ_250m_2019.zip", mode="wb")
unzip("TZ_250m_2019.zip", overwrite=T)

# download latest GeoSurvey prediction grids
download("https://www.dropbox.com/s/g4a80sddyo4o3i2/TZ_GS_preds.zip?raw=1", "TZ_GS_preds.zip", mode="wb")
unzip("TZ_GS_preds.zip", overwrite=T)

# stack grids
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Data setup ---------------------------------------------------------------
# project reference data coords to grid CRS
sdat.proj <- as.data.frame(project(cbind(sdat$lon, sdat$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(sdat.proj) <- c("x","y")
sdat <- cbind(sdat, sdat.proj)
coordinates(sdat) <- ~x+y
projection(sdat) <- projection(grids)

# extract gridded variables at site locations
sgrid <- extract(grids, sdat)
sdat <- as.data.frame(cbind(sdat, sgrid))

# Write output file --------------------------------------------------------
dir.create("Results", showWarnings=F)
write.csv(sdat, "./Results/TZ_sdat.csv", row.names = FALSE)
