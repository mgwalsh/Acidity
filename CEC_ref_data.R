# AfSIS cation exchange reference data setup
# M. Walsh, October 2017

# Required packages
# install.packages(c("downloader","rgdal","raster")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(lattice)
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
cecref <- merge(prof, samp, by="SSN")

