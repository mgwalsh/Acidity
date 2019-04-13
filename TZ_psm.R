# Stacked prediction maps of Tanzania topsoil measurements
# M. Walsh, April 2019

# Required packages
# install.packages(c("devtools","caret","MASS","randomForest","gbm","nnet","plyr","doParallel")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(MASS)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(plyr)
  require(doParallel)
})

# Data setup --------------------------------------------------------------
# SourceURL <- "https://github.com/mgwalsh/Acidity/blob/master/TZ_pH_data.R"
# source_url(SourceURL)
rm(list=setdiff(ls(), c("sdat","grids"))) ## scrubs extraneous objects in memory

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
sIndex <- createDataPartition(sdat$pH, p = 4/5, list = F, times = 1)
s_cal <- sdat[ sIndex,]
s_val <- sdat[-sIndex,]

# Soil property calibration labels
cp_cal <- s_cal$pH

# raster calibration features
gf_cal <- s_cal[,5:50]

# GLM <glmStepAIC> --------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# model training
gl <- train(gf_cal, cp_cal, 
            method = "glmStepAIC",
            family = "gaussian",
            preProc = c("center","scale"), 
            trControl = tc,
            metric ="RMSE")

# model outputs & predictions
gl
summary(gl)
gl.pred <- predict(grids, gl) ## spatial predictions
stopCluster(mc)
saveRDS(gl, "./Results/gl.rds")

# Random forest <randomForest> --------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)
tg <- expand.grid(mtry = seq(1,8, by=1)) ## model tuning steps

# model training
rf <- train(gf_cal, cp_cal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 501,
            tuneGrid = tg,
            trControl = tc)

# model outputs & predictions
print(rf) ## RMSEs accross tuning parameters
rf.pred <- predict(grids, rf) ## spatial predictions
stopCluster(mc)
saveRDS(rf, "./Results/rf.rds")

# Generalized boosting <gbm> ----------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

## for initial <gbm> tuning guidelines see @ https://stats.stackexchange.com/questions/25748/what-are-some-useful-guidelines-for-gbm-parameters
tg <- expand.grid(interaction.depth = seq(6,12, by=2), shrinkage = seq(0.02,0.1, by=0.02), n.trees = 501,
                  n.minobsinnode = 25) ## model tuning steps

# model training
gb <- train(gf_cal, cp_cal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg)

# model outputs & predictions
print(gb) ## RMSEs accross tuning parameters
gb.pred <- predict(grids, gb) ## spatial predictions
stopCluster(mc)
saveRDS(gb, "./Results/gb.rds")

# Neural network <nnet> ---------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)
tg <- expand.grid(size = seq(2,10, by=2), decay = c(0.0001,0.001,0.01,0.1)) ## model tuning steps

# model training
nn <- train(gf_cal, cp_cal, 
            method = "nnet",
            preProc = c("center","scale"), 
            tuneGrid = tg,
            trControl = tc,
            metric ="RMSE")

# model outputs & predictions
print(nn) ## RMSEs accross tuning parameters
nn.pred <- predict(grids, nn) ## spatial predictions
stopCluster(mc)
saveRDS(nn, "./Results/nn.rds")

# Model stacking setup ----------------------------------------------------
preds <- stack(gl.pred, rf.pred, gb.pred, nn.pred)
names(preds) <- c("gl","rf","gb","nn")
plot(preds, axes = F)

# extract model predictions
coordinates(s_val) <- ~x+y
projection(s_val) <- projection(preds)
spred <- extract(preds, s_val)
spred <- as.data.frame(cbind(s_val, spred))

# stacking model validation labels and features
cp_val <- spred$pH
gf_val <- spred[,51:54] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# model training
st <- train(gf_val, cp_val,
            method = "glm",
            trControl = tc)

# model outputs & predictions
st
summary(st)
plot(varImp(st))
st.pred <- predict(preds, st) ## spatial predictions
plot(st.pred, axes = F)

stopCluster(mc)
saveRDS(st, "./Results/st.rds")

# Write prediction grids --------------------------------------------------
gspreds <- stack(gl.pred, rf.pred, gb.pred, nn.pred, st.pred)
names(gspreds) <- c("gl","rf","gb","nn","st")
writeRaster(gspreds, filename="./Results/TZ_pH_2019.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)## ... change feature names here

# Write output data frame -------------------------------------------------
# coordinates(gsdat) <- ~x+y
# projection(gsdat) <- projection(grids)
# gspre <- extract(gspreds, gsdat)
# gsout <- as.data.frame(cbind(gsdat, gspre))
# write.csv(gsout, "./Results/TZ_pH_out.csv", row.names = F)
