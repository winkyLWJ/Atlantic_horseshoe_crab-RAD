library(sf)
library(terra)
library(sp)
library(raster)
library(dplyr)
library(dismo)
library(ecospat)
library(rgeos)
library(rJava)
library(knitr)
library(magrittr)
library(ENMeval)
library(rgdal)
setwd("D:/PhD/LP/LPSDM")
########read back species#######

sppts2 <- read.csv("LP_sampling_correct_removerep.csv")
sppdist2 <- dplyr::select(sppts2, longmid, latmid)
sp.occ <- sppdist2
sp.occ <-as.data.frame(sp.occ)

sp.occ

#####Read Bio tif files#######
setwd("D:/PhD/LP/LPSDM/environment/BIO_CLI")
bio_now <- list.files(pattern = ".tif$")
bio_now <- raster::stack(bio_now)
bio_now
plot(bio_now$bio11)
names(bio_now) <- c( "bio7","bio11","bio15","bio16","bio19")

######Read Ocean tif files#####
setwd("D:/PhD/LP/LPSDM/environment/Now_ocean")
ocean_now <- list.files(pattern = ".tif$")
ocean_now <- raster::stack(ocean_now)
ocean_now
names(ocean_now) <- c("ocean1", "ocean3", "ocean4", "ocean5","ocean6","ocean7")
ocean_now<-aggregate(ocean_now, fact = 10)

#######combine bioclimate and ocean####
predictors <- stack(bio_now,ocean_now) 

#########Current:make binary max ice sheet map ######
setwd("D:/PhD/LP/LPSDM/environment/nic_autoc2022056n_pl_a")

current_ice_sheet_shp <- shapefile("D:/PhD/LP/LPSDM/environment/nic_autoc2022056n_pl_a/nic_autoc2022056n_pl_a.shp")
plot(current_ice_sheet_shp)
current_ice_sheet <- rasterize(current_ice_sheet_shp, ice_sheet, 1)
plot(current_ice_sheet)
current_ice_sheet
summary(current_ice_sheet)
values(current_ice_sheet)
writeRaster(current_ice_sheet, "binary_current_ice_sheet_noreclassy.tif", overwrite=TRUE)

current_ice_sheet[is.na(current_ice_sheet)] <- 0
writeRaster(current_ice_sheet, "binary_current_ice_sheet.tif", overwrite=TRUE)


#########calculate blank grids distance to max ice sheet map######

allKmlLayers <- function(kmlfile){
  lyr <- ogrListLayers(kmlfile)
  mykml <- list()
  for (i in 1:length(lyr)) {
    mykml[i] <- readOGR(kmlfile,lyr[i])
  }
  names(mykml) <- lyr
  return(mykml)
}
#made change for Landscape google earth on 20190927

kmlfile <- "nic_autoc2022056n_pl_a.shp"
mykml <- allKmlLayers(kmlfile)
summary(mykml)

#From 'mykml' extract the 'woodded area' data
wa <- mykml$nic_autoc2022056n_pl_a
wa
#Convert "uni" to planar coordinates.
proj4string(wa) <- CRS("+init=epsg:4326")
CRS.new <- CRS("+proj=longlat +datum=WGS84 +no_defs")
wa <- spTransform(wa, CRS.new)
plot(wa)

#Frame the grids
rs <- raster(extent(-105.0001,-45.00014,-0.0001392488,54.99986), crs=projection(wa),nrows=660,ncols=720)
rs[] <- 1:ncell(rs)
rs
##make polygon to grids (as occupied)
x<-rasterize(wa,rs)
plot(x)
#then calculate blank grids distance to place of worship (occupied grids)
y<-distance(x)

setwd("D:/PhD/LP/LPSDM/environment/Current_icecap")
summary(y)
values(y)
writeRaster(y, "current_distance_icesheet.tif", overwrite=TRUE)
#Choose color scheme for visualization
colfunc <- colorRampPalette(c("yellow", "royalblue1"))
#plot the raster ("colfunc" to change the distance interval) 
plot(y, col=colfunc(300))

##STACK ice sheet######
setwd("D:/PhD/LP/LPSDM/environment/Current_icecap")
current_icecombine <- list.files(pattern = ".tif$")
current_icecombine <- stack(current_icecombine)

predictors2 <- stack(predictors, current_icecombine)



### ENMEval using all variables###
rm<-seq(from =1, to = 5, by = 0.5)
tune.args<- list(fc=c("L", "LQ", "LQH", "LQHP", "LQHPT", "H"),rm=rm)
set.seed(123)
options(java.parameters = "-Xmx32000m")
#if occurrence points are more than 25 (see ENMEval vignette)
enmeval_results <- ENMevaluate(sp.occ, predictors2, method="block", n.bg=10000,  algorithm='maxent.jar',tune.args =tune.args)
#enmeval_results@results 
#find the model that has the lowest delta AICc
write.csv(enmeval_results@results, "ENMeval_results.csv", row.names = F)                        
                         
####Maxent##
sp.modelt_bs <- maxent(
  x=predictors2,
  p=sp.occ,
  J = TRUE,
  path=paste0("D:/PhD/LP/LPSDM/Output/Current_maxent"),
  burnin=5000,
  nbg=10000,
  args=c(
    'betamultiplier=1',
    'noautofeature',
    'linear',
    'quadratic',
    'hinge',
    'product',
    'threshold',
    'writebackgroundpredictions=true',
    'doclamp',
    'threads=3',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=true',
    'removeduplicates=true',
    'replicatetype=bootstrap',
    'outputformat=Cloglog',
    'replicates=6',
    'maximumiterations=500000'
  )
)

### TRUE SKILL STATISTIC TEST FUNCTION #####
#to be able to run this script you need to have told the Maxent model to produce background predictions. If you are running MaxEnt in R this means putting the argument (after "args") "writebackgroundpredictions=true" as true not false. 

#FUNCTION CODE
TSS_calculations <- function (sample_clog, prediction_clog, n, th) {
  
  xx <- sum(sample_clog > th)
  yy <- sum(prediction_clog > th)
  xxx <- sum(sample_clog < th)
  yyy <- sum(prediction_clog < th)
  
  ncount <- sum(xx,yy,xxx,yyy)
  
  overallaccuracy <- (xx + yyy)/ncount 
  sensitivity <- xx / (xx + xxx)
  specificity <- yyy / (yy + yyy)
  tss <- sensitivity + specificity - 1
  
  #kappa calculations
  a <- xx + xxx
  b <- xx + yy
  c <- yy + yyy
  d <- xxx + yyy
  e <- a * b
  f <- c * d
  g <- e + f
  h <- g / (ncount * ncount)
  hup <- overallaccuracy - h
  hdown <- 1 - h
  
  kappa <- hup/hdown
  Po <- (xx + yyy) / ncount
  Pe <- ((b/ncount) * (a/ncount)) + ((d/ncount) * (c/ncount))
  Px1 <- Po - Pe
  Px2 <- 1 - Pe
  Px3 <- Px1/Px2
  
  tx1 <- xx + yyy
  tx2 <- 2 * a * c
  tx3 <- a - c
  tx4 <- xx - yyy
  tx5 <- ncount * (( 2 * a ) - tx4)
  tx6 <- ncount * tx1
  
  kappamax <- (tx6 - tx2 - (tx3 * tx4)) / ((tx5 - tx3) - (tx3 * tx4))
  
  return(tss) #returns the TSS value alone to an object file
  
  cat(" Maxent results for model with/n",a,"training sample predictions/n",c ,"background predictions/n/n TSS value:        ", tss,"/n Overall accuracy: ",overallaccuracy,"/n Sensitivity:      ",sensitivity,"/n Specificity:      ",specificity,"/n Kappa:            ",kappa,"/n Kappa max:        ",kappamax)
  
}


# Compile AUC and TSS results

setwd("D:/PhD/LP/LPSDM/Output/Current_maxent")
for (i in 0:5) { 
  
  #read in the background predictions
  
  backgroundPredictions <- read.csv(paste0("species_", i,"_backgroundPredictions.csv"))
  
  
  #we need the last column so will set the number as x
  e <- length(backgroundPredictions)
  
  #extract the cloglog/logistic results
  backgroundclog <-backgroundPredictions[,e]
  
  #now read in the sample predictions for testing
  samplepredictions <- read.csv(paste0("species_", i,"_samplePredictions.csv"))
  
  #we need the last column again of logistic or cloglog predictions so set a second x
  f <- length(samplepredictions)
  
  #extract the cloglog/logistic results for sample
  sampleclog <- samplepredictions[,f]
  
  #set n the number of pseudoabsences used for background predictions by MaxEnt
  n <- 10000
  
  #read in maxent results
  maxres <- read.csv(paste0("maxentResults.csv"))
  
  
  ### Set the threshold rule here
  th <- maxres[maxres[,1]=="species","Maximum.training.sensitivity.plus.specificity.Cloglog.threshold"]
  
  #run the function, the input values are the sampleclog values, then the background clog values, the sample number for the pseudo absences and then threshold value
  ts <- TSS_calculations(sampleclog,backgroundclog,n,th)
  
  if (i==0) {
    tss.df <- ts
  } else {
    tss.df <- rbind(tss.df, ts)
  }
  
}

tss.df <- as.data.frame(tss.df)
tss.df[7,1] <- summarise(tss.df, avg = mean(V1))

results.df <- cbind(maxres[,1], maxres[,6], tss.df)
colnames(results.df) <- c("Replicate", "Training.AUC", "TSS")
results.df

write.csv(results.df, file = "Maxent_AUC_TSS_results_MSS.csv", row.names = F)

# predict to entire dataset
setwd("D:/PhD/LP/LPSDM/Output")
sp.name <- "Current_MAP"
sp.pred.map <- predict(object=sp.modelt_bs,
                       x=predictors2
                       ,
                       filename=paste0(sp.name,"/",sp.name, "_R"),
                       na.rm=TRUE,
                       format='GTiff',
                       overwrite=TRUE,
                       doclamp= TRUE
)

sp.pred.map.ave <- calc(sp.pred.map,
                        fun = mean,
                        na.rm = T,
                        filename=paste0(sp.name,"/",sp.name,"_ave"),
                        format='GTiff',
                        overwrite=T)

#create binary maps
#extract the minimum training presence threshold
d <- as.data.frame(sp.modelt_bs@results)
a <- d["Maximum training sensitivity plus specificity Cloglog threshold", "species (average)"]

a <- d["Maximum.training.sensitivity.plus.specificity.Cloglog.threshold", "species (average)"]
#create matrix for reclassification
n <- c(0,a,0,a,1,1)
binmat <- matrix(n, ncol=3, byrow=T)
binmat
#reclassify current map

sp.pred.map.bin <- reclassify(x=sp.pred.map.ave, 
                              rcl=binmat, 
                              filename=paste0(sp.name,"/",sp.name, "_currentEnv_bin_ave"), 
                              format="GTiff",
                              overwrite=T)
values(sp.pred.map.bin)
#windows()
par(mfrow=c(1,2))
plot(sp.pred.map.ave, main="Predicted Suitability")
points(sp.occ$Longitude, sp.occ$Latitude, pch="+", cex=0.5)
plot(sp.pred.map.bin, main="Binary", legend=F)
#points(sp.occ$lon, sp.occ$lat, pch="+", cex=0.5)

### READ INPUTS ####

#read in predictor layers ##layer.3 is clay
#current Env

modelclimEnv <- predictors2
names(modelclimEnv)
ocean.ph <- dropLayer(modelclimEnv, c("bio7", "bio11", "bio15" , "bio16", "bio19","layer" ))
plot(ocean.ph)
predictors2
bio_now

##########################################################LGM Env layers
##STACK ice sheet######

###
######Read Ocean tif files#####
setwd("D:/PhD/LP/LPSDM/environment/LGM_ocean")
ocean_LGM <- list.files(pattern = ".tif$")
ocean_LGM <- raster::stack(ocean_LGM)
ocean_LGM
names(ocean_LGM) <- c("ocean1", "ocean3", "ocean4", "ocean5","ocean6","ocean7")
ocean_LGM<-aggregate(ocean_LGM, fact = 10)

#LGM Env layers
#########CCSM4_LGM#####

setwd("D:/PhD/LP/LPSDM/environment/LGM_BIO_Chelsa")
LGM_CCSM4 <- brick("LGM_CCSM4_envstack.tif")
names(LGM_CCSM4) <- c("CCSM4_bio07", "CCSM4_bio11", "CCSM4_bio15", "CCSM4_bio16", "CCSM4_bio19")

LGM_CCSM4 <- crop(LGM_CCSM4, bio_now )
LGM_CCSM4
LGM_CCSM4 <- mask(LGM_CCSM4, ocean_LGM$ocean1) 

LGM_CCSM4 <- stack(LGM_CCSM4, ocean_LGM)
LGM_CCSM4 

##STACK ice sheet######

setwd("D:/PhD/LP/LPSDM/environment/LGM_icecap")
icesheet <- list.files(pattern = ".tif$")
icesheet <- stack(icesheet)
plot(icesheet)

LGM_CCSM4 <- stack(LGM_CCSM4, icesheet)
LGM_CCSM4
names(LGM_CCSM4) <- c("bio7", "bio11", "bio15" , "bio16", "bio19","ocean1", "ocean3", "ocean4", "ocean5","ocean6","ocean7","distance") 
names(LGM_CCSM4) <- names(modelclimEnv)
plot()
plot(LGM_CCSM4$layer)
plot(bio_LGM$bio7)
LGM_CCSM4
predictors2

#########MIROC_LGM#####
setwd("D:/PhD/LP/LPSDM/environment/LGM_BIO_Chelsa")
LGM_MIROC <- brick("LGM_MIROC_envstack.tif")
names(LGM_MIROC) <- c("MIROC_bio07", "MIROC_bio11", "MIROC_bio15", "MIROC_bio16","MIROC_bio19")
LGM_MIROC
LGM_MIROC <- crop(LGM_MIROC, bio_now )

LGM_MIROC <- mask(LGM_MIROC, ocean_LGM$ocean1) 

LGM_MIROC <- stack(LGM_MIROC, ocean_LGM)

LGM_MIROC <- stack(LGM_MIROC, icesheet)
names(LGM_MIROC) <- names(modelclimEnv)



### 
#predict into LGM climate
setwd("D:/PhD/LP/LPSDM/Output2")
sp.name <- "LGM_MAP"
LGM.mod4a <- predict(object=sp.modelt_bs, 
                     x=LGM_CCSM4, 
                     filename=paste0(sp.name,"/",sp.name,"_MaxEnt_LGM_CCSM4"),
                     na.rm=TRUE,
                     format='GTiff',
                     overwrite=TRUE,
                     doclamp= TRUE)
LGM.mod4a.ave <- calc(LGM.mod4a,
                      fun = mean,
                      na.rm = T,
                      filename=paste0(sp.name,"/",sp.name,"_MaxEnt_LGM_CCSM4_ave"),
                      format='GTiff',
                      overwrite=T)

LGM.mod5a <- predict(object=sp.modelt_bs, 
                     x=LGM_MIROC,
                     filename=paste0(sp.name,"/",sp.name,"_MaxEnt_LGM_MIROC"),
                     na.rm=TRUE,
                     format='GTiff',
                     overwrite=TRUE,
                     doclamp= TRUE)
LGM.mod5a.ave <- calc(LGM.mod5a,
                      fun = mean,
                      na.rm = T,
                      filename=paste0(sp.name,"/",sp.name,"_MaxEnt_LGM_MIROC_ave"),
                      format='GTiff',
                      overwrite=T)

#LGM
LGM.ave <- (LGM.mod4a.ave + LGM.mod5a.ave) / 2

#write it out
#LGM
writeRaster(LGM.ave, filename=paste0(sp.name,"/",sp.name,"_LGM_ensemble"), format='GTiff', overwrite=TRUE)

#convert to binary
#averages

LGM.ave.bin <- reclassify(x=LGM.ave, 
                          rcl=binmat, 
                          filename=paste0(sp.name,"/",sp.name,"_LGM_ensem_bin"),                               format="GTiff",
                          overwrite=T)
#single models
#CCSM4
LGM.mod4a.bin <- reclassify(x=LGM.mod4a.ave, 
                            rcl=binmat, 
                            filename=paste0(sp.name,"/",sp.name,"_CCSM4_ave_bin"),                                 format="GTiff",
                            overwrite=T)

#MIROC
LGM.mod5a.bin <- reclassify(x=LGM.mod5a.ave, 
                            rcl=binmat, 
                            filename=paste0(sp.name,"/",sp.name,"_MIROC_ave_bin"),                                 format="GTiff",
                            overwrite=T)

sp.name <- " "

par(mfrow=c(2,2))
plot(sp.pred.map.ave, main = paste0(sp.name, " CURRENT"))
plot(LGM.ave, main = paste0(sp.name, " COMBINE_LGM"))
plot(LGM.mod4a.ave, main = paste0(sp.name, " CCSM4"))
plot(LGM.mod5a.ave, main = paste0(sp.name, " MIROC"))

par(mfrow=c(2,2))
plot(sp.pred.map.bin, main = paste0(sp.name, " CURRENT_Binary"))
plot(LGM.ave.bin, main = paste0(sp.name, " COMBINE_LGM_Binary"))
plot(LGM.mod4a.bin, main = paste0(sp.name, " CCSM4_Binary"))
plot(LGM.mod5a.bin, main = paste0(sp.name, " MIROC_Binary"))

plot(LGM.mod5a.bin, main = paste0(sp.name, " MIROC_Binary"))
plot(wa, add=TRUE)


