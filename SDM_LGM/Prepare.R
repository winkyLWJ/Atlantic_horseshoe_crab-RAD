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
setwd("/Users/LP_SDM")
############################################################Prepare species points##########
## Species Points ##
sppts <- read.csv("LP_sample.csv")
summary(sppts)

colnames(sppts) <- c("gibfID", "Lat", "Lon")

# Species points cleaning ####

#rarefy

#new centering code

gridsize <- 0.008333333 # value for half a degree = 0.5; value for 10 minutes: 0.167; 5 mins is 0.083; 30 secs is 0.0083

z <- sppts$Lon
sppts$longmid <- ifelse(z >= 0.0,
                        (z%/%1)+(((z-z%/%1)/gridsize)%/%1)*gridsize,
                        (((z*(-1.0))%/%1)+((((z*(-1.0))-(z*(-1.0))%/%1)/gridsize)%/%1)*gridsize)*(-1.0))


z <- sppts$Lat
sppts$latmid <- ifelse(z >= 0.0, (z%/%1)+(((z-z%/%1)/gridsize)%/%1)*gridsize,
                       (((z*(-1.0))%/%1)+((((z*(-1.0))-(z*(-1.0))%/%1)/gridsize)%/%1)*gridsize)*(-1.0))


#now remove repeats
sppdist <- sppts[c("gibfID","longmid","latmid")]

sppdist2 <- distinct(sppdist, gibfID, longmid, latmid, .keep_all = TRUE)
summary(sppdist2)

write.csv(sppdist2, file = "LP_rmrepeat.csv", row.names = F)


###Exclude those points in land using qgis####

########read back species#######

sppts2 <- read.csv("LP_onlysea.csv")
sppdist2 <- dplyr::select(sppts2, longmid, latmid)
sp.occ <- sppdist2
sp.occ <-as.data.frame(sp.occ)

sp.occ

#############################################################crop environments with same extent############
#########crop bioclim ########

setwd("/Users/WJ//LP_SDM/Environment/Current_BIO_Chelsa")
bio_now <- list.files(pattern = ".tif$")
bio_now <- raster::stack(bio_now)
ext<-extent(c(-105, -45, 0, 55))
bio_now<-crop(bio_now,ext)

plot(bio_now$CHELSA_bio1)

bio_now
names(bio_now) <- c("bio1", "bio2", "bio3", "bio4", "bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
plot(bio_now$bio12)
names(bio_now) 



#write it out
writeRaster(bio_now, filename = "RAS2_now_bioclim_chelsa30s_all_stack.tif", overwrite = T) 

#read back in
nowEnvRAS <- brick("RAS2_now_bioclim_chelsa30s_all_stack.tif")
nowEnvRAS

### Data Processing ####
### climate ###

#check for multicollinearity, climate

nowEnvRAS.df <- as.data.frame(nowEnvRAS, na.rm = T)
climcor <- cor(nowEnvRAS.df, method = "spearman")

climcor[upper.tri(climcor)] <- 0
diag(climcor) <- 0

climdata.new <- nowEnvRAS.df[,!apply(climcor,2,function(x) any(x > 0.75))] #drop variables with r > 0.75
head(climdata.new)
write.csv(climcor, file = "now_spearmans_0.75_chelsa30s.csv")

bio_now
#keep the following parameters: "bio7","bio11","bio15","bio16","bio19"
names(nowEnvRAS) <- c("bio1", "bio2", "bio3", "bio4", "bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
currentEnvRAS3 <- dropLayer(nowEnvRAS, c("bio1", "bio2","bio3", "bio4","bio5", "bio6", "bio8", "bio9","bio10", "bio13","bio14","bio17","bio18"))
names(currentEnvRAS3) 
currentEnvRAS3



########################Bathymetry#################
install.packages("ncdf4")
library(ncdf4)
ext<-extent(c(-105, -45, 0, 55))
plot(ext)

bio_now


########################Bathymetry#################

bath_LGM<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/21kya_Geophysical_Data/bathy_5m/w001001.adf")
ext<-extent(c(-105, -45, 0, 55))
bath_LGM<-crop(bath_LGM,ext)
bath_LGM<-disaggregate(bath_LGM, fact = 10)
plot(bath_LGM, col = topo.colors(100))
bath_LGM
bath_now<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/ETOPO1_Bed_g_geotiff/ETOPO1_Bed_g_geotiff.tif")
bath_now<-crop(bath_now,ext)
bath_now
crs(bath_now)<-crs(bath_LGM)
projection(bath_now)<-projection(bath_LGM)
bath_now<-disaggregate(bath_now,fact=2)
bath_now[bath_now > 0] <- NA
plot(bath_now, col = topo.colors(100))
values(bath_now)
bio03_cli
writeRaster(bath_now, "bath_now.tif", overwrite=TRUE)

bath_now <- resample(bath_now, bio_now)
bath_LGM <- resample(bath_LGM, bio_now)

#####################Aspect#######################

Asp_EW<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/21kya_Geophysical_Data/biogeo01_5m/w001001.adf")
Asp_EW<-crop(Asp_EW,ext)
Asp_EW<-disaggregate(Asp_EW, fact = 10)
plot(Asp_EW, col = topo.colors(100))
Asp_EW

Asp_NS<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/21kya_Geophysical_Data/biogeo02_5m/w001001.adf")
Asp_NS<-crop(Asp_NS,ext)
plot(Asp_NS, col = topo.colors(100))
Asp_NS
Asp_NS<-disaggregate(Asp_NS, fact = 10)

Asp_now<-terrain(bath_now, opt="aspect", unit="radians", neighbors=8)
Asp_EW_now <- sin(Asp_now)*100
plot(Asp_EW_now, col = topo.colors(100))
Asp_NS_now <- cos(Asp_now)*100
plot(Asp_NS_now, col = topo.colors(100))

writeRaster(Asp_EW, "Asp_EW_LGM.tif",  overwrite=TRUE)
writeRaster(Asp_EW_now, "Asp_EW_now.tif",  overwrite=TRUE)
writeRaster(Asp_NS, "Asp_NS_LGM.tif", overwrite=TRUE)
writeRaster(Asp_NS_now, "Asp_NS_now.tif", overwrite=TRUE)
Asp_EW_now
Asp_NS_now

Asp_EW_now <- resample(Asp_EW_now, bio_now)
Asp_NS_now <- resample(Asp_NS_now, bio_now)
Asp_EW_LGM <- resample(Asp_EW_LGM, bio_now)
Asp_NS_LGM <- resample(Asp_NS_LGM, bio_now)
#####################Distance to Shore#######################

Dist_LGM<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/21kya_Geophysical_Data/biogeo05_5m/w001001.adf")
Dist_LGM<-crop(Dist_LGM,ext)
plot(Dist_LGM, col = topo.colors(100))
Dist_LGM
Dist_LGM<-disaggregate(Dist_LGM, fact = 10)


Dist_now<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/GMT_intermediate_coast_distance_01d/GMT_intermediate_coast_distance_01d.tif")
Dist_now<-crop(Dist_now,ext)
Dist_now
Dist_now<-disaggregate(Dist_now,fact=30)
Dist_now<-aggregate(Dist_now, fact=25)
Dist_now[Dist_now < 0] <- NA
plot(Dist_now, col = topo.colors(100))
Dist_now <- resample(Dist_now, bio_now)
Dist_LGM <- resample(Dist_LGM, bio_now)

writeRaster(Dist_LGM, "Dist_LGM.tif", overwrite=TRUE)
writeRaster(Dist_now, "Dist_now.tif",overwrite=TRUE)


#####################Bathymetric Slope#######################

slp_LGM<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/21kya_Geophysical_Data/biogeo06_5m/w001001.adf")
slp_LGM<-crop(slp_LGM,ext)
slp_LGM
slp_LGM<-disaggregate(slp_LGM, fact = 10)

plot(slp_LGM, col = topo.colors(100))

slp_now<-terrain(bath_now, opt="slope", unit="degrees", neighbors=8)
slp_now<-slp_now*10
plot(slp_now, col = topo.colors(100))
slp_now <- resample(slp_now, bio_now)
slp_LGM <- resample(slp_LGM, bio_now)

writeRaster(slp_LGM, "slp_LGM.tif",overwrite=TRUE)
writeRaster(slp_now, "slp_now.tif",overwrite=TRUE)
slp_now

#####################Mean salinity#######################
sal_LGM<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/21kya_HadCM/biogeo08_5m/w001001.adf")
sal_LGM<-crop(sal_LGM,ext)
sal_LGM
sal_LGM<-disaggregate(sal_LGM, fact = 10)
plot(sal_LGM, col = topo.colors(100))

library(ncdf4)
sal_now<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/Ocean/woa18_decav_s00_04.nc",  varname = "s_an")
sal_now<-crop(sal_now,ext)
sal_now<-disaggregate(sal_now, fact = 30)
sal_now
sal_now<-sal_now*100
plot(sal_now, col = topo.colors(100))
writeRaster(sal_LGM, "sal_LGM.tif",overwrite=TRUE)
writeRaster(sal_now, "sal_now.tif",overwrite=TRUE)
plot(sal_now)
sal_now <- resample(sal_now, bio_now)
sal_LGM <- resample(sal_LGM, bio_now)

#####################Mean surface temperature#######################

tmp_LGM<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/21kya_HadCM/biogeo13_5m/w001001.adf")
tmp_LGM<-crop(tmp_LGM,ext)
tmp_LGM
tmp_LGM<-disaggregate(tmp_LGM, fact = 10)
plot(tmp_LGM, col = topo.colors(100))

tmp_now<-raster("/Users/WJ//LP_SDM/TQ_SDM/Environment/Ocean/woa18_decav_t00_04.nc",  varname = "t_an")
tmp_now<-crop(tmp_now,ext)
tmp_now
tmp_now<-disaggregate(tmp_now, fact = 30)
tmp_now<-tmp_now*100
plot(tmp_now, col = topo.colors(100))

tmp_now <- resample(tmp_now, bio_now)
tmp_LGM <- resample(tmp_LGM, bio_now)

writeRaster(tmp_LGM, "tmp_LGM.tif",overwrite=TRUE)
writeRaster(tmp_now, "tmp_now.tif",overwrite=TRUE)

######Rean Ocean tif files#####

setwd("/Users/WJ//LP_SDM/Environment/Now_ocean")
ocean_now <- list.files(pattern = ".tif$")
ocean_now <- stack(ocean_now)
ocean_now
names(ocean_now) <- c("ocean1", "ocean2", "ocean3", "ocean4", "ocean5","ocean6","ocean7")
### Data Processing ####
### climate ###

#check for multicollinearity, climate

ocean_now.df <- as.data.frame(ocean_now, na.rm = T)
ocean_climcor <- cor(ocean_now.df, method = "spearman")

ocean_climcor[upper.tri(ocean_climcor)] <- 0
diag(ocean_climcor) <- 0

ocean_climdatanew <- ocean_now.df[,!apply(ocean_climcor,2,function(x) any(x > 0.75))] #drop variables with r > 0.75
head(ocean_climdatanew)
write.csv(ocean_climcor, file = "ocean_spearmans_0.75_chelsa30s.csv")

