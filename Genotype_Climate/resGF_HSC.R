## Load library
library(vcfR)
library(resGF)
library(raster)
install.packages("resGF")
## Obtain some environmental data

setwd("D:/PhD/LP/LP_SNP_Cliamte/resGF")
data <-  load("LP_resGF.Rdata")
hsc.ext <- extent(-105.0001,-45.00014,-0.0001392488,54.99986)

files_0026 <- crop(hsc_2100_26.crop, hsc.ext)
files_0085 <- crop(hsc_2100_85.crop, hsc.ext)
files_5026 <- crop(hsc_2050_26.crop, hsc.ext)
files_5085 <- crop(hsc_2050_85.crop, hsc.ext)

names(files_5026) <- c("BO21_chlomean_ss_lonlat", "BO21_chlorange_ss_lonlat", "BO21_curvelmean_ss_lonlat", "BO21_curvelrange_ss_lonlat", "BO21_curvelmean_bdmin_lonlat", "BO21_curvelrange_bdmin_lonlat", "BO21_salinitymean_ss_lonlat" , "BO21_salinityrange_ss_lonlat", "BO21_salinitymean_bdmin_lonlat","BO21_salinityrange_bdmin_lonlat", "BO21_tempmean_ss_lonlat","BO21_temprange_ss_lonlat","BO21_tempmean_bdmin_lonlat" , "BO21_temprange_bdmin_lonlat")
names(files_0026) <- c("BO21_chlomean_ss_lonlat", "BO21_chlorange_ss_lonlat", "BO21_curvelmean_ss_lonlat", "BO21_curvelrange_ss_lonlat", "BO21_curvelmean_bdmin_lonlat", "BO21_curvelrange_bdmin_lonlat", "BO21_salinitymean_ss_lonlat" , "BO21_salinityrange_ss_lonlat", "BO21_salinitymean_bdmin_lonlat","BO21_salinityrange_bdmin_lonlat", "BO21_tempmean_ss_lonlat","BO21_temprange_ss_lonlat","BO21_tempmean_bdmin_lonlat" , "BO21_temprange_bdmin_lonlat")
names(files_5085) <- c("BO21_chlomean_ss_lonlat", "BO21_chlorange_ss_lonlat", "BO21_curvelmean_ss_lonlat", "BO21_curvelrange_ss_lonlat", "BO21_curvelmean_bdmin_lonlat", "BO21_curvelrange_bdmin_lonlat", "BO21_salinitymean_ss_lonlat" , "BO21_salinityrange_ss_lonlat", "BO21_salinitymean_bdmin_lonlat","BO21_salinityrange_bdmin_lonlat", "BO21_tempmean_ss_lonlat","BO21_temprange_ss_lonlat","BO21_tempmean_bdmin_lonlat" , "BO21_temprange_bdmin_lonlat")
names(files_0085) <- c("BO21_chlomean_ss_lonlat", "BO21_chlorange_ss_lonlat", "BO21_curvelmean_ss_lonlat", "BO21_curvelrange_ss_lonlat", "BO21_curvelmean_bdmin_lonlat", "BO21_curvelrange_bdmin_lonlat", "BO21_salinitymean_ss_lonlat" , "BO21_salinityrange_ss_lonlat", "BO21_salinitymean_bdmin_lonlat","BO21_salinityrange_bdmin_lonlat", "BO21_tempmean_ss_lonlat","BO21_temprange_ss_lonlat","BO21_tempmean_bdmin_lonlat" , "BO21_temprange_bdmin_lonlat")


nows <- stack(hsc.crop,ocean_now) 
e5026s <- stack(files_5026,ocean_now)
e0026s <- stack(files_0026,ocean_now)
e5085s <- stack(files_5085,ocean_now)
e0085s <- stack(files_0085,ocean_now)
head(e5026s)
plot(nows)
plot(e5026s)
library(RColorBrewer)
pal <- colorRampPalette(c("blue","red"))

##########################################################
## create a genind object
CRv <- read.vcfR("LP159_mis09LD.vcf")
CR = vcfR2genind(CRv)

## obtain spatial data
CR.coord  <-read.table("LP.coordinates.txt", header=T, stringsAsFactors=F, sep="\t", row.names=1)
colnames(CR.coord) <- c("x", "y")
## create a resistance surface using resistantGF
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
CR.coord.sp <- SpatialPointsDataFrame(CR.coord[,c('y','x')], proj4string=CRS(crs.wgs), data=CR.coord)
head(CR.coord.sp)


#install.packages("rSDM")
library(rSDM)
CR.coord.sp <- points2nearestcell(CR.coord.sp, nows, layer=19)


# Gradient Forest
library(gradientForest)
# extract points from
CR.clim.points <- raster::extract(nows, CR.coord.sp)
CR.clim.points[is.na(CR.clim.points)] <- 0 #change three NA values of bathymetry to 0
#generates the PCNMs
library(vegan)
CR.pcnm <- pcnm(dist(CR.coord.sp@data))
CR.keep <- round(length(which(CR.pcnm$value > 0))/2)
CR.pcnm.keep <- scores(CR.pcnm)[,1:CR.keep]  #keep half of positive ones as suggested by some authors

CR.nbvar <- length(colnames(CR.clim.points))
CR.env.gf <- cbind(CR.clim.points[ , c(seq(1, CR.nbvar, by=1)) ], CR.pcnm.keep)
CR.maxLevel <- log2(0.368*nrow(CR.env.gf)/2)
CR.env.gf <- as.data.frame(CR.env.gf)
CR.snp <- as.data.frame(CR$tab)

# run gradient forest
colnames(CR.snp) <- paste("A",1:length(colnames(CR.snp)), sep = "")
rownames(CR.snp) <- rownames(CR.env.gf)
CR.gf <- gradientForest(cbind(CR.env.gf, CR.snp), predictor.vars=colnames(CR.env.gf), response.vars=colnames(CR.snp), ntree=500, maxLevel=CR.maxLevel, trace=T, corr.threshold=0.50, compact = T, nbin = 201)

# generate resistance surface
CR.single_r <- resGF(CR.gf, nows, save.image = T)
dev.off()

pdf("LP_final.pdf")
plot(CR.single_r, col = pal(15))
dev.off()


CR.single_r <- resGF(CR.gf, e5026s, save.image = F)
dev.off()

pdf("LP_5026e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

CR.single_r <- resGF(CR.gf, e0026s, save.image = F)
dev.off()
pdf("LP_0026e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()


CR.single_r <- resGF(CR.gf, e5085s, save.image = F)
dev.off()
pdf("LP_5085e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()


CR.single_r <- resGF(CR.gf, e0085s, save.image = F)
dev.off()
pdf("LP_0085e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()


###45
files_5045 <- crop(hsc_2050_45.crop, hsc.ext)
files_0045 <- crop(hsc_2100_45.crop, hsc.ext)

names(files_5045) <- c("BO21_chlomean_ss_lonlat", "BO21_chlorange_ss_lonlat", "BO21_curvelmean_ss_lonlat", "BO21_curvelrange_ss_lonlat", "BO21_curvelmean_bdmin_lonlat", "BO21_curvelrange_bdmin_lonlat", "BO21_salinitymean_ss_lonlat" , "BO21_salinityrange_ss_lonlat", "BO21_salinitymean_bdmin_lonlat","BO21_salinityrange_bdmin_lonlat", "BO21_tempmean_ss_lonlat","BO21_temprange_ss_lonlat","BO21_tempmean_bdmin_lonlat" , "BO21_temprange_bdmin_lonlat")
names(files_0045) <- c("BO21_chlomean_ss_lonlat", "BO21_chlorange_ss_lonlat", "BO21_curvelmean_ss_lonlat", "BO21_curvelrange_ss_lonlat", "BO21_curvelmean_bdmin_lonlat", "BO21_curvelrange_bdmin_lonlat", "BO21_salinitymean_ss_lonlat" , "BO21_salinityrange_ss_lonlat", "BO21_salinitymean_bdmin_lonlat","BO21_salinityrange_bdmin_lonlat", "BO21_tempmean_ss_lonlat","BO21_temprange_ss_lonlat","BO21_tempmean_bdmin_lonlat" , "BO21_temprange_bdmin_lonlat")

e0045s <- stack(files_0045,ocean_now)
e5045s <- stack(files_5045,ocean_now)
head(CR.gf)
CR.single_r <- resGF(CR.gf, e5045s, save.image = F)
dev.off()
pdf("LP_5045e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

CR.single_r <- resGF(CR.gf, e0045s, save.image = F)
dev.off()
pdf("LP_0045e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()
