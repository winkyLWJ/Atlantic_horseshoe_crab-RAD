library(raster)
library(RColorBrewer)

### LP
setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Output/Current_MAP")
LP <- raster("Current_MAP_ave.tif")
plot(LP)

setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Future/Output/LP_2100.85")
LP_2100_RCP85 <- raster("LP_2100.85_fut.2100.85_ensemble.tif")
plot(LP_2100_RCP85)
LP_diff <- LP - LP_2100_RCP85

setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Diff")
writeRaster(LP_diff, "LP_diff_2100_RCP85.asc", "ascii")

pdf("LP_diff_2100_RCP85.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(LP_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()
plot(LP_diff)

#######RCP26
setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Future/Output/LP_2100.26")
LP_2100_RCP26 <- raster("LP_2100.26_fut.2100.26_ensemble.tif")
LP_diff <- LP - LP_2100_RCP26

setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Diff")
writeRaster(LP_diff, "LP_diff_2100_RCP26.asc", "ascii", overwrite=TRUE)

pdf("LP_diff_2100_RCP26.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(LP_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()
plot(LP_diff)


##RCP45
setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Future/Output/LP_2100.45")
LP_2100_RCP45 <- raster("LP_2100.45_fut.2100.45_ensemble.tif")
plot(LP_2100_RCP45)
LP_diff <- LP - LP_2100_RCP45

setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Diff")
writeRaster(LP_diff, "LP_diff_2100_RCP45.asc", "ascii")

pdf("LP_diff_2100_RCP45.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(LP_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()
plot(LP_diff)


#######RCP26
setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Future/Output/LP_2050.26")
LP_2050_RCP26 <- raster("LP_2050.26_fut.2050.26_ensemble.tif")
LP_diff <- LP - LP_2050_RCP26

setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Diff")
writeRaster(LP_diff, "LP_diff_2050_RCP26.asc", "ascii", overwrite=TRUE)

pdf("LP_diff_2050_RCP26.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(LP_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()
plot(LP_diff)

#######RCP45
setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Future/Output/LP_2050.45")
LP_2050_RCP45 <- raster("LP_2050.45_fut_ensemble.tif")
LP_diff <- LP - LP_2050_RCP45

setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Diff")
writeRaster(LP_diff, "LP_diff_2050_RCP45.asc", "ascii", overwrite=TRUE)

pdf("LP_diff_2050_RCP45.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(LP_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()
plot(LP_diff)


#######RCP85
setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Future/Output/LP_2050.85")
LP_2050_RCP85 <- raster("LP_2050.85_fut_ensemble.tif")
LP_diff <- LP - LP_2050_RCP85

setwd("D:/PhD/LP/LP_SNP_Cliamte/ENM/Diff")
writeRaster(LP_diff, "LP_diff_2050_RCP85.asc", "ascii", overwrite=TRUE)

pdf("LP_diff_2050_RCP85.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(LP_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()
plot(LP_diff)
