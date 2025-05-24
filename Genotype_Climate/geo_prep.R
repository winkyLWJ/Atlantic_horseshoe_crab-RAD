library(sdmpredictors)
library(raster)
getOption('timeout')
options(timeout = 600000000)


## download raster layers for present day
hsc <- load_layers(c(
  "BO21_chlomean_ss",
  "BO21_chlorange_ss",
  "BO21_curvelmean_ss",
  "BO21_curvelrange_ss",
  "BO21_curvelmean_bdmin",
  "BO21_curvelrange_bdmin",
  "BO21_salinitymean_ss",
  "BO21_salinityrange_ss",
  "BO21_salinitymean_bdmin",
  "BO21_salinityrange_bdmin",
  "BO21_tempmean_ss",
  "BO21_temprange_ss",
  "BO21_tempmean_bdmin",
  "BO21_temprange_bdmin"
))

hsc.ext <- extent(-105.0001,-45.00014,-0.0001392488,54.99986)
hsc.crop <- crop(hsc, hsc.ext)
hsc_current <- as.data.frame(hsc.crop, xy=TRUE)
write.csv(hsc_current, "current_bio.csv")

## extact env for LP
setwd("D:/PhD/LP/LPSDM")
LP_loc <- read.csv("LP_sampling_correct_removerep.csv")
xy <- LP_loc[, c(2,3)]
spdf <- SpatialPointsDataFrame(coords = xy, data = LP_loc, proj4string = crs(hsc))
temp <- extract(hsc.crop, spdf)

write.csv(temp, "LP_env.csv")


## download raster layers for 2050 RCP26
hsc_2050_26 <- load_layers(c(
  "BO21_RCP26_2050_chlomean_ss",
  "BO21_RCP26_2050_chlorange_ss",
  "BO21_RCP26_2050_curvelmean_ss",
  "BO21_RCP26_2050_curvelrange_ss",
  "BO21_RCP26_2050_curvelmean_bdmin",
  "BO21_RCP26_2050_curvelrange_bdmin",
  "BO21_RCP26_2050_salinitymean_ss",
  "BO21_RCP26_2050_salinityrange_ss",
  "BO21_RCP26_2050_salinitymean_bdmin",
  "BO21_RCP26_2050_salinityrange_bdmin",
  "BO21_RCP26_2050_tempmean_ss",
  "BO21_RCP26_2050_temprange_ss",
  "BO21_RCP26_2050_tempmean_bdmin",
  "BO21_RCP26_2050_temprange_bdmin"
))

hsc_2050_26.crop <- crop(hsc_2050_26, hsc.ext)
hsc_5026 <- as.data.frame(hsc_2050_26.crop, xy=TRUE)
write.csv(hsc_5026, "fut_5026.csv")

## download raster layers for 2050 RCP45
hsc_2050_45 <- load_layers(c(
  "BO21_RCP45_2050_chlomean_ss",
  "BO21_RCP45_2050_chlorange_ss",
  "BO21_RCP45_2050_curvelmean_ss",
  "BO21_RCP45_2050_curvelrange_ss",
  "BO21_RCP45_2050_curvelmean_bdmin",
  "BO21_RCP45_2050_curvelrange_bdmin",
  "BO21_RCP45_2050_salinitymean_ss",
  "BO21_RCP45_2050_salinityrange_ss",
  "BO21_RCP45_2050_salinitymean_bdmin",
  "BO21_RCP45_2050_salinityrange_bdmin",
  "BO21_RCP45_2050_tempmean_ss",
  "BO21_RCP45_2050_temprange_ss",
  "BO21_RCP45_2050_tempmean_bdmin",
  "BO21_RCP45_2050_temprange_bdmin"
))

hsc_2050_45.crop <- crop(hsc_2050_45, hsc.ext)
hsc_5045 <- as.data.frame(hsc_2050_45.crop, xy=TRUE)
write.csv(hsc_5045, "fut_5045.csv")

## download raster layers for 2050 RCP85
hsc_2050_85 <- load_layers(c(
  "BO21_RCP85_2050_chlomean_ss",
  "BO21_RCP85_2050_chlorange_ss",
  "BO21_RCP85_2050_curvelmean_ss",
  "BO21_RCP85_2050_curvelrange_ss",
  "BO21_RCP85_2050_curvelmean_bdmin",
  "BO21_RCP85_2050_curvelrange_bdmin",
  "BO21_RCP85_2050_salinitymean_ss",
  "BO21_RCP85_2050_salinityrange_ss",
  "BO21_RCP85_2050_salinitymean_bdmin",
  "BO21_RCP85_2050_salinityrange_bdmin",
  "BO21_RCP85_2050_tempmean_ss",
  "BO21_RCP85_2050_temprange_ss",
  "BO21_RCP85_2050_tempmean_bdmin",
  "BO21_RCP85_2050_temprange_bdmin"
))

hsc_2050_85.crop <- crop(hsc_2050_85, hsc.ext)
hsc_5085 <- as.data.frame(hsc_2050_85.crop, xy=TRUE)
write.csv(hsc_5085, "fut_5085.csv")


## download raster layers for 2100 RCP26
hsc_2100_26 <- load_layers(c(
  "BO21_RCP26_2100_chlomean_ss",
  "BO21_RCP26_2100_chlorange_ss",
  "BO21_RCP26_2100_curvelmean_ss",
  "BO21_RCP26_2100_curvelrange_ss",
  "BO21_RCP26_2100_curvelmean_bdmin",
  "BO21_RCP26_2100_curvelrange_bdmin",
  "BO21_RCP26_2100_salinitymean_ss",
  "BO21_RCP26_2100_salinityrange_ss",
  "BO21_RCP26_2100_salinitymean_bdmin",
  "BO21_RCP26_2100_salinityrange_bdmin",
  "BO21_RCP26_2100_tempmean_ss",
  "BO21_RCP26_2100_temprange_ss",
  "BO21_RCP26_2100_tempmean_bdmin",
  "BO21_RCP26_2100_temprange_bdmin"
))

hsc_2100_26.crop <- crop(hsc_2100_26, hsc.ext)
hsc_0026 <- as.data.frame(hsc_2100_26.crop, xy=TRUE)
write.csv(hsc_0026, "fut_0026.csv")

## download raster layers for 2100 RCP45
hsc_2100_45 <- load_layers(c(
  "BO21_RCP45_2100_chlomean_ss",
  "BO21_RCP45_2100_chlorange_ss",
  "BO21_RCP45_2100_curvelmean_ss",
  "BO21_RCP45_2100_curvelrange_ss",
  "BO21_RCP45_2100_curvelmean_bdmin",
  "BO21_RCP45_2100_curvelrange_bdmin",
  "BO21_RCP45_2100_salinitymean_ss",
  "BO21_RCP45_2100_salinityrange_ss",
  "BO21_RCP45_2100_salinitymean_bdmin",
  "BO21_RCP45_2100_salinityrange_bdmin",
  "BO21_RCP45_2100_tempmean_ss",
  "BO21_RCP45_2100_temprange_ss",
  "BO21_RCP45_2100_tempmean_bdmin",
  "BO21_RCP45_2100_temprange_bdmin"
))

hsc_2100_45.crop <- crop(hsc_2100_45, hsc.ext)
hsc_0045 <- as.data.frame(hsc_2100_45.crop, xy=TRUE)
write.csv(hsc_0045, "fut_0045.csv")

## download raster layers for 2100 RCP85
hsc_2100_85 <- load_layers(c(
  "BO21_RCP85_2100_chlomean_ss",
  "BO21_RCP85_2100_chlorange_ss",
  "BO21_RCP85_2100_curvelmean_ss",
  "BO21_RCP85_2100_curvelrange_ss",
  "BO21_RCP85_2100_curvelmean_bdmin",
  "BO21_RCP85_2100_curvelrange_bdmin",
  "BO21_RCP85_2100_salinitymean_ss",
  "BO21_RCP85_2100_salinityrange_ss",
  "BO21_RCP85_2100_salinitymean_bdmin",
  "BO21_RCP85_2100_salinityrange_bdmin",
  "BO21_RCP85_2100_tempmean_ss",
  "BO21_RCP85_2100_temprange_ss",
  "BO21_RCP85_2100_tempmean_bdmin",
  "BO21_RCP85_2100_temprange_bdmin"
))

hsc_2100_85.crop <- crop(hsc_2100_85, hsc.ext)
hsc_0085 <- as.data.frame(hsc_2100_85.crop, xy=TRUE)
write.csv(hsc_0085, "fut_0085.csv")
