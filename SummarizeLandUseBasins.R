

#Script to calculate percent land use for GLRI watersheds

library(tidyr)
library(dplyr)
library(sp)
library(rgdal, verbose=F)
library(SDMTools)
library(raster,verbose = T)
library(FedData)
library(raster)
library(sp)
library(rgdal)
library(sf)

#Watersheds
# basinPath <- "M:/QW Monitoring Team/GLRI_GIS/layers/20131299_GLRI_MASTER/~old/201401_HUC12Basins"
# basins_OGR <- readOGR(basinPath, 'GLRI_Basins_HUC12Linework', verbose=F)

# basinPath <- "M:/QW Monitoring Team/GLRI_GIS/layers/20190730_Pesticides"
# basins_OGR <- readOGR(basinPath, 'PesticideBasins', verbose=F)
# basins_OGR <- filter(basins_OGR, STAID %in% sites$site_no)

# basins = st_read(basinPath, 'PesticideBasins')
# basins <- filter(basins, STAID %in% sites$site_no)

#nlcd
# big_nlcd = get_nlcd(template=GL, label="GL")
big_nlcd = get_nlcd(template=basins_OGR, label="basins", year='2011')
# NLCDpath <- c('O:/BaseData/National/LandCover/NLCD16_2019Version/nlcd16lc.ovr')
# big_nlcd = raster(NLCDpath)
# 
# big_nlcd<-spTransform(big_nlcd,  crs(big_nlcd))

basins_OGR <-spTransform(basins_OGR, crs(big_nlcd))

#States
state_names <- c("minnesota","wisconsin", "michigan","ohio",
                 "illinois","pennsylvania","new york","indiana")

GL <- us_states(resolution = "high", states = state_names) %>%
  st_transform(crs =crs(big_nlcd))


nlcd_masked <- mask(big_nlcd, basins_OGR, inverse=F)

#Save first color of raster
firstcolor<-nlcd_masked@legend@colortable[1]
nlcd_masked@legend@colortable[1]<-firstcolor

#Make a short legned color table
short_colors<-c('#476BA0', "#BAD8EA", "#68AA63", '#ED0000', "#DBD83D", "#CCBA7C")
short_labels<-c('Open Water', 'Wetland', 'Forest', 'Developed', 'Agriculture', 'Other')
short_legend<-data.frame(short_labels, short_colors, stringsAsFactors = F)

#Change colors of raster
nlcd_masked@legend@colortable[1]<-"#ffffff"
nlcd_masked@legend@colortable[1]<-NA
nlcd_masked@legend@colortable[11+1]<-short_colors[1]
nlcd_masked@legend@colortable[21:24+1]<-short_colors[4]
nlcd_masked@legend@colortable[41:43+1]<-short_colors[3]
nlcd_masked@legend@colortable[c(31,52,71)+1]<-short_colors[6]
nlcd_masked@legend@colortable[81:82+1]<-short_colors[5]
nlcd_masked@legend@colortable[c(90,95)+1]<-short_colors[2]

# plot
plot(nlcd_masked)
plot(basins_OGR, add=T, col=NA, border='black', lwd=2)
plot(GL, add=T, border='grey', col=NA, lwd=1)





#Make legend for plotting
classes<-unique(as.character(nlcd_masked@data@attributes[[1]]$NLCD.2011.Land.Cover.Class))
classes<-classes[-which(classes %in% c('Unclassified', '', "Perennial Snow/Ice"))]

IDs<-nlcd_masked@data@attributes[[1]]$ID[match(classes, nlcd_masked@data@attributes[[1]]$NLCD.2011.Land.Cover.Class)]

colors<-nlcd_masked@legend@colortable[IDs+1]

legendtable<-data.frame(IDs, classes, colors, stringsAsFactors=F)



#Generate percent of each land cover
basins_OGR@data$perWater <- rep(NA, nrow(basins_OGR@data))
basins_OGR@data$perForest <- rep(NA, nrow(basins_OGR@data))
basins_OGR@data$perWetland <- rep(NA, nrow(basins_OGR@data))
basins_OGR@data$perUrban <- rep(NA, nrow(basins_OGR@data))
basins_OGR@data$perAg <- rep(NA, nrow(basins_OGR@data))
basins_OGR@data$perOther <- rep(NA, nrow(basins_OGR@data))

watershed_i <- 1
for (watershed_i in 1:length(basins_OGR)){
  nlcd_masked_i <-  mask(nlcd_masked, basins_OGR[watershed_i,], inverse=F)
  
  Vals_nlcdmasked = as.data.frame(freq(nlcd_masked_i))
  
  # Total counts
  total = Vals_nlcdmasked %>% dplyr::filter(value >= 1) %>%
    summarise_at('count',sum)
  
  # Total forest counts
  forest = Vals_nlcdmasked %>% dplyr::filter(value >= 41 & value <= 43) %>%
    summarise_at('count',sum)
  
  # Total wetland counts
  wetland = Vals_nlcdmasked %>% dplyr::filter(value >= 90 & value <= 95) %>%
    summarise_at('count',sum)
  
  # Total ag counts
  ag = Vals_nlcdmasked %>% dplyr::filter(value >= 81 & value <= 82) %>%
    summarise_at('count',sum)
  
  # Total urban counts
  urban = Vals_nlcdmasked %>% dplyr::filter(value >= 21 & value <= 24) %>%
    summarise_at('count',sum)
  
  # Total water counts
  water = Vals_nlcdmasked %>% dplyr::filter(value == 11) %>%
    summarise_at('count',sum)
  
  basins_OGR@data$perWater[watershed_i] = water/total
  basins_OGR@data$perForest[watershed_i] = forest/total
  basins_OGR@data$perWetland[watershed_i] = wetland/total
  basins_OGR@data$perAg[watershed_i] = ag/total
  basins_OGR@data$perUrban[watershed_i] = urban/total
  basins_OGR@data$perOther[watershed_i] = 1-c(water+forest+wetland+urban+ag)/total
  
}

basins_OGR@data$perWater <- unlist(basins_OGR@data$perWater)
basins_OGR@data$perForest <- unlist(basins_OGR@data$perForest)
basins_OGR@data$perWetlan <- unlist(basins_OGR@data$perWetland)
basins_OGR@data$perAg <- unlist(basins_OGR@data$perAg)
basins_OGR@data$perUrban <- unlist(basins_OGR@data$perUrban)
basins_OGR@data$perOther <- unlist(basins_OGR@data$perOther)


saveRDS(basins_OGR, file_out(file.path(path_to_data, "Rdata/watershedmetrics.rds")))
