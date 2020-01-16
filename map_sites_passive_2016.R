# devtools::install_github("tidyverse/ggplot2")
# devtools::install_github("r-spatial/sf")


library(httr)
library(ggplot2)
library(dplyr)
library(remake)
library(sf)
library(USAboundaries)

sites <- sites_2016 %>%
  rename(site_no = SiteID)

# get info about number of chemicals by site, merge with sites
nchems <- site_detection %>%
  group_by(site) %>%
  summarize(chemicals = sum(value))

graphData <- site_detection %>%
  filter(group=='AboveEAR0.001') %>%
  group_by(site) %>%
  select(site, value) %>%
  rename(EAR = value) %>%
  left_join(nchems) %>%
  data.frame()


get_flowlines <- function(streamorder, mapRange){
  postURL <- "https://cida.usgs.gov/nwc/geoserver/nhdplus/ows"
  
  filterXML <- paste0('<?xml version="1.0"?>',
                      '<wfs:GetFeature xmlns:wfs="http://www.opengis.net/wfs" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:gml="http://www.opengis.net/gml" service="WFS" version="1.1.0" outputFormat="shape-zip" xsi:schemaLocation="http://www.opengis.net/wfs http://schemas.opengis.net/wfs/1.1.0/wfs.xsd">',
                      '<wfs:Query xmlns:feature="https://gov.usgs.cida/nhdplus" typeName="feature:nhdflowline_network" srsName="EPSG:4326">',
                      '<ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">',
                      '<ogc:And>',
                      '<ogc:PropertyIsGreaterThan>',
                      '<ogc:PropertyName>streamorde</ogc:PropertyName>',
                      '<ogc:Literal>',streamorder-1,'</ogc:Literal>',
                      '</ogc:PropertyIsGreaterThan>',
                      '<ogc:BBOX>',
                      '<ogc:PropertyName>the_geom</ogc:PropertyName>',
                      '<gml:Envelope>',
                      '<gml:lowerCorner>',mapRange[3]," ",mapRange[1],'</gml:lowerCorner>',
                      '<gml:upperCorner>',mapRange[4]," ",mapRange[2],'</gml:upperCorner>',
                      '</gml:Envelope>',
                      '</ogc:BBOX>',
                      '</ogc:And>',
                      '</ogc:Filter>',
                      '</wfs:Query>',
                      '</wfs:GetFeature>')
  
  destination = file.path(tempdir(),"nhdflowline_network.zip")
  file <- POST(postURL, body = filterXML, write_disk(destination, overwrite=T))
  
  filePath <- tempdir()
  print("unzipping...")
  unzip(destination, exdir = filePath)
  
  flowLines <- st_read(filePath, layer = 'nhdflowline_network')

  return(flowLines)
}


getBasin <- function(sites, filePath = NA){
  
  postURL <- "https://cida.usgs.gov/nwc/geoserver/NWC/ows"
  
  filterXML <- paste0('<?xml version="1.0"?>',
                      '<wfs:GetFeature xmlns:wfs="http://www.opengis.net/wfs" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:gml="http://www.opengis.net/gml" service="WFS" version="1.1.0" outputFormat="shape-zip" xsi:schemaLocation="http://www.opengis.net/wfs http://schemas.opengis.net/wfs/1.1.0/wfs.xsd">',
                      '<wfs:Query xmlns:feature="https://owi.usgs.gov/NWC" typeName="feature:epa_basins" srsName="EPSG:4326">')
  
  
  if(length(sites) > 1){
    siteText <- ""
    for(i in sites){
      siteText <- paste0(siteText,'<ogc:PropertyIsEqualTo  matchCase="true">',
                         '<ogc:PropertyName>site_no</ogc:PropertyName>',
                         '<ogc:Literal>',i,'</ogc:Literal>',
                         '</ogc:PropertyIsEqualTo>')
    }
    
    filterXML <- paste0(filterXML,'<ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">',
                        '<ogc:Or>',siteText,'</ogc:Or>',
                        '</ogc:Filter>')
    
  } else {
    filterXML <- paste0(filterXML,
                        '<ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">',
                        '<ogc:PropertyIsEqualTo matchCase="true">',
                        '<ogc:PropertyName>site_no</ogc:PropertyName>',
                        '<ogc:Literal>',sites,'</ogc:Literal>',
                        '</ogc:PropertyIsEqualTo>',
                        '</ogc:Filter>')
  }
  
  filterXML <- paste0(filterXML,'</wfs:Query>',
                      '</wfs:GetFeature>')
  
  destination = tempfile(pattern = 'basins_shape', fileext='.zip')
  
  file <- POST(postURL, body = filterXML, write_disk(destination, overwrite=T))
  if(is.na(filePath)){
    filePath <- tempdir()
  }
  
  unzip(destination, exdir = filePath)
  basins = st_read(filePath, layer='epa_basins')
  return(basins)
}


getLakes <- function(){
  shapefile_loc <- "http://geo.glin.net/gis/shps/glin_gl_mainlakes.zip"
  
  destination = file.path(tempdir(),"glin_gl_mainlakes.zip")
  file <- GET(shapefile_loc, write_disk(destination, overwrite=T))
  filePath <- tempdir()
  unzip(destination, exdir = filePath)
  
  lakes <- st_read(filePath, layer = "gl_mainlakes")
  return(lakes)
}

mapRange <- c(-93.298145, -74.816895,40.937378, 48.058570 )
streamorder <- 5
crsLONGLAT <- 4326
crs_plot <- st_crs(102003)

#basins <- getBasin(sites$site_no)
filePath <- "M:/QW Monitoring Team/GLRI_GIS/layers/20131299_GLRI_MASTER/~old/201401_HUC12Basins"

basins = st_read(filePath, layer='GLRI_Basins_HUC12Linework')
basins <- st_transform(basins, crs = crsLONGLAT) %>%
  filter(STAID %in% sites$site_no)
basins <- basins %>%
  left_join(sites, by = c("STAID" = "site_no"))


flowlines <- get_flowlines(streamorder, mapRange)
lakes <- getLakes()

# set the state and county names of interest
state_names <- c("minnesota","wisconsin", "michigan","ohio",
                 "illinois","pennsylvania","new york","indiana")

# get STATE data
GL <- us_states(resolution = "high", states = state_names) %>%
  st_transform(crs = crsLONGLAT)

bb <- st_sfc(
  st_point(mapRange[c(1,3)]),
  st_point(mapRange[c(2,4)]),
  crs = crsLONGLAT) 

bb_proj <- st_transform(bb, crs = crs_plot)
b <- st_bbox(bb_proj)


# 
# metolachlor <- chemicalSummary %>%
#   filter(chnm == 'Metolachlor') %>%
#   group_by(site, date) %>%
#   summarize(sumEAR = sum(EAR)) %>%
#   group_by(site) %>%
#   summarize(maxsumEAR_metolachlor = max(sumEAR))
# 
# fipronil <- chemicalSummary %>%
#   filter(chnm == 'Fipronil') %>%
#   group_by(site, date) %>%
#   summarize(sumEAR = sum(EAR)) %>%
#   group_by(site) %>%
#   summarize(maxsumEAR_fipronil = max(sumEAR))
  
# NR1I2 <- get_AOP(gene = "NR1I2")
# names(NR1I2) <- c('site', 'NR_top_chem', 'NR_maxsumEAR', 'NR_max_nchems')
# PTGS2 <- get_AOP(gene = "PTGS2")
# names(PTGS2) <- c('site', 'PT_top_chem', 'PT_maxsumEAR','PT_max_nchems')
# TP53 <- get_AOP(gene = "TP53")
# names(TP53) <- c('site', 'TP_top_chem', 'TP_maxsumEAR', 'TP_max_nchems')
# 

sites <- left_join(sites, graphData, by = c('site_no' = 'site'))
  
  
sites_df <- st_as_sf(sites[, c("dec_long","dec_lat", 'chemicals', 'EAR')],
                     coords = c("dec_long","dec_lat"),
                     crs = crsLONGLAT)

minneapolis <- data.frame(longitude = -93.273882, 
                          latitude = 44.969226)
minneapolis <- st_as_sf(minneapolis, coords = c('longitude', 'latitude'), 
                        crs = crsLONGLAT)

sites_proj <- st_transform(sites_df, crs = crs_plot)
minneapolis <- st_transform(minneapolis, crs = crs_plot)

# add column if fipronil was detected at any given site based on NA vals
# sites_proj <- mutate(sites_proj, fipronil_detected = ifelse(is.na(maxsumEAR_fipronil), 'NO', 'YES'))
# now set NA values to zero
# sites_proj <- mutate(sites_proj, maxsumEAR_fipronil = ifelse(is.na(maxsumEAR_fipronil), 0, maxsumEAR_fipronil))

# Make a quick Map:

#base_map <- 
  
first_pass <- ggplot() + 
  # geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
  geom_sf(data = basins, alpha = 0.5, aes(fill = 'brown')) +
  # scale_fill_gradient(low = '#f5f5f5', high = '#543005', name = "% Agriculture") +
  geom_sf(data = flowlines, color = "lightblue") +
  geom_sf(data = GL, color = "gray50", fill=NA) +
  geom_sf(data = sites_proj, alpha = 0.9, shape = 16, aes(size = chemicals, color = EAR), show.legend = FALSE) +

  scale_size(range = c(2,14), breaks = c(5,10,15,20, 25, 30, 35), guide = 'legend') +
  scale_color_gradient(low = "#ffeda0", high = "#f03b20", breaks = c(0,5,10,15),
                       guide = 'colourbar') +

  geom_sf(data = minneapolis, pch = "\u2605", size = 8) +
  geom_text(aes(x = 214731.983109861, y = 838589.951769598, label = "Minneapolis", vjust = 3, hjust = 0.1), fontface = 'bold') +
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"],b["xmax"]), 
           ylim = c(b["ymin"],b["ymax"])) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = 'transparent'), #Bug that's apparently getting fixed
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(), 
        #legend.position = c(.75, .75),
        legend.direction = 'horizontal') +
  guides(fill = guide_colourbar(title.position = 'top')) +
  geom_point(alpha = 0, shape = 16,
             aes(x = rep(214731.983109861, 16), y = rep(838589.951769598, 16), size = sites_proj$nUniqueChems, color = sites_proj$nChem)) +
  guides(size = guide_legend(label.position = 'bottom', label.hjust = 0.5,  
                             title = '# Chemicals Detected', title.position = 'top', 
                             override.aes = list(alpha = 1, stroke = 2), order = 2), 
         fill = guide_colourbar(title.position = 'top', order = 1), 
         color = guide_colorbar(title.position = 'top', order = 3, 
                                title = '# Chemicals with\nEAR > 0.001')) +
  theme(legend.box = "vertical")

first_pass

#Working above here
#Need to get lakes and watershed %ag data


base_p <- ggplot() + 
  geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
  geom_sf(data = flowlines, color = "lightblue") +
  geom_sf(data = GL, color = "gray50", fill=NA) +
  
  geom_sf(data = minneapolis, pch = "\u2605", size = 8) +
  geom_text(aes(x = 214731.983109861, y = 838589.951769598, label = "Minneapolis", vjust = 3, hjust = 0.1), fontface = 'bold') +

  theme_minimal() +
  theme(panel.grid.major = element_line(colour = 'transparent'), #Bug that's apparently getting fixed
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(), 
        legend.position = c(.75, .75),
        legend.direction = 'horizontal')
# add stuff
metolachlor_p <- base_p + 
  geom_sf(data = basins, alpha = 0.5, aes(fill = `Ag..total`)) +
  scale_fill_gradient(low = '#f5f5f5', high = '#543005', name = "% Agriculture") +
  geom_sf(data = sites_proj, alpha = 0.8, shape = 21, stroke = 2, aes(size = maxsumEAR_metolachlor), col = 'red', show.legend = FALSE) +
  
  scale_size(range = c(2,20), guide = 'legend', breaks = c(0.001, 0.01, 0.1, 1), 
             labels = c(0.001, 0.01, 0.1, 1)) +
  #scale_color_gradient(low = "#ffeda0", high = "#f03b20", breaks = c(0,5,10,15),
  #                    guide = 'colourbar') +
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"],b["xmax"]), 
           ylim = c(b["ymin"],b["ymax"])) +
  geom_point(alpha = 0, shape = 21, col = 'red', 
             aes(x = rep(214731.983109861, 16), y = rep(838589.951769598, 16), size = sites_proj$maxsumEAR_metolachlor)) +
  guides(size = guide_legend(label.position = 'bottom', label.hjust = 0.5,  
                             title = 'max EAR', title.position = 'top', 
                             override.aes = list(alpha = 1, stroke = 2)), 
         fill = guide_colourbar(title.position = 'top')) +
  theme(legend.box = "horizontal")

  


fipronil_p <- base_p +
  geom_sf(data = basins, alpha = 0.8, aes(fill = Urban)) +
  scale_fill_gradient(low = '#f7f7f7', high = "#40004b", name = "% Urban") +
  geom_sf(data = sites_proj, alpha = 0.8, shape = 21, stroke = 2, aes(size = maxsumEAR_fipronil, col = fipronil_detected), show.legend = FALSE) +
  scale_size(range = c(2,20), breaks = c(0, 0.001, 0.01), labels = c(0, 0.001, 0.01)) +
  scale_color_manual(values = c('black', 'red')) +
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"],b["xmax"]), 
           ylim = c(b["ymin"],b["ymax"]))  +
  geom_point(alpha = 0, shape = 21, col = 'red', 
             aes(x = rep(214731.983109861, 16), y = rep(838589.951769598, 16), size = sites_proj$maxsumEAR_fipronil)) +
  guides(size = guide_legend(label.position = 'bottom', label.hjust = 0.5, order = 2,  
                             title = 'max EAR', title.position = 'top', override.aes = list(alpha = 1, stroke = 2, color = c('black', 'red', 'red'))), 
         fill = guide_colorbar(order = 1, title.position = 'top')) +
  theme(legend.box = "horizontal")

NR1I2_p <- base_p +
  geom_sf(data = basins, alpha = 0.8, fill = 'lightgray') +
  geom_sf(data = sites_proj, alpha = 0.8, shape = 16, 
          aes(size = NR_maxsumEAR, color = NR_top_chem), show.legend = FALSE) +
  scale_size(range = c(2,20), breaks = c(0.0001, 0.001, 0.01, 0.1), labels = c(0.0001, 0.001, 0.01, 0.1)) +
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"],b["xmax"]), 
           ylim = c(b["ymin"],b["ymax"]))  +
  geom_point(alpha = 0, shape = 16, 
             aes(x = rep(214731.983109861, 16), y = rep(838589.951769598, 16), 
                 size = sites_proj$NR_maxsumEAR, color = sites_proj$NR_top_chem)) +
  guides(size = guide_legend(label.position = 'right', label.hjust = 0.5, order = 1,  
                             title = 'max EAR', title.position = 'top', override.aes = list(alpha = 1)), 
         color = guide_legend(order = 2, title.position = 'top', 
                              override.aes = list(alpha = 1, size=8),
                              title = 'Top Contributor', ncol = 2)) +
  theme(legend.box = "horizontal", 
        legend.direction = 'vertical')

# add column if chem leading to PTGS2 was detected at any given site based on NA vals
sites_proj <- mutate(sites_proj, chemto_PTGS2_detected = ifelse(is.na(PT_maxsumEAR), 'NO', 'YES'))
# now set NA values to zero
sites_proj <- mutate(sites_proj, PT_maxsumEAR = ifelse(is.na(PT_maxsumEAR), 0, PT_maxsumEAR))
sites_proj$PT_top_chem[is.na(sites_proj$PT_top_chem)] <- "No Chemicals"

PTGS2_p <- base_p +
  geom_sf(data = basins, alpha = 0.8, fill = 'lightgray') +
  geom_sf(data = sites_proj, alpha = 0.8, shape = 16, 
          aes(size = PT_maxsumEAR, color = PT_top_chem), show.legend = FALSE) +
  scale_size(range = c(2,20), breaks = c(0.0001, 0.001, 0.01, 0.1), labels = c(0.0001, 0.001, 0.01, 0.1)) +
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"],b["xmax"]), 
           ylim = c(b["ymin"],b["ymax"]))  +
  geom_point(alpha = 0, shape = 16, 
             aes(x = rep(214731.983109861, 16), y = rep(838589.951769598, 16), 
                 size = sites_proj$PT_maxsumEAR, color = sites_proj$PT_top_chem)) +
  guides(size = guide_legend(label.position = 'right', label.hjust = 0.5, order = 1,  
                             title = 'max EAR', title.position = 'top', override.aes = list(alpha = 1)), 
         color = guide_legend(order = 2, title.position = 'top', 
                              override.aes = list(alpha = 1, size=8),
                              title = 'Top Contributor', ncol = 1)) +
  theme(legend.box = "horizontal", 
        legend.direction = 'vertical')

# add column if chem leading to PTGS2 was detected at any given site based on NA vals
# sites_proj <- mutate(sites_proj, chemto_PTGS2_detected = ifelse(is.na(PT_maxsumEAR), 'NO', 'YES'))
# now set NA values to zero
sites_proj <- mutate(sites_proj, TP_maxsumEAR = ifelse(is.na(TP_maxsumEAR), 0, TP_maxsumEAR))
sites_proj$TP_top_chem[is.na(sites_proj$TP_top_chem)] <- "No Chemicals"

TP53_p <- base_p +
  geom_sf(data = basins, alpha = 0.8, fill = 'lightgray') +
  geom_sf(data = sites_proj, alpha = 0.8, shape = 16, 
          aes(size = TP_maxsumEAR, color = TP_top_chem), show.legend = FALSE) +
  scale_size(range = c(2,20), breaks = c(0.001, 0.01, 0.1, 1), labels = c(0.001, 0.01, 0.1, 1)) +
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"],b["xmax"]), 
           ylim = c(b["ymin"],b["ymax"]))  +
  geom_point(alpha = 0, shape = 16, 
             aes(x = rep(214731.983109861, 16), y = rep(838589.951769598, 16), 
                 size = sites_proj$TP_maxsumEAR, color = sites_proj$TP_top_chem)) +
  guides(size = guide_legend(label.position = 'right', label.hjust = 0.5, order = 1,  
                             title = 'max EAR', title.position = 'top', override.aes = list(alpha = 1)), 
         color = guide_legend(order = 2, title.position = 'top', 
                              override.aes = list(alpha = 1, size=8),
                              title = 'Top Contributor', ncol = 1)) +
  theme(legend.box = "horizontal", 
        legend.direction = 'vertical')
  

  
  


base_map
  
ggsave(first_pass, filename = "site_map_nchems.png", height = 6, width = 8)  
ggsave(metolachlor_p, filename = "site_map_metolachlor_nolog.png", height = 6, width = 8)  
ggsave(fipronil_p, filename = "site_map_fipronil_nolog.png", height = 6, width = 8)  
ggsave(NR1I2_p, filename = "site_map_NR1I2_nolog.png", height = 8, width = 10.5)  
ggsave(PTGS2_p, filename = "site_map_PTGS2_nolog.png", height = 8, width = 10.5)  
ggsave(TP53_p, filename = "site_map_TP53_nolog.png", height = 8, width = 10.5)  

