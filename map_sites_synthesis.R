# devtools::install_github("tidyverse/ggplot2")
# devtools::install_github("r-spatial/sf")

# install.packages("nhdR")
library(nhdR)
library(nhdplusTools)
library(httr)
library(ggplot2)
library(dplyr)
# library(remake)
library(sf)
library(USAboundaries)
library(ggsn)

EAR_mix_max <- EAR_mix %>%
  group_by(site) %>%
  summarize(EARexceed = length(which(EARsum > 0.01)),
            EARsum = max(EARsum, na.rm = TRUE))
            

tox_file_max$chem_site$dec_lat[which(tox_file_max$chem_site$SiteID == "04232007")] <- 43.230583
tox_file_max$chem_site$dec_lon[which(tox_file_max$chem_site$SiteID == "04232007")] <-  -77.616528

sites <- tox_file_max$chem_site %>%
  rename(SiteID = SiteID) %>%
  group_by()

# get info about number of chemicals by site, merge with sites
nchems <- tox_file_max$chem_data %>%
  filter(Value>0) %>%
  group_by(SiteID) %>%
  dplyr::summarize(chemicals = n())


graphData <- nchems %>% 
  left_join(sites, by = "SiteID") %>%
  filter(!is.na(dec_lon) & !is.na(dec_lat)) %>%
  left_join(EAR_mix_max, by = c("SiteID" = "site"))

mapRange <- c(-93.298145, -74.016895,40.937378, 48.058570 )
streamorder <- 5
crsLONGLAT <- 4326
# crs_plot <- st_crs(102003)
# crs_plot <- st_crs(crsLONGLAT)

crs_plot <- st_crs('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

#basins <- getBasin(sites$site_no)
# filePath <- "M:/QW Monitoring Team/GLRI_GIS/layers/20131299_GLRI_MASTER/~old/201401_HUC12Basins"

# basins = st_read(filePath, layer='GLRI_Basins_HUC12Linework')
# basins <- st_transform(basins, crs = crsLONGLAT) %>%
#   dplyr::filter(STAID %in% sites$site_no)
# basins <- basins %>%
#   left_join(sites, by = c("STAID" = "site_no"))

# basins_OGR <- readRDS( file_in(file.path(path_to_data, "Rdata/watershedmetrics.rds")))
# basins_OGR <- basins_OGR[basins_OGR$STAID %in% sites$site_no,]
# basins <- st_as_sf(basins_OGR) %>%
#   st_transform(crs=crsLONGLAT) %>%
#   dplyr::filter(STAID %in% sites$site_no)
# basins <- basins %>%
#   left_join(sites, by = c("STAID" = "site_no"))


lakespath <- "C:/Users/lloken/OneDrive - DOI/GIS/Great_Lakes"
lakes <- st_read(lakespath, layer="Great_Lakes")

# set the state and county names of interest
state_names <- c("minnesota","wisconsin", "michigan","ohio",
                 "illinois","pennsylvania","new york","indiana", "west virginia")

# get STATE data
GL <- us_states(resolution = "high", states = state_names) %>%
  st_transform(crs = crsLONGLAT)

bb <- st_sfc(
  st_point(mapRange[c(1,3)]),
  st_point(mapRange[c(2,4)]),
  crs = crsLONGLAT) 

bb_proj <- st_transform(bb, crs = crs_plot)
b <- st_bbox(bb_proj)



# nhd_plus_get(vpu = 4, "NHDSnapshot")
# 
# 
# nhd_get(state = c("MN", "WI", "IL", "IN", "MI", "OH", "PA", "NY"))
# 
# state_flowlines <- nhd_load(state = c("MN", "WI", "IL", "IN", "MI", "OH", "PA", "NY"), 
#                             dsn = "NHDFlowline")
# 
# small_flowline <- filter(state_flowlines, VisibilityFilter > 3000000)
# 
# indexes <- get_flowline_index(small_flowline,
#                               sites_df,
#                               search_radius = 0.1,
#                               max_matches = 1)
# 
# plot(st_geometry(small_flowline))
# 
# # start_point <- st_sfc(st_point(c(-89.362239, 43.090266)), crs = 4269)
# start_point <- st_sfc(st_point(c(-89.6162, 48.01211)), crs = 4269)
# start_comid <- rep(NA, nrow(sites_df))
# flowline <- list()
# for (i in 1:length(start_comid)){
#   start_point <- st_geometry(sites_df$geometry[i,])
#   start_comid[i] <- discover_nhdplus_id(start_point)
#   flowline[[i]] <- navigate_nldi(list(featureSource = "comid",
#                                  featureID = start_comid[i]),
#                             mode = "upstreamTributaries",
#                             distance_km = 9000)
#   subset_file <- tempfile(fileext = ".gpkg")
#   subset <- subset_nhdplus(comids = flowline[[i]]$nhdplus_comid,
#                            output_file = subset_file,
#                            nhdplus_data = "download",
#                            flowline_only = TRUE,
#                            return_data = TRUE)
#   if (i == 1){
#     
#     flowline_out <- subset$NHDFlowline_Network
#     
#   } else if (i > 1){
#     
#     flowline_out <- rbind(flowline_out, subset$NHDFlowline_Network)
#   }
#   print(i)
# }


# saveRDS(flowline_out, file.path(path_to_data, "Rdata", "NHDPlus_Flowlines_Synthesis.rds"))
flowline_out <- readRDS(file.path(path_to_data, "Rdata", "NHDPlus_Flowlines_Synthesis.rds"))

# plot(st_geometry(flowline_out))

#Test plot. look at NHD
# ggplot() + 
#   geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
#   geom_sf(data = flowlines_all, color = "orange") +
#   geom_sf(data = GL, color = "gray50", fill=NA)

indexes <- get_flowline_index(flowline_out,
                              sites_df,
                              search_radius = 0.1,
                              max_matches = 1)

indexes <- left_join(sf::st_sf(id = c(1:nrow(sites_df)),
                               geom = sf::st_geometry(sites_df)),
                   indexes, by = "id")

flowlines_small <- flowline_out %>%
  filter(reachcode %in% indexes$REACHCODE)

flowlines_med <- flowline_out %>%
  filter(terminalpa %in% unique(flowlines_small$terminalpa) &
           # gnis_name %in% unique(flowlines_small$gnis_name) &
           streamorde >= 4)


sites_df <- st_as_sf(graphData[, c("dec_lon","dec_lat", 'chemicals', 'EARsum', 'EARexceed')],
                     coords = c("dec_lon","dec_lat"),
                     crs = crsLONGLAT)

sites_proj <- st_transform(sites_df, crs = crs_plot)

colors_EARmap <- brewer.pal(n = 9, name = "YlOrRd")[c(2,5,8,9)]

  
first_pass <- ggplot() + 
  geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
  geom_sf(data = flowlines_med, color = "lightblue", size = 1) +
  geom_sf(data = GL, color = "gray50", fill=NA) +

  geom_sf(data = sites_df, alpha = 0.7, shape = 21, aes(fill = chemicals), show.legend = TRUE, size = 5) +

  scale_fill_gradientn(colours=colors_EARmap, breaks = c(0,50,100, 150, 200, 250),
                       guide = 'colourbar') + 

  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"]+1000,b["xmax"]+10000),
           ylim = c(b["ymin"],b["ymax"]-100000)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        # panel.background = element_blank(),
        axis.text = element_text(),
        axis.title = element_blank(),
        axis.ticks = element_line(), 
        legend.position = c(.8, .9),
        legend.direction = 'horizontal') +

  guides(fill = guide_colorbar(title.position = 'top', 
                                title = expression('# Chemicals Detected'))) 
  # theme(legend.box = "horizontal") 

first_pass

ggsave(first_pass, filename = file_out(file.path(path_to_data, "Figures/synthesis_map2_nchems.png")), height = 5, width = 8)  


second_pass <- ggplot() + 
  geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
  geom_sf(data = flowlines_med, color = "lightblue", size = 1) +
  geom_sf(data = GL, color = "gray50", fill=NA) +
  
  geom_sf(data = sites_df, alpha = 0.7, shape = 21, aes(fill = EARexceed), show.legend = TRUE, size = 5) +
  
  scale_fill_gradientn(colours=colors_EARmap, breaks = c(0,50, 100, 150),
                       guide = 'colourbar') + 
  
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"]+1000,b["xmax"]+10000),
           ylim = c(b["ymin"],b["ymax"]-100000)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        # panel.background = element_blank(),
        axis.text = element_text(),
        axis.title = element_blank(),
        axis.ticks = element_line(), 
        legend.position = c(.8, .9),
        legend.direction = 'horizontal') +
  
  guides(fill = guide_colorbar(title.position = 'top', 
                               title = expression(paste('Number of ', EAR[mix], ' > 0.01')))) 
# theme(legend.box = "horizontal") 

second_pass

ggsave(second_pass, filename = file_out(file.path(path_to_data, "Figures/synthesis_map2_numberexceedEAR.png")), height = 5, width = 8)  




third_pass <- ggplot() + 
  geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
  geom_sf(data = flowlines_med, color = "lightblue", size = 1) +
  geom_sf(data = GL, color = "gray50", fill=NA) +
  
  geom_sf(data = sites_df, alpha = 0.7, shape = 21, aes(fill = EARsum), show.legend = TRUE, size = 5) +
  
  scale_fill_gradientn(colours=colors_EARmap, breaks = c(0,1,2),
                       guide = 'colourbar') + 
  
  coord_sf(crs = crs_plot,
           xlim = c(b["xmin"]+1000,b["xmax"]+10000),
           ylim = c(b["ymin"],b["ymax"]-100000)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        # panel.background = element_blank(),
        axis.text = element_text(),
        axis.title = element_blank(),
        axis.ticks = element_line(), 
        legend.position = c(.8, .9),
        legend.direction = 'horizontal') +
  
  guides(fill = guide_colorbar(title.position = 'top', 
                               title = expression(paste('Max ', EAR[mix])))) 
# theme(legend.box = "horizontal") 

third_pass

ggsave(third_pass, filename = file_out(file.path(path_to_data, "Figures/synthesis_map2_maxEARmix.png")), height = 5, width = 8)  
