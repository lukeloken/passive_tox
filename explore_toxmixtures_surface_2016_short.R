
# #########
# Mixtures
# #########

#### Setup ####
library(toxEval)
library(ToxMixtures)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

#Need to run (ToxEval_Passive2016.R and ToxEval_WaterSamples2016.R)

# Take maximum EAR for water samples with more than one sample
cs_surf <- chemicalSummary_surface %>%
  # mutate(shortName = factor(shortName, site_order)) %>%
  # filter(`date`>=date_filter[1], `date`<=date_filter[2],
  #        shortName %in% site_order) %>%
  # mutate(shortName = factor(shortName, site_order)) %>%
  group_by(site, endPoint, shortName, Lake, site_grouping, CAS, chnm, Bio_category, Class) %>%
  summarise(EAR = max(EAR, na.rm=T), .groups = 'drop') 

# Prepare passive data
# cs_pass <- chemicalSummary %>%
#   mutate(shortName = factor(shortName, site_order))

# thresholds for mixtures analysis
group_by_this <- "endPoint"
ear_threshold <- .001 #Ear threshold
site_threshold <- 1.9 #How many sites

#Surface samples

cs_surf <- chemicalSummary_surface #%>%
# mutate(shortName = factor(shortName, site_order)) %>%
# filter(`date`>=date_filter[1], `date`<=date_filter[2],
#        shortName %in% site_order) %>%
# mutate(shortName = factor(shortName, site_order))

#Site table to add land use
site_table <- unique(cs_surf[c('shortName', 'site')])
site_order_surface <- site_table$shortName[match(surface_ID_order, site_table$site)]

site_table_v2 <- site_table %>%
  mutate(LU = c("Natural", "Natural", "Crops", "AgMix", 
                "Urban", "Urban", "Crops", "AgMix", 
                "Crops", "Urban", "Urban", "Crops", 
                "Crops", "Urban", "AgMix", "AgMix"),
         Disturbance = c(1,2,1,4,3,4,4,3,2,2,5,5,3,1,2,1)) %>%
  mutate(LU = factor(LU, c("Urban", "Crops", "AgMix", "Natural"))) %>%
  arrange(LU, desc(Disturbance))


allcombos_surf <- top_mixes(cs_surf, group_by_this, 
                            ear_threshold,
                            .9) 

prioritycombos_surf <- filter(allcombos_surf, n_samples>10)
prioritymixtures_pass <- overall_mixtures(prioritycombos_surf, "max")

class_key <- class_key_fnx(cs_surf)

prioritymixtures_fancy <- overall_df_format(prioritymixtures_pass, class_key)


allmixtures_surf <- overall_mixtures(allcombos_surf, "max")


# top combos EAR > 0.01 and site_treshold > 1.9
topcombos_surf <- top_mixes(cs_surf,
                            group_by_this, 
                            0.01,
                            site_threshold) 


topmixtures_surf <- overall_mixtures(topcombos_surf, "max")

genetable_surf <- gene_mixtures(topmixtures_surf)  


genefunctions_surf <- gene_functions(topmixtures_surf, 
                                     species = c("Homo sapiens", 
                                                 "Danio rerio", 
                                                 "Xenopus tropicalis"))

geneentrez_surf <- gene_functions(topmixtures_surf, 
                                  columns = "entrez", 
                                  species = "Homo sapiens")

genesummary_surf <- gene_summary(genefunctions_surf)

#Calculate EARmix for each sample and endpoint (note includes all EARs)
EARmix_surf <- EAR_mixtures(cs_surf, group_by_this)

#Add site information and grouping variables
EARmix_surf <- EARmix_surf %>% 
  left_join(unique(cs_surf[c("site", "shortName", "Lake")])) %>%
  left_join(site_table_v2)  


EARsite_surf <- site_mixtures(EARmix_surf, ear_cutoff = 0.01) %>%
  full_join(site_table_v2) %>%
  select(site, shortName, LU, everything())

data.frame(EARsite_surf)

#plotting
surf_gene_fig <- plot_genebar(EARmix_surf, ear_cutoff = 0.01, 
                              site_threshold, type="gene", 
                              facet_col = "LU", fill = LU) + 
  theme(legend.position = "none") +
  scale_x_log10nice(name = expression(paste("max ", EAR[mixture], " by site")))


surf_endpoint_fig <- plot_genebar(EARmix_surf, ear_cutoff = 0.01, 
                                  site_threshold, type="endPoint", 
                                  facet_col = "LU", fill = LU) + 
  theme(legend.position = "none") +
  scale_x_log10nice(name = expression(paste("max ", EAR[mixture], " by site")))


ggsave(file.path(path_to_data, "Figures", "SurfaceGeneByLanduse.png"),  
       plot = surf_gene_fig, height=4, width=8)


ggsave(file.path(path_to_data, "Figures", "SurfaceEndpointByLanduse.png"),  
       plot = surf_endpoint_fig, height=4, width=8)



       ggsave(file_out(file.path(path_to_data, "Figures/StackBox_EARBySite_watersamples_allyear.png")),
       
       