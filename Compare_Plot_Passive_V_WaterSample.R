

#Compare pocix and water sample directly
#Use chemical summaries from both ToxEval prep scripts

#R objects
chemicalSummary #Select data from POCIS 
chemicalSummary_allpocis #all data from POCIS
chemicalSummary_bench # Select data using custom benchmarks

chemicalSummary_surface #Sam's data from surface
chemicalSummary_bench_surface # Sam's surface data using custom benchmarks



chemicalSummary_surface_annual <- chemicalSummary_surface %>%
  # filter(EAR>0) %>%
  filter(CAS %in% unique(tox_list$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list$chem_site$SiteID)) %>%
  mutate(Class = as.character(Class),
         chnm = as.character(chnm)) %>%
  # mutate(chnm = factor(chnm, levels(chem_detection$chnm))) %>%
  dplyr::group_by(site, chnm, date, Class, CAS) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  dplyr::group_by(site, chnm, Class) %>%
  mutate(shortName = site_order[match(site, site_ID_order)])

# chemicalSummary_surface_annual$EAR[which(chemicalSummary_surface_annual$EAR==0)] <- 10^-7

chemicalSummary_surface_annual$Detected<- "Detected"
chemicalSummary_surface_annual$Detected[which(chemicalSummary_surface_annual$EAR==0)] <- "belowMDL"

# chemicalSummary_surface_annual$EAR[which(chemicalSummary_surface_annual$EAR==0)] <- 10^-7


chemicalSummary_surface_annual$chnm[which(chemicalSummary_surface_annual$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))] <- 
  paste0(chemicalSummary_surface_annual$chnm[which(chemicalSummary_surface_annual$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))], "*")
chemicalSummary_surface_annual$Class <- gsub("Deg - ", "", chemicalSummary_surface_annual$Class)



chemicalSummary_surface_prepped <- chemicalSummary_surface %>%
  # filter(EAR>0) %>%
  filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(CAS %in% unique(tox_list$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list$chem_site$SiteID)) %>%
  # mutate(chnm = factor(chnm, levels(chem_detection$chnm))) %>%
  drop_na(chnm) %>%
  dplyr::group_by(site, chnm, Class, CAS) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  dplyr::group_by(site, chnm, Class) %>%
  mutate(method = "water")

chemicalSummary_surface_prepped$WaterDetected<- "Detected"
chemicalSummary_surface_prepped$WaterDetected[which(chemicalSummary_surface_prepped$EAR==0)] <- "belowMDL"

# chemicalSummary_surface_prepped$EAR[which(chemicalSummary_surface_prepped$EAR==0)] <- 10^-7


chemicalSummary_passive_prepped <- chemicalSummary %>%
  # filter(EAR>0) %>%
  filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  # mutate(chnm = factor(chnm, levels(chem_detection$chnm))) %>%
  drop_na(chnm) %>%
  dplyr::group_by(site, chnm, Class, CAS) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  dplyr::group_by(site, chnm, Class) %>%
  mutate(method = "passive")

chemicalSummary_passive_prepped$PassiveDetected<- "Detected"
chemicalSummary_passive_prepped$PassiveDetected[which(chemicalSummary_passive_prepped$EAR==0)] <- "belowMDL"


chemicalSummary_passive_prepped_spread <- chemicalSummary_passive_prepped %>%
  rename(passive = EAR) %>%
  select(-method, -CAS)

chemicalSummary_passive_surface_spread <- chemicalSummary_surface_prepped %>%
  rename(surface = EAR) %>%
  select(-method, -CAS)

#Spread each table before joining
chemicalSummary_merged <- full_join(chemicalSummary_passive_surface_spread, chemicalSummary_passive_prepped_spread) 

chemicalSummary_merged$shortName = site_order[match(chemicalSummary_merged$site, site_ID_order)]

AnyMDL <- which(chemicalSummary_merged$WaterDetected == "belowMDL" | chemicalSummary_merged$PassiveDetected == "belowMDL")

chemicalSummary_merged$Detected <- "Detected"
chemicalSummary_merged$Detected[AnyMDL] <- "belowMDL"

# chemicalSummary_merged$Detected = "Detected"
# chemicalSummary_merged$Detected[which(chemicalSummary_merged$passive==10^-7 | chemicalSummary_merged$water==10^-7)] <- 'belowMDL'

chemicalSummary_merged$chnm[which(chemicalSummary_merged$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))] <- 
  paste0(chemicalSummary_merged$chnm[which(chemicalSummary_merged$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))], "*")
chemicalSummary_merged$Class <- gsub("Deg - ", "", chemicalSummary_merged$Class)

chemicalSummary_merged <- chemicalSummary_merged %>%
  group_by() %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide')),
         chnm = factor(chnm, levels(chem_detection$chnm))) %>%
  drop_na(chnm)


mainchems <- c('Diuron', 'Metolachlor', 'Atrazine', 'Deisopropylatrazine*', 'Acetochlor', 'Simazine', 'Prometon', 'Dimethenamid', 'Propazine', 'Metalaxyl', 'Azoxystrobin', 'Imidacloprid')

chemicalSummary_mainchems <- chemicalSummary_surface_annual %>% 
  group_by() %>%
  filter(chnm %in% mainchems) %>%
  mutate(chnm = factor(chnm, mainchems),
         Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide')),
         shortName = factor(shortName, site_order)) 

# chemicalSummary_mainchems$Detected <- 'Detected'
# chemicalSummary_mainchems$Detected[which(chemicalSummary_mainchems$EAR==10^-7)] <- 'belowMDL'


chemicalSummary_merged_mainchems <- chemicalSummary_merged %>%
  group_by() %>%
  filter(chnm %in% mainchems) %>%
  mutate(chnm = factor(chnm, mainchems),
         Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide')),
         shortName = factor(shortName, site_order)) 

  
ggplot(chemicalSummary_mainchems, aes(x=EAR)) +
  scale_x_log10nice() + 
  geom_vline(data=chemicalSummary_merged_mainchems, aes(xintercept=passive), col='blue', size=1.5) + 
  geom_vline(data=chemicalSummary_merged_mainchems, aes(xintercept=surface), col='red', size=1.5, linetype='dotted') + 
  # geom_vline(xintercept = .001) +
  # geom_histogram(colour='black', fill='white') +
  geom_density(alpha=.2, fill= '#FF6666') + 
  facet_grid(shortName~chnm, scales = 'free') +
  theme_bw()


# ggplot(chemicalSummary_mainchems, aes(x=as.Date(date), y=EAR, col='#FF6666')) +
#   scale_y_log10nice() + 
#   geom_segment(data=chemicalSummary_merged_mainchems, aes(x=date_filter[1], xend = date_filter[2], y=passive, yend=passive), col='blue', size=1) +
#   # geom_vline(xintercept = .001) +
#   geom_point() +
#   geom_path() + 
#   facet_grid(shortName~chnm, scales = 'free_y') +
#   theme_bw()


mainchem_timeseries_bysite <- ggplot(chemicalSummary_mainchems, aes(x=as.Date(date), y=EAR, col=Class)) +
  scale_y_log10nice(name='EAR') + 
  scale_x_date(date_breaks='4 months', date_labels = '%b') + 
  geom_segment(data=chemicalSummary_merged_mainchems, aes(x=date_filter[1], xend = date_filter[2], y=passive, yend=passive), col='black', size=1.5) +
  # geom_vline(xintercept = .001) +
  geom_point(aes(shape=Detected)) +
  geom_path() + 
  # facet_grid(chnm~shortName) +
  facet_grid(chnm~shortName, scales = 'free_y') +
  theme_bw() +
  scale_color_brewer(palette = 'Dark2') +
  scale_shape_manual(values=c(1,16)) +
  labs(x='Date') +
  theme(legend.position='bottom', legend.title=element_blank())

print(mainchem_timeseries_bysite)

ggsave(file_out(file.path(path_to_data, "Figures/Timeseries_EAR_Bychem.png")), mainchem_timeseries_bysite, height=12, width=15, units='in')



Scatterplot_passive_versus_water <- ggplot(chemicalSummary_merged, aes(x=passive, y=surface, col=Class)) +
  scale_y_log10nice(name="water sample EAR") + 
  scale_x_log10nice(name = "passive EAR") +
  geom_hline(yintercept = 10^-3) +
  geom_vline(xintercept = 10^-3) +
  geom_abline(linetype='dashed') +
  geom_point(aes(shape=Detected), size=2,stroke=1.5, fill='white') +
  scale_shape_manual(values=c(21,16)) +
  scale_color_brewer(palette = 'Dark2') +
  facet_wrap(~chnm, ncol=4) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position='bottom', legend.title = element_blank()) 

print(Scatterplot_passive_versus_water)

ggsave(file_out(file.path(path_to_data, "Figures/Scatterplot_EAR_Passive_V_WaterSample.png")), Scatterplot_passive_versus_water, height=10, width=7, units='in')



chemicalSummary_merged_hits <- chemicalSummary_merged %>%
  select(chnm, passive, water) %>%
  gather(key='method', value='EAR', 2:3) %>%
  filter(EAR>=0.001) %>%
  group_by(chnm, method) %>%
  tally() %>%
  spread(method, n)

ggplot(chemicalSummary_merged_hits) +
  geom_histogram(aes(x=water-passive), binwidth=.5, color='black', fill='white') +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())






#using figures from "Plot_SurfaceTox_2016" and "Plot_PassiveTox_2016" make a combined EAR and TQ figure


EARbox_passive_v2 <- EARbox_passive +
  # scale_x_discrete(position = 'top') 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ggtitle('Passive Samples') +
  theme(plot.title = element_text(hjust = 0.5))
  
EARbox_surface_v2 <- EARbox_surface +
  theme(axis.text.y = element_text(hjust=0.5), axis.ticks.y = element_blank()) +
  ggtitle('Surface Water Samples') +
  theme(plot.title = element_text(hjust = 0.5))


box_by_box_EAR <- grid.arrange(grobs=list(EARbox_passive_v2, EARbox_surface_v2), nrow=1)



png(file.path(path_to_data, "Figures/EAR_Passive_V_WaterSample_Boxplot.png"), height=6, width=5, units='in', res=400)

grid.newpage()
boxes_EAR<-grid.draw(cbind(ggplotGrob(EARbox_passive_v2), ggplotGrob(EARbox_surface_v2), size = "first"))

dev.off()
