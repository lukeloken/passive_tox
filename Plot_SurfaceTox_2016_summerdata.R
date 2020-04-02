

#Calculate chemical summary for surface water pesticide data
# tox_list_surface <- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides.xlsx")))

#R objects
chemicalSummary_surface #Sam's data from surface
chemicalSummary_bench_surface # Sam's surface data using custom benchmarks

#Check out the standard boxplots
plot_tox_boxplots(chemicalSummary_bench_surface, 
                  category = "Chemical", 
                  sum_logic = FALSE,
                  x_label = "Toxicity Quotient")

plot_tox_boxplots(chemicalSummary_surface, 
                  category = "Chemical", 
                  sum_logic = FALSE)

date_filter <-as.Date(c("2016-06-01", "2016-07-21"))

colors_EAR2 <- c('grey90', colors_EAR)


# create tables by chemical and by sites
# these are detections above the minimum detection
# some data do not have benchmarks so first identify if they were detected

#Identify how many detections by chemical total
chem_freq_all<-tox_list_surface$chem_data %>%
  filter(`Sample Date`>=date_filter[1], `Sample Date`<=date_filter[2]) %>%
  filter(Value>0) %>%
  # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  filter(SiteID %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  group_by(CAS, SiteID) %>%
  summarize(Value_mean = mean(Value, na.rm=T)) %>%
  tally(name = "LabDetected") %>%
  left_join(tox_list_surface$chem_info)

chem_freq_all$compound[which(is.na(chem_freq_all$compound))] <- chem_freq_all$`Chemical Name`[which(is.na(chem_freq_all$compound))] 

chem_freq_all<-rename(chem_freq_all, chnm = compound) %>%
  select(CAS, LabDetected, chnm, Class)

#Use all chemicals in surface EAR dataset
chemicalSummary2_surface <- chemicalSummary_surface %>%
  filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(EAR>0) %>%
  # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  dplyr::group_by(site, CAS, chnm, date) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  dplyr::group_by(site, CAS,chnm) %>%
  summarize(EAR=mean(EAR, na.rm=T))


#Use all chemicals in surface TQ dataset
chemicalSummary_bench2 <- chemicalSummary_bench_surface %>%
  filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(EAR>0) %>%
  # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  filter(date >= date_filter[1], date < date_filter[2]) %>%
  dplyr::group_by(site, chnm, CAS) %>% 
  summarize(EAR = max(EAR, na.rm=T))



chnm_list <- tox_list_surface$chem_data %>%
  filter(`Sample Date`>=date_filter[1], `Sample Date`<=date_filter[2]) %>%
  filter (Value > 0) %>%
  filter(SiteID %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  group_by (SiteID, CAS, pCode) %>%
  summarize(Value = max(Value, na.rm=T)) %>%
  group_by(SiteID) %>%
  summarize(n=n()) %>%
  left_join(tox_list_surface$chem_site) %>%
  rename(shortName = `Short Name`,
         site = SiteID)

#only look at data analyzed with POCIS

# chemicalSummary2_surface <- chemicalSummary_surface %>%
#   filter(date>=date_filter[1], date<=date_filter[2]) %>%
#   filter(EAR>0) %>%
#   filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
#   filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
#   dplyr::group_by(site, chnm, date) %>% 
#   summarize(EAR = max(EAR, na.rm=T)) %>%
#   dplyr::group_by(site, chnm) %>%
#   summarize(EAR=mean(EAR, na.rm=T))
# 
# 
# chem_freq_all<-chemicalSummary_surface %>%
#   filter(date>=date_filter[1], date<=date_filter[2]) %>%
#   filter(EAR>0) %>%
#   # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
#   filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
#   filter(date >= date_filter[1], date < date_filter[2]) %>%
#   group_by(chnm, site) %>%
#   summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
#   tally(name = "EARDetected")
# 
# chemicalSummary_bench2 <- chemicalSummary_bench_surface %>%
#   filter(date>=date_filter[1], date<=date_filter[2]) %>%
#   filter(EAR>0) %>%
#   # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
#   filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
#   filter(date >= date_filter[1], date < date_filter[2]) %>%
#   dplyr::group_by(site, chnm) %>% 
#   summarize(EAR = max(EAR, na.rm=T))




# Calculate chemical frequencies
chem_freq_EAR0.01<-chemicalSummary2_surface %>%
  group_by(CAS, site) %>%
  filter(EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.01")

chem_freq_EAR0.001<-chemicalSummary2_surface %>%
  group_by(CAS, site) %>%
  filter(EAR >= 0.001 & EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

chem_freq_EAR0.0001<-chemicalSummary2_surface %>%
  group_by(CAS, site) %>%
  filter(EAR < 0.001 & EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")

chem_freq_Detected<-chemicalSummary2_surface %>%
  group_by(CAS, site) %>%
  filter(EAR < 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "EARDetected")


chem_detection <- full_join(chem_freq_EAR0.01, chem_freq_EAR0.001) %>%
  full_join(chem_freq_EAR0.0001) %>%
  full_join(chem_freq_Detected) %>%
  left_join(unique(chemicalSummary_surface[c("CAS", "chnm", "Class")])) %>%
  rowwise() %>%
  mutate(chnm = as.character(chnm))


chem_detection <- full_join(chem_detection, chem_freq_all) %>%
  mutate(UnknownTox = LabDetected - 
           sum(AboveEAR0.01, AboveEAR0.001, AboveEAR0.0001, EARDetected, na.rm=T)) %>%
  select(-LabDetected)
  


chem_detection$chnm[which(chem_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))] <- 
  paste0(chem_detection$chnm[which(chem_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))], "*")
chem_detection$Class <- gsub("Deg - ", "", chem_detection$Class)

chem_detection <- chem_detection %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  arrange(Class, desc(AboveEAR0.01), desc(AboveEAR0.001), desc(AboveEAR0.0001), desc(EARDetected), desc(UnknownTox)) %>%
  tidyr::gather(key=group, value=value, -chnm, -Class, -CAS) %>%
  mutate(chnm = factor(chnm, levels=unique(chnm)), 
         group = factor(group, levels=rev(c("AboveEAR0.01", "AboveEAR0.001", "AboveEAR0.0001", "EARDetected", "UnknownTox"))))

chem_detection$value[which(is.na(chem_detection$value))] <- 0



##########################################

#Make similar table for TQ

chem_freq_TQ1<-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR >= 1) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ1")


chem_freq_TQ0.1<-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR >= 0.1 & EAR<1) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ0.1")

chem_freq_TQ0.01<-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR < 0.1 & EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ0.01")

chem_freq_TQDetected<-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "Detected")

# TQ_detection <- full_join(chem_freq_TQ1, chem_freq_TQ0.1) %>%
#   full_join(chem_freq_TQ0.01) %>%
#   full_join(chem_freq_TQDetected) %>%
#   full_join(chem_freq_all) %>% 
#   # left_join(unique(chemicalSummary_surface[c("chnm", "Class")])) %>%
#   rowwise() %>%
#   mutate(UnknownTox = LabDetected - 
#            sum(AboveTQ1, AboveTQ0.1, AboveTQ0.01, Detected, na.rm=T)) %>%
#   select(-LabDetected)



TQ_detection <- full_join(chem_freq_TQ1, chem_freq_TQ0.1) %>%
    full_join(chem_freq_TQ0.01) %>%
    full_join(chem_freq_TQDetected) %>%
  full_join(unique(chem_detection[,c('CAS', 'chnm', 'Class')])) %>%
  # left_join(unique(chemicalSummary_bench_surface[c("CAS", "chnm", "Class")])) %>%
  rowwise() %>%
  mutate(chnm = as.character(chnm))


TQ_detection <- full_join(TQ_detection, chem_freq_all[,c('CAS', 'LabDetected')]) %>%
  mutate(UnknownTox = LabDetected - 
           sum(AboveTQ1, AboveTQ0.1, AboveTQ0.01, Detected, na.rm=T)) %>%
  select(-LabDetected)



TQ_detection$chnm[which(TQ_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))] <- 
  paste0(TQ_detection$chnm[which(TQ_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide", "Deg - Insecticide"))], "*")
TQ_detection$Class <- gsub("Deg - ", "", TQ_detection$Class)

TQ_detection <- TQ_detection %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  tidyr::gather(key=group, value=value, -chnm, -Class, -CAS) %>%
  mutate(chnm = factor(chnm, levels=levels(chem_detection$chnm)), 
         group = factor(group, levels=rev(c("AboveTQ1", "AboveTQ0.1", "AboveTQ0.01", "Detected", "UnknownTox"))))

TQ_detection$value[which(is.na(TQ_detection$value))] <- 0


# ########################################
# Plot frequency of chemical and EAR/TQ by chemical
# ########################################

#Plot barplot by chemical 
chemicalbyEAR2 <- ggplot(data=chem_detection, aes(x=chnm, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Class, space="free", scales="free") +
  labs(x='Chemical', y='Number of rivers', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        axis.title.x=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR2, labels = c("Unknown EAR", expression(paste("EAR < 10"^"-4")), expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  ggtitle(paste0("Surface water samples between ", date_filter[1], " and ", date_filter[2]))

print(chemicalbyEAR2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByEAR_Bychem_AllSummersamples.png")), plot = chemicalbyEAR2, height=6, width=8)


# 
# #Plot barplot by chemical horiztonal
# chemicalbyEAR2_horiztonal <- ggplot(data=chem_detection, aes(x=chnm, y=value, fill=group)) + 
#   geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
#   coord_flip() +
#   facet_grid(Class~., space="free", scales="free") +
#   labs(x='Chemical', y='Number of streams', fill = 'EAR') + 
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=.5),
#         axis.title.x=element_blank(),
#         legend.title = element_blank(),
#         axis.text.x = element_text(size=8),
#         axis.text.y = element_text(size=8)) + 
#   scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
#   scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
#   theme(legend.position = 'bottom') +
#   guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
#   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
# ggtitle("Surface Water Samples Jun 1 to Jul 31, 2016")
# 
# print(chemicalbyEAR2_horiztonal)
# 
# ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByEAR_Bychem2_horiztonal_surface.png")), plot = chemicalbyEAR2_horiztonal, height=6, width=4)
# 


#Plot barplot by chemical TQ
TQ_detection2 <- ggplot(data=TQ_detection , aes(x=chnm, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Class, space="free", scales="free") +
  labs(x='Chemical', y='Number of rivers', fill = 'TQ') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        axis.title.x=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR2, labels = c("Unknown TQ",  expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-1")), expression(paste("TQ > 1")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-1")))) +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  ggtitle(paste0("Surface water samples between ", date_filter[1], " and ", date_filter[2]))

print(TQ_detection2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByTQ_Bychem_AllSummersamples.png")), plot = TQ_detection2, height=6, width=8)





# #######################
# Summarize EAR by site
# #######################

site_freq_EAR0.01<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.01")

site_freq_EAR0.001<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.001 & EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

site_freq_EAR0.0001<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR < 0.001, EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")

site_freq_EARDetected<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR > 0 & EAR < 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "EARDetected")

site_table <- unique(chemicalSummary_surface[c('shortName', 'site')])
site_order_surface <- site_table$shortName[match(site_ID_order, site_table$site)]

site_detection <- full_join(site_freq_EAR0.01, site_freq_EAR0.001) %>%
  full_join(site_freq_EAR0.0001) %>%
  full_join(site_freq_EARDetected) %>%
  full_join(chnm_list[c('site', 'n')]) %>%
  left_join(unique(chemicalSummary_surface[c("site", "shortName", "Lake")])) %>%
  rowwise() %>%
  mutate(Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")),
         shortName = factor(shortName, site_order_surface),
         UnknownTox = n - sum(AboveEAR0.01, AboveEAR0.001, AboveEAR0.0001, EARDetected, na.rm=T)) %>%
  select(-n)

site_detection <- site_detection %>%
  arrange(Lake, desc(AboveEAR0.01),desc(AboveEAR0.001), desc(AboveEAR0.0001), desc(EARDetected), desc(UnknownTox)) %>%
  tidyr::gather(key=group, value=value, -site, -shortName, -Lake) %>%
  mutate(group = factor(group, levels=rev(c("AboveEAR0.01", "AboveEAR0.001", "AboveEAR0.0001", "EARDetected", 'UnknownTox'))))

site_detection$value[which(is.na(site_detection$value))] <- 0



#Plot barplot by site 
sitebyEAR2 <- ggplot(data=site_detection, aes(x=shortName, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Lake, space="free", scales="free") +
  labs(x='River', y='Number of chemicals', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR2, labels = c("Unknown EAR", expression(paste("EAR < 10"^"-4")), expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,70.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  ggtitle(paste0("Surface water between ", date_filter[1], " and ", date_filter[2]))

print(sitebyEAR2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_BySite_AllSummersamples.png")), plot = sitebyEAR2, height=4, width=6)




# #######################
# Summarize EAR by site. All chemicals. Not just those with POCIS
# #######################


chemicalSummary2_surface_all <- chemicalSummary_surface %>%
  filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(EAR>0) %>%
  # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  dplyr::group_by(site, chnm, date) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  dplyr::group_by(site, chnm) %>%
  summarize(EAR=mean(EAR, na.rm=T))



site_freq_EAR0.01_all<-chemicalSummary2_surface_all %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.01")

site_freq_EAR0.001_all<-chemicalSummary2_surface_all %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.001 & EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

site_freq_EAR0.0001_all<-chemicalSummary2_surface_all %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR < 0.001, EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")

site_freq_EARDetected_all<-chemicalSummary2_surface_all %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR > 0 & EAR < 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "EARDetected")

site_table <- unique(chemicalSummary_surface[c('shortName', 'site')])
site_order_surface <- site_table$shortName[match(site_ID_order, site_table$site)]

site_detection_all <- full_join(site_freq_EAR0.01_all, site_freq_EAR0.001_all) %>%
  full_join(site_freq_EAR0.0001_all) %>%
  full_join(site_freq_EARDetected_all) %>%
  left_join(unique(chemicalSummary_surface[c("site", "shortName", "Lake")])) %>%
  rowwise() %>%
  mutate(Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")),
         shortName = factor(shortName, site_order_surface))



site_detection_all <- site_detection_all %>%
  arrange(Lake, desc(AboveEAR0.01),desc(AboveEAR0.001), desc(AboveEAR0.0001), desc(EARDetected)) %>%
  tidyr::gather(key=group, value=value, -site, -shortName, -Lake) %>%
  mutate(group = factor(group, levels=rev(c("AboveEAR0.01", "AboveEAR0.001", "AboveEAR0.0001", "EARDetected"))))

site_detection_all$value[which(is.na(site_detection_all$value))] <- 0



#Plot barplot by site 
sitebyEAR2_all <- ggplot(data=site_detection_all, aes(x=shortName, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Lake, space="free", scales="free") +
  labs(x='Stream', y='Number of chemicals', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR, labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,45), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  ggtitle(paste0("Surface water (all chemicals) between ", date_filter[1], " and ", date_filter[2]))

print(sitebyEAR2_all)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_BySite2_allchemicals_watersamples.png")), plot = sitebyEAR2_all, height=4, width=6)






#########################################################################
#summarize and organize EAR and TQ for side by side horizontal boxplots
#some of this code should be eliminated as it is repeated above

chemicalSummary_surface #Sam's data from surface
chemicalSummary_bench_surface # Sam's surface data using custom benchmarks


chemicalSummary3_surface <- chemicalSummary_surface %>%
  filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  mutate(chnm = as.character(chnm))
chemicalSummary3_surface$chnm[which(chemicalSummary3_surface$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(chemicalSummary3_surface$chnm[which(chemicalSummary3_surface$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
chemicalSummary3_surface$Class <- gsub("Deg - ", "", chemicalSummary3_surface$Class)

chemicalSummary3_surface <- chemicalSummary3_surface %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  arrange((Class), desc(EAR)) 

chemicalSummary3_surface_maxbySite <- chemicalSummary3_surface %>%
  group_by(chnm, site, Class) %>%
  summarize(EAR = max(EAR)) %>%
  arrange(Class, chnm, (EAR)) 

chemicalSummary3_surface_medianAcrossSites <- chemicalSummary3_surface_maxbySite %>%
  group_by(chnm, Class) %>%
  summarize(EAR = median(EAR[which(EAR>0)])) 

chemicalSummary3_surface_medianAcrossSites$EAR[which(is.na(chemicalSummary3_surface_medianAcrossSites$EAR))] <-0
chemicalSummary3_surface_medianAcrossSites <- chemicalSummary3_surface_medianAcrossSites  %>% 
  arrange(desc(Class), (EAR))

chemicalSummary3_surface_maxbySite <- chemicalSummary3_surface_maxbySite %>%
  group_by() %>%
  mutate(chnm = factor(chnm, chemicalSummary3_surface_medianAcrossSites$chnm)) %>%
  arrange(Class, chnm, desc(EAR))

#summarize and organize TQ
chemicalSummary_bench3_surface <- chemicalSummary_bench_surface %>%
  filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  mutate(chnm = as.character(chnm))
chemicalSummary_bench3_surface $chnm[which(chemicalSummary_bench3_surface$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(chemicalSummary_bench3_surface $chnm[which(chemicalSummary_bench3_surface$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
chemicalSummary_bench3_surface$Class <- gsub("Deg - ", "", chemicalSummary_bench3_surface$Class)

chemicalSummary_bench3_surface  <- chemicalSummary_bench3_surface  %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  arrange((Class), desc(EAR)) 

chemicalSummary_bench3_surface_maxbySite <- chemicalSummary_bench3_surface  %>%
  group_by(chnm, site, Class) %>%
  summarize(EAR = max(EAR)) %>%
  arrange(Class, chnm, (EAR)) 

chemicalSummary_bench3_surface_medianAcrossSites <- chemicalSummary_bench3_surface_maxbySite %>%
  group_by(chnm, Class) %>%
  summarize(EAR = median(EAR[which(EAR>0)])) 

chemicalSummary_bench3_surface_medianAcrossSites$EAR[which(is.na(chemicalSummary_bench3_surface_medianAcrossSites$EAR))] <-0
chemicalSummary_bench3_surface_medianAcrossSites <- chemicalSummary_bench3_surface_medianAcrossSites  %>% 
  arrange(desc(Class), (EAR))

chemicalSummary_bench3_surface_maxbySite <- chemicalSummary_bench3_surface_maxbySite %>%
  group_by() %>%
  mutate(chnm = factor(chnm, chemicalSummary_bench3_surface_medianAcrossSites$chnm)) %>%
  arrange(Class, chnm, desc(EAR))


#combine chemical order
#use chemical order from passive data, exclude other chemicals

# chemorder2 <- rev(intersect(chemicalSummary3_medianAcrossSites$chnm, chemicalSummary_bench3_medianAcrossSites$chnm))
# 
# chemorder3 <- c(chemorder2, chemicalSummary3_medianAcrossSites$chnm[-which(chemicalSummary3_medianAcrossSites$chnm %in% chemorder2)], chemicalSummary_bench3_medianAcrossSites$chnm[-which(chemicalSummary_bench3_medianAcrossSites$chnm %in% chemorder2)])

chemicalSummary3_surface_maxbySite <- chemicalSummary3_surface_maxbySite %>%
  filter(chnm %in% chemorder3) %>%
  mutate(chnm = factor(chnm, rev(chemorder3))) %>%
  arrange(Class, chnm, desc(EAR))

chemicalSummary_bench3_surface_maxbySite <- chemicalSummary_bench3_surface_maxbySite %>%
  filter(chnm %in% chemorder3) %>%
  mutate(chnm = factor(chnm, rev(chemorder3))) %>%
  arrange(Class, chnm, desc(EAR))


EARbox_surface <- ggplot(chemicalSummary3_surface_maxbySite, aes(x=chnm, y=EAR)) +
  # geom_vline(xintercept = 4.5) +
  geom_hline(yintercept = .001, linetype=2) +
  geom_boxplot(aes(fill=Class)) +
  theme_bw() +
  scale_y_log10nice(name = 'max (EAR)') +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme(legend.position='none', axis.title.y=element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=8), axis.text.y = element_text(size=8)) +
  scale_fill_brewer(palette = 'Dark2')


TQbox_surface <- ggplot(chemicalSummary_bench3_surface_maxbySite, aes(x=chnm, y=EAR)) +
  # geom_vline(xintercept = c(4.5)) +
  geom_hline(yintercept = .1, linetype=2) + 
  geom_boxplot(aes(fill=Class)) +
  theme_bw() +
  scale_y_log10nice(name = 'max (TQ)') +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme(legend.position='none', axis.title.y=element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=8), axis.text.y = element_text(size=8)) +
  scale_fill_brewer(palette = 'Dark2')

  


box_by_box_surface <- grid.arrange(grobs=list(EARbox_surface, TQbox_surface), nrow=1, top='Surface water samples')


TQ_box_surface_withLegend <-TQbox_surface + 
  theme(legend.position='bottom', legend.title = element_blank()) + 
  guides(color = guide_legend(nrow = 1, title.position='top', title.hjust=0.5)) 

mylegend<-g_legend(TQ_box_surface_withLegend)

rm(TQ_box_surface_withLegend)


boxes_withLegend_surface<-grid.arrange(box_by_box_surface, mylegend, nrow=2, heights=c(15,1))


ggsave(file_out(file.path(path_to_data, "Figures/SideBoxes_EARandTQ_ByChemical_watersamples.png")), plot = boxes_withLegend_surface, height=6, width=7)


