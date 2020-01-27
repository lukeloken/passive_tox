

#Calculate chemical summary for passive pesticide data
# tox_list<- create_toxEval(file_in(file.path(path_to_data, "Data/PassiveForToxEval.xlsx")))
chemicalSummary
chemicalSummary_allpocis

# create tables by chemical and by sites
# these are detections above the POCIS minimum detection
# allpocis data only useful for presence/absence
chem_freq_allpocis<-chemicalSummary_allpocis %>%
  filter(EAR > 0) %>%
  group_by(chnm, site) %>%
  summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "Detected")

chemicalSummary2 <- chemicalSummary %>%
  filter(EAR>0) %>%
  dplyr::group_by(site, chnm) %>% 
  summarize(EAR = max(EAR, na.rm=T))

chem_freq_EAR0.001<-chemicalSummary2 %>%
  group_by(chnm, site) %>%
  filter(EAR >= 0.001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

chem_freq_EAR0.0001<-chemicalSummary2 %>%
  group_by(chnm, site) %>%
  filter(EAR < 0.001 & EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")

chem_detection <- full_join(chem_freq_EAR0.001, chem_freq_EAR0.0001) %>%
  full_join(chem_freq_allpocis) %>%
  left_join(unique(chemicalSummary_allpocis[c("chnm", "Class")])) %>%
  rowwise() %>%
  mutate(EAR_total = sum(AboveEAR0.001, AboveEAR0.0001, na.rm=T),
         OnlyDetected = Detected - EAR_total,
         chnm = as.character(chnm)) %>%
  select(-Detected, -EAR_total )


chem_detection$chnm[which(chem_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(chem_detection$chnm[which(chem_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
chem_detection$Class <- gsub("Deg - ", "", chem_detection$Class)

chem_detection <- chem_detection %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  arrange(Class, desc(AboveEAR0.001), desc(AboveEAR0.0001), desc(OnlyDetected)) %>%
  tidyr::gather(key=group, value=value, -chnm, -Class) %>%
  mutate(chnm = factor(chnm, levels=unique(chnm)), 
         group = factor(group, levels=rev(c("AboveEAR0.001", "AboveEAR0.0001", "OnlyDetected"))))

chem_detection$value[which(is.na(chem_detection$value))] <- 0

##########################################


#Plot barplot by chemical 
chemicalbyEAR2 <- ggplot(data=chem_detection, aes(x=chnm, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Class, space="free", scales="free") +
  labs(x='Chemical', y='Number of streams', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        axis.title.x=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

print(chemicalbyEAR2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByEAR_Bychem2.png")), plot = chemicalbyEAR2, height=4, width=6)



#Plot barplot by chemical horiztonal
chemicalbyEAR2_horiztonal <- ggplot(data=chem_detection, aes(x=chnm, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  coord_flip() +
  facet_grid(Class~., space="free", scales="free") +
  labs(x='Chemical', y='Number of streams', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        axis.title.x=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

print(chemicalbyEAR2_horiztonal)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByEAR_Bychem2_horiztonal.png")), plot = chemicalbyEAR2_horiztonal, height=6, width=4)





#Summarize by site

site_freq_allpocis<-chemicalSummary_allpocis %>%
  filter(EAR > 0) %>%
  group_by(site, chnm) %>%
  dplyr::summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "Detected")


site_freq_EAR0.001<-chemicalSummary2 %>%
  group_by(site, chnm) %>%
  filter(EAR >= 0.001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

site_freq_EAR0.0001<-chemicalSummary2 %>%
  group_by(site, chnm) %>%
  filter(EAR < 0.001 & EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")

site_order <- unique(chemicalSummary$shortName)


site_detection <- full_join(site_freq_EAR0.001, site_freq_EAR0.0001) %>%
  full_join(site_freq_allpocis) %>%
  left_join(unique(chemicalSummary[c("site", "shortName", "Lake")])) %>%
  rowwise() %>%
  mutate(EAR_total = sum(AboveEAR0.001, AboveEAR0.0001, na.rm=T),
         OnlyDetected = Detected - EAR_total,
         Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")),
         shortName = factor(shortName, site_order)) %>%
  select(-Detected, -EAR_total )



site_detection <- site_detection %>%
  arrange(Lake, desc(AboveEAR0.001), desc(AboveEAR0.0001), desc(OnlyDetected)) %>%
  tidyr::gather(key=group, value=value, -site, -shortName, -Lake) %>%
  mutate(group = factor(group, levels=rev(c("AboveEAR0.001", "AboveEAR0.0001", "OnlyDetected"))))

site_detection$value[which(is.na(site_detection$value))] <- 0



#Plot barplot by site 
sitebyEAR2 <- ggplot(data=site_detection, aes(x=shortName, y=value, fill=group)) + 
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
  scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,35.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

print(sitebyEAR2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_BySite2.png")), plot = sitebyEAR2, height=4, width=6)






