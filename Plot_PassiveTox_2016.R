

#Calculate chemical summary for passive pesticide data
# tox_list<- create_toxEval(file_in(file.path(path_to_data, "Data/PassiveForToxEval.xlsx")))

#R objects
chemicalSummary #Select data from POCIS 
chemicalSummary_allpocis #all data from POCIS
chemicalSummary_bench # Select data using custom benchmarks

#Check out the standard boxplots
plot_tox_boxplots(chemicalSummary_bench, 
                  category = "Chemical", 
                  sum_logic = FALSE,
                  x_label = "Toxicity Quotient")

plot_tox_boxplots(chemicalSummary, 
                  category = "Chemical", 
                  sum_logic = FALSE)



# create tables by chemical and by sites
# these are detections above the POCIS minimum detection
# allpocis data only useful for presence/absence
# chem_freq_allpocis<-chemicalSummary_allpocis %>%
#   filter(EAR > 0) %>%
#   group_by(chnm, site) %>%
#   summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
#   tally(name = "Detected")

chem_freq_allpocis <- tox_list_allpocis$chem_data %>%
  group_by(CAS, chnm) %>%
  select(chnm, CAS, Value) %>%
  filter(Value > 0) %>%
  tally(name = 'Detected') %>%
  left_join(unique(tox_list_allpocis$chem_data[,c('chnm', 'CAS')])) %>%
  distinct() %>%
  arrange(chnm) %>%
  left_join(unique(tox_list_allpocis$chem_info[,c('CAS', 'Class')])) %>%
  # filter(Class !="") %>%
  arrange(Class, desc(Detected))


chemicalSummary2 <- chemicalSummary %>%
  filter(EAR>0) %>%
  dplyr::group_by(site, CAS) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  arrange(desc(EAR)) %>%
  filter(CAS %in% chem_freq_allpocis$CAS)

chem_freq_EAR0.01<-chemicalSummary2 %>%
  group_by(CAS, site) %>%
  filter(EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.01")

chem_freq_EAR0.001<-chemicalSummary2 %>%
  group_by(CAS, site) %>%
  filter(EAR >= 0.001 & EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

chem_freq_EAR0.0001<-chemicalSummary2 %>%
  group_by(CAS, site) %>%
  filter(EAR < 0.001 & EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")

chem_freq_EARAboveZero<-chemicalSummary2 %>%
  group_by(CAS, site) %>%
  filter(EAR < 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "BelowEAR0.0001")

chem_detection <- full_join(chem_freq_EAR0.01, chem_freq_EAR0.001) %>%
  full_join(chem_freq_EAR0.0001) %>%
  full_join(chem_freq_allpocis) %>%
  full_join(chem_freq_EARAboveZero) %>%
  left_join(unique(chemicalSummary_allpocis[c("chnm", "Class")])) %>%
  rowwise() %>%
  mutate(EAR_total = sum(AboveEAR0.01, AboveEAR0.001, AboveEAR0.0001, BelowEAR0.0001, na.rm=T),
         OnlyDetected = Detected - EAR_total,
         chnm = as.character(chnm)) %>%
  select(-Detected, -EAR_total )


chem_detection$chnm[which(chem_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(chem_detection$chnm[which(chem_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
chem_detection$Class <- gsub("Deg - ", "", chem_detection$Class)

chem_detection <- chem_detection %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  arrange(Class, desc(AboveEAR0.01), desc(AboveEAR0.001), desc(AboveEAR0.0001), desc(BelowEAR0.0001), desc(OnlyDetected)) %>%
  tidyr::gather(key=group, value=value, -chnm, -Class) %>%
  mutate(chnm = factor(chnm, levels=unique(chnm)), 
         group = factor(group, levels=rev(c("AboveEAR0.01", "AboveEAR0.001", "AboveEAR0.0001", "BelowEAR0.0001", "OnlyDetected"))))

chem_detection$value[which(is.na(chem_detection$value))] <- 0

#Reverse chemical order for horizontal barplot
chem_detection$chnm2 = factor(chem_detection$chnm, rev(levels(chem_detection$chnm)))


##########################################

#Make similar table for TQ
chemicalSummary_bench2 <- chemicalSummary_bench %>%
  filter(EAR>0) %>%
  dplyr::group_by(site, CAS) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  arrange(desc(EAR)) %>%
  filter(CAS %in% chem_freq_allpocis$CAS)

chem_freq_TQ1<-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR >= 1) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ1")

chem_freq_TQ0.1<-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR >= 0.1 & EAR < 1) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ0.1")

chem_freq_TQ0.01<-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR < 0.1 & EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ0.01")

chem_freq_TQAbove0 <-chemicalSummary_bench2 %>%
  group_by(CAS, site) %>%
  filter(EAR > 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "BelowTQ0.01")



TQ_detection <- full_join(chem_freq_TQ1, chem_freq_TQ0.1) %>%
  full_join(chem_freq_TQ0.01) %>%
  full_join(chem_freq_TQAbove0) %>%
  full_join(chem_freq_allpocis) %>%
  left_join(unique(chemicalSummary_allpocis[c("chnm", "Class")])) %>%
  rowwise() %>%
  mutate(EAR_total = sum(AboveTQ1, AboveTQ0.01, AboveTQ0.1, na.rm=T),
         OnlyDetected = Detected - EAR_total,
         chnm = as.character(chnm)) %>%
  select(-Detected, -EAR_total )


TQ_detection$chnm[which(TQ_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(TQ_detection$chnm[which(TQ_detection$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
TQ_detection$Class <- gsub("Deg - ", "", TQ_detection$Class)

TQ_detection <- TQ_detection %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  tidyr::gather(key=group, value=value, -chnm, -Class) %>%
  mutate(chnm = factor(chnm, levels=levels(chem_detection$chnm)), 
         group = factor(group, levels=rev(c("AboveTQ1", "AboveTQ0.1", "AboveTQ0.01", "OnlyDetected"))))

TQ_detection$value[which(is.na(TQ_detection$value))] <- 0

TQ_detection$chnm2 = factor(TQ_detection$chnm, rev(levels(TQ_detection$chnm)))


# ########################################
# Plot frequency of chemical and EAR/TQ by chemical
# ########################################

#Make color ramp
colors_EAR <- brewer.pal(n = 9, name = "YlOrRd")[c(2,4,7,9)]

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
  scale_fill_manual(values = colors_EAR, labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

print(chemicalbyEAR2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByEAR_Bychem2.png")), plot = chemicalbyEAR2, height=4, width=6)


#Plot barplot by chemical horiztonal
chemicalbyEAR2_horiztonal <- ggplot(data=chem_detection, aes(x=chnm2, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  coord_flip() +
  facet_grid(Class~., space="free", scales="free") +
  labs(x='Chemical', y='Number of streams', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        # axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        # axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR, labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(size=8, angle=0, vjust=1, hjust=.5))

print(chemicalbyEAR2_horiztonal)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByEAR_Bychem2_horiztonal.png")), plot = chemicalbyEAR2_horiztonal, height=6, width=4)



#Plot barplot by chemical TQ
TQ_detection2 <- ggplot(data=TQ_detection , aes(x=chnm, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Class, space="free", scales="free") +
  labs(x='Chemical', y='Number of streams', fill = 'TQ') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        axis.title.x=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR, labels = c("Detected", expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-1")), expression(paste("TQ > 1")))) +
  #   scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-1")), expression(paste("TQ > 1")))) +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

print(TQ_detection2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByTQ_Bychem2.png")), plot = TQ_detection2, height=4, width=6)



#Plot barplot by chemical horiztonal
TQ_detection2_horiztonal <- ggplot(data=TQ_detection, aes(x=chnm2, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  coord_flip() +
  facet_grid(Class~., space="free", scales="free") +
  labs(x='Chemical', y='Number of streams', fill = 'TQ') + 
  theme_bw() +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) +
  # scale_y_reverse(limits=c(15.5,0), expand=c(0,0)) + 
  # scale_x_discrete(position='top') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        # axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        # axis.text.x = element_text(size=8, angle=0, hjust=1),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR, labels = c("Detected", expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-1")), expression(paste("TQ > 1")))) +
  #   scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-1")), expression(paste("TQ > 1")))) +
  # scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(size=8, angle=0, vjust=1, hjust=.5))

print(TQ_detection2_horiztonal)

ggsave(file_out(file.path(path_to_data, "Figures/StackBar_ByTQ_Bychem2_horiztonal.png")), plot = chemicalbyEAR2_horiztonal, height=6, width=4)


#Combine EAR and TQ horizontal plots

EAR_box1 <- chemicalbyEAR2_horiztonal +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.title=element_text(size=12, hjust=.5, vjust=0),
        legend.text = element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.key.size = unit(.5, 'line')) +
  guides(fill=guide_legend(ncol=1, reverse=T, label.hjust=0)) + 
  ggtitle('EAR')

EAR_box1

TQ_box2 <- TQ_detection2_horiztonal +
  theme(axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=.5, vjust=0),
        legend.text = element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.key.size = unit(.5, 'line')) +
  guides(fill=guide_legend(ncol=1, reverse=T, label.hjust=0)) + 
  ggtitle('TQ')
        

TQ_box2

png(file.path(path_to_data, "Figures/SideBar_EAR_TQ_byChem.png"), height=6, width=5, units='in', res=400)

grid.newpage()
boxes_EAR<-grid.draw(cbind(ggplotGrob(EAR_box1), ggplotGrob(TQ_box2), size = "first"))

dev.off()

# #######################
# Summarize EAR by site
# #######################

site_freq_allpocis<-chemicalSummary_allpocis %>%
  filter(EAR > 0) %>%
  group_by(site, chnm) %>%
  dplyr::summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "Detected")

site_freq_EAR0.01<-chemicalSummary2 %>%
  group_by(site, chnm) %>%
  filter(EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.01")


site_freq_EAR0.001<-chemicalSummary2 %>%
  group_by(site, chnm) %>%
  filter(EAR >= 0.001 & EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

site_freq_EAR0.0001<-chemicalSummary2 %>%
  group_by(site, chnm) %>%
  filter(EAR < 0.001 & EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")



site_detection <- full_join(site_freq_EAR0.01, site_freq_EAR0.001) %>%
  full_join(site_freq_EAR0.0001) %>%
  full_join(site_freq_allpocis) %>%
  left_join(unique(chemicalSummary[c("site", "shortName", "Lake")])) %>%
  rowwise() %>%
  mutate(EAR_total = sum(AboveEAR0.01, AboveEAR0.001, AboveEAR0.0001, na.rm=T),
         OnlyDetected = Detected - EAR_total,
         Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")),
         shortName = factor(shortName, site_order)) %>%
  dplyr::select(-Detected, -EAR_total )



site_detection <- site_detection %>%
  arrange(Lake, desc(AboveEAR0.01), desc(AboveEAR0.001), desc(AboveEAR0.0001), desc(OnlyDetected)) %>%
  tidyr::gather(key=group, value=value, -site, -shortName, -Lake) %>%
  mutate(group = factor(group, levels=rev(c("AboveEAR0.01", "AboveEAR0.001", "AboveEAR0.0001", "OnlyDetected"))))

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
  scale_fill_manual(values = colors_EAR, labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(limits=c(0,35.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

print(sitebyEAR2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_BySite2.png")), plot = sitebyEAR2, height=4, width=6)








#########################################################################
#summarize and organize EAR and TQ for side by side horizontal boxplots
#some of this code should be eliminated as it is repeated above
chemicalSummary3_passive <- chemicalSummary %>%
  mutate(chnm = as.character(chnm))
chemicalSummary3_passive$chnm[which(chemicalSummary3_passive$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(chemicalSummary3_passive$chnm[which(chemicalSummary3_passive$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
chemicalSummary3_passive$Class <- gsub("Deg - ", "", chemicalSummary3_passive$Class)

chemicalSummary3_passive <- chemicalSummary3_passive %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  arrange((Class), desc(EAR)) 

chemicalSummary3_passive_maxbySite <- chemicalSummary3_passive %>%
  group_by(chnm, site, Class) %>%
  summarize(EAR = max(EAR)) %>%
  arrange(Class, chnm, (EAR)) 

chemicalSummary3_passive_medianAcrossSites <- chemicalSummary3_passive_maxbySite %>%
  group_by(chnm, Class) %>%
  summarize(EAR = median(EAR[which(EAR>0)])) 

chemicalSummary3_passive_medianAcrossSites$EAR[which(is.na(chemicalSummary3_passive_medianAcrossSites$EAR))] <-0
chemicalSummary3_passive_medianAcrossSites <- chemicalSummary3_passive_medianAcrossSites  %>% 
  arrange(desc(Class), (EAR))

chemicalSummary3_passive_maxbySite <- chemicalSummary3_passive_maxbySite %>%
  group_by() %>%
  mutate(chnm = factor(chnm, chemicalSummary3_passive_medianAcrossSites$chnm)) %>%
  arrange(Class, chnm, desc(EAR))

#summarize and organize TQ
chemicalSummary_bench3_passive <- chemicalSummary_bench %>%
  mutate(chnm = as.character(chnm))
chemicalSummary_bench3_passive $chnm[which(chemicalSummary_bench3_passive$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(chemicalSummary_bench3_passive $chnm[which(chemicalSummary_bench3_passive$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
chemicalSummary_bench3_passive$Class <- gsub("Deg - ", "", chemicalSummary_bench3_passive$Class)

chemicalSummary_bench3_passive  <- chemicalSummary_bench3_passive  %>%
  mutate(Class = factor(Class, c('Herbicide', 'Fungicide', 'Insecticide'))) %>%
  arrange((Class), desc(EAR)) 

chemicalSummary_bench3_passive_maxbySite <- chemicalSummary_bench3_passive  %>%
  group_by(chnm, site, Class) %>%
  summarize(EAR = max(EAR)) %>%
  arrange(Class, chnm, (EAR)) 

chemicalSummary_bench3_passive_medianAcrossSites <- chemicalSummary_bench3_passive_maxbySite %>%
  group_by(chnm, Class) %>%
  summarize(EAR = median(EAR[which(EAR>0)])) 

chemicalSummary_bench3_passive_medianAcrossSites$EAR[which(is.na(chemicalSummary_bench3_passive_medianAcrossSites$EAR))] <-0
chemicalSummary_bench3_passive_medianAcrossSites <- chemicalSummary_bench3_passive_medianAcrossSites  %>% 
  arrange(desc(Class), (EAR))

chemicalSummary_bench3_passive_maxbySite <- chemicalSummary_bench3_passive_maxbySite %>%
  group_by() %>%
  mutate(chnm = factor(chnm, chemicalSummary_bench3_passive_medianAcrossSites$chnm)) %>%
  arrange(Class, chnm, desc(EAR))


#combine chemical order
chemorder1 <- rev(chemicalSummary3_passive_medianAcrossSites$chnm)

chemorder2 <- rev(intersect(chemicalSummary3_passive_medianAcrossSites$chnm, chemicalSummary_bench3_passive_medianAcrossSites$chnm))

chemorder3 <- c(chemorder2, chemicalSummary3_passive_medianAcrossSites$chnm[-which(chemicalSummary3_passive_medianAcrossSites$chnm %in% chemorder2)], chemicalSummary_bench3_passive_medianAcrossSites$chnm[-which(chemicalSummary_bench3_passive_medianAcrossSites$chnm %in% chemorder2)])

chemicalSummary3_passive_maxbySite <- chemicalSummary3_passive_maxbySite %>%
  mutate(chnm = factor(chnm, rev(chemorder3))) %>%
  arrange(Class, chnm, desc(EAR))

chemicalSummary_bench3_passive_maxbySite <- chemicalSummary_bench3_passive_maxbySite %>%
  mutate(chnm = factor(chnm, rev(chemorder3))) %>%
  arrange(Class, chnm, desc(EAR))


EARbox_passive <- ggplot(chemicalSummary3_passive_maxbySite, aes(x=chnm, y=EAR)) +
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


TQbox_passive <- ggplot(chemicalSummary_bench3_passive_maxbySite, aes(x=chnm, y=EAR)) +
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

  


box_by_box_passive <- grid.arrange(grobs=list(EARbox_passive, TQbox_passive), nrow=1, top='POCIS passive samples')


TQ_box_passive_withLegend <-TQbox_passive + 
  theme(legend.position='bottom', legend.title = element_blank()) + 
  guides(color = guide_legend(nrow = 1, title.position='top', title.hjust=0.5)) 

mylegend<-g_legend(TQ_box_passive_withLegend)

rm(TQ_box_withLegend_passive)


boxes_withLegend_passive<-grid.arrange(box_by_box_passive, mylegend, nrow=2, heights=c(15,1))


ggsave(file_out(file.path(path_to_data, "Figures/SideBoxes_Passive_EARandTQ_ByChemical.png")), plot = boxes_withLegend_passive, height=6, width=7)


