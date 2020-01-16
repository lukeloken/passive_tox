

#Calculate chemical summary for passive pesticide data
# tox_list<- create_toxEval(file_in(file.path(path_to_data, "Data/PassiveForToxEval.xlsx")))
chemicalSummary
chemicalSummary_allpocis

#create table by chemical and by site
#these are detections above the POCIS minimum detection
#Only useful for presence/absence
chem_freq_allpocis<-chemicalSummary_allpocis %>%
  filter(EAR > 0) %>%
  group_by(chnm, site) %>%
  summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "Detected")

site_freq_allpocis<-chemicalSummary_allpocis %>%
  filter(EAR > 0) %>%
  group_by(site, chnm) %>%
  dplyr::summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
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
  left_join(unique(chemicalSummary[c("chnm", "Class")])) %>%
  rowwise() %>%
  mutate(EAR_total = sum(AboveEAR0.001, AboveEAR0.0001, na.rm=T),
         OnlyDetected = Detected - EAR_total) %>%
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


chem_freq_EAR<-chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(chnm, site) %>%
  summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "Detected")

chem_freq_EAR0.001<-chemicalSummary %>%
  filter(EAR > 0.001) %>%
  group_by(chnm, site) %>%
  summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "AboveEARThreshold")

chem_detection <- full_join(chem_freq_EAR, chem_freq_EAR0.001) %>%
  left_join(unique(chemicalSummary[c("chnm", "Class")])) %>%
  arrange(Class, chnm)






#Using the'select' data evaluate EAR levels and summarize
#Calculate EAR max for each site and chemical
EAR_max<-chemicalSummary %>%
  filter(EAR > 0) %>%
  mutate(chnm = as.character(chnm)) %>%
  group_by(chnm, site, shortName, Class) %>%
  summarize (EAR = max(EAR)) #%>%
  #left_join(unique(chemicalSummary[c("chnm", "Class")]))

#Add * to deg compounds and change class
EAR_max$chnm[which(EAR_max$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(EAR_max$chnm[which(EAR_max$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
EAR_max$Class <- gsub("Deg - ", "", EAR_max$Class)

#Categorize by EAR level. These will translate into plotting colors 
EAR_max$Group <- rep(NA,(nrow(EAR_max)))
EAR_max$Group[which(EAR_max$EAR<10^-4)] <- "mdl to 0.0001"
EAR_max$Group[which(EAR_max$EAR>=10^-4 & EAR_max$EAR<10^-3)] <- "0.0001 to 0.001"
EAR_max$Group[which(EAR_max$EAR>=10^-3)] <- "> 0.001" 
# EAR_max$Group[which(EAR_max$EAR>=10^-2)] <- "> 0.01" 


EAR_max <- EAR_max %>%
  mutate(Group = factor(Group, c( 'mdl to 0.0001', "0.0001 to 0.001",  "> 0.001")),
         Class = factor(Class, c("Herbicide", "Fungicide", "Insecticide"))) %>%
  arrange(chnm, Class, Group)

#Sorting by Class then group so data plot in order of toxicity
EAR_table <- EAR_max %>%
  group_by(Group, chnm) %>%
  summarize (n = n()) %>%
  left_join(unique(chemicalSummary[c("chnm", "Class")])) %>%
  dplyr::arrange(Class, desc(Group), -n)


chem_order <- unique(EAR_table$chnm)

EAR_max <- EAR_max %>%
  group_by() %>%
  mutate(chnm = factor(chnm, chem_order)) %>%
  dplyr::arrange(chnm, Group)

#Ordering for sites. This sort order will be used to plot chemicals by site
Site_table <- EAR_max %>%
  group_by(site, Group) %>%
  summarize(n=n()) %>%
  dplyr::arrange(desc(Group), desc(n))

# site_order <- unique(Site_table$site)
site_order <- unique(chemicalSummary$shortName)

EAR_site <- EAR_max %>%
  left_join(unique(chemicalSummary[c("site", "shortName", "Lake")])) %>%
  group_by() %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  dplyr::arrange(site) %>%
  mutate(Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")))

Site_summary <- EAR_site %>%
  group_by(shortName, Group) %>%
  summarize (n=n()) %>%
  spread(key=Group, value = n) %>%
  mutate( Total = sum(`mdl to 0.0001`, `0.0001 to 0.001`, `> 0.001`, na.rm=T))



chem_freq_EAR<-chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(chnm, site) %>%
  summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "Detected")
             

chem_freq_EAR0.001<-chemicalSummary %>%
  filter(EAR > 0.001) %>%
  group_by(chnm, site) %>%
  summarize(EAR_mean = mean(EAR, na.rm=T)) %>%
  tally(name = "AboveEARThreshold")

chem_detection <- full_join(chem_freq_EAR, chem_freq_EAR0.001) %>%
  left_join(unique(chemicalSummary[c("chnm", "Class")])) %>%
  arrange(Class, chnm)


ggplot(chem_detection, (aes(x=chnm, y=Detected, fill=Class))) + 
  geom_col() +
  coord_flip() +
  labs(x='Chemical', y='Number of sites detected') +
  theme_bw() +
  scale_y_continuous(limits=c(0,15), expand=c(0,0)) #+
  # theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) + 
  # theme(axis.text.y = element_text(angle=90, vjust=0.5, hjust=1))

ggplot(chem_detection, (aes(x=chnm, y=AboveEARThreshold, fill=Class))) + 
  geom_col() +
  coord_flip() +
  labs(x='Chemical', y='Number of sites with EAR > 0.001') + 
  theme_bw() + 
  scale_y_continuous(limits=c(0,15), expand=c(0,0)) #+
  # theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))


# Stacked barplot by EAR
chemicalbyEAR <- ggplot(data=EAR_max, aes(x=chnm, fill=Group)) + 
  geom_bar(color = 'grey', width=.8, size=.1) + 
  coord_flip() +
  labs(x='Chemical', y='Number of sites', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5)) + 
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits=c(0,15), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='top', title.hjust=0.5, reverse=T))

# ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR.png")), plot = chemicalbyEAR, height=5, width=5)



# Stacked barplot by EAR
chemicalbyEAR_facet <- ggplot(data=EAR_max, aes(x=chnm, fill=Group)) + 
  geom_bar(color = 'grey', width=.8, size=.1) + 
  coord_flip() +
  labs(x='Chemical', y='Number of sites', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5)) + 
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits=c(0,15), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='top', title.hjust=0.5, reverse=T)) + 
  facet_grid(~Class)

# ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_facet.png")), plot = chemicalbyEAR_facet, height=5, width=5)


#Playing with grouping by classes. 
chemicalbyEAR_byClass <- ggplot(data=EAR_max, aes(x=chnm, fill=Group)) + 
  geom_bar(color = 'grey', width=.8, size=.1) + 
  # geom_bar(position = "dodge") +
  facet_grid(~Class, space="free", scales="free") +
  labs(x='Chemical', y='Number of sites', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        axis.title.x=element_blank()) + 
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits=c(0,15.5), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))


ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_3Classes.png")), plot = chemicalbyEAR_byClass, height=4, width=6)


EAR_max_sortedbyChem <- EAR_max %>%
  group_by() %>%
  arrange()

#Playing with grouping by site. 
chemicalbyEAR_bysite <- ggplot(data=EAR_site, aes(x=shortName, fill=Group)) + 
  geom_bar(color = 'grey', width=.8, size=.1) + 
  # geom_bar(position = "dodge") +
  facet_grid(~Lake, space="free", scales="free") +
  labs(x='Chemical', y='Number of chemicals', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        axis.title.x=element_blank()) + 
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits=c(0,28), expand=c(0,0)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

print(chemicalbyEAR_bysite)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_BySite.png")), plot = chemicalbyEAR_bysite, height=4, width=5)




#Plot using identity rather than count. 
# Work in progress, but this will be more powerful




chem_freq_allpocis$EAR2 <- chem_freq_allpocis$Detected-4

chem_freq_allpocis$EAR2[which(chem_freq_allpocis$EAR2<0)]<-0
gathertest<-tidyr::gather(chem_freq_allpocis, key=group, value=value, -chnm)

#Playing with using stat=identify to plot height. 
chemicalbyEAR_bysite2 <- ggplot(data=chem_detection, aes(x=chnm, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Class, space="free", scales="free") +
  labs(x='Chemical', y='Number of chemicals', fill = 'EAR') + 
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

print(chemicalbyEAR_bysite2)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_BySite2.png")), plot = chemicalbyEAR_bysite2, height=4, width=6)

#Map figure
# make_tox_map(chemicalSummary, 
#              chem_site = tox_list$chem_site,
#              category = 'Chemical Class',
#              mean_logic = FALSE)
#   
