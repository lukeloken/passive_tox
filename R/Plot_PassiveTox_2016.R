

#Calculate chemical summary for passive pesticide data

tox_list<- create_toxEval(file_in(file.path(path_to_data, "Data/PassiveForToxEval.xlsx")))

chemicalSummary


EAR_max<-chemicalSummary %>%
  filter(EAR > 0) %>%
  mutate(chnm = as.character(chnm)) %>%
  group_by(chnm, site) %>%
  summarize (EAR = max(EAR)) %>%
  left_join(unique(chemicalSummary[c("chnm", "Class")]))

#Add * to deg compounds and change class
EAR_max$chnm[which(EAR_max$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))] <- 
  paste0(EAR_max$chnm[which(EAR_max$Class %in% c("Deg - Fungicide", "Deg - Herbicide"))], "*")
EAR_max$Class <- gsub("Deg - ", "", EAR_max$Class)

EAR_max$Group <- rep(NA,(nrow(EAR_max)))
EAR_max$Group[which(EAR_max$EAR<10^-4)] <- "mdl to 0.0001"
EAR_max$Group[which(EAR_max$EAR>=10^-4 & EAR_max$EAR<10^-3)] <- "0.0001 to 0.001"
EAR_max$Group[which(EAR_max$EAR>=10^-3)] <- "> 0.001" 
# EAR_max$Group[which(EAR_max$EAR>=10^-2)] <- "> 0.01" 

EAR_max <- EAR_max %>%
  mutate(Group = factor(Group, c( 'mdl to 0.0001', "0.0001 to 0.001",  "> 0.001")),
         Class = factor(Class, c("Herbicide", "Fungicide", "Insecticide"))) %>%
  arrange(chnm, Class, Group)

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

#Ordering for sites
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

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR.png")), plot = chemicalbyEAR, height=5, width=5)



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

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_ByEAR_facet.png")), plot = chemicalbyEAR_facet, height=5, width=5)


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

#Map figure


make_tox_map(chemicalSummary, 
             chem_site = tox_list$chem_site,
             category = 'Chemical Class',
             mean_logic = FALSE)
  
