

library(ggrepel)

#passive mdl tox_list
tox_list_passive_mdl <- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016MDL_ToxEval.xlsx")))

#Combine exclusions
exclusions <- tox_list_passive_mdl$exclusions %>%
  bind_rows(tox_list_surface$exclusions) %>%
  distinct()

tox_list_passive_mdl$exclusions <- exclusions

#Rename to match passive
tox_list_passive_mdl$chem_site <- tox_list_passive_mdl$chem_site %>%
  left_join(tox_list_passive_mdl$chem_site[,c("SiteID", 'site_grouping')]) 

#Manually change chemical classes and State IDs
tox_list_passive_mdl$chem_info$Class[tox_list_passive_mdl$chem_info$CAS == '78-48-8'] = 'Herbicide'

#Flag some ACCs
ACClong <- get_ACC(tox_list_surface$chem_info$CAS)
# ACClong <- remove_flags(ACClong, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical","ACCLessThan"))
ACClong <- remove_flags(ACClong)


cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))

chemicalSummary_passive_mdl <- get_chemical_summary(tox_list_passive_mdl, 
                                                ACClong, 
                                                filtered_ep)




#surface mdl tox_list
tox_list_surface_mdl <- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides_MDL.xlsx")))

#Combine exclusions
exclusions <- tox_list_surface_mdl$exclusions %>%
  bind_rows(tox_list_surface_mdl$exclusions) %>%
  distinct()

tox_list_surface_mdl$exclusions <- exclusions

#Rename to match passive
tox_list_surface_mdl$chem_site <- tox_list_surface_mdl$chem_site %>%
  left_join(tox_list_surface_mdl$chem_site[,c("SiteID", 'site_grouping')]) 

#Manually change chemical classes and State IDs
tox_list_surface_mdl$chem_info$Class[tox_list_surface_mdl$chem_info$CAS == '78-48-8'] = 'Herbicide'

#Flag some ACCs
ACClong <- get_ACC(tox_list_surface_mdl$chem_info$CAS)
# ACClong <- remove_flags(ACClong, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical","ACCLessThan"))
ACClong <- remove_flags(ACClong)


cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))

chemicalSummary_surface_mdl <- get_chemical_summary(tox_list_surface_mdl, 
                                                    ACClong, 
                                                    filtered_ep)


#Summarize mdls

#Add groups to summary
surface_mdl_summary <- chemicalSummary_surface_mdl %>%
  group_by(CAS, chnm) %>%
  select(CAS, chnm, EAR) %>%
  summarize(EAR = max(EAR, na.rm=T))
  
#Add groups to summary
passive_mdl_summary <- chemicalSummary_passive_mdl %>%
  group_by(CAS, chnm) %>%
  select(CAS, chnm, EAR) %>%
  summarize(EAR = max(EAR, na.rm=T))



#Plot mdls against each other

mdl_EAR_join <- surface_mdl_summary %>%
  rename(EAR_surface = EAR) %>%
  full_join(passive_mdl_summary) %>%
  rename(EAR_passive = EAR) %>%
  drop_na(EAR_passive, EAR_surface) %>%
  left_join(unique(chemicalSummary_passive_mdl[,c('CAS', 'Class')])) %>%
  drop_na(Class)

lmmodel <- lm(log10(mdl_EAR_join$EAR_passive) ~ log10(mdl_EAR_join$EAR_surface))

mdlcompare <- ggplot(mdl_EAR_join, aes(x=EAR_surface, y=EAR_passive, label=chnm)) +
  geom_point(col='red', size=2) +
  geom_abline(linetype='dashed') +
  geom_vline(xintercept = .001) + 
  geom_hline(yintercept = .001) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # xlim(range(c(join_mdl$Surface_mdl, join_mdl$Passive_mdl))) +
  # ylim(range(c(join_mdl$Surface_mdl, join_mdl$Passive_mdl))) +
  # xlim(c(0,.02)) +
  # ylim(c(0,.02)) +
  labs(x=expression(paste('Surface mdl (', mu, 'g L'^'-1', ')')), 
       y=expression(paste('Passive mdl (', mu, 'g L'^'-1', ')'))) +
  geom_text_repel(size=3, color='black', segment.color='black', segment.alpha=.5, segment.size=.5, alpha=.7) +
  scale_x_log10nice(name = 'Surface EAR mdl', limits=c(0.000000005, .02), expand=c(0,0)) +
  scale_y_log10nice(name = 'Passive EAR mdl', limits=c(0.000000005, .02), expand=c(0,0)) 

print(mdlcompare)

ggsave(file_out(file.path(path_to_data, "Figures/MDLComparison.png")), plot = mdlcompare, height=6, width=6)


mdl_EAR_join2 <- mdl_EAR_join %>%
  filter(chnm %in% unique(chemicalSummary_merged$chnm))

mdlcompare2 <- ggplot(mdl_EAR_join2, aes(x=EAR_surface, y=EAR_passive, label=chnm)) +
  geom_point(col='red', size=2)+
  geom_abline(linetype='dashed') +
  geom_vline(xintercept = .001) + 
  geom_hline(yintercept = .001) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # xlim(range(c(join_mdl$Surface_mdl, join_mdl$Passive_mdl))) +
  # ylim(range(c(join_mdl$Surface_mdl, join_mdl$Passive_mdl))) +
  # xlim(c(0,.02)) +
  # ylim(c(0,.02)) +
  labs(x=expression(paste('Surface mdl (', mu, 'g L'^'-1', ')')), 
       y=expression(paste('Passive mdl (', mu, 'g L'^'-1', ')'))) +
  geom_text_repel(size=3, color='black', segment.color='black', alpha=.5, segment.size=.5) +
  scale_x_log10nice(name = 'Surface EAR mdl', limits=c(0.000000005, .006), expand=c(0,0)) +
  scale_y_log10nice(name = 'Passive EAR mdl', limits=c(0.000000005, .006), expand=c(0,0)) 

print(mdlcompare2)

ggsave(file_out(file.path(path_to_data, "Figures/MDLComparison2.png")), plot = mdlcompare2, height=6, width=6)






surface_mdl <- read_excel(file_in(file.path(path_to_data, 'Data/pesticides_dls.xlsx')), sheet='Data') %>%
  select(CAS, Value) %>%
  rename(Surface_mdl = Value)


passive_mdl <- read_excel(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016MDL_ToxEval.xlsx")), sheet='Data') %>%
  select(CAS, Value, chnm) %>%
  rename(Passive_mdl = Value) 



join_mdl <- full_join(surface_mdl, passive_mdl) %>%
  drop_na(Surface_mdl, Passive_mdl) %>%
  mutate(Diff = Surface_mdl - Passive_mdl,
         PerDiff = (Surface_mdl - Passive_mdl)/(Surface_mdl + Passive_mdl)*200,
         PerGreater = Surface_mdl/Passive_mdl) %>%
  arrange(desc(PerDiff))

ggplot(join_mdl, aes(x=Surface_mdl, y=Passive_mdl, label=chnm)) +
  geom_point()+
  geom_abline() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(range(c(join_mdl$Surface_mdl, join_mdl$Passive_mdl))) +
  ylim(range(c(join_mdl$Surface_mdl, join_mdl$Passive_mdl))) +
  # xlim(c(0,.02)) +
  # ylim(c(0,.02)) +
  labs(x=expression(paste('Surface mdl (', mu, 'g L'^'-1', ')')), 
       y=expression(paste('Passive mdl (', mu, 'g L'^'-1', ')'))) +
  geom_text() 
  # scale_x_log10nice() +
  # scale_y_log10nice()




