


tox_list<- create_toxEval(file_in(file.path(path_to_data, "Data/PassiveForToxEval.xlsx")))
tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))

chemicalSummary <- get_chemical_summary(tox_list, 
                                        ACClong, 
                                        filtered_ep)

chemicalSummary <- tox_list$chem_site %>%
  rename(site = SiteID,
         shortName = `Short Name`) %>% 
  right_join(chemicalSummary)

chemicalSummary






chem_class_plot <- plot_tox_boxplots(chemicalSummary,
                                     category = 'Chemical Class')
chem_class_plot
ggsave(file_out(file.path(path_to_data, "Figures/PassiveBoxplot_chem_class.png")), height=5, width=5)

plot_stacks <- plot_tox_stacks(chemicalSummary, 
                               tox_list$chem_site, 
                               category = "Chemical Class")
plot_stacks
ggsave(file_out(file.path(path_to_data, "Figures/PassiveStacks_chem_class.png")), height=5, width=7)


plot_heat <- plot_tox_heatmap(chemicalSummary, 
                              tox_list$chem_site, 
                              category = "Chemical Class",
                              font_size = 7)
plot_heat

make_tox_map(chemicalSummary, 
             chem_site = tox_list$chem_site,
             category = 'Chemical Class',
             mean_logic = FALSE)



bio_plot <- plot_tox_boxplots(chemicalSummary, 
                              category = 'Biological',
                              mean_logic = FALSE,
                              hit_threshold = NA,
                              title = 'Summing EARs for chemicals within a grouping of a sample, taking the max of each site',
                              plot_ND = TRUE)
bio_plot
ggsave(file_out(file.path(path_to_data, "Figures/Passive_by_Biogroup.png")), height=5, width=7)
