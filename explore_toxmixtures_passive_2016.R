

#### Setup ####
library(toxEval)
library(ToxMixtures)
library(dplyr)
library(RColorBrewer)
library(ggplot2)


chemicalSummary2_surf <- chemicalSummary_surface %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  filter(`date`>=date_filter[1], `date`<=date_filter[2],
         shortName %in% site_order) %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  group_by(site, endPoint, shortName, Lake, site_grouping, CAS, chnm, Bio_category, Class) %>%
  summarise(EAR = max(EAR, na.rm=T)) %>%
  ungroup()
  

chemicalSummary2 <- chemicalSummary %>%
  mutate(shortName = factor(shortName, site_order))

plot_stackbar(chemicalSummary2, x='shortName', x_label='Site', y_label="Number of chemicals",
              fill='EAR', stack='chnm', 
              breaks=c(.0001, .001, .01)) + 
  theme(axis.text.x = element_text(angle=90, hjust=1))


plot_stackbar(chemicalSummary2, x='chnm', x_label='Chemical', y_label="Number of rivers",
              fill='EAR', stack='site', 
              breaks=c(.0001, .001, .01)) + 
  theme(axis.text.x = element_text(angle=90, hjust=1))


plot_stackbar(chemicalSummary_surface, x='shortName', x_label='Site', y_label="Number of chemicals",
              fill='EAR', stack='chnm', 
              breaks=c(.0001, .001, .01)) + 
  theme(axis.text.x = element_text(angle=90, hjust=1))


plot_stackbar(chemicalSummary_surface, x='chnm', x_label='Chemical', y_label="Number of rivers",
              fill='EAR', stack='site', 
              breaks=c(.0001, .001, .01)) + 
  theme(axis.text.x = element_text(angle=90, hjust=1))


# cs <- chemical_summary
# group_by_this <- "endPoint"
# ear_threshold <- 0.001
# site_threshold_percent <- 10
# n_sites <- length(unique(chemical_summary$site))
# site_threshold <- ceiling(n_sites / site_threshold_percent)


#Priority endpoints for passive data

priorityEndpoints <- filter(chemicalSummary2, EAR>0.001) %>%
  dplyr::select(chnm, endPoint) %>%
  unique() %>%
  left_join(chemicalSummary2) %>%
  arrange(desc(Bio_category), EAR)

name_order <- unique(priorityEndpoints$endPoint)

priorityEndpoints <- priorityEndpoints %>%
  mutate( endPoint = factor(endPoint,name_order),
          Bio_category = factor(Bio_category)) %>%
  left_join(unique(select(join_criteria(), geneSymbol,endPoint)), by = "endPoint")
  # left_join(unique(join_criteria()), by = "endPoint") 


endpoint_bychem_boxplot <- ggplot(unique(priorityEndpoints[,c('endPoint', 'EAR', 'site', 'chnm', 'Bio_category')]), aes(y=endPoint, x=EAR)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, aes(color=Bio_category), alpha=0.5, width=0, height=.1, size=1.5) + 
  geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = 'EAR') +
  labs(y='ToxCast assay name', color = 'Bio Category', fill = 'Bio Category') + 
  # facet_grid(Bio_category~chnm, scales = "free") + 
  facet_grid(~chnm) + 
  theme_tox() +
  theme(legend.position='bottom', legend.title.align=0.5,
        axis.text = element_text(size=8), strip.text =  element_text(size=8), 
        strip.background=element_rect(fill='white', color=NA)) +
  guides(col = guide_legend(nrow = 2))

print(endpoint_bychem_boxplot)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_byChemical.png")), endpoint_bychem_boxplot, height=6, width=8.5, units='in')

endpoint_bychemgene_boxplot <- ggplot(priorityEndpoints, aes(y=geneSymbol, x=EAR)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='red', alpha=0.5, width=0, height=.2, size=1.5) + 
  geom_boxplot(outlier.shape=NA, fill='red', alpha=.2, width=.5) + 
  scale_x_log10nice(name = 'EAR') +
  labs(y='Gene Target') + 
  # facet_grid(Bio_category~chnm, scales = "free") + 
  facet_grid(~chnm) + 
  theme_tox() +
  # theme(legend.position='bottom', legend.title.align=0.5) +
  theme(axis.text = element_text(size=8), strip.text =  element_text(size=8), 
        strip.background=element_rect(fill='white', color=NA)) +
  guides(col = guide_legend(nrow = 2)) +
  scale_y_discrete(limits = unique(priorityEndpoints$geneSymbol)[order(unique(priorityEndpoints$geneSymbol), decreasing = T)])

print(endpoint_bychemgene_boxplot)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_byChemical.png")), endpoint_bychemgene_boxplot, height=6, width=8.5, units='in')

genes_out <- priorityEndpoints %>%
  select(chnm, endPoint, Bio_category, geneSymbol) %>%
  distinct() %>%
  group_by(geneSymbol, Bio_category) %>% 
  summarize(chnm = paste(unique(chnm), collapse = ","),
            endPoint = paste(unique(endPoint), collapse = ",")) %>%
  ungroup()  


write.csv(genes_out, file = file_out(file.path(path_to_data, 'Data', 'priority_genes.csv')), row.names = F)






#Priority endpoints for water data

priorityEndpoints_surf <- filter(chemicalSummary2_surf, EAR>0.001) %>%
  dplyr::select(chnm, endPoint) %>%
  unique() %>%
  left_join(chemicalSummary2_surf) %>%
  arrange(desc(Bio_category), EAR)

name_order_surf <- unique(priorityEndpoints_surf$endPoint)

priorityEndpoints_surf <- priorityEndpoints_surf %>%
  mutate( endPoint = factor(endPoint,name_order_surf),
          Bio_category = factor(Bio_category)) %>%
  left_join(unique(select(join_criteria(), geneSymbol,endPoint)), by = "endPoint")
  # filter(chnm %in% c('Atrazine'))
# left_join(unique(join_criteria()), by = "endPoint") 


endpoint_bychem_boxplot_surf <- ggplot(unique(priorityEndpoints_surf[,c('endPoint', 'EAR', 'site', 'chnm', 'Bio_category')]), aes(y=endPoint, x=EAR)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, aes(color=Bio_category), alpha=0.5, width=0, height=.1, size=1.5) + 
  geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = 'EAR') +
  labs(y='ToxCast assay name', color = 'Bio Category', fill = 'Bio Category') + 
  # facet_grid(Bio_category~chnm, scales = "free") + 
  facet_grid(~chnm) + 
  theme_tox() +
  theme(legend.position='bottom', legend.title.align=0.5,
        axis.text = element_text(size=8), strip.text =  element_text(size=8), 
        strip.background=element_rect(fill='white', color=NA)) +
  guides(col = guide_legend(nrow = 2))

print(endpoint_bychem_boxplot_surf)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_byChemical_surface.png")), endpoint_bychem_boxplot_surf, height=6, width=8.5, units='in')

endpoint_bychemgene_boxplot_surf <- ggplot(priorityEndpoints_surf, aes(y=geneSymbol, x=EAR)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='blue', alpha=0.5, width=0, height=.2, size=1.5) + 
  geom_boxplot(outlier.shape=NA, fill='blue', alpha=.2, width=.5) + 
  scale_x_log10nice(name = 'EAR') +
  labs(y='Gene Target') + 
  # facet_grid(Bio_category~chnm, scales = "free") + 
  facet_grid(~chnm) + 
  theme_tox() +
  # theme(legend.position='bottom', legend.title.align=0.5) +
  theme(axis.text = element_text(size=8), strip.text =  element_text(size=8), 
        strip.background=element_rect(fill='white', color=NA)) +
  guides(col = guide_legend(nrow = 2)) +
  scale_y_discrete(limits = unique(priorityEndpoints_surf$geneSymbol)[order(unique(priorityEndpoints_surf$geneSymbol), decreasing = T)])

print(endpoint_bychemgene_boxplot_surf)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_byChemical_surf.png")), endpoint_bychemgene_boxplot_surf, height=6, width=8.5, units='in')

genes_out_surf <- priorityEndpoints_surf %>%
  select(chnm, endPoint, Bio_category, geneSymbol) %>%
  distinct() %>%
  group_by(geneSymbol, Bio_category) %>% 
  summarize(chnm = paste(unique(chnm), collapse = ","),
            endPoint = paste(unique(endPoint), collapse = ",")) %>%
  ungroup()  


write.csv(genes_out_surf, file = file_out(file.path(path_to_data, 'Data', 'priority_genes_surface.csv')), row.names = F)





#Mixtures


group_by_this <- "endPoint"
ear_threshold <- .001
site_threshold <- 2

#Passive
all_combos <- top_mixes(chemicalSummary2, group_by_this, 
                        ear_threshold,
                        1) 

overall_max_n_all <- overall_mixtures(all_combos, "max")


endpoint_review <- overall_max_n_all %>%
  dplyr::select(endPoint, CASs) %>%
  tidyr::separate(CASs, sep='\\|', into = letters[1:10], fill='right') %>%
  select_if(~sum(!is.na(.)) > 0) %>%
  gather(key = chem_mix_nu, value=CAS, -1) %>%
  drop_na(CAS) %>%
  dplyr::select(-chem_mix_nu)
  # left_join(dplyr::select(ACClong, endPoint, ACC, ACC_value, chnm))
  # left_join(unique(dplyr::select(ACClong, endPoint, chnm)))
  
ACCs <- dplyr::select(ACClong, endPoint, chnm, CAS, ACC, ACC_value) %>%
  filter(endPoint %in% endpoint_review$endPoint, 
         CAS %in% endpoint_review$CAS)

endpoint_out <- endpoint_review %>%
  left_join(ACCs) %>%
  arrange(endPoint, ACC) %>%
  distinct()

  data.frame(endpoint_out)
  

top_combos <- top_mixes(chemicalSummary2, group_by_this, 
                        ear_threshold,
                        site_threshold) 


overall_max_n_chem <- overall_mixtures(top_combos, "max")


summed_EARs <- chemicalSummary2 %>%
  group_by(site, date, !!sym(group_by_this)) %>%
  summarize(sum_ear_endpoint = sum(EAR)) %>%
  filter(!!sym(group_by_this) %in% top_combos$endPoint) %>%
  ungroup() %>%
  distinct() %>%
  left_join(unique(select(chemicalSummary2, endPoint, Bio_category)))

endpoint_rank <- summed_EARs %>%
  group_by(endPoint) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  arrange((sum_ear_median))

summed_EARs <- summed_EARs %>%
  left_join(overall_max_n_chem) %>%
  left_join(unique(select(join_criteria(), geneSymbol,endPoint)), by = "endPoint") %>%
  mutate(endPoint = factor(endPoint, endpoint_rank$endPoint))

top_mix_box <- ggplot(unique(summed_EARs[,c('endPoint', 'sum_ear_endpoint', 'site', 'Bio_category')]), aes(y=endPoint, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='red', alpha=0.5, width=0, height=.1, size=1.5) +
  geom_boxplot(outlier.shape=NA, fill='red', alpha=.2, width=.5) +
  # geom_jitter(shape=16, aes(color=Bio_category), alpha=0.5, width=0, height=.1, size=1.5) + 
  # geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='ToxCast assay name', fill='Bio Category', color='Bio Category') + 
  theme_tox() +
  theme(legend.position='right', legend.title.align=0.5) +
  guides(col = guide_legend(ncol = 1))

print(top_mix_box)
  
ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot.png")), top_mix_box, height=4, width=6, units='in')


gene_rank <- summed_EARs %>%
  group_by(geneSymbol) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  arrange((sum_ear_median))

summed_EARs <- summed_EARs %>%
  mutate(geneSymbol = factor(geneSymbol, gene_rank$geneSymbol))

top_mixgene_box <- ggplot(summed_EARs, aes(y=geneSymbol, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='red', alpha=0.5, width=0, height=.1, size=1.5) +
  geom_boxplot(outlier.shape=NA, fill='red', alpha=.2, width=.5) +
  # geom_jitter(shape=16, aes(color=Bio_category), alpha=0.5, width=0, height=.1, size=1.5) + 
  # geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='Gene Target', fill='Bio Category', color='Bio Category') + 
  theme_tox() +
  theme(legend.position='right', legend.title.align=0.5) +
  guides(col = guide_legend(ncol = 1))

print(top_mixgene_box)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot.png")), 
       top_mixgene_box, height=4, width=4, units='in')




#Surface samples
all_combos_surf <- top_mixes(chemicalSummary_surface, group_by_this, 
                        ear_threshold,
                        1) 

overall_max_n_all_surf <- overall_mixtures(all_combos_surf, "max")


endpoint_review_surf <- overall_max_n_all_surf %>%
  dplyr::select(endPoint, CASs) %>%
  tidyr::separate(CASs, sep='\\|', into = letters[1:10], fill='right') %>%
  select_if(~sum(!is.na(.)) > 0) %>%
  gather(key = chem_mix_nu, value=CAS, -1) %>%
  drop_na(CAS) %>%
  dplyr::select(-chem_mix_nu)
# left_join(dplyr::select(ACClong, endPoint, ACC, ACC_value, chnm))
# left_join(unique(dplyr::select(ACClong, endPoint, chnm)))

ACCs_surf <- dplyr::select(ACClong, endPoint, chnm, CAS, ACC, ACC_value) %>%
  filter(endPoint %in% endpoint_review_surf$endPoint, 
         CAS %in% endpoint_review_surf$CAS)

endpoint_out_surf <- endpoint_review_surf %>%
  left_join(ACCs_surf) %>%
  arrange(endPoint, ACC) %>%
  distinct()

data.frame(endpoint_out_surf)


top_combos_surf <- top_mixes(chemicalSummary_surface, group_by_this, 
                        ear_threshold,
                        site_threshold) 


overall_max_n_chem_surf <- overall_mixtures(top_combos_surf, "max")


summed_EARs_surf <- chemicalSummary_surface %>%
  group_by(site, date, !!sym(group_by_this)) %>%
  summarize(sum_ear_endpoint = sum(EAR)) %>%
  filter(!!sym(group_by_this) %in% top_combos_surf$endPoint) %>%
  ungroup() %>%
  distinct() %>%
  # left_join(unique(select(priorityEndpoints_surf, endPoint))) %>%
  group_by(site, endPoint) %>%
  summarize(sum_ear_endpoint = max(sum_ear_endpoint, na.rm=T )) %>%
  ungroup() %>%
  left_join(unique(select(chemicalSummary_surface, endPoint, Bio_category)))

endpoint_rank_surf <- summed_EARs_surf %>%
  group_by(endPoint) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  arrange((sum_ear_median))

summed_EARs_surf <- summed_EARs_surf %>%
  left_join(overall_max_n_chem_surf) %>%
  left_join(unique(select(join_criteria(), geneSymbol,endPoint)), by = "endPoint") %>%
  mutate(endPoint = factor(endPoint, endpoint_rank_surf$endPoint)) %>%
  select(-chems, -CASs, -n_chems) %>%
  distinct()

top_mix_box_surf <- ggplot(unique(summed_EARs_surf[c('endPoint', 'sum_ear_endpoint', 
                                                     'Bio_category')]), aes(y=endPoint, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='blue', alpha=0.5, width=0, height=.1, size=1.5) +
  geom_boxplot(outlier.shape=NA, fill='blue', alpha=.2, width=.5) +
  # geom_jitter(shape=16, aes(color=Bio_category), alpha=0.5, width=0, height=.1, size=1.5) + 
  # geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='ToxCast assay name', fill='Bio Category', color='Bio Category') + 
  theme_tox() +
  theme(legend.position='right', legend.title.align=0.5) +
  guides(col = guide_legend(ncol = 1))

print(top_mix_box_surf)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_surf.png")), top_mix_box, height=4, width=6, units='in')


gene_rank_surf <- summed_EARs_surf %>%
  group_by(geneSymbol) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  arrange((sum_ear_median))

summed_EARs_surf <- summed_EARs_surf %>%
  mutate(geneSymbol = factor(geneSymbol, gene_rank_surf$geneSymbol))

top_mixgene_box_surf <- ggplot(summed_EARs_surf, aes(y=geneSymbol, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='blue', alpha=0.5, width=0, height=.1, size=1.5) +
  geom_boxplot(outlier.shape=NA, fill='blue', alpha=.2, width=.5) +
  # geom_jitter(shape=16, aes(color=Bio_category), alpha=0.5, width=0, height=.1, size=1.5) + 
  # geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='Gene Target', fill='Bio Category', color='Bio Category') + 
  theme_tox() +
  theme(legend.position='right', legend.title.align=0.5) +
  guides(col = guide_legend(ncol = 1))

print(top_mixgene_box_surf)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot.png")), 
       top_mixgene_box, height=4, width=4, units='in')



gene_combine <- summed_EARs_surf %>%
  mutate(method = 'Water') %>%
  full_join(select(summed_EARs, -date))

gene_combine$method[is.na(gene_combine$method)] <- 'Passive'
gene_combine$geneSymbol[is.na(gene_combine$geneSymbol)] <- paste0(gene_combine$Bio_category[is.na(gene_combine$geneSymbol)], "*")


gene_rank_combine <- gene_combine %>%
  group_by(geneSymbol) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  arrange((sum_ear_median))

gene_combine <- gene_combine %>%
  mutate(geneSymbol = factor(geneSymbol, gene_rank_combine$geneSymbol),
         method = factor(method, c('Passive', 'Water')))



top_mixgene_box_twomethod <- ggplot(gene_combine, aes(y=geneSymbol, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  scale_fill_manual(values = c('red', 'blue')) + 
  scale_color_manual(values = c('red', 'blue')) + 
  # geom_point(position = position_dodge(width=0.75), aes(group=method, color=method), shape=16,  alpha=0.5, size=1.5) +
  geom_jitter(shape=16, aes(color=method), alpha=0.5, width=0, height=.1, size=1) +
  geom_boxplot(position = position_dodge(preserve = "single"), aes(fill=method), outlier.shape=NA, alpha=.2, width=0.6) +
  # geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='Gene Target', fill='Method') + 
  theme_tox() +
  theme(legend.position='none', legend.title.align=0.5) +
  guides(fill = guide_legend(nrow = 1)) +
  facet_grid(~method) +
  theme(strip.background = element_rect(fill=NA, color=NA))

print(top_mixgene_box_twomethod)


ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_twoMethods.png")), 
       top_mixgene_box_twomethod, height=5, width=5, units='in')




#Pieces of code from package development





overall_max_n_chem <- overall_mixtures(top_combos, "max")  %>%
  mutate(endPoint = factor(endPoint, endpoint_rank$endPoint)) %>%
  arrange(desc(endPoint))

class_key <- class_key_fnx(chemicalSummary2)

overall_df_fancy <- overall_df_format(overall_max_n_chem,
                                      class_key)


data.frame(top_combos)
data.frame(overall_max_n_chem)
data.frame(overall_df_fancy)


knitr::kable(head(overall_df_fancy))
