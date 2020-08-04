



# #########
# Mixture analysis for surface pesticide samples 12 months
# #########

library(ToxMixtures)

group_by_this <- "endPoint"
ear_threshold <- .001
site_threshold <- 1.9 #Number of sites minimum (10%)


all_combos_surf <- top_mixes(chemicalSummary_surface, 
                             group_by_this, 
                             ear_threshold,
                             0.9) 

priority_combos_surf <- filter(all_combos_surf, n_samples>10)
priority_max <- overall_mixtures(priority_combos_surf, "max")

class_key <- class_key_fnx(chemicalSummary2_summer)

priority_max_fancy <- overall_df_format(priority_max, class_key)

write.csv(priority_max, file = file_out(file.path(path_to_data, 'Data', 'priority_mixtures_surface_12months.csv')), row.names = F)

overall_max_n_all_surf <- overall_mixtures(all_combos_surf, "max")


endpoint_review_surf <- overall_max_n_all_surf %>%
  dplyr::select(endPoint, CASs) %>%
  tidyr::separate(CASs, sep='\\|', into = letters[1:10], fill='right') %>%
  # select_if(~sum(!is.na(.)) > 0) %>%
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


top_combos_surf <- top_mixes(chemicalSummary_surface,
                             group_by_this, 
                             ear_threshold,
                             site_threshold) 


overall_max_n_chem_surf <- overall_mixtures(top_combos_surf, "max")

gene_summary_all <- gene_summary(overall_max_n_chem_surf)

write.csv(gene_summary_all, file = file_out(file.path(path_to_data, 'Data', 'gene_summary_mixtures_surface_12months.csv')), row.names = F)

gene_summary_entrez <- gene_summary(overall_max_n_chem_surf, columns = 'entrez') %>%
  filter(nchar(ENTREZ_GENE_SUMMARY)>0)

write.csv(gene_summary_entrez, file = file_out(file.path(path_to_data, 'Data', 'gene_entrez_mixtures_surface_12months.csv')), row.names = F)



summed_EARs_surf <- chemicalSummary_surface %>%
  group_by(site, date, !!sym(group_by_this)) %>%
  summarize(sum_ear_endpoint = sum(EAR, na.rm = TRUE)) %>%
  filter(!!sym(group_by_this) %in% top_combos_surf$endPoint) %>%
  ungroup() %>%
  distinct() %>%
  # left_join(unique(select(priorityEndpoints_surf, endPoint))) %>%
  group_by(site, endPoint) %>%
  summarize(sum_ear_endpoint = max(sum_ear_endpoint, na.rm = TRUE), 
            n=n()) %>%
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

ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_surf_12months.png")), top_mix_box_surf, height=5, width=6, units='in')

summed_EARs_surf$geneSymbol[is.na(summed_EARs_surf$geneSymbol)] <- paste0(summed_EARs_surf$Bio_category[is.na(summed_EARs_surf$geneSymbol)], "*")

gene_rank_surf <- summed_EARs_surf %>%
  group_by(geneSymbol) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  mutate(order = ifelse(grepl("\\*", geneSymbol), 2, 1)) %>%
  arrange(desc(order), sum_ear_median)

summed_EARs_surf <- summed_EARs_surf %>%
  left_join(select(gene_rank_surf, geneSymbol, order)) %>%
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

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_surf_12months.png")), 
       top_mixgene_box_surf, height=5, width=4, units='in')





