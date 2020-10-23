

#### Setup ####
library(toxEval)
library(ToxMixtures)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

#Need to run (ToxEval_Passive2016.R and ToxEval_WaterSamples2016.R)

chemicalSummary2_surf <- chemicalSummary_surface %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  filter(`date`>=date_filter[1], `date`<=date_filter[2],
         shortName %in% site_order) %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  group_by(site, endPoint, shortName, Lake, site_grouping, CAS, chnm, Bio_category, Class) %>%
  summarise(EAR = max(EAR, na.rm=T), .groups = 'drop') 
  

chemicalSummary2 <- chemicalSummary %>%
  mutate(shortName = factor(shortName, site_order))


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

priorityEndpoints$geneSymbol[is.na(priorityEndpoints$geneSymbol)] <- paste0(priorityEndpoints$Bio_category[is.na(priorityEndpoints$geneSymbol)], "*")


priorityGenes <- priorityEndpoints %>%
  group_by(chnm, geneSymbol, site) %>%
  summarize(EAR = max(EAR, na.rm=T))

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

endpoint_bychemgene_boxplot <- ggplot(priorityGenes, aes(y=geneSymbol, x=EAR)) +
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
            endPoint = paste(unique(endPoint), collapse = ","),
            .groups = 'drop') 


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
  # filter(chnm %in% c('2,4-Dichlorophenoxyacetic acid'))
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
            endPoint = paste(unique(endPoint), collapse = ","), .groups = 'drop')


write.csv(genes_out_surf, file = file_out(file.path(path_to_data, 'Data', 'priority_genes_surface.csv')), row.names = F)




# #########
# Mixtures
# #########
# group_by_this <- "AOP"
group_by_this <- "endPoint"
ear_threshold <- .001
site_threshold <- 1.9

#Passive
all_combos <- top_mixes(chemicalSummary2, group_by_this, 
                        ear_threshold,
                        1) 

overall_max_n_all <- overall_mixtures(all_combos, "max")

endpoint_review <- overall_max_n_all %>%
  dplyr::select(endPoint, CASs) %>%
  tidyr::separate(CASs, sep='\\|', into = letters[1:10], fill='right') %>%
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


gene_table <- gene_mixtures(overall_max_n_chem)  


full_genesummary <- gene_functions(overall_max_n_chem, species = c("Homo sapiens", 
                                                                 "Danio rerio", 
                                                                 "Xenopus tropicalis"))

entrez_genesummary <- gene_functions(overall_max_n_chem, 
                                   columns = "entrez", 
                                   species = "Homo sapiens") 


int <- full_genesummary %>%
  select(geneSymbol, Species, gene_name, pathway_name, ENTREZ_GENE_SUMMARY, 
         GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, 
         KEGG_PATHWAY, OMIM_DISEASE) %>%
  distinct()

write.csv(int, file = file_out(file.path(path_to_data, 'Data', 'passive_prioritygene_interesting.csv')), row.names = F)


int2 <- int %>%
  group_by(geneSymbol) %>%
  summarize(orthologs = paste(unique(Species[!is.na(Species)]), 
                              collapse= "|"),
            pathways = paste(unique(pathway_name[!is.na(pathway_name)]), 
                             collapse= ","),
            KEGG_PATHWAY = paste(unique(KEGG_PATHWAY[!is.na(KEGG_PATHWAY)]),
                                 collapse= ","),
            OMIM_DISEASE = paste(unique(OMIM_DISEASE[!is.na(OMIM_DISEASE)]), 
                                 collapse= ","),
            ENTREZ_GENE_SUMMARY = paste(unique(ENTREZ_GENE_SUMMARY[!is.na(ENTREZ_GENE_SUMMARY)]),
                                        collapse= "|"),
            GOTERM_BP_DIRECT = paste(unique(GOTERM_BP_DIRECT[!is.na(GOTERM_BP_DIRECT)]), 
                                     collapse= ","),
            GOTERM_CC_DIRECT = paste(unique(GOTERM_CC_DIRECT[!is.na(GOTERM_CC_DIRECT)]), 
                                     collapse= ","),
            GOTERM_MF_DIRECT = paste(unique(GOTERM_MF_DIRECT[!is.na(GOTERM_MF_DIRECT)]), 
                                     collapse= ","),
            .groups = "drop") %>%
  mutate(pathways = unlist(lapply(strsplit(pathways, ","), 
                function(l) paste(l, collapse = "|"))),
         KEGG_PATHWAY = unlist(lapply(strsplit(KEGG_PATHWAY, ","), 
                                      function(l) paste(l, collapse = "|"))),
         OMIM_DISEASE = unlist(lapply(strsplit(OMIM_DISEASE, ","), 
                                      function(l) paste(l, collapse = "|"))),
         GOTERM_BP_DIRECT = unlist(lapply(strsplit(GOTERM_BP_DIRECT, ","), 
                                      function(l) paste(l, collapse = "|"))),
         GOTERM_CC_DIRECT = unlist(lapply(strsplit(GOTERM_CC_DIRECT, ","), 
                                      function(l) paste(l, collapse = "|"))),
         GOTERM_MF_DIRECT = unlist(lapply(strsplit(GOTERM_MF_DIRECT, ","), 
                                      function(l) paste(l, collapse = "|")))) %>%
  left_join(unique(select(overall_max_n_chem, genes, AOP_ids, AOP_names)), by = c("geneSymbol" = "genes"))


data.frame(int2)

all_kegg <- paste(unique(int2$KEGG_PATHWAY), collapse= "|")
all_kegg <- unique(unlist(strsplit(all_kegg, "\\|")))
all_kegg <- unique(unlist(strsplit(all_kegg, ":")))
length(all_kegg[grepl('hsa', all_kegg)])
length(all_kegg[grepl('dre', all_kegg)])
length(all_kegg[grepl('xtr', all_kegg)])

all_kegg <- all_kegg[!grepl('hsa', all_kegg)]
all_kegg <- all_kegg[!grepl('dre', all_kegg)]
all_kegg <- all_kegg[!grepl('xtr', all_kegg)]

int2 <- full_join(gene_table, int2, by = c("geneSymbol"))

int2[is.na(int2)] <- ""


write.csv(int2, file = file_out(file.path(path_to_data, 'Data', 'passive_prioritygene_interesting2.csv')), row.names = F)



entrez_genesummary[,1:2]

sample_EARmix <- EAR_mixtures(chemicalSummary2, group_by_this)

AOP_table <- left_join(sample_EARmix, AOP_crosswalk, by = "endPoint") %>%
  filter(!is.na(ID)) %>%
  group_by(ID, AOP_title, site) %>%
  summarize(AOPsum = sum(EARsum), 
            genes = paste(unique(geneSymbol[!is.na(geneSymbol)]), collapse = "|"),
            endPoints = paste(unique(endPoint[!is.na(endPoint)]), collapse = "|")) %>%
  filter(AOPsum > ear_threshold) %>%
  group_by(ID, AOP_title, genes, endPoints) %>%
  summarize(AOP_sum_mean = round(mean(AOPsum, na.rm=T), 5),
            AOP_sum_max = round(max(AOPsum, na.rm=T), 5), 
            sites =  paste(unique(site[!is.na(site)]), collapse = "|")) %>%
  select(ID, AOP_sum_max, AOP_sum_mean, everything())
  
data.frame(AOP_table)

write.csv(AOP_table, file = file.path(path_to_data, "Data", "PassiveAOP_results.csv"), row.names = FALSE)

site_EARmix <- site_mixtures(sample_EARmix, ear_cutoff = 0.001)


#plotting

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



summed_EARs_maxGene <- sample_EARmix %>%
  group_by(site, geneSymbol) %>%
  summarize(max_ear_gene = max(EARsum, na.rm = TRUE),
            .groups = "drop") 

# summed_EARs_maxGene <- summed_EARs %>%
#   group_by(site, geneSymbol) %>%
#   summarize(sum_ear_endpoint = max(sum_ear_endpoint))

priority_genes <- summed_EARs_maxGene %>%
  group_by(geneSymbol) %>%
  mutate(n_exceed = length(which(max_ear_gene > ear_threshold))) %>%
  filter(n_exceed > site_threshold) %>%
  select(-n_exceed) 

gene_rank <- priority_genes %>%
  group_by(geneSymbol) %>%
  summarize(sum_ear_median = median(max_ear_gene, na.rm=T)) %>%
  mutate(order = ifelse(grepl("\\*", geneSymbol), 2, 1)) %>%
  arrange(desc(order), sum_ear_median)

priority_genes <- priority_genes %>%
  left_join(select(gene_rank, geneSymbol, order)) %>%
  mutate(geneSymbol = factor(geneSymbol, gene_rank$geneSymbol))


top_mixgene_box <- ggplot(priority_genes, aes(y=geneSymbol, x=max_ear_gene)) +
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
chemicalSummary2_summer <- chemicalSummary_surface %>%  
  mutate(shortName = factor(shortName, site_order)) %>%
  filter(`date`>=date_filter[1], `date`<=date_filter[2],
         shortName %in% site_order) %>%
  mutate(shortName = factor(shortName, site_order))
  
all_combos_surf <- top_mixes(chemicalSummary2_summer, group_by_this, 
                        ear_threshold,
                        .9) 

priority_combos_surf <- filter(all_combos_surf, n_samples>10)
priority_max <- overall_mixtures(priority_combos_surf, "max")

class_key <- class_key_fnx(chemicalSummary2_summer)

priority_max_fancy <- overall_df_format(priority_max, class_key)

write.csv(priority_max, file = file_out(file.path(path_to_data, 'Data', 'priority_mixtures_surface.csv')), row.names = F)

overall_max_n_all_surf <- overall_mixtures(all_combos_surf, "max")

sample_EARmix_surf <- EAR_mixtures(chemicalSummary2_summer, group_by_this)



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


top_combos_surf <- top_mixes(chemicalSummary2_summer,
                             group_by_this, 
                             ear_threshold,
                             site_threshold) 


overall_max_n_chem_surf <- overall_mixtures(top_combos_surf, "max")

gene_table_surf <- gene_mixtures(overall_max_n_chem_surf)  


summed_EARs_surf <- chemicalSummary_surface %>%
  filter(`date`>=date_filter[1], `date`<=date_filter[2],
         shortName %in% site_order) %>%
  group_by(site, date, !!sym(group_by_this)) %>%
  summarize(sum_ear_endpoint = sum(EAR)) %>%
  filter(!!sym(group_by_this) %in% top_combos_surf$endPoint) %>%
  ungroup() %>%
  distinct() %>%
  # left_join(unique(select(priorityEndpoints_surf, endPoint))) %>%
  group_by(site, endPoint) %>%
  summarize(sum_ear_endpoint = max(sum_ear_endpoint, na.rm=T ), 
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

ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_surf.png")), top_mix_box_surf, height=4, width=6, units='in')

# summed_EARs_surf$geneSymbol[is.na(summed_EARs_surf$geneSymbol)] <- paste0(summed_EARs_surf$Bio_category[is.na(summed_EARs_surf$geneSymbol)], "*")

sample_EARmix_surf$geneSymbol[is.na(sample_EARmix_surf$geneSymbol)] <- 
  paste0(sample_EARmix_surf$endPoint[is.na(sample_EARmix_surf$geneSymbol)], "*")
  
sample_EARmix_surf$geneSymbol[grepl("Tanguay", sample_EARmix_surf$endPoint)] <- 
  paste0("Zebrafish*")


summed_EARs_maxGene_surf <- sample_EARmix_surf %>%
  group_by(site, geneSymbol) %>%
  summarize(max_ear_gene = max(EARsum, na.rm = TRUE),
            .groups = "drop") 

priority_genes_surf <- summed_EARs_maxGene_surf %>%
  group_by(geneSymbol) %>%
  mutate(n_exceed = length(which(max_ear_gene > ear_threshold))) %>%
  filter(n_exceed > site_threshold) %>%
  select(-n_exceed) 

gene_rank_surf <- priority_genes_surf %>%
  group_by(geneSymbol) %>%
  summarize(sum_ear_median = median(max_ear_gene, na.rm=T)) %>%
  mutate(order = ifelse(grepl("\\*", geneSymbol), 2, 1)) %>%
  arrange(desc(order), sum_ear_median)

priority_genes_surf <- priority_genes_surf %>%
  left_join(select(gene_rank_surf, geneSymbol, order)) %>%
  mutate(geneSymbol = factor(geneSymbol, gene_rank_surf$geneSymbol))


# 
# 
# summed_EARs_maxGene_surf <- summed_EARs_surf %>%
#   group_by(site, geneSymbol) %>%
#   summarize(sum_ear_endpoint = max(sum_ear_endpoint))
# 
# gene_rank_surf <- summed_EARs_maxGene_surf %>%
#   group_by(geneSymbol) %>%
#   summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
#   mutate(order = ifelse(grepl("\\*", geneSymbol), 2, 1)) %>%
#   arrange(desc(order), sum_ear_median)
# 
# summed_EARs_maxGene_surf <- summed_EARs_maxGene_surf %>%
#   left_join(select(gene_rank_surf, geneSymbol, order)) %>%
#   mutate(geneSymbol = factor(geneSymbol, gene_rank_surf$geneSymbol))
# 


top_mixgene_box_surf <- ggplot(priority_genes_surf, aes(y=geneSymbol, x=max_ear_gene)) +
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

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_surf.png")), 
       top_mixgene_box_surf, height=4, width=4, units='in')



gene_combine <- priority_genes_surf %>%
  mutate(method = 'Water') %>%
  bind_rows(mutate(priority_genes, method = 'Passive'))

gene_rank_combine <- gene_combine %>%
  group_by(geneSymbol, order) %>%
  summarize(max_ear_gene  = median(max_ear_gene , na.rm=T)) %>%
  arrange(desc(order), max_ear_gene)

gene_combine <- gene_combine %>%
  mutate(geneSymbol = factor(geneSymbol, gene_rank_combine$geneSymbol),
         method = factor(method, c('Passive', 'Water')))



top_mixgene_box_twomethod <- ggplot(gene_combine, aes(y=geneSymbol, x=max_ear_gene)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  scale_fill_manual(values = c('red', 'blue')) + 
  scale_color_manual(values = c('red', 'blue')) + 
  # geom_point(position = position_dodge(width=0.75), aes(group=method, color=method), shape=16,  alpha=0.5, size=1.5) +
  geom_jitter(shape=16, aes(color=method), alpha=0.5, width=0, height=.2, size=1) +
  geom_boxplot(position = position_dodge(preserve = "single"), 
               aes(fill=method), outlier.shape=NA, alpha=.2, width=0.6, color='grey20') +
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





endpoint_combine <- summed_EARs_surf %>%
  mutate(method = 'Water') %>%
  bind_rows(mutate(summed_EARs, method = 'Passive')) %>%
  select(site, endPoint, sum_ear_endpoint, n, method) %>%
  distinct()

endpoint_rank_combine <- endpoint_combine %>%
  group_by(endPoint) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  arrange(sum_ear_median)

endpoint_combine <- endpoint_combine %>%
  mutate(endPoint = factor(endPoint, endpoint_rank_combine$endPoint),
         method = factor(method, c('Passive', 'Water')))



top_mix_endpoint_box_twomethod <- ggplot(endpoint_combine,
                                         aes(y=endPoint, x=sum_ear_endpoint)) +
  scale_fill_manual(values = c('red', 'blue')) + 
  scale_color_manual(values = c('red', 'blue')) + 
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(aes(color=method), shape=16, alpha=0.5, width=0, height=.2, size=1) +
  geom_boxplot(aes(fill=method), outlier.shape=NA, alpha=.2, width=.6, col='grey20') +
  # geom_jitter(shape=16, aes(color=Bio_category), alpha=0.5, width=0, height=.1, size=1.5) + 
  # geom_boxplot(outlier.shape=NA, aes(fill=Bio_category), alpha=.2, width=.5) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='ToxCast assay name') + 
  theme_tox() +
  theme(legend.position='none', legend.title.align=0.5) +
  guides(col = guide_legend(ncol = 1)) +
  facet_wrap(~method) +
  theme(strip.background = element_rect(fill=NA, color=NA))

print(top_mix_endpoint_box_twomethod)


ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_twoMethods.png")), 
       top_mix_endpoint_box_twomethod, height=5, width=6, units='in')






#Get biological relevance for geneSymbols

head(panther)
head(join_criteria())
head(david_full)

#All endpoints with EARs > 10-3
priorityEndpoints 
priorityEndpoints_surf

#List of endpoints, genes, and pathways
overall_max_n_chem
overall_max_n_all_surf



full_panther_david <- separate_rows(overall_max_n_chem, genes) %>%
  left_join(panther, by=c('genes' = 'gene_abbr')) %>%
  left_join(david_full, by=c('genes' = 'GeneSymbol', 'species' = 'Species')) 

select_panther_david <- full_panther_david %>% 
  select( endPoint, chems, CASs, genes, species, gene_name,
          ENTREZ_GENE_SUMMARY, OMIM_DISEASE, KEGG_PATHWAY, pathway)



data.frame(select_panther_david) 

full_panther_david[,1:10]

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
