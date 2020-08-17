

# ########################################################
# Mixture analysis for surface pesticide samples 12 months
# ########################################################

library(ToxMixtures)
library(ggplot2)

cs <- readRDS(file.path(path_to_data, "Rdata", "chemicalSummary_degs_estimated.rds"))

metolachlor_CAS <- filter(tox_list_surface$chem_info, parent_pesticide == "Metolachlor")$CAS

cs <- cs[-which(cs$endPoint == "TOX21_p53_BLA_p5_ratio" & cs$CAS %in% metolachlor_CAS),]

# cs <- chemicalSummary_surface
group_by_this <- "endPoint"
ear_threshold <- .01
site_threshold <- 0.9 #Number of sites minimum (e.g., 10%)

sample_n <- cs %>%
  select(date, site) %>%
  distinct() %>%
  group_by(site) %>%
  summarize(sample_n = n(), .groups = "drop")

sample_EARmix <- EAR_mixtures(cs, group_by_this, ear_cutoff = 0.001)

saveRDS(sample_EARmix, file.path(path_to_data, "Rdata", "EAR_Mix_PesticideSurfaceWater.rds"))

plot1 <- ggplot(filter(sample_EARmix, EARsum > 0.01),
                aes(x=date, y=EARsum)) +
  geom_point(size = 2, col='magenta4') +
  facet_grid(site~endPoint) +
  theme_bw() +
  scale_y_log10() +
  # theme(axis.text.x = element_blank()) +
  theme(strip.text = element_text(size=6)) +
  ggtitle("EARmix > 0.01") + 
  scale_x_datetime(date_labels ="%b", date_breaks = "3 months")

print(plot1)

ggsave(file.path(path_to_data, "Figures", "endPoint_Excedances_TS_surface.png"), plot1, height=12, width=20)


EARmix_freq <- sample_EARmix %>%
  select(site, endPoint, EARsum) %>%
  distinct() %>%
  group_by(site, endPoint) %>%
  filter(EARsum >= 0.01) %>%
  group_by(site, endPoint) %>%
  summarize(exceed_n = n(), .groups = "drop") %>%
  full_join(sample_n, by= "site") %>%
  mutate(exceed_percent = exceed_n / sample_n) %>%
  filter(!is.na(endPoint))

plot2 <- ggplot(EARmix_freq, aes(x=site, y=endPoint, fill=exceed_percent)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90), axis.text = element_text(size=8)) +
  ggtitle("EARmix > 0.01")
  
print(plot2)

ggsave(file.path(path_to_data, "Figures", "endPoint_Excedances_Tile_surface.png"), plot2, height=10, width=8)

EARmix_gene_freq <- EARmix_freq %>%
  left_join(unique(select(sample_EARmix, endPoint, geneSymbol)), by = 'endPoint') %>%
  distinct() %>%
  group_by(site, geneSymbol) %>%
  summarize(exceed_percent = max(exceed_percent), 
            n_endpoints=n(),
            .groups = "drop") %>%
  filter(!is.na(geneSymbol))

EARmix_gene_freq %>%
  select(geneSymbol, n_endpoints) %>%
  group_by(geneSymbol) %>%
  summarize(n_endpoint = max(n_endpoints)) %>%
  arrange(desc(n_endpoint))

plot3 <- ggplot(EARmix_gene_freq, aes(x=site, y=geneSymbol, fill=exceed_percent)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90), axis.text = element_text(size=8)) +
  ggtitle("EARmix > 0.01")

print(plot3)

ggsave(file.path(path_to_data, "Figures", "Gene_Excedances_Tile_surface.png"), plot3, height=10, width=8)

sample_table_out <- sample_EARmix %>%
  # filter(site %in% c("04249000", "04119400")) %>%
  group_by(site) %>%
  summarise(all_chems = paste(chems, collapse= "|"),
            all_CASs = paste(CASs, collapse= "|"),
            all_genes = paste(unique(geneSymbol[!is.na(geneSymbol)]), collapse = "|"),
            all_gene_names = paste(unique(geneSymbol[!is.na(geneSymbol)]), collapse = "|"),
            all_endPoints = paste(unique(endPoint[!is.na(endPoint)]), collapse = "|"),
            max_n_chems_per_endPoint = max(n_chems),
            .groups = "drop") %>%
  mutate(all_chems =   unlist(lapply(strsplit(all_chems, "\\|"), 
                                     function(l) paste(unique(l), collapse = "|"))),
         all_CASs =   unlist(lapply(strsplit(all_CASs, "\\|"), 
                                     function(l) paste(unique(l), collapse = "|")))) %>%
  select(site, max_n_chems_per_endPoint, all_genes, everything())
  

data.frame(sample_table_out)

write.csv(sample_table_out, file.path(path_to_data, "Data", "GeneSummaryBySite_SurfaceWater.csv"), row.names=F)

top_combos_surf <- top_mixes(cs, 
                             group_by_this, 
                             ear_threshold,
                             site_threshold) 

#List all chemicals and endpoints
all_max <- overall_mixtures(top_combos_surf, "max") %>%
  select(genes, endPoint, n_chems, chems, everything()) %>%
  arrange(genes)

priority_max_fancy <- overall_df_format(all_max, class_key) 

write.csv(all_max, file = file_out(file.path(path_to_data, 'Data', 'priority_endpoints_fullchemlist_surface_pesticides.csv')), row.names = F)

knitr::kable(head(priority_max_fancy))


# filter using number of samples
priority_combos_surf <- filter(top_combos_surf, n_samples>=20)



#Connect AOPs and first try at panther pathways
priority_max <- overall_mixtures(priority_combos_surf, "max")

class_key <- class_key_fnx(cs)

priority_max_fancy <- overall_df_format(priority_max, class_key)


write.csv(priority_max, file = file_out(file.path(path_to_data, 'Data', 'priority_mixtures_surface_withDegs_12months.csv')), row.names = F)

# Gene summaries and details annotations
gene_summary_all <- gene_summary(priority_max)
gene_summary_entrez <- gene_summary(priority_max, columns = 'entrez') %>%
  filter(nchar(ENTREZ_GENE_SUMMARY)>0)


# write.csv(gene_summary_all, file = file_out(file.path(path_to_data, 'Data', 'gene_summary_mixtures_surface_12months.csv')), row.names = F)

# write.csv(gene_summary_entrez, file = file_out(file.path(path_to_data, 'Data', 'gene_entrez_mixtures_surface_12months.csv')), row.names = F)



# ################################
# No guarantee things work below. 
# Need ACClong from toxeval file creation
# ################################

endpoint_review_surf <- priority_max %>%
  dplyr::select(endPoint, CASs) %>%
  tidyr::separate(CASs, sep='\\|', into = letters[1:10], fill='right') %>%
  gather(key = chem_mix_nu, value=CAS, -1) %>%
  drop_na(CAS) %>%
  dplyr::select(-chem_mix_nu)


ACCs_surf <- dplyr::select(ACClong, endPoint, chnm, CAS, ACC, ACC_value) %>%
  filter(endPoint %in% endpoint_review_surf$endPoint, 
         CAS %in% endpoint_review_surf$CAS)

endpoint_out_surf <- endpoint_review_surf %>%
  left_join(ACCs_surf) %>%
  arrange(endPoint, ACC) %>%
  distinct()

data.frame(endpoint_out_surf)






summed_EARs_surf <- cs %>%
  group_by(site, date, !!sym(group_by_this)) %>%
  summarize(sum_ear_endpoint = sum(EAR, na.rm = TRUE), .groups="drop") %>%
  filter(!!sym(group_by_this) %in% priority_combos_surf$endPoint) %>%
  distinct() %>%
  # left_join(unique(select(priorityEndpoints_surf, endPoint))) %>%
  group_by(site, endPoint) %>%
  summarize(sum_ear_endpoint = max(sum_ear_endpoint, na.rm = TRUE),  
            n=n(), .groups="drop") %>%
  left_join(unique(select(cs, endPoint, Bio_category)))

endpoint_rank_surf <- summed_EARs_surf %>%
  group_by(endPoint) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T), 
            .groups = "drop") %>%
  arrange((sum_ear_median))

summed_EARs_surf <- summed_EARs_surf %>%
  left_join(priority_max) %>%
  left_join(unique(select(join_criteria(), geneSymbol,endPoint)), by = "endPoint") %>%
  mutate(endPoint = factor(endPoint, endpoint_rank_surf$endPoint)) %>%
  select(-chems, -CASs, -n_chems) %>%
  distinct() 


top_mix_box_surf <- ggplot(unique(summed_EARs_surf[c('endPoint', 'sum_ear_endpoint', 
                                                     'Bio_category')]), aes(y=endPoint, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='blue', alpha=0.5, width=0, height=.1, size=1.5) +
  geom_boxplot(outlier.shape=NA, fill='blue', alpha=.2, width=.5) +
  scale_x_log10(name = expression(paste(EAR[mixture]))) +
  # scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='ToxCast assay name', fill='Bio Category', color='Bio Category') + 
  theme_tox() +
  theme(legend.position='right', legend.title.align=0.5) +
  guides(col = guide_legend(ncol = 1))

print(top_mix_box_surf)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_surf_withDegs_12months.png")), top_mix_box_surf, height=9, width=7, units='in')

summed_EARs_surf$geneSymbol[is.na(summed_EARs_surf$geneSymbol)] <- paste0(summed_EARs_surf$Bio_category[is.na(summed_EARs_surf$geneSymbol)], "*")

median_EARs_byGene_surf <- summed_EARs_surf %>%
  group_by(geneSymbol, site) %>%
  summarize(sum_ear_endpoint = max(sum_ear_endpoint, na.rm=T), 
            .groups = "drop")

gene_rank_surf <- median_EARs_byGene_surf %>%
  group_by(geneSymbol) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T),
            .groups = "drop") %>%
  mutate(order = ifelse(grepl("\\*", geneSymbol), 2, 1)) %>%
  arrange(desc(order), sum_ear_median)

median_EARs_byGene_surf <- median_EARs_byGene_surf %>%
  left_join(select(gene_rank_surf, geneSymbol, order)) %>%
  mutate(geneSymbol = factor(geneSymbol, gene_rank_surf$geneSymbol))

top_mixgene_box_surf <- ggplot(median_EARs_byGene_surf, 
                               aes(y=geneSymbol, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='blue', alpha=0.5, width=0, height=.1, size=1.2) +
  geom_boxplot(outlier.shape=NA, fill='blue', alpha=.2, width=.5) +
  scale_x_log10(name = expression(paste(EAR[mixture]))) +
  # scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='Gene Target', fill='Bio Category', color='Bio Category') + 
  theme_tox() +
  theme(legend.position='right', legend.title.align=0.5) +
  guides(col = guide_legend(ncol = 1))

print(top_mixgene_box_surf)

ggsave(file_out(file.path(path_to_data, "Figures", "PriorityGenesBoxplot_surf_withDegs_12months.png")), top_mixgene_box_surf, height=5, width=4, units='in')





