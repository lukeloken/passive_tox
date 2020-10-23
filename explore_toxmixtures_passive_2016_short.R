
# #########
# Mixtures
# #########

#### Setup ####
library(toxEval)
library(ToxMixtures)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

#Need to run (ToxEval_Passive2016.R and ToxEval_WaterSamples2016.R)

# Take maximum EAR for water samples with more than one sample
# cs_surf <- chemicalSummary_surface %>%
#   mutate(shortName = factor(shortName, site_order)) %>%
#   filter(`date`>=date_filter[1], `date`<=date_filter[2],
#          shortName %in% site_order) %>%
#   mutate(shortName = factor(shortName, site_order)) %>%
#   group_by(site, endPoint, shortName, Lake, site_grouping, CAS, chnm, Bio_category, Class) %>%
#   summarise(EAR = max(EAR, na.rm=T), .groups = 'drop') 

# Prepare passive data
cs_pass <- chemicalSummary %>%
  mutate(shortName = factor(shortName, site_order))

# thresholds for mixtures analysis
group_by_this <- "endPoint"
ear_threshold <- .001 #Ear threshold
site_threshold <- 1.9 #How many sites

# ################
# Passive analysis
# ################

# Calculate EARMix and EARAOP for each sample and endpoint/AOP
EAR_mix_pass <- EAR_mixtures(cs_pass, group_by_this)
AOP_mix_pass <- EAR_mixtures(cs_pass, "AOP")

# Summarize endpoints, chemicals, and genes by site that exceed EAR thresholds
EARsite_pass <- site_mixtures(EAR_mix_pass, ear_cutoff = 0.001)

AOPsite_pass <- site_mixtures(AOP_mix_pass, ear_threshold) %>%
  rename(all_chems_AOP = all_chems,
         all_CASs_AOP = all_CASs)

site_mix_all_pass <- EARsite_pass %>%
  rename(all_chems_endPoint = all_chems,
         all_CASs_endPoint = all_CASs) %>%
  full_join(AOPsite_pass, by = "site")

write.csv(site_mix_all_pass, file.path(path_to_data, "Data", "PassiveSiteMixtures_Table1.csv"), row.names = FALSE)


#Get top endpoint combos
topcombos_pass <- top_mixes(cs_pass, group_by_this, 
                            ear_threshold,
                            site_threshold) 


topmixtures_pass <- overall_mixtures(topcombos_pass, "max")


topaopcombos_pass <- top_mixes(cs_pass, 
                             group_by_what = "AOP",
                             ear_threshold,
                             site_threshold)

topAOPmixtures_pass <- overall_mixtures(topaopcombos_pass, "max")


genetable_pass <- gene_mixtures(topmixtures_pass)  


genefunction_pass <- gene_functions(topmixtures_pass, species = c("Homo sapiens", 
                                                               "Danio rerio", 
                                                               "Xenopus tropicalis"))

geneentrez_pass <- gene_functions(topmixtures_pass, 
                                columns = "entrez", 
                                species = "Homo sapiens") %>%
  filter(Species == "Homo sapiens")

gene_summarypass <- gene_summary(genefunction_pass)

gene_mixespass <- gene_mixtures(topmixtures_pass)

gene_outpass <- full_join(gene_mixespass, gene_summarypass, 
                      by = "geneSymbol")


write.csv(gene_outpass, file = file_out(file.path(path_to_data, 'Data', 'PassiveGeneMixtures_Table2.csv')), row.names = F)


gene_wb <- create_Excel_wb_gene(gene_summarypass, gene_mixespass)

saveWorkbook(gene_wb, file = file.path(path_to_data, "Data", "Gene_Workbook_passive.xlsx"), overwrite = TRUE)

# priority endpoints
priorityEndpoints <- filter(cs_pass, EAR>0.001) %>%
  dplyr::select(chnm, endPoint) %>%
  unique() %>%
  left_join(cs_pass) %>%
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

#plotting

box_gene_pass <- plot_genebar(EAR_mix_pass, ear_threshold, site_threshold, type="geneSymbol", fill = 'nothing') +
  scale_fill_manual(values = c("darkgreen")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none')

box_endpoint_pass <- plot_genebar(EAR_mix_pass, ear_threshold, site_threshold, type="endPoint", fill = 'nothing') +
  scale_fill_manual(values = c("darkred")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none')

box_AOP_pass <- plot_genebar(AOP_mix_pass, ear_threshold, site_threshold, type="AOP", fill = 'nothing') +
  scale_fill_manual(values = c("darkblue")) + 
  scale_x_log10nice(name = expression(paste(EAR[AOP]))) +
  theme(legend.position = 'none')

png(file.path(path_to_data, "Figures", "PriorityBoxplot_3panel_v2.png"), height = 4, width = 9, units = "in", res = 300)
grid.newpage()
boxes_pass <- grid.draw(cbind(ggplotGrob(box_endpoint_pass), 
                              ggplotGrob(box_gene_pass), 
                              ggplotGrob(box_AOP_pass)))
dev.off()

# Custom plot to look at individual chemicals on endpoint/genes
priority_combos <- cs_pass %>%
  filter(EAR > ear_threshold) %>%
  unite(chnm, endPoint, col = "chnm_endPoint", sep = "|", remove = FALSE)
  
cs_pass_modified <- cs_pass %>%
  unite(chnm, endPoint, col = "chnm_endPoint", sep = "|", remove = FALSE) %>%
  filter(chnm_endPoint %in% priority_combos$chnm_endPoint) %>%
  rename(EARsum = EAR) %>%
  left_join(join_criteria())

# plot_genebar(cs_pass_modified, ear_threshold, site_threshold, type = "endPoint", facet_col = "chnm")

gene_bychem_boxplot <- plot_genebar(cs_pass_modified, ear_threshold, site_thres = 0.9, 
                                    type = "geneSymbol", facet_col = "chnm", fill = "nothing") +
  scale_x_log10nice(name = expression(paste("EAR"))) +
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6))

print(gene_bychem_boxplot)  

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_byChemical_v2.png")), gene_bychem_boxplot, height=4, width=6.7, units='in')


# ################
# Surface samples
# ################

cs_surf <- chemicalSummary_surface %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  filter(`date`>=date_filter[1], `date`<=date_filter[2],
         shortName %in% site_order) %>%
  mutate(shortName = factor(shortName, site_order))


# Calculate EARMix and EARAOP for each sample and endpoint/AOP
EAR_mix_surf <- EAR_mixtures(cs_surf, group_by_this)
AOP_mix_surf <- EAR_mixtures(cs_surf, "AOP")

# Summarize endpoints, chemicals, and genes by site that exceed EAR thresholds
EARsite_surf <- site_mixtures(EAR_mix_surf, ear_cutoff = 0.001)

AOPsite_surf <- site_mixtures(AOP_mix_surf, ear_threshold) %>%
  rename(all_chems_AOP = all_chems,
         all_CASs_AOP = all_CASs)

site_mix_all_surf <- EARsite_surf %>%
  rename(all_chems_endPoint = all_chems,
         all_CASs_endPoint = all_CASs) %>%
  full_join(AOPsite_surf, by = "site")



#Get top endpoint combos
topcombos_surf <- top_mixes(cs_surf, group_by_this, 
                            ear_threshold,
                            site_threshold) 


topmixtures_surf <- overall_mixtures(topcombos_surf, "max")

topaopcombos_surf <- top_mixes(cs_surf, 
                               group_by_what = "AOP",
                               ear_threshold,
                               site_threshold)

topAOPmixtures_surf <- overall_mixtures(topaopcombos_surf, "max")


genetable_surf <- gene_mixtures(topmixtures_surf)  


genefunction_surf <- gene_functions(topmixtures_surf, species = c("Homo sapiens", 
                                                                  "Danio rerio", 
                                                                  "Xenopus tropicalis"))

geneentrez_surf <- gene_functions(topmixtures_surf, 
                                  columns = "entrez", 
                                  species = "Homo sapiens") %>%
  filter(Species == "Homo sapiens")

gene_summarysurf <- gene_summary(genefunction_surf)

gene_mixessurf <- gene_mixtures(topmixtures_surf)

gene_outsurf <- full_join(gene_mixessurf, gene_summarysurf, 
                          by = "geneSymbol")


write.csv(gene_outsurf, file = file_out(file.path(path_to_data, 'Data', 'SurfaceGeneMixtures_Table2.csv')), row.names = F)


#plotting

box_gene_surf <- plot_genebar(EAR_mix_surf, ear_threshold, site_threshold, type="geneSymbol", fill = 'nothing') +
  scale_fill_manual(values = c("darkgreen")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none')

box_endpoint_surf <- plot_genebar(EAR_mix_surf, ear_threshold, site_threshold, type="endPoint", fill = 'nothing') +
  scale_fill_manual(values = c("darkred")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none')

box_AOP_surf <- plot_genebar(AOP_mix_surf, ear_threshold, site_threshold, type="AOP", fill = 'nothing') +
  scale_fill_manual(values = c("darkblue")) + 
  scale_x_log10nice(name = expression(paste(EAR[AOP]))) +
  theme(legend.position = 'none')

# Custom plot to look at individual chemicals on endpoint/genes
priority_combos_surf <- cs_surf %>%
  filter(EAR > ear_threshold) %>%
  unite(chnm, endPoint, col = "chnm_endPoint", sep = "|", remove = FALSE)

cs_surf_modified <- cs_surf %>%
  unite(chnm, endPoint, col = "chnm_endPoint", sep = "|", remove = FALSE) %>%
  filter(chnm_endPoint %in% priority_combos_surf$chnm_endPoint) %>%
  rename(EARsum = EAR) %>%
  left_join(join_criteria())

# plot_genebar(cs_pass_modified, ear_threshold, site_threshold, type = "endPoint", facet_col = "chnm")

gene_bychem_boxplot_surf <- plot_genebar(cs_surf_modified, ear_threshold, site_thres = 1.9, 
                                    type = "geneSymbol", facet_col = "chnm", fill = "nothing") +
  scale_x_log10nice(name = expression(paste("EAR"))) +
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6))

print(gene_bychem_boxplot_surf)  

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_byChemical_surface_v2.png")), gene_bychem_boxplot_surf, height=5, width=12, units='in')




#Combined plot
EAR_mix_surf <- EAR_mix_surf %>% 
  mutate(method = "Surface water")

EAR_mix_combined <- EAR_mix_pass %>% 
  mutate(method = "Passive") %>%
  mutate(date = as.Date(paste0(date, "-07-01"))) %>%
  bind_rows(EAR_mix_surf)

EAR_mix_combined$geneSymbol[grepl("Tanguay_ZF", EAR_mix_combined$endPoint)] <- "Zebrafish*"

gene_combined_boxplot <- plot_genebar(EAR_mix_combined, ear_threshold, site_threshold, 
             type="geneSymbol", fill = method, facet_col = "method") +
  # scale_fill_manual(values = c("darkgreen")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture])), expand = c(0.03,0.03)) +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(fill = NA, color = NA))

print(gene_combined_boxplot)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityGenesBoxplot_twoMethods_v2.png")), 
       gene_combined_boxplot, height=4.5, width=5, units='in')


endPoint_combined_boxplot <- plot_genebar(EAR_mix_combined, ear_threshold, site_threshold, 
                                      type="endPoint", fill = method, facet_col = "method") +
  # scale_fill_manual(values = c("darkgreen")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(fill = NA, color = NA))

print(endPoint_combined_boxplot)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_twoMethods_v2.png")), 
       endPoint_combined_boxplot, height=5, width=6, units='in')

#Combined AOP
AOP_mix_surf <- AOP_mix_surf %>% 
  mutate(method = "Surface water")

AOP_mix_combined <- AOP_mix_pass %>% 
  mutate(method = "Passive") %>%
  mutate(date = as.Date(paste0(date, "-07-01"))) %>%
  bind_rows(AOP_mix_surf)

AOP_combined_boxplot <- plot_genebar(AOP_mix_combined, ear_threshold, site_threshold, 
                                      type="AOP", fill = method, facet_col = "method") +
  # scale_fill_manual(values = c("darkgreen")) + 
  scale_x_log10nice(name = expression(paste(EAR[AOP])), expand = c(0.03,0.03)) +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(fill = NA, color = NA))

print(AOP_combined_boxplot)

ggsave(file_out(file.path(path_to_data, "Figures/PriorityAOPBoxplot_twoMethods_v2.png")), 
       AOP_combined_boxplot, height=4.5, width=4, units='in')




# #############
# Tables for SI
# #############
table_SI7 <- EAR_mix_pass %>%
  select(endPoint) %>%
  distinct() %>%
  left_join(AOP_crosswalk, by = "endPoint") %>%
  filter(!is.na(ID)) %>%
  select(AOP_number = ID, AOP_title, 
         key_event_number = KE, key_event_name = `Key Event Name`,
         key_event_type = `KeyEvent Type`) %>%
  arrange(AOP_number, key_event_number) %>%
  distinct()


KE_table <- filter(table_SI7, !is.na(key_event_name)) %>%
  select(key_event_number, key_event_name, key_event_type) %>%
  distinct() %>%
  arrange(key_event_number)

table_SI7$key_event_name <- KE_table$key_event_name[match(table_SI7$key_event_number, KE_table$key_event_number)]


table_SI7 <- table_SI7 %>%
  arrange(AOP_number, key_event_number) %>%
  distinct()

write.csv(table_SI7, file.path(path_to_data, "SI tables", "AOP_table_SI7.csv"), 
          row.names = F)


table_SI8 <- cs_pass %>%
  select(chnm, CAS, endPoint, EAR) %>%
  full_join(select(cs_surf, chnm, CAS, endPoint, EAR)) %>%
  # filter(EAR > 0) %>%
  distinct() %>%
  left_join(select(end_point_info, assay_component_endpoint_name, intended_target_gene_symbol), 
            by = c("endPoint" = "assay_component_endpoint_name")) %>%
  rename(geneSymbol = intended_target_gene_symbol) %>%
  mutate(geneSymbol = toupper(geneSymbol)) %>%
  distinct()


table_SI8$geneSymbol[which(table_SI8$geneSymbol == "CYP3A23/3A1")] <- "CYP3A23|CYP3A1"

table_SI8$geneSymbol[which(table_SI8$geneSymbol == "THRA|THRB|THRB|THRB |THRA")] <- "THRA|THRB"

# unique(table_SI8$geneSymbol)

table_SI8 <- table_SI8 %>%
  mutate(all_genes = c(strsplit(geneSymbol, split = "\\|")),
         all_genes = sapply(all_genes, function(x) x[!(x %in% "")]),
         all_genes = sapply(all_genes, function(x) paste(unique(x), collapse = "|")),
         geneSymbol = ifelse(all_genes == "NA", "", all_genes)) %>% 
  select(-all_genes) %>%
  separate_rows(geneSymbol, sep = "\\|") %>%
  distinct()



AOP_tojoin <- AOP_crosswalk %>%
  filter(endPoint %in% table_SI8$endPoint) %>%
  mutate(AOPs = "yes") %>%
  select(endPoint, AOPs) %>%
  distinct()

DAVID_tojoin <- david_full %>%
  filter(geneSymbol %in% c(table_SI8$geneSymbol)) %>%
  mutate(DAVID = "yes") %>%
  select(geneSymbol, DAVID) %>%
  distinct()

Panther_tojoin <- panther %>%
  filter(geneSymbol %in% c(table_SI8$geneSymbol)) %>%
  mutate(PANTHER = "yes") %>%
  select(geneSymbol, PANTHER) %>%
  distinct()

table_SI8 <- table_SI8 %>%
  left_join(AOP_tojoin, by = "endPoint") %>%
  left_join(DAVID_tojoin, by = "geneSymbol") %>%
  left_join(Panther_tojoin, by = "geneSymbol") %>%
  mutate(AOPs = ifelse(is.na(AOPs), "no", AOPs),
         DAVID = ifelse(is.na(DAVID), "no", DAVID),
         PANTHER = ifelse(is.na(PANTHER), "no", PANTHER)) %>%
  select(geneSymbol, chemical_name = chnm, CAS, ToxCast_assay_endpoint = endPoint, AOPs, DAVID, PANTHER) %>%
  mutate(geneSymbol = ifelse(geneSymbol == "", NA, geneSymbol)) %>%
  arrange(geneSymbol, chemical_name, ToxCast_assay_endpoint) %>%
  distinct()

table_SI8

write.csv(table_SI8, file.path(path_to_data, "SI tables", "Gene_table_SI8.csv"), 
          row.names = F)

# ####
# End
# ####

