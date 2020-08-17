# ##########################
# Mixture analysis for OWCs 
# ##########################

library(toxEval)
library(ToxMixtures)
library(ggplot2)
library(dplyr)

path_to_data <- c('C:/Users/lloken/DOI/Corsi, Steven R - GLRI CECs/2016/Manuscripts/Pesticides_passive')

toxeval_OWC <- create_toxEval((file.path("C:/Users/lloken/OneDrive - DOI/Loken_Personal", "toxEval_InputFile20200304.xlsx")))

ACClong <- get_ACC(toxeval_OWC$chem_info$CAS)
ACClong <- remove_flags(ACClong)
cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))


chemicalSummary_OWC <- get_chemical_summary(toxeval_OWC, 
                                            ACClong, 
                                            filtered_ep) %>%
  mutate(date = as.Date(date, tz = 'America/Chicago'))

# The Chemical Summary (cs) table is the gateway into the ToxMixtures package
# Note on the toxeval file, the date column contains date_times and their are multiple date_times per site. I'm wondering if this is correct
# For example, the following code shows two different times associated with this site
table(filter(toxeval_OWC$chem_data, SiteID == "4157000")[,4])

cs <- chemicalSummary_OWC

#Terms for grouping and filtering
group_by_this <- "endPoint"
ear_threshold <- .1
site_threshold <- floor(length(unique(toxeval_OWC$chem_site$SiteID))*.1) #Number of sites minimum (e.g., 10%)

#Calculate EARmix for all samples. 
EAR_mix <- EAR_mixtures(cs, group_by_this)

#Summary table by site. Endpoints and genes included for all Mixtures that exceed the threshold
site_mix <- site_mixtures(EAR_mix, ear_threshold)

site_mix

write.csv(site_mix, file.path(path_to_data, "Data", "GeneSummaryBySite_OWC.csv"), row.names=F)

#Look at all mixtures that have EARmix exceedances in enough sites
top_combos <- top_mixes(cs, 
                        group_by_this, 
                        ear_threshold,
                        site_threshold) 

#Filter combos to the ones with the most chemicals 
all_max <- overall_mixtures(top_combos, "max") %>%
  select(genes, endPoint, n_chems, chems, everything()) %>%
  arrange(genes)

class_key <- class_key_fnx(cs)

all_max_fancy <- overall_df_format(all_max, class_key) 

write.csv(all_max, file = file_out(file.path(path_to_data, 'Data', 'priority_endpoints_OWC.csv')), row.names = F)

#List of any chemical that contributed to an EAR exceedance
unique(unlist(strsplit(paste(all_max$chems, collapse = "|"), "\\|")))

chem_table <- as.data.frame(table(unlist(strsplit(paste(all_max$chems, collapse = "|"), "\\|"))))
names(chem_table) <-  c('chnm', 'frequency')

chem_table %>% arrange(desc(frequency))

#List of all genes with an EAR exceedance
unique(unlist(strsplit(all_max$genes, ",")))


# Gene summaries and details annotations
gene_summary_all <- gene_summary(all_max)
gene_summary_entrez <- gene_summary(all_max, columns = "entrez", species = "Homo sapiens") %>%
  filter(nchar(ENTREZ_GENE_SUMMARY)>0)

gene_summary_interesting <- gene_summary_all %>%
  select(gene_abbr, species, gene_name, pathway_name, ENTREZ_GENE_SUMMARY, 
         GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, 
         KEGG_PATHWAY, OMIM_DISEASE)

write.csv(gene_summary_all, file = file_out(file.path(path_to_data, 'Data', 'gene_summary_full_OWC.csv')), row.names = F)

write.csv(gene_summary_interesting, file = file_out(file.path(path_to_data, 'Data', 'gene_summary_slim_OWC.csv')), row.names = F)

write.csv(gene_summary_entrez, file = file_out(file.path(path_to_data, 'Data', 'gene_entrez_OWC.csv')), row.names = F)

#Filter EARmix to those in top_combos. Need to have enough EAR exceedances
EAR_mix_priority <- EAR_mix %>%
  filter(endPoint %in% top_combos$endPoint)

gene_priority_fig <- plot_genebar(EAR_mix_priority, ear_threshold, type = "gene") +
  ggtitle(paste("EAR > ", ear_threshold, " in ", site_threshold, "sites or more"))

gene_priority_fig

endPoint_priority_fig <- plot_genebar(EAR_mix_priority, ear_threshold, type = "endPoint") +
  ggtitle(paste("EAR > ", ear_threshold, " in ", site_threshold, "sites or more"))

endPoint_priority_fig


ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot_OWC.png")), endPoint_priority_fig, height=12, width=7, units='in')


ggsave(file_out(file.path(path_to_data, "Figures", "PriorityGenesBoxplot_OWC.png")), gene_priority_fig, height=7, width=7, units='in')





