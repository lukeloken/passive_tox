
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
cs_surf <- chemicalSummary_surface %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  filter(`date`>=date_filter[1], `date`<=date_filter[2],
         shortName %in% site_order) %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  group_by(site, endPoint, shortName, Lake, site_grouping, CAS, chnm, Bio_category, Class) %>%
  summarise(EAR = max(EAR, na.rm=T), .groups = 'drop') 

# Prepare passive data
cs_pass <- chemicalSummary %>%
  mutate(shortName = factor(shortName, site_order))

# thresholds for mixtures analysis
group_by_this <- "endPoint"
ear_threshold <- .001 #Ear threshold
site_threshold <- 1.9 #How many sites

#Passive
allcombos_pass <- top_mixes(cs_pass, group_by_this, 
                            ear_threshold,
                            1) 

allmixtures_pass <- overall_mixtures(allcombos_pass, "max")

endpoints_pass <- allmixtures_pass %>%
  dplyr::select(endPoint, CASs) %>%
  tidyr::separate(CASs, sep='\\|', into = letters[1:10], fill='right') %>%
  gather(key = chem_mix_nu, value=CAS, -1) %>%
  drop_na(CAS) %>%
  dplyr::select(-chem_mix_nu)
# left_join(dplyr::select(ACClong, endPoint, ACC, ACC_value, chnm))
# left_join(unique(dplyr::select(ACClong, endPoint, chnm)))

ACCs <- dplyr::select(ACClong, endPoint, chnm, CAS, ACC, ACC_value) %>%
  filter(endPoint %in% endpoints_pass$endPoint, 
         CAS %in% endpoints_pass$CAS)

endpoints_pass <- endpoints_pass %>%
  left_join(ACCs) %>%
  arrange(endPoint, ACC) %>%
  distinct()

data.frame(endpoints_pass)


topcombos_pass <- top_mixes(cs_pass, group_by_this, 
                            ear_threshold,
                            site_threshold) 


topmixtures_pass <- overall_mixtures(topcombos_pass, "max")


genetable_pass <- gene_mixtures(topmixtures_pass)  


genesummary_pass <- gene_summary(topmixtures_pass, species = c("Homo sapiens", 
                                                               "Danio rerio", 
                                                               "Xenopus tropicalis"))

geneentrez_pass <- gene_summary(topmixtures_pass, 
                                columns = "entrez", 
                                species = "Homo sapiens") %>%
  filter(species == "Homo sapiens")


geneinteresting_pass <- genesummary_pass %>%
  select(gene_abbr, species, gene_name, pathway_name, ENTREZ_GENE_SUMMARY, 
         GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, 
         KEGG_PATHWAY, OMIM_DISEASE) %>%
  distinct()

write.csv(geneinteresting_pass, 
          file = file_out(file.path(path_to_data, 'Data',
                                    'passive_prioritygene_interesting.csv')), 
          row.names = F)


geneoutput_pass <- geneinteresting_pass %>%
  group_by(gene_abbr) %>%
  summarize(orthologs = paste(unique(species[!is.na(species)]), 
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
  left_join(unique(select(topmixtures_pass, genes, AOP_ids, AOP_names)), by = c("gene_abbr" = "genes"))


data.frame(geneoutput_pass)

all_kegg <- paste(unique(geneoutput_pass$KEGG_PATHWAY), collapse= "|")
all_kegg <- unique(unlist(strsplit(all_kegg, "\\|")))
all_kegg <- unique(unlist(strsplit(all_kegg, ":")))
length(all_kegg[grepl('hsa', all_kegg)])
length(all_kegg[grepl('dre', all_kegg)])
length(all_kegg[grepl('xtr', all_kegg)])

all_kegg <- all_kegg[!grepl('hsa', all_kegg)]
all_kegg <- all_kegg[!grepl('dre', all_kegg)]
all_kegg <- all_kegg[!grepl('xtr', all_kegg)]

geneoutput_pass <- full_join(genetable_pass, geneoutput_pass, 
                             by = c("genes" = "gene_abbr")) %>%
  rename(GeneSymbol = genes) 

geneoutput_pass[is.na(geneoutput_pass)] <- ""


write.csv(geneoutput_pass, file = file_out(file.path(path_to_data, 'Data', 'passive_prioritygene_interesting2.csv')), row.names = F)




EARmix_pass <- EAR_mixtures(cs_pass, group_by_this)

EARsite_pass <- site_mixtures(EARmix_pass, ear_cutoff = 0.001)

data.frame(EARsite_pass)

#plotting

plot_genebar(EARmix_pass, ear_threshold, site_threshold, type="gene")

plot_genebar(EARmix_pass, ear_threshold, site_threshold, type="endPoint")





#Surface samples


cs_surf <- chemicalSummary_surface %>%
  mutate(shortName = factor(shortName, site_order)) %>%
  filter(`date`>=date_filter[1], `date`<=date_filter[2],
         shortName %in% site_order) %>%
  mutate(shortName = factor(shortName, site_order))

allcombos_surf <- top_mixes(cs_surf, group_by_this, 
                             ear_threshold,
                             .9) 

prioritycombos_surf <- filter(allcombos_surf, n_samples>10)
prioritymixtures_pass <- overall_mixtures(prioritycombos_surf, "max")

class_key <- class_key_fnx(cs_surf)

prioritymixtures_fancy <- overall_df_format(prioritymixtures_pass, class_key)

# write.csv(prioritymixtures_fancy, file = file_out(file.path(path_to_data, 'Data', 'priority_mixtures_surface.csv')), row.names = F)

allmixtures_surf <- overall_mixtures(allcombos_surf, "max")



endpoints_surf <- allmixtures_surf %>%
  dplyr::select(endPoint, CASs) %>%
  tidyr::separate(CASs, sep='\\|', into = letters[1:10], fill='right') %>%
  # select_if(~sum(!is.na(.)) > 0) %>%
  gather(key = chem_mix_nu, value=CAS, -1) %>%
  drop_na(CAS) %>%
  dplyr::select(-chem_mix_nu)
# left_join(dplyr::select(ACClong, endPoint, ACC, ACC_value, chnm))
# left_join(unique(dplyr::select(ACClong, endPoint, chnm)))

ACCs_surf <- dplyr::select(ACClong, endPoint, chnm, CAS, ACC, ACC_value) %>%
  filter(endPoint %in% endpoints_surf$endPoint, 
         CAS %in% endpoints_surf$CAS)

endpoints_surf <- endpoints_surf %>%
  left_join(ACCs_surf) %>%
  arrange(endPoint, ACC) %>%
  distinct()

data.frame(endpoints_surf)


topcombos_surf <- top_mixes(cs_surf,
                             group_by_this, 
                             ear_threshold,
                             site_threshold) 


topmixtures_surf <- overall_mixtures(topcombos_surf, "max")

genetable_surf <- gene_mixtures(topmixtures_surf)  


genesummary_surf <- gene_summary(topmixtures_surf, species = c("Homo sapiens", 
                                                               "Danio rerio", 
                                                               "Xenopus tropicalis"))

geneentrez_surf <- gene_summary(topmixtures_surf, 
                                columns = "entrez", 
                                species = "Homo sapiens") %>%
  filter(species == "Homo sapiens")


geneinteresting_surf <- genesummary_surf %>%
  select(gene_abbr, species, gene_name, pathway_name, ENTREZ_GENE_SUMMARY, 
         GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, 
         KEGG_PATHWAY, OMIM_DISEASE) %>%
  distinct()

write.csv(geneinteresting_surf, 
          file = file_out(file.path(path_to_data, 'Data',
                                    'surfacesummer_prioritygene_interesting.csv')), 
          row.names = F)


geneoutput_surf <- geneinteresting_surf %>%
  group_by(gene_abbr) %>%
  summarize(orthologs = paste(unique(species[!is.na(species)]), 
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
  left_join(unique(select(topmixtures_surf, genes, AOP_ids, AOP_names)), by = c("gene_abbr" = "genes"))


data.frame(geneoutput_surf)


geneoutput_surf <- full_join(genetable_surf, geneoutput_surf, 
                             by = c("genes" = "gene_abbr")) %>%
  rename(GeneSymbol = genes) 

geneoutput_surf[is.na(geneoutput_surf)] <- ""


write.csv(geneoutput_surf, file = file_out(file.path(path_to_data, 'Data', 'surfacesummer_prioritygene_interesting2.csv')), row.names = F)


EARmix_surf <- EAR_mixtures(cs_surf, group_by_this)

EARsite_surf <- site_mixtures(EARmix_surf, ear_cutoff = 0.001)

data.frame(EARsite_surf)

#plotting

plot_genebar(EARmix_surf, ear_threshold, site_threshold, type="gene")

plot_genebar(EARmix_surf, ear_threshold, site_threshold, type="endPoint")


