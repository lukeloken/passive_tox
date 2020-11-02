
# Analyze synthesis dataset
# Loken Oct 2020


#### Setup ####
path_to_data <- c('C:/Users/lloken/DOI/Corsi, Steven R - GLRI CECs/2016/Manuscripts/Pesticides_passive')

library(toxEval)
library(ToxMixtures)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(openxlsx)


tox_file <- create_toxEval(file.path(path_to_data, 
                                     "ToxEvalFiles",
                                     "USGS_Water_data_synthesis 1 with exclusions.xlsx"))

#Add some grouping factors
tox_file$chem_site <- tox_file$chem_site %>%
  mutate(State = factor(State, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY')),
         Lake  = factor(Lake, c('Superior', 'Michigan', 'Huron', 'Erie', 'Ontario')))

class_table <- table(tox_file$chem_info$Class) 
study_table <- table(tox_file$chem_info$type) 


#Create a new tox_eval file using the maximum concentration from each site/chemical
tox_file_max <- tox_file

tox_file_max$chem_info$`Chemical Name` <- gsub(", .*", "", tox_file_max$chem_info$`Chemical Name`)

tox_file_max$chem_info <- tox_file_max$chem_info %>%
  group_by(CAS) %>%
  summarize(across(c(`Chemical Name`, Class), function (x) x[x != "" & !is.na(x)][1]), 
            .groups = "drop")

tox_file_max$chem_data <- tox_file_max$chem_data %>%
  filter(type %in% c("Passive", "Water"),
         !is.na(Value)) %>%
  mutate(Value_new = ifelse(files %in% c("Phase2_Pharm.xlsx", 
                                          "Phase2_Pesticides.xlsx"), 
                            Value/1000,
                            Value)) %>%
  group_by(SiteID, CAS) %>%
  summarize(Value = max(Value_new),
            `Sample Start Date` = min(`Sample Date`),
            `Sample End Date` = max(`Sample Date`),
            .groups = "drop") %>%
  mutate(`Sample Date` = Sys.Date())

#Process using toxEval
ACClong <- get_ACC(tox_file_max$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))

cs_synth <- get_chemical_summary(tox_file_max, 
                                 ACClong, 
                                 filtered_ep)

#Check to ensure exclusions are excluded

exclude_test <- filter(cs_synth, CAS == tox_file$exclusions$CAS[4] & 
         endPoint == tox_file$exclusions$endPoint[4] )

if(nrow(exclude_test) > 0){
  warning("Exclusions are not excluded")
}

# Add grouping variables
cs_synth <- tox_file_max$chem_site %>%
  rename(site = SiteID, shortName = `Short Name`) %>% 
  right_join(cs_synth) 


#Quick look at data, does it make sense? 
ggplot(tox_file_max$chem_data[tox_file_max$chem_data$Value>0,], aes(x=Value)) +
  geom_histogram(bins=30) +
  scale_x_log10()


cs_summary <- cs_synth %>%
  filter(EAR > 0.00) %>%
  group_by(chnm, CAS, site) %>%
  summarize(EAR = sum(EAR), .groups = "drop") %>%
  group_by(chnm, CAS) %>%
  summarize(n = n(), n_exceed = length(which(EAR > .1)),
            n_exceed2 = length(which(EAR > 0.01)),
            MaxEAR = max(EAR)) %>%
  arrange(desc(n_exceed), desc(n_exceed2), desc(n), desc(MaxEAR), chnm)  %>%
  filter(n_exceed > 1 | n_exceed2 > 10 | n > 25)

cs_plot <- cs_synth %>%
  filter(EAR > 0) %>%
  filter(chnm %in% cs_summary$chnm) %>%
  mutate(chnm = factor(chnm, rev(cs_summary$chnm))) %>%
  mutate(Class2 = case_when(grepl("icide", Class) ~ "Pesticides",
                            Class == "PAHs" ~ "PAHs",
                            Class == "Solvent" ~ "Other (some are pharms)",
                            Class == "Flavor/Fragrance" ~ "Other (some are pharms)",
                            Class == "Food Additive/Plasticizer" ~ "Other (some are pharms)",
                            Class == "Dye/Pigment" ~ "Other (some are pharms)",
                            Class == "Fire retardant" ~ "Other (some are pharms)",
                            Class == "Antioxidants" ~ "Other (some are pharms)",
                            Class == "Other" ~ "Other (some are pharms)",
                            Class == "Pharmaceuticals" ~ "Pharmaceuticals"))
  

cs_plot$Class2[which(is.na(cs_plot$Class2))] <- "Other (some are pharms)"

cs_plot <- cs_plot %>%
  mutate(Class2 = factor(Class2, c("Pesticides", "Pharmaceuticals",
                                 "PAHs", "Other (some are pharms)")))

Atrazine_test <- cs_plot %>%
  filter(chnm == "Atrazine" & site == "04024000")

Atrazine_raw <- tox_file$chem_data %>%
  filter(CAS == "1912-24-9", SiteID  == "04024000")

summary(data.frame(Atrazine_test))


#Check to ensure exclusions are excluded
exclude_test2 <- filter(cs_plot, CAS == tox_file$exclusions$CAS[4] & 
         endPoint == tox_file$exclusions$endPoint[4] )


if(nrow(exclude_test2) > 0){
  warning("Exclusions are not excluded")
}

synthesis_chnm_stackbar <- plot_stackbar(cs_plot, x = "chnm", stack = "site", 
                                    fill = "EAR", 
                                    facet = "Class2", scales = "free_y", ncol = 1,
                                    breaks = c(0.001, 0.01,  0.1, 1),
                                    palette = c("black", brewer.pal(9, "YlOrRd")[c(9,6,3, 1)])) +
  facet_grid(rows = vars(Class2), scales = "free_y", space="free") +
  # theme(axis.text.x = element_blank()) +
  labs(x = "Chemical name", y = "Number of sites",
       fill = expression(paste("max ",EAR[chem]))) + 
  scale_y_continuous(expand = c(0,0,0,1)) +
  theme(strip.background = element_rect(fill = NA, color = NA), axis.ticks.y = element_blank()) + 
  theme(legend.position = c(0.98, 0.15), legend.justification = "right") +
  coord_flip()

print(synthesis_chnm_stackbar)

ggsave(file.path(path_to_data, "Figures", "Synthesis_chnm_Stackbar.png"),
       synthesis_chnm_stackbar, height = 9, width = 6)


synthesis_site_stackbar <- plot_stackbar(cs_synth, x = "shortName", stack = "chnm", 
                                    fill = "EAR", 
                                    facet = "Lake", scales = "free_y", ncol = 1,
                                    breaks = c(0.001, 0.01,  0.1, 1),
                                    palette = c("black", brewer.pal(9, "YlOrRd")[c(9,6,3, 1)])) +
  facet_grid(rows = vars(Lake), scales = "free_y", space="free") +
  # theme(axis.text.x = element_blank()) +
  labs(x = "Site", y = "Number of chemicals",
       fill = expression(paste("max ",EAR[chem]))) + 
  scale_y_continuous(expand = c(0,0,0,5)) +
  theme(strip.background = element_rect(fill = NA, color = NA), axis.ticks.y = element_blank()) + 
  theme(legend.position = c(0.99, 0.8), legend.justification = "right") + 
  # theme(legend.position = "bottom", legend.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(override.aes = list(size = 0.3))) +
  coord_flip()

print(synthesis_site_stackbar)

ggsave(file.path(path_to_data, "Figures", "Synthesis_site_Stackbar.png"),
       synthesis_site_stackbar, height = 9, width = 5)


# #################
# Mixtures analysis
# #################



# thresholds for mixtures analysis
group_by_this <- "endPoint"
ear_threshold <- 1 #Ear threshold
site_threshold <- length(unique(tox_file_max$chem_data$SiteID))/10 #How many sites


# Calculate EARMix and EARAOP for each sample and endpoint/AOP
EAR_mix <- EAR_mixtures(cs_synth, group_by_this)
AOP_mix <- EAR_mixtures(cs_synth, "AOP")

# Summarize endpoints, chemicals, and genes by site that exceed EAR thresholds
EARsite <- site_mixtures(EAR_mix, ear_threshold)

AOPsite <- site_mixtures(AOP_mix, ear_threshold) %>%
  rename(all_chems_AOP = all_chems,
         all_CASs_AOP = all_CASs)

site_mix_all <- EARsite %>%
  rename(all_chems_endPoint = all_chems,
         all_CASs_endPoint = all_CASs) %>%
  full_join(AOPsite, by = "site")

write.csv(site_mix_all, file.path(path_to_data, "Data", 
                                  "SynthesisSiteMixtures_Table1.csv"), 
          row.names = FALSE)


#Get top endpoint combos
topcombos <- top_mixes(cs_synth, group_by_this, 
                       ear_threshold,
                       site_threshold) 


topmixtures <- overall_mixtures(topcombos, "max")


topaopcombos <- top_mixes(cs_synth, 
                          group_by_what = "AOP",
                          ear_threshold,
                          site_threshold)

topAOPmixtures <- overall_mixtures(topaopcombos, "max")


genetable <- gene_mixtures(topmixtures)  


genefunction <- gene_functions(topmixtures, species = c("Homo sapiens", 
                                                        "Danio rerio", 
                                                        "Xenopus tropicalis"))

geneentrez <- gene_functions(topmixtures, 
                             columns = "entrez", 
                             species = "Homo sapiens") %>%
  filter(Species == "Homo sapiens")

gene_summary <- gene_summary(genefunction)

gene_mixes <- gene_mixtures(topmixtures) %>%
  filter(geneSymbol != "")

gene_out <- full_join(gene_mixes, gene_summary, 
                      by = "geneSymbol")


write.csv(gene_out, file = file.path(path_to_data, 'Data', 
                                     'SynthesisGeneMixtures_Table2.csv'), 
          row.names = F)


gene_wb <- create_Excel_wb_gene(gene_summary, gene_mixes)

saveWorkbook(gene_wb, file = file.path(path_to_data, "Data", 
                                       "Gene_Workbook_synthesis.xlsx"), 
             overwrite = TRUE)

#Change name for plotting
EAR_mix$geneSymbol[grepl("Tanguay_ZF", EAR_mix$endPoint)] <- "zebrafish*"
EAR_mix$geneSymbol[grepl("TOX21_DT40", EAR_mix$endPoint)] <- "cytotoxicity*"
EAR_mix$geneSymbol[EAR_mix$geneSymbol == "NA"] <- NA

#Identify other geneless endpoints and rename manually
odd_endpoints <- unique(EAR_mix$endPoint[is.na(EAR_mix$geneSymbol)])
odd_table <- filter(end_point_info, assay_component_endpoint_name %in% odd_endpoints) %>%
  select(assay_component_endpoint_name, biological_process_target) 


odd_table$biological_process_target[grepl("detection of DNA", 
                                          odd_table$biological_process_target)] <- "DNA detection*"
odd_table$biological_process_target[grepl("regulation", 
                                          odd_table$biological_process_target)] <- "regulation*"
odd_table$biological_process_target[grepl("apoptosis", 
                                          odd_table$biological_process_target)] <- "cell death"
odd_table$biological_process_target[grepl("cell death", 
                                          odd_table$biological_process_target)] <- "cell death*"
odd_table$biological_process_target[grepl("cell cycle", 
                                          odd_table$biological_process_target)] <- "cell cycle*"
odd_table$biological_process_target[grepl("cell proliferation", 
                                          odd_table$biological_process_target)] <- "cell proliferation*"
odd_table$biological_process_target[grepl("protein stablization", 
                                          odd_table$biological_process_target)] <- "other*"
odd_table$biological_process_target[grepl("depolarization", 
                                          odd_table$biological_process_target)] <- "other*"

EAR_mix$geneSymbol[is.na(EAR_mix$geneSymbol)] <- 
  odd_table$biological_process_target[match(EAR_mix$endPoint[is.na(EAR_mix$geneSymbol)], 
                                            odd_table$assay_component_endpoint_name)]


# Plotting
box_gene <- plot_genebar(EAR_mix, ear_threshold, site_threshold, type="geneSymbol", fill = "date") +
  scale_fill_manual(values = c("darkgreen")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none')

box_endpoint <- plot_genebar(EAR_mix, ear_threshold, site_threshold, type="endPoint", fill = 'date') +
  scale_fill_manual(values = c("darkred")) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none')

box_AOP <- plot_genebar(AOP_mix, ear_threshold, site_threshold, type="AOP", fill = 'date') +
  scale_fill_manual(values = c("darkblue")) + 
  scale_x_log10nice(name = expression(paste(EAR[AOP]))) +
  theme(legend.position = 'none')

png(file.path(path_to_data, "Figures", "Synthesis_PriorityBoxplot_3panel_v2.png"), height = 8, width = 12, units = "in", res = 300)
grid.newpage()
boxes <- grid.draw(cbind(ggplotGrob(box_endpoint), 
                              ggplotGrob(box_gene), 
                              ggplotGrob(box_AOP)))
dev.off()


#Subset to only look at genes, endpoints, and AOPs associated with Estrogen and Androgen receptors

EAR_mix_sub <- EAR_mix %>%
  filter(geneSymbol %in% c("AR", "ESR1", "ESR2")) %>%
  mutate(group = case_when(geneSymbol == "AR" ~ "AR",
                           grepl("ESR", geneSymbol) ~ "ESR")) %>%
  mutate(group = factor(group, c("AR", "ESR")))

AOP_mix_sub <- cs_synth %>%
  filter(endPoint %in% EAR_mix_sub$endPoint) %>%
  EAR_mixtures("AOP") %>%
  mutate(group = case_when(grepl("Androgen", AOP_title, ignore.case = TRUE) ~ "AR",
                           grepl("Estrogen", AOP_title, ignore.case = TRUE) ~ "ESR",
                           grepl("Estradiaol", AOP_title) ~ "ESR",
                           grepl("testosterone", AOP_title) ~ "AR",
                           grepl("ER agonism", AOP_title) ~ "ESR")) %>%
  mutate(group = factor(group, c("AR", "ESR")))



#plotting
box_gene_sub <- plot_genebar(EAR_mix_sub, ear_cutoff = 0.1, site_threshold, type="geneSymbol", fill = "group") +
  scale_fill_manual(values = c("darkblue", "darkred", "grey")) +
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none', axis.title.y = element_blank()) +
  ggtitle("Gene symbol")

box_endpoint_sub <- plot_genebar(EAR_mix_sub, ear_cutoff = 0.1, site_threshold, type="endPoint", fill = 'group') +
  scale_fill_manual(values = c("darkblue", "darkred", "grey")) +
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  theme(legend.position = 'none', axis.title.y = element_blank()) +
  ggtitle("ToxCast Assay")

box_AOP_sub <- plot_genebar(AOP_mix_sub, ear_cutoff = 0.1, site_threshold, type="AOP", fill = 'group') +
  scale_fill_manual(values = c("darkblue", "darkred", "grey")) +
  scale_x_log10nice(name = expression(paste(EAR[AOP]))) +
  theme(legend.position = 'none', axis.title.y = element_blank()) +
  ggtitle("AOP-wiki#")

png(file.path(path_to_data, "Figures", "Synthesis_AR_PriorityBoxplot_3panel_v2.png"), height = 4, width = 9, units = "in", res = 300)
grid.newpage()
boxes <- grid.draw(cbind(ggplotGrob(box_endpoint_sub), 
                         ggplotGrob(box_gene_sub), 
                         ggplotGrob(box_AOP_sub)))
dev.off()

