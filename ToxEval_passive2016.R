
#Load all toxEval files and generate chemical summary for each

############################################
# Select data 
# include all chemicals that have known POCIS kinetics

tox_list<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016_ToxEval.xlsx")))
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

#Add chemical class for chemicals that have new names from toxEval
chemicalSummary$Class[is.na(chemicalSummary$Class)] <- tox_list$chem_info$Class[match(chemicalSummary$CAS[is.na(chemicalSummary$Class)], tox_list$chem_info$CAS)]


###############################################################
# all pocis data
# These data were above the minimum detection in the POCIS
# EAR is irrelevant, but used for presence/absence

tox_list_allpocis<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016AllPOCIS_ToxEval.xlsx")))
tox_list_allpocis$chem_site$site_grouping <- factor(tox_list_allpocis$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))
ACClong_allpocis <- get_ACC(tox_list_allpocis$chem_info$CAS)
ACClong_allpocis <- remove_flags(ACClong_allpocis)

chemicalSummary_allpocis <- get_chemical_summary(tox_list_allpocis, 
                                        ACClong, 
                                        filtered_ep)

chemicalSummary_allpocis <- tox_list_allpocis$chem_site %>%
  rename(site = SiteID,
         shortName = `Short Name`) %>% 
  right_join(chemicalSummary_allpocis)

chemicalSummary_allpocis

#Add chemical class for chemicals that have new names from toxEval
chemicalSummary_allpocis$Class[is.na(chemicalSummary_allpocis$Class)] <- tox_list$chem_info$Class[match(chemicalSummary_allpocis$CAS[is.na(chemicalSummary_allpocis$Class)], tox_list$chem_info$CAS)]



##################################
# custom benchmarks
# same file as select data, but use custom benchmarks
# use this for toxicity quotient (TQ)

tox_bench_list<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016Benchmarks_ToxEval.xlsx")))
tox_bench_list$chem_site$site_grouping <- factor(tox_bench_list$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))


chemicalSummary_bench <- get_chemical_summary(tox_bench_list)

chemicalSummary_bench <- tox_bench_list$chem_site %>%
  rename(site = SiteID,
         shortName = `Short Name`) %>% 
  right_join(chemicalSummary_bench) 

#Add chemical class for chemicals that have new names from toxEval
chemicalSummary_bench$Class[is.na(chemicalSummary_bench$Class)] <- tox_bench_list$chem_info$Class[match(chemicalSummary_bench$CAS[is.na(chemicalSummary_bench$Class)], tox_bench_list$chem_info$CAS)]


rm(ACClong, ACClong_allpocis, cleaned_ep, filtered_ep )
