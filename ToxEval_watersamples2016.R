
#Load all toxEval files and generate chemical summary for each

############################################
# Select data 
# include all chemicals that have known POCIS kinetics

tox_list_surface <- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides.xlsx")))
# tox_list_surface$chem_site$site_grouping <- factor(tox_list_surface$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))
ACClong <- get_ACC(tox_list_surface$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))

chemicalSummary_surface <- get_chemical_summary(tox_list_surface, 
                                        ACClong, 
                                        filtered_ep)

chemicalSummary_surface <- tox_list_surface$chem_site %>%
  rename(site = SiteID,
         shortName = `Short Name`) %>% 
  right_join(chemicalSummary_surface) 

#Add chemical class for chemicals that have new names from toxEval
# chemicalSummary_surface$Class[is.na(chemicalSummary_surface$Class)] <- tox_list_surface$chem_info$Class[match(chemicalSummary_surface$CAS[is.na(chemicalSummary_surface$Class)], tox_list_surface$chem_info$CAS)]



##################################
# custom benchmarks
# same file as select data, but use custom benchmarks
# use this for toxicity quotient (TQ)

tox_bench_list_surface<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides_Bench.xlsx")))
# tox_bench_list_surface$chem_site$site_grouping <- factor(tox_bench_list_surface$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))


chemicalSummary_bench_surface <- get_chemical_summary(tox_bench_list_surface)

#Error is happening here
chemicalSummary_bench_surface <- tox_bench_list_surface$chem_site %>%
  rename(site = SiteID,
         shortName = `Short Name`) %>% 
  right_join(chemicalSummary_bench_surface) 

#Add chemical class for chemicals that have new names from toxEval
# chemicalSummary_bench$Class[is.na(chemicalSummary_bench$Class)] <- tox_bench_list_surface$chem_info$Class[match(chemicalSummary_bench$CAS[is.na(chemicalSummary_bench$Class)], tox_bench_list_surface$chem_info$CAS)]

