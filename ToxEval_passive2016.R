
#Load all toxEval files and generate chemical summary for each

############################################
# Select data 
# include all chemicals that have known POCIS kinetics

tox_list<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016_ToxEval.xlsx")))
tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical","ACCLessThan"))

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

#Site order for plotting
site_order <- unique(chemicalSummary$shortName)
site_ID_order <- unique(chemicalSummary$site)


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




#Create csv for publication
#Chemicals included in analysis
chemicals_in <- tox_input_list$Data %>%
  group_by(CAS, chnm) %>%
  select(chnm, CAS, Value) %>%
  filter(Value > 0) %>%
  tally(name = 'Detections') %>%
  full_join(tox_input_list$Data[,c('chnm', 'CAS')]) %>%
  distinct() %>%
  arrange(chnm) %>%
  left_join(tox_input_list$Chemicals[,c('CAS', 'Class')]) %>%
  filter(Class !="")

# chemicals_in$Detections[is.na(chemicals_in$Detections)] <- 0

chemicals_out <- chemicalSummary %>%
  group_by(Class, chnm, CAS, site) %>%
  select(Class, chnm, CAS, EAR, site) %>%
  summarize(EAR = max(EAR, na.rm=T)) %>%
  filter(EAR > 0) %>%
  tally(name = 'Detections') %>%
  full_join(unique(chemicalSummary[,c('chnm', 'CAS', 'Class')])) %>%
  distinct() %>%
  arrange(Class, desc(Detections), as.character(chnm)) %>%
  filter(Class !="")

# chemicals_out$Detections[is.na(chemicals_out$Detections)] <- 0

missing_ToxCast <- setdiff(chemicals_in$CAS, chemicals_out$CAS)

chemicals_in[chemicals_in$CAS %in% missing_ToxCast,]

chemicals_combine <- chemicals_in %>%
  filter(CAS %in% missing_ToxCast) %>%
  bind_rows(chemicals_out) %>%
  mutate(Class = factor(Class, c('Herbicide', 'Deg - Herbicide', 'Fungicide', 'Deg - Fungicide', 'Insecticide', 'Deg - Insecticide'))) %>%
  select(Class, chnm, CAS, Detections) %>%
  arrange(Class, desc(Detections), chnm) 


#Find number of EAR 'hits'
chemicalSummary_EARHits <- chemicalSummary %>%
  filter(EAR>=.001) %>%
  dplyr::group_by(chnm, site) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  tally(name = "EAR_Hits") %>%
  group_by() %>%
  mutate(chnm = as.character(chnm))

#Find number of TQ 'hits'
chemicalSummary_TQHits <- chemicalSummary_bench %>%
  filter(EAR>=.1) %>%
  dplyr::group_by(chnm, site) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>%
  tally(name = "TQ_Hits") %>%
  group_by() %>%
  mutate(chnm = as.character(chnm))

missing_TQ <- setdiff(chemicals_combine$CAS, chemicalSummary_bench$CAS)

chemicals_combine[chemicals_combine$CAS %in% missing_TQ,]


chemicals_combine<- full_join(chemicals_combine, chemicalSummary_EARHits) %>%
  full_join(chemicalSummary_TQHits) %>%
  arrange(Class, desc(Detections), desc(EAR_Hits), desc(TQ_Hits)) %>%
  rename('EAR hits' = EAR_Hits,
         'TQ hits' = TQ_Hits,
         'Chemical Name' = chnm)

chemicals_combine$`EAR hits`[chemicals_combine$CAS %in% missing_ToxCast] <- "*"
chemicals_combine$`TQ hits`[chemicals_combine$CAS %in% missing_TQ] <- "*"

chemicals_combine$`EAR hits`[is.na(chemicals_combine$`EAR hits`)] <- ""
chemicals_combine$`TQ hits`[is.na(chemicals_combine$`TQ hits`)] <- ""

write.csv(chemicals_combine, file = file_out(file.path(path_to_data, "Data/Chemicals_characterisitcs.csv")), row.names=F)

# write.table(chemicals_combine, file = file_out(file.path(path_to_data, "Data/Chemicals_characterisitcs.txt")), row.names=F, sep=',')




rm(ACClong, ACClong_allpocis, cleaned_ep, filtered_ep )
