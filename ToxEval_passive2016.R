
#Load all toxEval files and generate chemical summary for each

############################################
# Select data 
# include all chemicals that have known POCIS kinetics

tox_list<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016_ToxEval.xlsx")))
tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))

#load surface exclusions to combine with passive exclusions
tox_list_surface <- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides.xlsx")))


#Combine exclusions
exclusions <- tox_list$exclusions %>%
  # dplyr::select(-chnm) %>%
  bind_rows(tox_list_surface$exclusions) %>%
  distinct()

tox_list$exclusions <- exclusions


ACClong <- get_ACC(tox_list$chem_info$CAS)
# ACClong <- remove_flags(ACClong, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical"))
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

#Site order for plotting
site_order <- unique(chemicalSummary$shortName)
site_ID_order <- unique(chemicalSummary$site)

conc_table <- tox_list$chem_data %>%
  filter(Value > 0)

summary(conc_table) 

missingTox <- conc_table %>%
  filter(!CAS %in% chemicalSummary$CAS) %>%
  group_by(chnm) %>%
  summarize(min = min(Value),
            mean = mean(Value), 
            max = max(Value), 
            n = n()) %>%
  arrange(desc(n))

summary(missingTox)

ggplot(tox_list$chem_data[tox_list$chem_data$Value>0,], aes(x=Value)) +
  geom_histogram(bins=30) +
  scale_x_log10()

###############################################################
# all pocis data
# These data were above the minimum detection in the POCIS
# EAR is irrelevant, but used for presence/absence

tox_list_allpocis<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016AllPOCIS_ToxEval.xlsx")))
tox_list_allpocis$chem_site$site_grouping <- factor(tox_list_allpocis$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))
ACClong_allpocis <- get_ACC(tox_list_allpocis$chem_info$CAS)
# ACClong_allpocis <- remove_flags(ACClong_allpocis, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical","ACCLessThan"))
ACClong_allpocis <- remove_flags(ACClong_allpocis)


cleaned_ep_allpocis <- clean_endPoint_info(end_point_info)
filtered_ep_allpocis <- filter_groups(cleaned_ep_allpocis, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))


chemicalSummary_allpocis <- get_chemical_summary(tox_list_allpocis, 
                                                 ACClong_allpocis, 
                                                 filtered_ep_allpocis)

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



#Figure out how many chemicals are in each classification
# Conc versus no conc
# Detect versus no detect
# EAR available yes/no

#Create csv for publication

#Chemicals detected in pocis
chemicals_allpocis <- tox_list_allpocis$chem_data %>%
  group_by(CAS, chnm) %>%
  select(chnm, CAS, Value) %>%
  filter(Value > 0) %>%
  tally(name = 'POCIS_Detect') %>%
  left_join(unique(tox_list_allpocis$chem_data[,c('chnm', 'CAS')])) %>%
  distinct() %>%
  arrange(chnm) %>%
  left_join(unique(tox_list_allpocis$chem_info[,c('CAS', 'Class')])) %>%
  # filter(Class !="") %>%
  arrange(Class, desc(POCIS_Detect))

chemicals_notdetected <- tox_list_allpocis$chem_data %>%
  group_by(CAS, chnm) %>%
  select(chnm, CAS, Value) %>%
  summarize(maxValue = max(Value, na.rm=T),
            n=n())  %>%
  filter(maxValue == 0)

pocischemsintoxcast <- chemicalSummary_allpocis %>%
  filter(EAR>0) %>%
  select(chnm, CAS, site) %>%
  distinct() %>%
  select(-site) %>%
  group_by(chnm, CAS) %>%
  summarize(n = n())
  
pocischemsinbenchmarks <- tox_list_allpocis$chem_data %>%
  group_by(CAS, chnm) %>%
  select(chnm, CAS, Value) %>%
  summarize(maxValue = max(Value, na.rm=T),
            n=n())  %>%
  filter(maxValue > 0) %>%
  filter(CAS %in% tox_bench_list$benchmarks$CAS)

#Chemicals included in concentration analysis
chemicals_in <- tox_list$chem_data %>%
  group_by(CAS, chnm) %>%
  select(chnm, CAS, Value) %>%
  filter(Value > 0) %>%
  tally(name = 'Detections') %>%
  left_join(unique(tox_list$chem_data[,c('chnm', 'CAS')])) %>%
  distinct() %>%
  arrange(chnm) %>%
  left_join(unique(tox_list$chem_info[,c('CAS', 'Class')])) %>%
  # filter(Class !="") %>%
  filter(CAS %in% chemicals_allpocis$CAS)

#Chemicals included in concentration analysis
#Note that Tebupirimos had a concentration for one site, yet was below detection limit for all POCIS data
chemicals_zeroconc <- tox_list$chem_data %>%
  group_by(CAS, chnm) %>%
  select(chnm, CAS, Value) %>%
  summarize(maxValue = max(Value, na.rm=T),
            n=n())  %>%
  filter(CAS %in% chemicals_notdetected$CAS & 
           chnm %in% chemicals_notdetected$chnm)

# chemicals_in$Detections[is.na(chemicals_in$Detections)] <- 0

chemicals_out <- chemicalSummary %>%
  group_by(Class, chnm, CAS, site) %>%
  select(Class, chnm, CAS, EAR, site) %>%
  summarize(EAR = max(EAR, na.rm=T)) %>%
  filter(EAR > 0) %>%
  tally(name = 'Detections') %>%
  full_join(unique(chemicalSummary[,c('chnm', 'CAS', 'Class')])) %>%
  distinct() %>%
  arrange(Class, desc(Detections), as.character(chnm)) 
  # filter(Class !="")

# chemicals_out$Detections[is.na(chemicals_out$Detections)] <- 0

missing_ToxCast <- setdiff(chemicals_in$CAS, chemicals_out$CAS)

chemicals_in[chemicals_in$CAS %in% missing_ToxCast,]

chemicals_combine <-   chemicals_in %>%
  full_join(chemicals_allpocis) %>%
  arrange(chnm) %>%
  # filter(CAS %in% missing_ToxCast) %>%
  # bind_rows(chemicals_out) %>%
  mutate(Class = factor(Class, c('Herbicide', 'Deg - Herbicide', 'Fungicide', 'Deg - Fungicide', 'Insecticide', 'Deg - Insecticide'))) %>%
  select(Class, chnm, CAS, Detections, POCIS_Detect) %>%
  arrange(Class, desc(Detections), chnm) 

chemicals_out2 <- rename(chemicals_out, chnm_ne = chnm)

chemicals_combine_rename <- chemicals_combine[which(chemicals_combine$CAS %in% chemicals_out2$CAS),] %>%
  left_join(chemicals_out2) %>%
  select(-chnm) %>%
  rename(chnm = chnm_ne) %>%
  full_join(chemicals_combine[-which(chemicals_combine$CAS %in% chemicals_out2$CAS),])

chemicals_combine <- chemicals_combine_rename %>%
  mutate(Class = factor(Class, c('Herbicide', 'Deg - Herbicide', 'Fungicide', 'Deg - Fungicide', 'Insecticide', 'Deg - Insecticide')))
  

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

missing_TQ <- setdiff(chemicals_in$CAS, unique(chemicalSummary_bench$CAS))

chemicals_combine[chemicals_combine$CAS %in% missing_TQ,]


chemicals_combine<- full_join(chemicals_combine, chemicalSummary_EARHits) %>%
  full_join(chemicalSummary_TQHits) %>%
  select(Class, CAS, chnm, POCIS_Detect, Detections, EAR_Hits, TQ_Hits) %>%
  arrange(Class, desc(Detections), desc(EAR_Hits), desc(TQ_Hits), desc(POCIS_Detect)) %>%
  rename('EAR hits' = EAR_Hits,
         'TQ hits' = TQ_Hits,
         'Chemical Name' = chnm)

chemicals_combine$`EAR hits`[chemicals_combine$CAS %in% missing_ToxCast] <- "*"
chemicals_combine$`TQ hits`[chemicals_combine$CAS %in% missing_TQ] <- "*"

chemicals_combine$`EAR hits`[is.na(chemicals_combine$`EAR hits`)] <- ""
chemicals_combine$`TQ hits`[is.na(chemicals_combine$`TQ hits`)] <- ""

chemicals_combine$Detections[which(chemicals_combine$Detections>0)] <- "Yes"
chemicals_combine$Detections[is.na(chemicals_combine$Detections)] <- "No"

chemicals_combine$`ToxCast Available` <- 'No'
chemicals_combine$`ToxCast Available`[chemicals_combine$CAS %in% pocischemsintoxcast$CAS] <- "Yes"

chemicals_combine$`Benchmark Available` <- 'No'
chemicals_combine$`Benchmark Available`[chemicals_combine$CAS %in% pocischemsinbenchmarks$CAS] <- "Yes"

chemicals_combine <- chemicals_combine %>%
  rename(`Concentration Available` = Detections) %>%
  rename(Detections = POCIS_Detect)

write.csv(chemicals_combine, file = file_out(file.path(path_to_data, "Data/Chemicals_characterisitcs.csv")), row.names=F)

# write.table(chemicals_combine, file = file_out(file.path(path_to_data, "Data/Chemicals_characterisitcs.txt")), row.names=F, sep=',')


chemicals_notavailable <- chemicals_combine %>%
  select(CAS, `Chemical Name`, `Concentration Available`, `ToxCast Available`, `Benchmark Available`,  Detections) %>%
  arrange(desc(`Concentration Available`), desc(`ToxCast Available`), desc(`Benchmark Available`), desc(Detections)) 

write.csv(chemicals_notavailable, file = file_out(file.path(path_to_data, "Data/Chemicals Lacking Evaluations.csv")), row.names=F)


rm(ACClong_allpocis, cleaned_ep, filtered_ep )

