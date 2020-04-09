

pocisPrep <- chem_freq_allpocis %>%
  rename(chnm_passive = chnm,
         detect_passive = Detected)


waterPrep <- chem_freq_all %>%
  rename(chnm_water = chnm,
         detect_water = LabDetected)


chem_combine <- full_join(pocisPrep, waterPrep) %>%
  select(CAS, chnm_passive,  chnm_water, detect_passive, detect_water, Class) %>%
  mutate(both_detect = case_when(!is.na(detect_passive) & !is.na(detect_water) ~ 'Both',
                                 !is.na(detect_passive) & is.na(detect_water) ~ 'Passive',
                                 is.na(detect_passive) & !is.na(detect_water) ~ 'Water')) %>%
  arrange(both_detect, chnm_water, chnm_passive)
  

data.frame(chem_combine)

 
# chem_combine$chnm_water

