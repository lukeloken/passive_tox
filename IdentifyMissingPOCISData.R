library(dataRetrieval)

Pcodes_POCIS <- readNWISpCode(51404:51635)

# Pcodes_POCIS <- filter(Pcodes, grepl("POCIS", parameter_units)) %>%
#   rename(PARMCODE = parameter_cd) %>%
#   mutate(PARMCODE = as.character(PARMCODE))

sample_ids <- read_excel(file.path(path_to_data, "RawData", "PassiveSampleIDs.xlsx"))

full_data <- read_excel(file.path(path_to_data, "RawData", "CSR14408_CSQA_Alvarez.xlsx"))

#From data_reader script
AllPOCIS_2016
cas_df
Alvarez_pocis_orig <- tox_list_allpocis$chem_data %>%
  left_join(Pcodes_POCIS, by = c("CAS" = "casrn"))

GLRI_data <- full_data %>%
  filter(LABID %in% c(sample_ids$`Lab ID`))

unique(GLRI_data$REMARKCODE)
unique(GLRI_data$QUALIFIER)
unique(GLRI_data$PARMCODE)

POCIS_data <- filter(GLRI_data, PARMCODE %in% 51404:51635) %>%
  mutate(parameter_cd = as.character(PARMCODE)) %>%
  left_join(Pcodes_POCIS, by = ("parameter_cd"))

unique(POCIS_data$REMARKCODE)
unique(POCIS_data$QUALIFIER)
unique(POCIS_data$PARMCODE)

unique(Alvarez_pocis_orig$CAS)

CAS_inPcodes <- setdiff(unique(POCIS_data$casrn), unique(Alvarez_pocis_orig$CAS))
CAS_inAlvarez <- setdiff(unique(Alvarez_pocis_orig$CAS), unique(POCIS_data$casrn))

data.frame(filter(Alvarez_pocis_orig, CAS %in% CAS_inAlvarez))

data.frame(filter(POCIS_data, casrn %in% CAS_inPcodes))


missing_data <- filter(Pcodes_POCIS, casrn %in% CAS_inPcodes) 

missing_data$includedpreviously <- "No"

missing_data2 <- full_join(missing_data, Pcodes_POCIS)

missing_data2$includedpreviously[is.na(missing_data2$includedpreviously)] <- "Yes"

missing_data2$chnm <- unlist(lapply(str_split(missing_data2$parameter_nm, ", "), function(l) l[1]))
missing_data2$casrn <- paste0("'", missing_data2$casrn)

  
write.csv(missing_data2, file.path(path_to_data, "Data", "MissingPOCISData.csv"), row.names = F)
