library(drake)
library(tidyverse)
library(googledrive)
library(readxl)
library(data.table)
library(toxEval)
library(openxlsx)
library(dataRetrieval)
library(dplyr)

# path_to_data <- Sys.getenv("PASSIVE_PATH")

options(drake_make_menu = FALSE)

dir.create("data", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create(file.path("data","raw"), showWarnings = FALSE)
dir.create(file.path("data","clean"), showWarnings = FALSE)

# Go from raw files to R objects:
source("R/analyze/data_reader_2016_Loken.R")
# source("R/analyze/get_sites_ready.R")
source("R/analyze/get_chem_info.R")
# source("R/analyze/create_tox_file.R")

#From toxeval tutorial
path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)
tox_list_example <- create_toxEval(full_path)

#Load chemical lists
cas_df <- read.csv(file.path(path_to_data, 
                             "Data", 
                             "POCIS_pesticides_samplingrates_AllPcodes.csv")) %>%
  rename(CAS = casrn) # %>%
  # mutate(chnm = tolower(chnm)) 

cas_df$CAS <- gsub("'", '', cas_df$CAS)
cas_df$CAS[which(cas_df$UpdatedCAS != "")] <- cas_df$UpdatedCAS[which(cas_df$UpdatedCAS != "")]

cas_df <- cas_df %>%
  select(parameter_cd, parameter_nm, CAS, chnm, POCIS_sampling_rate, Class) %>%
  mutate(parameter_cd = as.character(parameter_cd), 
         POCIS_sampling_rate = as.numeric(POCIS_sampling_rate))

head(cas_df)


#Benchmarks and chemicals from Sam's paper
#from Sam's surface water tox eval file
benchmarks_df <- read_excel(file_in(file.path(path_to_data, 'Data/WQ_pesticides_Bench.xlsx')), sheet='Benchmarks') 
cas_df_surf <- read_excel(file_in(file.path(path_to_data, 'Data/WQ_pesticides.xlsx')), sheet='Chemicals') 

# Load the sites
sites_2016 = readxl::read_excel(file_in(file.path(path_to_data, "RawData/GLRI 2016 POCIS pesticide data report.xlsx")), sheet = "Site List", skip = 2) %>%
  rename(site_no = `USGS Station ID`,
         Date_in = Date...6,
         Time_in = Time...7,
         Date_out = Date...8,
         Time_out = Time...9) %>%
  mutate(Time_in = str_pad(Time_in, 4, pad = "0"),
         Time_out = str_pad(Time_out, 4, pad = "0"),
         Date_in = as.Date(Date_in),
         Date_out = as.Date(Date_out)) %>%
  mutate(Datetime_in = as.POSIXct(paste(Date_in, Time_in), format="%Y-%m-%d %H%M"),
         Datetime_out = as.POSIXct(paste(Date_out, Time_out), format="%Y-%m-%d %H%M")) %>%
  # mutate(`Short Name` = substr(`Site Name`, start = 1, stop = 8)) %>%
  select(`Short Name`, `site_no`, Lake, Date_in, Date_out, Datetime_in, Datetime_out, Deployed)

locations <- readNWISsite(sites_2016$site_no)[c("site_no", "dec_lat_va", "dec_long_va", 'station_nm')]
names(locations) <- c("site_no", 'dec_lat', 'dec_long', 'site_grouping')
locations$site_grouping <- str_sub(locations$site_grouping,-2,-1)


#Data from lab
full_data <- read_excel(file.path(path_to_data, "RawData", "CSR14408_CSQA_Alvarez.xlsx"))
sample_ids <- read_excel(file.path(path_to_data, "RawData", "PassiveSampleIDs.xlsx")) %>%
  mutate(LABID = as.character(`Lab ID`),
         SiteID = paste0("0", `Station ID`)) %>%
  select(-`Station ID`, -Received)

blank_data <- full_data %>% 
  filter(PARMCODE %in% 51404:51635, 
         LABID %in% 20180370102:20180370104)

blank_data$Value <- NA
blank_data$comment <- ""
blank_data$comment[grep("<",blank_data$FINAL)] <- "<"
blank_data$comment[grep("E",blank_data$FINAL)] <- "Est"
blank_data$comment[grep("DELETED",blank_data$FINAL)] <- "Deleted"

blank_data$Value_ngPOCIS <- blank_data$FINAL
blank_data$Value_ngPOCIS[grepl("<", blank_data$Value_ngPOCIS)] <- 0
blank_data$Value_ngPOCIS[grepl("U-DELETED", blank_data$Value_ngPOCIS)] <- NA
blank_data$Value_ngPOCIS[grepl("E", blank_data$Value_ngPOCIS)] <- 
  gsub("E","", blank_data$Value_ngPOCIS[grepl("E", blank_data$Value_ngPOCIS)])
blank_data$Value_ngPOCIS <- gsub(" ","",blank_data$Value_ngPOCIS)

blank_data$Value_ngPOCIS  <- as.numeric(blank_data$Value_ngPOCIS)

blank_hits <- filter(blank_data, Value_ngPOCIS>0) %>%
  select(LABID, PARMCODE, Value_ngPOCIS, FINAL) 

blank_hits$LABID_water = ifelse(blank_hits$LABID == 20180370103, 
                                "20180370120", 
                                ifelse(blank_hits$LABID == 20180370104,
                                      NA, NA)) 
                                  
                                  
sites_2016 <- left_join(sites_2016, locations) %>%
  rename(SiteID = site_no) %>%
  filter(SiteID %in% unique(sample_ids$SiteID)) %>%
  distinct() %>%
  group_by(`Short Name`, SiteID, Lake, site_grouping) %>%
  summarize_all(.funs=mean)

fake_sites <- sites_2016[1,]

GLRI_data <- left_join(sample_ids, full_data, by = "LABID") %>%
  left_join(cas_df, by = c("PARMCODE" = "parameter_cd")) %>%
  filter(PARMCODE %in% 51404:51635)

# data.frame(filter(GLRI_data, PARMCODE == "51409"))
# data.frame(filter(GLRI_data, chnm == "Propoxur"))

POCIS_MDLs <- GLRI_data %>%
  group_by(PARMCODE, CAS, chnm) %>%
  summarize(MDL_ngPOCIS = min(as.numeric(gsub("<", "", unique(FINAL)[grepl("<", unique(FINAL))])), 
                              na.rm = TRUE)) %>% 
  ungroup() %>%
  left_join(cas_df) 

POCIS_MDLs[which(!is.finite(POCIS_MDLs$MDL_ngPOCIS)),]

#Replace two MDLs with data
POCIS_MDLs$MDL_ngPOCIS[which(POCIS_MDLs$PARMCODE == 51554)] <- 0.4
POCIS_MDLs$MDL_ngPOCIS[which(POCIS_MDLs$PARMCODE == 51528)] <- 7.1

POCIS_MDLs <- POCIS_MDLs %>%
  mutate(MDL_ngL_30day = MDL_ngPOCIS / POCIS_sampling_rate / 30)

data.frame(POCIS_MDLs)
POCIS_MDLs[which(!is.finite(POCIS_MDLs$MDL_ngPOCIS)),]

GLRI_data$comment <- ""
GLRI_data$comment[grep("<",GLRI_data$FINAL)] <- "<"
GLRI_data$comment[grep("E",GLRI_data$FINAL)] <- "Est"
GLRI_data$comment[grep("DELETED",GLRI_data$FINAL)] <- "Deleted"

GLRI_data$Value_ngPOCIS <- GLRI_data$FINAL
GLRI_data$Value_ngPOCIS[grepl("<", GLRI_data$Value_ngPOCIS)] <- 0
GLRI_data$Value_ngPOCIS[grepl("U-DELETED", GLRI_data$Value_ngPOCIS)] <- NA
GLRI_data$Value_ngPOCIS[grepl("E", GLRI_data$Value_ngPOCIS)] <- gsub("E","",GLRI_data$Value_ngPOCIS[grepl("E", GLRI_data$Value_ngPOCIS)])
GLRI_data$Value_ngPOCIS <- gsub(" ","",GLRI_data$Value_ngPOCIS)

GLRI_data$Value_ngPOCIS  <- as.numeric(GLRI_data$Value_ngPOCIS)

# data.frame(GLRI_data[which(is.na(GLRI_data$Value)),])

GLRI_data <- full_join(GLRI_data, select(POCIS_MDLs, PARMCODE, MDL_ngPOCIS, MDL_ngL_30day))


Site_data <- GLRI_data %>%
  filter(LABID %in% 20180370105:20180370120) %>%
  left_join(sites_2016, by = "SiteID") %>%
  mutate(Deploy_length_d = as.numeric(difftime(Datetime_out, Datetime_in, units = "days"))) %>%
  mutate(Value = 0.001 *signif(Value_ngPOCIS / (POCIS_sampling_rate * Deploy_length_d), 3),
         `Sample Date` = year(Date_in)) %>%
  arrange(`Sample Name`, SiteID, PARMCODE) %>%
  select(all_of(names(sites_2016)), everything())
   
Site_data$Value[which(Site_data$Value_ngPOCIS == 0)] <- 0

head(data.frame(Site_data))
unique(Site_data$PARMCODE)
# unique(GLRI_data$REMARKCODE)
# unique(GLRI_data$QUALIFIER)
# unique(GLRI_data$PARMCODE)

blank_comparison <- Site_data %>%
  filter(PARMCODE %in% blank_hits$PARMCODE) %>%
  group_by(PARMCODE, chnm) %>%
  summarize(Value_mean_POCIS = mean(Value_ngPOCIS), 
            Value_min_POCIS = min(Value_ngPOCIS), 
            Value_min_noNAPOCIS = min(Value_ngPOCIS, na.rm = TRUE), 
            Value_median_POCIS = median(Value_ngPOCIS), 
  ) %>%
  right_join(blank_hits) %>%
  mutate(percent_diff = 100 * Value_ngPOCIS / Value_mean_POCIS)

100*blank_comparison$Value_ngPOCIS / blank_comparison$Value_median_POCIS
100*blank_comparison$Value_ngPOCIS / blank_comparison$Value_min_POCIS


#Load pesticide data
# WW_2016 = generic_file_opener(file_in(file.path(path_to_data, "RawData/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
#                               sheet = "select data",
#                               n_max=60,
#                               site_sheet = "Site List",
#                               year = 2016,
#                               skip_site = 2)
# 
# WW_2016$Value[which(WW_2016$chnm == "Tebupirimfos")] <- 0

# AllPOCIS_2016 = generic_file_opener(file_in(file.path(path_to_data, "Data/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
                              # sheet = "all data",
                              # n_max=300,
                              # site_sheet = "Site List",
                              # year = 2016,
                              # skip_site = 2)

#testing extra CAS files
# AllPOCIS_2016 = generic_file_opener(file_in(file.path(path_to_data, "RawData/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
#                                     sheet = "all data",
#                                     n_max=300,
#                                     site_sheet = "Site List",
#                                     year = 2016,
#                                     skip_site = 2)



# Calculate new pocis rates based on deployment length
# 
# POCIS_uptake_rates <- read_excel(file_in(file.path(path_to_data, 
#                                                    "Data",
#                                                    "additional GLRI pesticide Rs values.xlsx")), 
#                                  skip = 1)
# names(POCIS_uptake_rates)[1:2] <- c("chnm", "Rs")
# 
# POCIS_uptake_rates <- POCIS_uptake_rates %>%
#   select(chnm, Rs) %>%
#   mutate(Rs = as.numeric(Rs)) %>%
#   filter(!is.na(Rs)) %>%
#   mutate(chnm2 = tolower(chnm)) %>%
#   select(-chnm)
# 
# POCIS_uptake_rates$chnm2[grepl("desamino metribuzin", POCIS_uptake_rates$chnm2)] <- "metribuzin da"

# Cw = N / (Rs * t) 
# Cw (Time weighted concentration (ng/L)
# N pocis concentration (ng/pocis)
# Rs sampling rate (L/D)
# t time (d)

# new_POCIS <- AllPOCIS_2016 %>%
#   mutate(chnm2 = tolower(chnm),
#          MDL = NA) %>%
#   left_join(POCIS_uptake_rates, by = "chnm2") %>%
#   filter(!is.na(Rs)) %>%
#   left_join(select(ungroup(sites_2016), SiteID, Datetime_in, Datetime_out), by = "SiteID") %>%
#   mutate(Deploy_length_d = as.numeric(difftime(Datetime_out, Datetime_in, units = "days"))) %>%
#   mutate(Value = 0.001 *signif(Value / (Rs * Deploy_length_d), 3),
#          `Sample Date` = year(`Date Deployed`)) %>%
#   select(all_of(names(WW_2016))) %>%
#   arrange(chnm)
# 
# new_POCIS_test <- filter(new_POCIS, Value > 0)
# unique(new_POCIS_test$chnm)
# 
# #Merge with other time weighted concentrations
# WW_2016 <- WW_2016 %>%
#   filter(chnm != "Hydroxysimazine") %>%
#   full_join(new_POCIS)
# 
# 



# average the two replicates for Saginow River Station
# WW_2016_forToxEval <- WW_2016 %>%
#   select(CAS, chnm, SiteID, Value, `Sample Date`) %>%
#   group_by(SiteID,  chnm, CAS) %>%
#   summarize(Value = mean(Value), `Sample Date` = min(`Sample Date`))
# 
# WW_2016_forGLRIDB <- WW_2016 %>%
#   select(CAS, chnm, Value, `Sample Date`, `Date Deployed`, `Date Retrieved`,
#          SiteID, comment, MDL, `Rep number`) 
  

WW_2016_forToxEval <- Site_data %>%
  ungroup() %>%
  filter(SiteID %in% sites_2016$SiteID, 
         # SiteID != "04085427",
         PARMCODE %in% filter(cas_df, !is.na(POCIS_sampling_rate))$parameter_cd) %>% 
  select(CAS, chnm, SiteID, Value, `Sample Date`, PARMCODE, Class) %>%
  group_by(SiteID,  chnm, CAS, PARMCODE, Class) %>%
  summarize(Value = mean(Value), `Sample Date` = min(`Sample Date`), .groups = "drop")

WW_2016_forGLRIDB <- Site_data %>%
  ungroup() %>%
  filter(SiteID %in% sites_2016$SiteID, 
         PARMCODE %in% filter(cas_df, !is.na(POCIS_sampling_rate))$parameter_cd,
         `Lab ID` %in% 20180370105:20180370120) %>% 
  mutate(MDL = MDL_ngL_30day/1000) %>%
  select(SiteID, CAS, chnm,  PARMCODE, Value, `Sample Date`, Class, 
        `Date Deployed` = `Date_in`, `Date Retrieved` = `Date_out`, 
         comment, MDL, `Sample Name`)

WW_2016_forGLRIDB$`Rep number` <- NA
WW_2016_forGLRIDB$`Rep number`[grepl("Site", WW_2016_forGLRIDB$`Sample Name`)] <- "1"
WW_2016_forGLRIDB$`Rep number`[WW_2016_forGLRIDB$`Sample Name` == "Replicate 1"] <- "2"

WW_2016_forGLRIDB <- select(WW_2016_forGLRIDB, -`Sample Name`) %>%
  arrange(SiteID, PARMCODE, `Rep number`)

AllPOCIS_forToxEval <- Site_data %>%
  ungroup() %>%
  filter(`Lab ID` %in% 20180370105:20180370120)  %>%
  select(CAS, chnm, SiteID, Value_ngPOCIS, `Sample Date`, PARMCODE, Class) %>%
  group_by(SiteID,  chnm, CAS, PARMCODE, Class) %>%
  summarize(Value = mean(Value_ngPOCIS), `Sample Date` = min(`Sample Date`))

  
mdl_2016_forToxEval <- Site_data %>%
  ungroup() %>%
  select(CAS, chnm, PARMCODE, Class, MDL_ngL_30day, `Sample Date`) %>%
  group_by(chnm, CAS, PARMCODE, Class) %>%
  summarize(Value = 0.001*mean(MDL_ngL_30day)) %>%
  mutate(`Sample Date` = as.Date("1986-11-13"))

surface_mdl <- read_excel(file_in(file.path(path_to_data, 'Data/pesticides_dls.xlsx')), sheet='Data')

exclude <- read.csv(file.path(path_to_data, "Data/exclude32.csv"), stringsAsFactors = FALSE) %>%
  left_join(select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by = "CAS") %>%
  select(CAS, endPoint)

if (any(grepl("\\?", exclude$CAS))){
  stop('One of the exclusions CAS contains a "?". Check exclusions spreadsheet.')
}

  
# tox_list <- list("Data" = WW_2016_forToxEval, 
#                  "Chemicals" = chem_info,
#                  "Sites" = sites_2016)
                 
tox_input_list <- list("Data" = WW_2016_forToxEval, 
                 "Chemicals" = cas_df,
                 "Sites" = sites_2016, 
                 "Exclude" = exclude)

tox_benchmark_input_list <- list("Data" = WW_2016_forToxEval, 
                       "Chemicals" = cas_df,
                       "Benchmarks" = benchmarks_df,
                       "Sites" = sites_2016)

allpocis_input_list <- list("Data" = AllPOCIS_forToxEval, 
                       "Chemicals" = cas_df,
                       "Sites" = sites_2016, 
                       "Exclude" = exclude)

glridb_input_list <- list("Data" =  WW_2016_forGLRIDB,
                            "Chemicals" = cas_df,
                            "Sites" = sites_2016, 
                            "Exclude" = exclude)


mdl_2016_forToxEval$SiteID <- sites_2016$SiteID[1]
surface_mdl$SiteID <- sites_2016$SiteID[1]


tox_input_mdl_list <- list("Data" = mdl_2016_forToxEval, 
                       "Chemicals" = cas_df,
                       "Sites" = fake_sites, 
                       "Exclude" = exclude)

tox_input_surf_mdl_list <- list("Data" = surface_mdl, 
                           "Chemicals" = cas_df_surf,
                           "Sites" = fake_sites, 
                           "Exclude" = exclude)

saveOutput = openxlsx::write.xlsx(tox_input_list, file = file_out(file.path(path_to_data, "ToxEvalFiles/Passive2016_ToxEval.xlsx")))

saveOutput = openxlsx::write.xlsx(tox_benchmark_input_list, file = file_out(file.path(path_to_data, "ToxEvalFiles/Passive2016Benchmarks_ToxEval.xlsx")))

saveOutput = openxlsx::write.xlsx(allpocis_input_list, file = file_out(file.path(path_to_data, "ToxEvalFiles/Passive2016AllPOCIS_ToxEval.xlsx")))

saveOutput = openxlsx::write.xlsx(tox_input_mdl_list, file = file_out(file.path(path_to_data, "ToxEvalFiles/Passive2016MDL_ToxEval.xlsx")))

saveOutput = openxlsx::write.xlsx(tox_input_surf_mdl_list, file = file_out(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides_MDL.xlsx")))

saveOutput = openxlsx::write.xlsx(glridb_input_list, file = file_out(file.path(path_to_data, "ToxEvalFiles/Passive2016_ToxEval_forGLRIDB.xlsx")))


POCIS_conc_table <- WW_2016_forToxEval %>%
  group_by() %>%
  select(chnm, CAS) %>%
  unique()

POCIS_detect_table <- AllPOCIS_forToxEval %>%
  group_by() %>%
  select(chnm, CAS) %>%
  unique()

#Table of mdl for pocis data
POCIS_mdl <- Site_data %>%
  ungroup() %>%
  select(chnm, CAS, PARMCODE, MDL_ngPOCIS) %>%
  distinct()

#Calculate percent diff for replicate
# SaginawSummary <- WW_2016_forGLRIDB %>%

SaginawSummary <- Site_data %>%
  ungroup() %>%
  filter(`Lab ID` %in% 20180370105:20180370120)  %>%
  select(CAS, chnm, SiteID, Value_ngPOCIS, `Sample Date`, PARMCODE, Class) %>%
  arrange(PARMCODE) %>%
  filter(SiteID == "04157005") %>%
  group_by(chnm, PARMCODE) %>%
  summarize(across(Value_ngPOCIS, 
                   .fns = list(mean = mean, min = min, max = max, sd = sd), 
                   na.rm = TRUE), 
            n = length(which(Value_ngPOCIS>0))) %>%
  mutate(Value_diff = Value_ngPOCIS_max - Value_ngPOCIS_min,
         Value_perdiff = Value_diff/Value_ngPOCIS_mean*100,
         onedetect = ifelse(Value_ngPOCIS_mean*1000 > 0 & Value_ngPOCIS_min*1000 == 0, "yes", "no"))

  
summary(filter(SaginawSummary, Value_ngPOCIS_min != 0))
data.frame(filter(SaginawSummary, Value_ngPOCIS_min > 0 ))
data.frame(filter(SaginawSummary,  n > 0 ))


table_SI2 <- WW_2016_forToxEval %>%
  filter(Value > 0) %>%
  group_by(CAS, chnm, PARMCODE) %>%
  summarize(across(Value, .fns = list(min = min, max = max, 
                                      median = median, mean = mean, 
                                      n_detects = length)),
            .groups = "drop") %>%
  full_join(unique(select(ungroup(WW_2016_forToxEval), chnm, CAS, PARMCODE))) %>%
  left_join(dplyr::select(cas_df, Class, CAS, PARMCODE = parameter_cd, POCIS_sampling_rate)) %>%
  distinct() %>%
  mutate(Class = factor(Class, c("Herbicide", "Fungicide", "Insecticide", "Deg - Herbicide", "Deg - Fungicide", "Deg - Insecticide", "Other"))) %>%
  arrange(Class, desc(Value_n_detects, Value_mean))

names(table_SI2) <- gsub("Value_", "", names(table_SI2))

no_pocis <- AllPOCIS_forToxEval %>%
  group_by(CAS, chnm, PARMCODE) %>%
  summarize(detects = length(which(Value>0)), .groups = "drop") %>%
  # mutate(CAS = paste0("'", CAS)) %>%
  filter(!PARMCODE %in% table_SI2$PARMCODE)  %>%
  left_join(distinct(select(cas_df, CAS, PARMCODE = parameter_cd, Class, POCIS_sampling_rate))) %>%
  # mutate(detects = 0) %>%
  rename(chemical_name = chnm,
         class = Class) 

table_SI2 <- table_SI2 %>%
  # mutate(chnm2 = tolower(chnm)) %>%
  # full_join(POCIS_uptake_rates, by = c("chnm2")) %>%
  # select(-chnm2) %>%
  select(class = Class, chemical_name = chnm, CAS, PARMCODE, POCIS_sampling_rate, min_value = min,
         max_value = max, median_value = median, mean_value = mean,
         detects = n_detects) %>%
  full_join(no_pocis) %>%
  left_join(unique(select(WW_2016_forGLRIDB, MDL_ngL_30day, CAS, PARMCODE))) %>%
  mutate(detects = ifelse(is.na(detects), 0, detects),
         CAS = paste0("'", CAS),
         class = factor(class, c("Herbicide", "Fungicide", "Insecticide",
                                 "Deg - Herbicide", "Deg - Fungicide", "Deg - Insecticide", "Other"))) %>%
  mutate(MDL = 0.001 * MDL_ngL_30day) %>%
  select(-MDL_ngL_30day, Rs = POCIS_sampling_rate) %>%
  select(class, chemical_name, CAS, Rs, MDL, detects, everything())

table_SI2$class[is.na(table_SI2$class)] <- cas_df$Class[match(tolower(table_SI2$chemical_name[is.na(table_SI2$class)]), cas_df$chnm)]

table_SI2 <- table_SI2 %>%
  arrange(class, desc(detects), desc(mean_value))

print(data.frame(table_SI2))

write.csv(table_SI2, file.path(path_to_data, "SI tables", "chemical_summary_SI2.csv"), 
          row.names = F)

rm(WW_2016_forToxEval, AllPOCIS_forToxEval, locations, exclude, cas_df, benchmarks_df)

