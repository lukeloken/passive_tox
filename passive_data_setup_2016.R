library(drake)
library(tidyverse)
library(googledrive)
library(readxl)
library(data.table)
library(toxEval)
library(openxlsx)
library(dataRetrieval)

path_to_data <- Sys.getenv("PASSIVE_PATH")

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

#Load chemical lists
cas_df <- read.csv(file_in(file.path(path_to_data, 'Data/Pesticides_2016_monitoring_list_GLRI_USGS.csv'))) %>%
  rename(chnm = Chemical.Name) %>%
  mutate(chnm = tolower(chnm))

cas_df_new <- read_excel(file_in(file.path(path_to_data, 'Data/pesticide CAS numbers.xlsx'))) %>%
  rename(chnm = `Parameter Name`,
         CAS = `CAS Number`, 
         ParaCode = `Parameter Code`) %>%
  mutate(chnm = tolower(chnm))

#Currently not usin this one. 
#Find single list of chemicals to use
chem_info = read.csv(file_in(file.path(path_to_data, "Data/chemical_classes.csv")), stringsAsFactors = FALSE)
chem_info = rename(chem_info, Chemical = Chemical.Name)


setdiff(cas_df$chnm, cas_df_new$chnm)
setdiff( cas_df_new$chnm, cas_df$chnm)


#Load pesticide data
WW_2016 = generic_file_opener(file_in(file.path(path_to_data, "Data/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
                              sheet = "select data",
                              n_max=60,
                              site_sheet = "Site List",
                              year = 2016,
                              skip_site = 2)

AllPOCIS_2016 = generic_file_opener(file_in(file.path(path_to_data, "Data/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
                              sheet = "all data",
                              n_max=300,
                              site_sheet = "Site List",
                              year = 2016,
                              skip_site = 2)

# average the two replicates for Saginow River Station
WW_2016_forToxEval <- WW_2016 %>%
  select(CAS, chnm, SiteID, Value, `Sample Date`) %>%
  group_by(SiteID,  chnm, CAS) %>%
  summarize(Value = mean(Value), `Sample Date` = min(`Sample Date`))


# AOP_crosswalk = read.csv(file_in(file.path(path_to_data, "Data/AOP_crosswalk.csv")))

sites_2016 = readxl::read_excel(file_in(file.path(path_to_data, "Data/GLRI 2016 POCIS pesticide data report.xlsx")), sheet = "Site List", skip = 2) %>%
  rename(site_no = `USGS Station ID`) %>%
  # mutate(`Short Name` = substr(`Site Name`, start = 1, stop = 8)) %>%
  select(`Short Name`, `site_no`, Lake)
  
locations <- readNWISsite(sites_2016$site_no)[c("site_no", "dec_lat_va", "dec_long_va", 'station_nm')]
names(locations) <- c("site_no", 'dec_lat', 'dec_long', 'site_grouping')
locations$site_grouping <- str_sub(locations$site_grouping,-2,-1)


sites_2016<-left_join(sites_2016, locations) %>%
  rename(SiteID = site_no) %>%
  filter(SiteID %in% unique(WW_2016_forToxEval$SiteID)) %>%
  distinct()


exclude = get_exclude(file.path(path_to_data, "Data/exclude.csv"))

# tox_list <- list("Data" = WW_2016_forToxEval, 
#                  "Chemicals" = chem_info,
#                  "Sites" = sites_2016)
                 
tox_input_list <- list("Data" = WW_2016_forToxEval, 
                 "Chemicals" = cas_df,
                 "Sites" = sites_2016, 
                 "Exclude" = exclude)

saveOutput = openxlsx::write.xlsx(tox_input_list, file = file_out(file.path(path_to_data, "Data/PassiveForToxEval.xlsx")))

