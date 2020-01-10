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
# source("R/analyze/get_chem_info.R")
# source("R/analyze/create_tox_file.R")


cas_df <- read.csv(file_in(file.path(path_to_data, 'Data/Pesticides_2016_monitoring_list_GLRI_USGS.csv'))) %>%
  rename(chnm = Chemical.Name) %>%
  mutate(chnm = tolower(chnm))

WW_2016 = generic_file_opener(file_in(file.path(path_to_data, "Data/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
                              sheet = "select data",
                              n_max=60,
                              site_sheet = "Site List",
                              year = 2016,
                              skip_site = 2)

WW_2016_forToxEval <- WW_2016 %>%
  select(CAS, chnm, SiteID, Value, `Sample Date`)

chem_info = read.csv(file_in(file.path(path_to_data, "Data/chemical_classes.csv")), stringsAsFactors = FALSE)
chem_info = rename(chem_info, Chemical = Chemical.Name)

# AOP_crosswalk = read.csv(file_in(file.path(path_to_data, "Data/AOP_crosswalk.csv")))

sites_2016 = readxl::read_excel(file_in(file.path(path_to_data, "Data/GLRI 2016 POCIS pesticide data report.xlsx")), sheet = "Site List", skip = 2) %>%
  rename(site_no = `USGS Station ID`) %>%
  mutate(`Short Name` = substr(`Site Name`, start = 1, stop = 8)) %>%
  select(`Short Name`, `site_no`)
  
locations <- readNWISsite(sites_2016$site_no)[c("site_no", "dec_lat_va", "dec_long_va")]
names(locations) <- c("site_no", 'dec_lat', 'dec_long')

sites_2016<-left_join(sites_2016, locations) %>%
  rename(SiteID = site_no)


# tox_list <- list("Data" = WW_2016_forToxEval, 
#                  "Chemicals" = chem_info,
#                  "Sites" = sites_2016)
                 
tox_input_list <- list("Data" = WW_2016_forToxEval, 
                 "Chemicals" = cas_df,
                 "Sites" = sites_2016)

saveOutput = openxlsx::write.xlsx(tox_input_list, file = file_out(file.path(path_to_data, "Data/PassiveForToxEval.xlsx")))

