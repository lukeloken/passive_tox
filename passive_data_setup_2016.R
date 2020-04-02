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
cas_df <- read.csv(file_in(file.path(path_to_data, 'Data/Pesticides_2016_monitoring_list_GLRI_USGS_Loken.csv'))) %>%
  rename(chnm = Chemical.Name) %>%
  mutate(chnm = tolower(chnm)) %>%
  select(-CAS, -CAS_wrongformat) %>%
  rename(CAS = CAS2)

cas_df$CAS <- gsub(' ', '', cas_df$CAS)

# cas_df_new <- read_excel(file_in(file.path(path_to_data, 'Data/pesticide CAS numbers.xlsx'))) %>%
  # rename(chnm = `Parameter Name`,
  #        CAS = `CAS Number`, 
  #        ParaCode = `Parameter Code`) %>%
  # mutate(chnm = tolower(chnm)) 

# match_yes<-length(intersect(cas_df$chnm, cas_df_new$chnm))
# match_no <- length(unique(setdiff( cas_df_new$chnm, cas_df$chnm), setdiff( cas_df$chnm,  cas_df_new$chnm)))


#Find single list of chemicals to use
# chem_info = read.csv(file_in(file.path(path_to_data, "Data/chemical_classes.csv")), stringsAsFactors = FALSE)
# chem_info = rename(chem_info, chnm = Chemical.Name)
# 
# chem_info2 = read.csv(file_in(file.path(path_to_data, "Data/chemical_classes2.csv")), stringsAsFactors = FALSE)
# chem_info = rename(chem_info, chnm = Chemical.Name)
# 
# cas_df_combine <- full_join(cas_df, cas_df_new) %>%
#   select(-ParaCode, -M) 
  # full_join(chem_info) %>%
  # full_join(chem_info2)



#Benchmarks
#from Sam's surface water tox eval file
benchmarks_df <- read_excel(file_in(file.path(path_to_data, 'Data/WQ_pesticides_Bench.xlsx')), sheet='Benchmarks') 


#Load pesticide data
WW_2016 = generic_file_opener(file_in(file.path(path_to_data, "RawData/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
                              sheet = "select data",
                              n_max=60,
                              site_sheet = "Site List",
                              year = 2016,
                              skip_site = 2)

# AllPOCIS_2016 = generic_file_opener(file_in(file.path(path_to_data, "Data/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
                              # sheet = "all data",
                              # n_max=300,
                              # site_sheet = "Site List",
                              # year = 2016,
                              # skip_site = 2)

#testing extra CAS files
AllPOCIS_2016 = generic_file_opener(file_in(file.path(path_to_data, "RawData/GLRI 2016 POCIS pesticide data report.xlsx")), cas_df, 
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


AllPOCIS_forToxEval <- AllPOCIS_2016 %>%
  select(CAS, chnm, SiteID, Value, `Sample Date`) %>%
  group_by(SiteID,  chnm, CAS) %>%
  summarize(Value = mean(Value), `Sample Date` = min(`Sample Date`))


mdl_2016_forToxEval <- WW_2016 %>%
  select(CAS, chnm, MDL, `Sample Date`) %>%
  group_by(chnm, CAS) %>%
  summarize(Value = mean(MDL)) %>%
  mutate(`Sample Date` = as.Date("1986-11-13"))


surface_mdl <- read_excel(file_in(file.path(path_to_data, 'Data/pesticides_dls.xlsx')), sheet='Data')

cas_df_surf <- read_excel(file_in(file.path(path_to_data, 'Data/pesticides_dls.xlsx')), sheet='Chemicals')

# AOP_crosswalk = read.csv(file_in(file.path(path_to_data, "Data/AOP_crosswalk.csv")))

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
  select(`Short Name`, `site_no`, Lake, Date_in, Date_out, Datetime_in, Datetime_out)

locations <- readNWISsite(sites_2016$site_no)[c("site_no", "dec_lat_va", "dec_long_va", 'station_nm')]
names(locations) <- c("site_no", 'dec_lat', 'dec_long', 'site_grouping')
locations$site_grouping <- str_sub(locations$site_grouping,-2,-1)


sites_2016<-left_join(sites_2016, locations) %>%
  rename(SiteID = site_no) %>%
  filter(SiteID %in% unique(WW_2016_forToxEval$SiteID)) %>%
  distinct() %>%
  group_by(`Short Name`, SiteID, Lake, site_grouping) %>%
  summarize_all(.funs=mean)

fake_sites <- sites_2016[1,]


exclude = get_exclude(file.path(path_to_data, "Data/exclude.csv")) %>%
  dplyr::select(-chnm) %>%
  full_join(tox_list_example$exclusions)

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


POCIS_conc_table <- WW_2016_forToxEval %>%
  group_by() %>%
  select(chnm, CAS) %>%
  unique()

POCIS_detect_table <- AllPOCIS_forToxEval %>%
  group_by() %>%
  select(chnm, CAS) %>%
  unique()

#Table of mdl for pocis data
POCIS_mdl <- WW_2016 %>%
  select(chnm, MDL, CAS, Class) %>%
  distinct()

rm(WW_2016_forToxEval, WW_2016, AllPOCIS_forToxEval, locations, exclude, cas_df, benchmarks_df)

