
# Workflow for analyzing passive sampler pesticide data
# Data are from 15 streams flowing into the Great Lakes
# Merge data with Sam Oliver's paper

# Luke Loken
# January 2020

path_to_data <- c('C:/Users/lloken/OneDrive - DOI/GLRI_Pesticides')

#Load passive sampler data and generate toxEval file
#Saves excel file for both select and all data
source('passive_data_setup_2016.R')

#Evaluate toxicity using ToxEval for passive data
#Loads two objects into environment (chemicalSummary and chemicalSummary_allpocis)
source('ToxEval_passive2016.R')

#Plot number of chemicals by site, and number of sites by chemical barplots
source('Plot_PassiveTox_2016.R')

#Make map figure of EAR, number of chemicals detected with watersheds colored by % ag
source('R/map_sites_passive_2016.R')