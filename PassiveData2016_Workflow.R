
# Workflow for analyzing passive sampler pesticide data
# Data are from 15 streams flowing into the Great Lakes
# Merge data with Sam Oliver's paper

# Luke Loken
# January 2020

# path_to_data <- c('C:/Users/lloken/OneDrive - DOI/GLRI_Pesticides')
path_to_data <- c('C:/Users/lloken/DOI/Corsi, Steven R - GLRI CECs/2016/Manuscripts/Pesticides_passive')


#load libraries
library(gridExtra)
library(RColorBrewer)

#load custom functions
source('R/functions/g_legend.R')
source('R/functions/ScaleYLog10Nice.R')


#Load passive sampler data and generate toxEval file
#Saves toxeval file (excel) for select data, all data, and TQ benchmarks data
source('passive_data_setup_2016.R')

#Evaluate toxicity using ToxEval for passive data
#Loads two objects into environment (chemicalSummary and chemicalSummary_allpocis)
source('ToxEval_passive2016.R')

#Evaluate toxicity using ToxEval for water samples data
source('ToxEval_watersamples2016.R')


#Plot number of chemicals by site, and number of sites by chemical barplots
source('Plot_PassiveTox_2016.R')

#Plot number of chemicals by site, and number of sites by chemical barplots for water samples
source('Plot_SurfaceTox_2016.R')

#Combine water and passive samples and compare. 
source('Compare_Plot_Passive_V_WaterSample.R')

#Make map figure of EAR, number of chemicals detected with watersheds colored by % ag
source('map_sites_passive_2016.R')

#Scatterplots of land use versus chemical detections
source('R/Plot_scatterplots_perAg_chemicals.R')
