# Passive sampler analysis

## Fork by Luke Loken

Rather than usin `drake`, I have been mapping data and outputs using a single character directory. At the start of the `PassiveData2016_Workflow` script I point to the project folder

```
path_to_data <- c('C:/Users/lloken/OneDrive - DOI/GLRI_Pesticides')

```

Then I call a handful of smaller scripts usin `source` to process the data and generate figures. 

```
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
```

These scripts will generate folders within the `path_to_data` location.

# Previous code and description


We're trying the `drake` R package for orchestrating the workflow. 

To initially get the datasets and clean up the data, run the `passive_data_setup.R` script. If there are any issues, you can try to diagnose the problems using `vis_drake_graph`.

## Converting to OneDrive

Halfway through this project, we switched from GoogleDrive to OneDrive. We're trying to take advantage of the automatic snycing of the data folders. To do this, each collaborator has to:

1. Go to the canonical shared drive location
2. Click the "Sync" button and open in Microsoft OneDrive
3. In Windows Explore* right-click on the shared folder and choose "Alway Keep on this device". (Note: later when you want to remove this, you'll need to right-click -> Settings -> Stop Syncing before you can delete it)
4. Create a system variable with the path to the local copy of the canonical shared drive:

```
rprofile_path = file.path(Sys.getenv("HOME"), ".Rprofile")

write('Sys.setenv(PASSIVE_PATH = "C:/Users/ldecicco/DOI/Corsi, Steven R - Manuscript")',
      rprofile_path, 
      append =  TRUE)

cat('Your Rprofile has been updated to include PASSIVE_PATH
    Please restart R for changes to take effect.')
```
You MUST restart your R session (restart RStudio)!

Then, you can call the chemicalSummary via:
```
chemicalSummary <- readRDS(file = file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","clean","chemical_summary.rds"))
```

If you've run the drake plans, this will also work:
```
loadd(chemicalSummary)
```


## Old GoogleDrive stuff:

The next drake workplan is in "passive_analysis.R".

Once that file is run sucessfully, there will be an Excel file: "data/clean/passive.xlsx". This file can be used for regular `toxEval` workflows.

Updating the data on the Google Drive must follow this pattern to work correctly:
https://www.labnol.org/internet/update-files-in-google-drive/28928/

## Disclaimer

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.