

#### Setup ####
library(toxEval)
# library(ToxMixtures)
library(dplyr)
library(RColorBrewer)
library(ggplot2)


# path_to_data <- c('C:/Users/lloken/DOI/Corsi, Steven R - GLRI CECs/2016/Manuscripts/Pesticides_passive')
# tox_list<- create_toxEval(file.path(path_to_data, "ToxEvalFiles/Passive2016_ToxEval.xlsx"))
# 
# # tox_list <- create_toxEval(path_to_file)
# ACC <- get_ACC(tox_list$chem_info$CAS)
# ACC <- remove_flags(ACC = ACC)
# 
# cleaned_ep <- clean_endPoint_info(end_point_info)
# filtered_ep <- filter_groups(cleaned_ep, 
#                              remove_groups = c('Background Measurement', 'Undefined', 'NA', 'Cell Cycle'),
#                              assays = c('ACEA','APR','ATG','NVS','OT','TOX21','CEETOX','CLD','TANGUAY','NHEERL_PADILLA','NCCT','NHEERL_HUNTER','NHEERL_NIS','NHEERL_MED','UPITT'))
# 
# chemicalSummary <- get_chemical_summary(tox_list, ACC, filtered_ep)
######################################

# plot_stackbar(chemical_summary, x='shortName', x_label='Site', y_label="Number of chemicals",
#               fill='EAR', stack='chnm', 
#               breaks=c(.0001, .001, .01)) + 
#   theme(axis.text.x = element_text(angle=90, hjust=1))
# 

chemicalSummary2 <- chemicalSummary %>%
  mutate(shortName = factor(shortName, site_order))

plot_stackbar(chemicalSummary2, x='shortName', x_label='Site', y_label="Number of chemicals",
              fill='EAR', stack='chnm', 
              breaks=c(.0001, .001, .01)) + 
  theme(axis.text.x = element_text(angle=90, hjust=1))


plot_stackbar(chemicalSummary2, x='chnm', x_label='Chemical', y_label="Number of rivers",
              fill='EAR', stack='site', 
              breaks=c(.0001, .001, .01)) + 
  theme(axis.text.x = element_text(angle=90, hjust=1))


# cs <- chemical_summary
# group_by_this <- "endPoint"
# ear_threshold <- 0.001
# site_threshold_percent <- 10
# n_sites <- length(unique(chemical_summary$site))
# site_threshold <- ceiling(n_sites / site_threshold_percent)

group_by_this <- "endPoint"
ear_threshold <- .001
site_threshold <- 2

top_combos <- top_mixes(chemicalSummary2, group_by_this, 
                        ear_threshold,
                        site_threshold) 


overall_max_n_chem <- overall_mixtures(top_combos, "max")


summed_EARs <- chemicalSummary2 %>%
  group_by(site, date, !!sym(group_by_this)) %>%
  summarize(sum_ear_endpoint = sum(EAR)) %>%
  filter(!!sym(group_by_this) %in% top_combos$endPoint) %>%
  ungroup() %>%
  distinct()

endpoint_rank <- summed_EARs %>%
  group_by(endPoint) %>%
  summarize(sum_ear_median = median(sum_ear_endpoint, na.rm=T)) %>%
  arrange((sum_ear_median))

summed_EARs <- summed_EARs %>%
  left_join(overall_max_n_chem) %>%
  mutate(endPoint = factor(endPoint, endpoint_rank$endPoint))

top_mix_box <- ggplot(summed_EARs, aes(y=endPoint, x=sum_ear_endpoint)) +
  geom_vline(xintercept=0.001, linetype='dashed') +
  geom_jitter(shape=16, color='red', alpha=0.5, width=0, height=.1, size=1.5) + 
  geom_boxplot(outlier.shape=NA, fill='red', alpha=.2, width=.5) + 
  scale_x_log10nice(name = expression(paste(EAR[mixture]))) +
  labs(y='ToxCast assay name') + 
  theme_tox()
  
ggsave(file_out(file.path(path_to_data, "Figures/PriorityEndpointsBoxplot.png")), top_mix_box, height=4, width=5, units='in')


overall_max_n_chem <- overall_mixtures(top_combos, "max")  %>%
  mutate(endPoint = factor(endPoint, endpoint_rank$endPoint)) %>%
  arrange(desc(endPoint))

class_key <- class_key_fnx(chemicalSummary2)

overall_df_fancy <- overall_df_format(overall_max_n_chem,
                                      class_key)


data.frame(top_combos)
data.frame(overall_max_n_chem)
data.frame(overall_df_fancy)


knitr::kable(head(overall_df_fancy))
