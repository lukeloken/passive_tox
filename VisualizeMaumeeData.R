#Look at Maumee data for 2016

maumee_data <- read_excel(file_in(file.path(path_to_data, 'Data/Maumee_All_data_FromFiles.xlsx')), sheet = 'Data') 

maumee_sites <- read_excel(file_in(file.path(path_to_data, 'Data/Maumee_All_data_FromFiles.xlsx')), sheet = 'Sites') 

maumee_chems <- read_excel(file_in(file.path(path_to_data, 'Data/Maumee_All_data_FromFiles.xlsx')), sheet = 'Chemicals') 


head(maumee_df)

summary(maumee_data_merge)

maumee_data_merge <- maumee_data %>%
  left_join(maumee_chems) %>%
  left_join(maumee_sites) %>%
  mutate(Date = as.Date(`Sample Date`),
         ShortName = `Short Name`) %>%
  filter(year(Date) == 2016) %>%
  group_by(SiteID, ShortName, Date) 

maumee_data_test <- maumee_data_merge %>%
# filter(grepl('<', remark_cd)==FALSE) %>%
  summarize(n_chems_tested = length(CAS))

maumee_data_detected <- maumee_data_merge %>%
  filter(grepl('<', remark_cd)==FALSE) %>%
  summarize(n_chems_detected = length(CAS)) %>%
  full_join(maumee_data_test) %>%
  mutate(label = paste0(n_chems_detected, ' of ', n_chems_tested))


timeseries_plot <- ggplot(maumee_data_detected, aes(x=Date, y=ShortName)) + 
  geom_point(shape=4, size=2, stroke=2) +
  geom_label_repel(aes(label=label)) + 
  theme_bw() +
  ggtitle('Maumee sites 2016: Water and sediment contaminants detected and analyzed')

print(timeseries_plot)

ggsave(file.path(path_to_data, 'Figures', 'Maumee', 'Maumee_sites_bytime.png'), timeseries_plot)


