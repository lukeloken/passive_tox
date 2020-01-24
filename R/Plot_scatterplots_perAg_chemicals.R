


sites <- sites_2016 %>%
  rename(site_no = SiteID)

# get info about number of chemicals by site, merge with sites
nchems <- site_detection %>%
  group_by(site) %>%
  summarize(chemicals = sum(value))

graphData <- site_detection %>%
  filter(group=='AboveEAR0.001') %>%
  dplyr::group_by(site) %>%
  dplyr::select(site, value) %>%
  rename(EAR = value) %>%
  left_join(nchems) %>%
  data.frame()



basins_OGR <- readRDS( file_in(file.path(path_to_data, "Rdata/watershedmetrics.rds")))


sites <- left_join(sites, graphData, by = c('site_no' = 'site')) %>% 
  left_join(basins_OGR@data, by = c('site_no' = 'STAID' ))


par(mfrow=c(1,2))
plot(sites$chemicals ~ sites$perAg)
AGmodel <-lm(sites$chemicals ~ sites$perAg)
abline(AGmodel)
summary(AGmodel)
residuals(AGmodel)
plot(residuals(AGmodel) ~ sites$perUrban)
abline(lm(residuals(AGmodel) ~ sites$perUrban))             


AGUrbanmodel <-lm(sites$chemicals ~ (sites$perAg+sites$perUrban))
summary(AGUrbanmodel)
residuals(AGUrbanmodel)

AGUrbanmodel2 <-lm(sites$EAR ~ (sites$perAg+sites$perUrban))
summary(AGUrbanmodel2)
residuals(AGUrbanmodel)


par(mfrow=c(1,1))
chems <- ggplot(data=sites, aes(x=(perAg+perUrban)*100, y=chemicals)) +
  geom_point(size=2) +
  theme_bw() +
  labs(x='Percent Agriculture plus Urban', y='Number of chemicals detected')

chems_AGonly <- ggplot(data=sites, aes(x=(perAg)*100, y=chemicals)) +
  geom_point(size=2) +
  theme_bw() +
  labs(x='Percent Agriculture', y='Number of chemicals detected')


ears<- ggplot(data=sites, aes(x=(perAg+perUrban)*100, y=EAR)) +
  geom_point(size=2) +
  theme_bw() +
  labs(x='Percent Agriculture plus Urban', y='Number of chemicals above EAR 0.001')

ears_AGonly<- ggplot(data=sites, aes(x=(perAg)*100, y=EAR)) +
  geom_point(size=2) +
  theme_bw() +
  labs(x='Percent Agriculture', y='Number of chemicals above EAR 0.001')



plotsbyAgUrban <- grid.arrange(grobs = list(chems,ears), nrow=1)


ggsave(file_out(file.path(path_to_data, "Figures/Scatter_chemsears_bylanduse.png")), plot = plotsbyAgUrban, height=4, width=8)


