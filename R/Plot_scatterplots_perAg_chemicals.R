
library(easypackages)
from_import("plyr", "ddply")
# library(plyr)

basins_OGR <- readRDS( file_in(file.path(path_to_data, "Rdata/watershedmetrics.rds")))

sites <- sites_2016 %>%
  rename(site_no = SiteID)

# get info about number of chemicals by site, merge with sites
nchems <- site_detection %>%
  group_by(site) %>%
  summarize(chemicals = sum(value))

graphData <- site_detection %>%
  filter(group %in% c('AboveEAR0.001', 'AboveEAR0.01')) %>%
  dplyr::group_by(site) %>%
  dplyr::select(site, value) %>%
  rename(EAR = value) %>%
  dplyr::summarize(EAR=sum(EAR)) %>%
  left_join(nchems) %>%
  data.frame()

sites <- left_join(sites, graphData, by = c('site_no' = 'site')) %>% 
  left_join(basins_OGR@data, by = c('site_no' = 'STAID' ))


graphData2 <- site_detection %>%
  dplyr::group_by(site) %>%
  dplyr::select(site, value, group) %>%
  spread(key=group, value=value) %>%
  mutate(detect = OnlyDetected + AboveEAR0.0001 + AboveEAR0.001,
         low = AboveEAR0.0001 + AboveEAR0.001,
         high = AboveEAR0.001) %>%
  rename(STAID = site) %>%
  left_join(basins_OGR@data) %>%
  mutate(perAgUrban = (perAg + perUrban)*100) %>%
  select(STAID, perAgUrban, detect, low, high) %>%
  gather(key='group', value='value', 3:5) %>%
  mutate(group=factor(group,c( 'detect', 'low', 'high')))



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







#plot percent Ag/Ur versus each EAR category

labels_EAR <- c('Detected', expression(paste("EAR > ", "10"^"-4")), expression(paste("EAR > ", "10"^"-3")))

find_hull <- function(df) {
  df[chull(df$perAgUrban, df$value), ]}

  hulls <- ddply(graphData2, "group", find_hull) 
# 
# hull <- graphData2 %>%
#   filter(group == 'detect') %>%
#   slice(chull(perAgUrban, value))

scatter_v2 <- ggplot(data=graphData2, aes(x=perAgUrban, y=value, group=group, fill=group, color=group)) +
  # geom_smooth(method='lm', se=F, size=2) + 
  geom_polygon(data=hulls, alpha=.3) +
  geom_point(size=2, pch=21, color='black', stroke=1) +
  theme_bw() +
  labs(x='Percent agriculture + urban', y='Number of chemicals') +
  scale_fill_brewer(palette = "YlOrRd", labels=labels_EAR) +
  scale_color_brewer(palette = "YlOrRd", labels=labels_EAR) +
  theme(legend.position='bottom', legend.title = element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_alpha(guide = 'none')

  
print(scatter_v2)

ggsave(file_out(file.path(path_to_data, "Figures/Scatter_Detections_bylanduse_v2.png")), plot = scatter_v2, height=4, width=4)


scatter_v3 <- ggplot(data=graphData2, aes(x=perAgUrban, y=value, shape=group, group=group, fill=group, color=group)) +
  geom_smooth(method='lm', se=F, size=2) +
  # geom_polygon(data=hulls, alpha=.3) +
  geom_point(size=2, color='black', stroke=1) +
  theme_bw() +
  labs(x='Percent agriculture + urban', y='Number of chemicals') +
  scale_fill_brewer(palette = "YlOrRd", labels=labels_EAR) +
  scale_color_brewer(palette = "YlOrRd", labels=labels_EAR) +
  scale_shape_manual(values = c(21,22,24), labels=labels_EAR) + 
  theme(legend.position='bottom', legend.title = element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_alpha(guide = 'none') +
  scale_size(guide='none') +
  guides(color = guide_legend(override.aes = list(linetype = 0)))

print(scatter_v3)

ggsave(file_out(file.path(path_to_data, "Figures/Scatter_Detections_bylanduse_v3.png")), plot = scatter_v3, height=4, width=4)

