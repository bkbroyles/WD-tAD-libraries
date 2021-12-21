##need to add WD repeats and DW repeats to fig 1
library(tidyverse)
library(ggpubr)

lib<- read_rds('Datasets/design_library.rds') %>% 
  mutate(
    g_count = str_count(aa_seq, 'G')
  )

lib <- lib %>% mutate(
  w = str_count(aa_seq, 'W'),
  d = str_count(aa_seq, 'D'),
  tots = g_count + w + d
)

hold <- lib %>% filter(w == 1, d == 1, tots == 20)
hold$index <- c(1,1,2:11)
hold$` ` <- c('DW','DW','DW', rep('WD',9))


plot_1w1d <- ggplot(hold, aes(index, slope, color = ` `))+
  geom_hline(yintercept = 0, size = 0.8)+
  geom_point(size = 2)+
  ylim(-2,2.5)+
  theme_pubr()+
  scale_x_continuous(breaks = 1:12)+ylab(' ')

#ggsave('1w1d.tiff', plot_1w1d, height = 4, width = 3, units = 'in', dpi = 800)


#saveRDS(plot_1w1d, 'fig1_1w1d.rds')
