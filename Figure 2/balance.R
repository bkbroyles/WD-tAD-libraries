##Figure 2 components
## live% of aro - acidic count
## live% of aa x pos
## live% of tetrapeptides x pos (summarise column)
## mixing feature
## table of ML performances

##packages
library(tidyverse)
library(caret)
library(ggpubr)
library(pROC)

#load in wd12 data
wd12 <- read_rds('Datasets/design_library.rds') %>% mutate(
  g_count = str_count(aa_seq, 'G'),
  d_count = str_count(aa_seq, 'D'),
  w_count = str_count(aa_seq, 'W')
) %>% filter(set == 'combinatorial', g_count ==8) %>% 
  mutate(
    balance = w_count - d_count
  )


#find w - d balance first ----
lib <- wd12 %>% filter(g_count == 8)
x <- table(lib$balance, lib$binary_stop)

y <- tibble(die = x[,'die'], live = x[,'live'], balance = row.names(x)) %>% 
  mutate(
    total = live + die,
    live_percent = live/total * 100
  )

y$balance <- as.numeric(y$balance)

p3 <- y %>% filter(total > 10) %>% 
  ggplot(., aes(balance, live_percent, fill = live_percent))+
  geom_hline(yintercept = 33.5, linetype = 'dashed')+
  geom_col(show.legend = F, color = 'black')+
  theme_pubr()+scale_x_continuous(breaks = seq(-8,8,2))+
  #annotate('text', label = 'Library Baseline', x = -8.8, y = 36, hjust = 0)+
  scale_fill_gradient2(midpoint = 33.5, high = '#F8766D', low = '#00BFC4', mid = 'white')+
  coord_cartesian(ylim = c(0,70), xlim = c(-8.5,8.5))+xlab('Balance (W - D)')+
  ylab('Functional tAD %')

##
#ggsave('balance.tiff',p3, height = 4, width = 3, dpi = 800, units = 'in')
