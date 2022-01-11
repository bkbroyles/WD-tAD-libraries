##load packages/functions
library(tidyverse)
library(ggpubr)
aavector <- str_split('AVILMWYFSTNQCGPRHKDE', '')[[1]]

#load design lib
dlib <- read_rds('Datasets/design_library.rds') %>% 
  mutate(g_count = str_count(aa_seq,'G'))


##fig 1 - B----
#get w_first and d_first subsets
lib <- dlib %>% filter(set %in% c('W_first', 'D_first'))

#two sequences were missing data, I will grab them from other sets
lib2 <- dlib %>% filter(set != 'W_first', set != 'D_first') %>% 
  filter(aa_seq %in% c('GGGGGGGGGGGGGGGGGGWD'))
lib2$set <- 'W_first'
#join the two libraries and add index column
lib <- rbind(lib, lib2)
lib$Index <- 20 - str_count(lib$aa_seq, 'G')

#change subset names to smaller names to help graph visiualization 
#lib <- mutate(lib, set = ifelse(grepl('WD', aa_seq),'b1','b2'))

lib$` ` <- lib$set

##fig1 - B - plot ----
p2 <- ggplot(lib, aes(Index, slope, color = ` `, fill = ` `))+
  #Below - these add horizontal control sequence lines
  geom_hline(yintercept = 0, color = 'black', linetype = 'solid', size = 0.8)+
  #Below two lines are to plot library data
  #geom_smooth(se = F, show.legend = F, size = 4)+
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), size = 3,
              se = F, show.legend = F)+
  geom_jitter(height = 0,size = 0.9, alpha = 1, width = 0.05,
              show.legend = F)+
  #These annotates add labels to h-lines
  #Below here is figure finalizing things
  theme_pubr()+scale_x_continuous(breaks = seq(2,20,2))+
  ylab('Growth Slope')+xlab('Length of WD-module')+
  coord_cartesian(clip = "off",ylim = c(-2, 1.8), xlim = c(2, 20))+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_color_manual(values = c("#F8766D","#00BFC4"))

#ggsave('d-wfirst.tiff',p2, height = 4, width = 3, dpi = 800, units = 'in')

#saveRDS(p2, 'fig1_D-Wfirst.rds')
