##load packages/functions
library(tidyverse)
library(ggpubr)
aavector <- str_split('AVILMWYFSTNQCGPRHKDE', '')[[1]]

#load design lib and control sequences
dlib <- read_rds('Datasets/design_library.rds') %>% 
  mutate(g_count = str_count(aa_seq,'G'))

known_sample <- read_rds('Datasets/known_sample.rds')

#Grab internal and terminal sets
#need internal and terminal wd sets, exclude sequences with 19 and 20 G's
lib <- dlib %>% filter(set %in% c('internal_WD', 'terminal_WD')) %>% 
  filter(g_count < 19)

#Long sequences were missing, and short sequences can be supplemented
lib2 <- dlib %>%
  filter(aa_seq %in% c('DWDWDWDWDWDWDWDWDWDW', 'WDWDWDWDWDWDWDWDWDWD',
                       'GGGGGGGGGGGGGGGGGGWD'))

#just need to fix set names since these sequences are from different subsets
lib2$set[c(1,3)] <- 'internal_WD'
lib2$set[c(2,4,5)] <- 'terminal_WD'

#join the two libraries together, and add index column to order sequences in the figure
lib <- rbind(lib, lib2)
lib$Index <- 20 - str_count(lib$aa_seq, 'G')

#copy set to empty column name to remove name from plot legend
lib$` ` <- lib$set

##fig1 - A - plot ----
p1 <- ggplot(lib, aes(Index, slope, color = ` `, fill = ` `))+
  #First - add horizontal control sequence lines
  geom_segment(x = 1, xend = 20, y = 0, yend = 0, 
               color = 'black', linetype = 'solid', size = 0.8)+
  geom_segment(x = 1, xend = 20, y = known_sample$avg[1], yend = known_sample$avg[1],
               color = '#62007a', linetype = 'dotdash', size = 0.6)+
  geom_segment(x = 1, xend = 20, y = known_sample$avg[2], yend = known_sample$avg[2],
               color = '#62007a', linetype = 'dotdash', size = 0.6)+
  geom_segment(x = 1, xend = 20, y = known_sample$avg[3], yend = known_sample$avg[3],
               color = '#007315', linetype = 'longdash', size = 0.6)+
  geom_segment(x = 1, xend = 20, y = known_sample$avg[4], yend = known_sample$avg[4],
               color = '#007315', linetype = 'longdash', size = 0.6)+
  geom_segment(x = 1, xend = 20, y = known_sample$avg[5], yend = known_sample$avg[5],
               color = '#826200', linetype = 'dashed', size = 0.6)+
  #Second - add smooth trend lines
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), size = 3,
              se = F, show.legend = F)+
  #Third - add individual points
  geom_jitter(height = 0,size = 0.9, alpha = 1, width = 0.05)+
  #Fourth - add control sequences labels
  annotate('text', label='VP16-minx2', x=20.1, y=known_sample$avg[4], color='#007315', hjust=0, size=4)+
  annotate('text', label='Gal4 (860-872)', x=20.1, y=known_sample$avg[2], color='#62007a', hjust=0, size=4)+
  annotate('text', label='VP16-min', x=20.1, y=known_sample$avg[3]+0.1, color='#007315', hjust=0, size=4)+
  annotate('text', label='Gal4 (840-857)', x=20.1, y=known_sample$avg[1]+0.04, color='#62007a', hjust=0, size=4)+
  annotate('text', label='Stop Codon', x=20.1, y=0, color='black', hjust=0, size=4)+
  annotate('text', label='20 G\'s', x=20.1, y=known_sample$avg[5], color='#826200', hjust=0, size=4)+
  #Fifth - finalize look of the plot
  theme_pubr(legend = 'none')+scale_x_continuous(breaks = seq(2,20,2))+
  ylab('Growth Slope')+xlab('Length of WD-module')+
  coord_cartesian(clip = "off",xlim = c(2, 25),ylim = c(-1.2,2.5))+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_color_manual(values = c("#F8766D","#00BFC4"))

p1

#ggsave('internal-terminal.tiff', p1, height = 4, width = 5, units = 'in', dpi = 800)

#saveRDS(p1, 'Figure 1/fig1_internal-terminalWD.rds')

