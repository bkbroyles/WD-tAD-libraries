##load packages/functions
library(tidyverse)
library(ggpubr)
aavector <- str_split('AVILMWYFSTNQCGPRHKDE', '')[[1]]

#load design lib
dlib <- read_rds('Datasets/design_library.rds') %>% 
  mutate(g_count = str_count(aa_seq,'G'))

#g in vs g out----

#grab two subsets
lib <- dlib %>% filter(set %in% c('g_in','g_out')) %>% mutate(
  set = ifelse(set == 'g_in', 'G_inwards','G_outwards')
)

#Find missing sequence in other libs, the sequence is in both subsets so add it to both
lib2 <- dlib %>% filter(aa_seq == 'WDWDWDWDWDWDWDWDWDWD')
lib2 <- rbind(lib2, lib2)
lib2$set[1:2] <- 'G_inwards'
lib2$set[3:4] <- 'G_outwards'

#join two libraries
lib <- rbind(lib,lib2)

#do not include 20 G or 19 G sequences
lib <- lib %>% filter(g_count < 19)

#add index for plotting
lib$Index <- 20 - str_count(lib$aa_seq, 'G')

lib$` ` <- lib$set

p4 <- ggplot(lib, aes(Index, slope, color = ` `, fill = ` `))+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0), 
            fill = 'grey85', color = 'grey85')+
  geom_hline(yintercept = 0, color = 'black', linetype = 'solid', size = 0.8)+
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), size = 2.2,
              se = F, show.legend = F)+
  geom_jitter(height = 0,size = 2, alpha = 1, width = 0.1, show.legend = F)+
  theme_pubr()+scale_x_continuous(breaks = seq(2,20,2))+
  ylab('Growth Slope')+xlab('Length of WD-module')+
  coord_cartesian(clip = "off",ylim = c(-1.2, 2.5), xlim = c(2, 20))+ylab(' ')+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_color_manual(values = c("#F8766D","#00BFC4"))

p4

  #ggsave('g-in_g-out.tiff',p4,height = 4,width = 3,units = 'in',dpi = 400)

#saveRDS(p4, 'Figure 1/fig1_Gin-out.rds')
