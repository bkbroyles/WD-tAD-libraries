##load packages/functions
library(tidyverse)
library(ggpubr)
aavector <- str_split('AVILMWYFSTNQCGPRHKDE', '')[[1]]

#load design lib
dlib <- read_rds('Datasets/design_library.rds') %>% 
  mutate(g_count = str_count(aa_seq,'G'))


##fig1 - C ----
#will include 1WD in chart form
## rest move to supp.
lib <- dlib %>% filter(set == 'set_g_spacing')

lib$Index <- 0

lib$Index[grep('WD',lib$aa_seq)] <- 1
lib$Index[grep('WGD',lib$aa_seq)] <- 2
lib$Index[grep('WGGD',lib$aa_seq)] <- 3
lib$Index[grep('WGGGD',lib$aa_seq)] <- 4
lib$Index[grep('WGGGGD',lib$aa_seq)] <- 5
lib$Index[grep('WGGGGGD',lib$aa_seq)] <- 6
lib$Index[grep('WGGGGGGD',lib$aa_seq)] <- 7
lib$Index[grep('WGGGGGGGD',lib$aa_seq)] <- 8
lib$Index[grep('WGGGGGGGGD',lib$aa_seq)] <- 9
lib$Index[grep('WGGGGGGGGGD',lib$aa_seq)] <- 10

lib$wd <- str_count(lib$aa_seq, 'W')
lib$wd2 <- lib$wd %>% as.factor
lib <- lib %>% filter(slope > -3)

#add max and min spacing indicator
lib$spacing <- 'other'
lib$spacing[1:5] <- 'min'
lib$spacing[c(13,20,27,31,38)] <- 'max'

lib$spacing <- as.factor(lib$spacing)

ro1 <- which(lib$spacing == 'other')
ro2 <- which(lib$spacing == 'max')
ro3 <- which(lib$spacing == 'min')

lib <- lib[c(ro1,ro2,ro3),]

lib$`G's between\n   W and D` <- lib$spacing

p3 <- ggplot(lib, aes(wd2, slope))+
  geom_hline(yintercept = 0, color = 'black', linetype = 'solid', size = 0.8)+
  geom_boxplot(fill = 'grey48', outlier.shape = '.', show.legend = F)+
  scale_fill_distiller(palette = 3)+ #i like 3,4,8,9,10, 16
  geom_jitter(aes(color = spacing), height = 0, width = 0.05, 
              size = 4, alpha = 1, show.legend = T)+scale_color_manual(values = c("#F8766D","#00BFC4",'grey'))+
  #annotate('text', label = 'Stop codon', hjust = 0, x = 0.45, y = 0, vjust = 1)+
  ylim(-2,2.5)+
  theme_pubr()+ylab('Growth Slope')+xlab('W & D count')+ylab(' ')

#ggsave('g-spacing.tiff',p3, height = 4, width = 3, units = 'in', dpi = 800)

#saveRDS(p3, 'Make figure1/fig1_gspacing.rds')


