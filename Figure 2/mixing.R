##Mixing graph for hahn and wd12
library(tidyverse)
library(ggpubr)

wd12 <- read_rds('Datasets/design_library.rds') %>% mutate(
  g_count = str_count(aa_seq, 'G')
) %>% filter(set == 'combinatorial', g_count == 8)

wd12 <- wd12 %>% mutate(
  acid = str_count(aa_seq,'[DE]'),
  aro = str_count(aa_seq,'[WYF]'),
  aro_acid = aro + acid
) %>% filter(aro_acid > 1)


str_locate_all(wd12$aa_seq[1:10],'[DE]')[[1]][,1]

acid_pos <- lapply(wd12$aa_seq,function(x){
  str_locate_all(x,'[DE]')[[1]][,1]
})

aro_pos <- lapply(wd12$aa_seq,function(x){
  str_locate_all(x,'[WYF]')[[1]][,1]
})

wd12$mixing <- 0
for(i in 1:nrow(wd12)){
  x <- acid_pos[[i]]
  y <- aro_pos[[i]]
  
  df <- data.frame(pos = c(x,y), aa = c(rep('D',length(x)),rep('W',length(y))), stringsAsFactors = F)
  df <- df[order(df$pos),]
  
  counter <- 0
  aa <- df$aa[1]
  
  for(j in 2:nrow(df)){
    aa2 <- df$aa[j]
    if(aa != aa2){
      aa <- aa2
      counter <- counter + 1
    }
  }
  wd12$mixing[i] <- counter
}

##add feature columns, W-D balance, mixing, and tetrapeptides
lib <- wd12 %>% filter(aro == 6, acid == 6)

plotme <- tibble(mixing = 1:11, 
                 live_n = 0, die_n = 0)

for(i in 1:nrow(plotme)){
  lib2 <- lib %>% filter(mixing == plotme$mixing[i])
  x <- table(lib2$binary_stop)
  plotme$live_n[i] <- x['live']
  plotme$die_n[i] <- x['die']
}

plotme <- plotme %>% mutate(
  total = live_n + die_n,
  live_percent = live_n/total * 100
)

##keep composition the same
p4 <- plotme %>% filter(total > 30) %>% 
  ggplot(., aes(mixing, live_percent, fill = live_percent))+
  geom_hline(yintercept = 58.1, linetype = 'dashed')+
  geom_col(show.legend = F,color = 'black')+
  theme_pubr()+scale_x_continuous(breaks = 1:10)+
  scale_fill_gradient2(midpoint = 58.1, high = '#00BFC4', low = '#F8766D', mid = 'white')+
  coord_cartesian(ylim = c(0,80))+xlab('Aromatic/Acidic Mixing')+
  ylab('Functional tAD %')

#ggsave('mixing-wd12.tiff',p4, height = 4, width = 3, dpi = 800, units = 'in')

#saveRDS(p4, 'mixing_wd12.rds')
