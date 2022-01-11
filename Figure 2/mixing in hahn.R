##mixing in hahn library
lib <- read_rds('Datasets/Gcn4_random_library.rds')

lib <- lib %>% mutate(
  acid = str_count(aa_seq,'[DE]'),
  aro = str_count(aa_seq,'[WYF]'),
  aro_acid = aro + acid
) %>% filter(aro_acid > 1)


str_locate_all(lib$aa_seq[1:10],'[DE]')[[1]][,1]

acid_pos <- lapply(lib$aa_seq,function(x){
  str_locate_all(x,'[DE]')[[1]][,1]
})

aro_pos <- lapply(lib$aa_seq,function(x){
  str_locate_all(x,'[WYF]')[[1]][,1]
})

lib$mixing <- 0
for(i in 1:nrow(lib)){
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
lib$mixing[i] <- counter
}


#saveRDS(lib, 'savemeplz.rds')
lib <- read_rds('Figure 2/savemeplz.rds')

table(lib$mixing, lib$binary_stop)

##6 acidic 6 aromatic----
df4 <- lib %>% filter(aro == 6 & acid == 6)
x <- table(df4$mixing,df4$binary_stop)
y <- tibble(die = x[,'die'],live = x[,'live'],
            mixing = row.names(x)) %>% mutate(
              total = live + die,
              live_percent = live/total * 100
            )
sum(y$live)
sum(y$total)
1148/3018 #live% 38.0

y$mixing <- as.numeric(y$mixing)

p3 <- y %>% filter(total > 30) %>% 
  ggplot(., aes(mixing, live_percent, fill = live_percent))+
  geom_hline(yintercept = 38, linetype = 'dashed')+
  geom_col(show.legend = F,color = 'black')+
  theme_pubr()+scale_x_continuous(breaks = 1:10)+
  scale_fill_gradient2(midpoint = 38, high = '#F8766D', low = '#00BFC4', mid = 'white')+
  coord_cartesian(ylim = c(0,80))+xlab('Aromatic/Acidic Mixing')+
  ylab('Functional tAD %')



#ggsave('mixing-hahn.tiff',p3, height = 4, width = 3, dpi = 800, units = 'in')

#saveRDS(p3, 'Figure 2/fig2_mixing_hahn.rds')
