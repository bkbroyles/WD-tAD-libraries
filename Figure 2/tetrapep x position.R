##tetrapeptides x position
##tetrapeptides overall average
library(tidyverse)
library(ggpubr)

#load in lib data----
wd12 <- read_rds('Datasets/design_library.rds') %>% mutate(
  g_count = str_count(aa_seq, 'G')
) %>% filter(set == 'combinatorial', g_count == 8) 

##make tetrapeptide combos columns----
x <- c('W','D')

hold <- expand.grid(x,x,x,x)
hold$pep <- paste(hold$Var1,hold$Var2,hold$Var3,hold$Var4,sep = '')

for(i in 1:9){
  for(j in hold$pep){
    col <- paste('pos_',i,'_',j, sep = '')
    wd12[,col] <- 0
  }
}

## fill columns
for(i in 1:nrow(wd12)){
  for(j in 1:9){
    x <- str_sub(wd12$aa_seq[i], j+8, j+11)
    col <- paste('pos_',j,'_',x,sep='')
    wd12[i,col] <- 1
  }
}














#split wd12 into live and die
live <- wd12 %>% filter(binary_stop == 'live')
die <- wd12 %>% filter(binary_stop == 'die')

for(i in 1:9){
  x <- table(str_sub(live$aa_seq,i+8,i+11))
  y <- table(str_sub(die$aa_seq,i+8,i+11))
  l_col <- paste('live_',i, sep = '')
  d_col <- paste('die_',i, sep = '')
  for(j in 1:nrow(hold)){
    hold[j,l_col] <- x[which(names(x) == hold$pep[j])]
    hold[j,d_col] <- y[which(names(y) == hold$pep[j])]
  }
}

##calculate totals and live_percent
for(i in 1:9){
  l_col <- paste('live_',i,sep='')
  d_col <- paste('die_',i,sep='')
  t_col <- paste('total_',i,sep='')
  lp_col <- paste('live_percent_',i,sep='')
  hold[,t_col] <- hold[,l_col] + hold[,d_col]
  hold[,lp_col] <- hold[,l_col]/hold[,t_col] * 100
}


##add live avg to plot_me

plot_me <- hold %>% select(pep,starts_with('live_per'))

plot_me <- pivot_longer(plot_me, cols = c(paste('live_percent_',1:9, sep = '')))

plot_me$position <- gsub('live_percent_','',plot_me$name) %>% as.numeric()

seqs <- c('DDWW','DWWD','WWDW','DWWW','WDDW',
          'WDWW','DWDW','WWWD','WDWD','WWDD',
          'WWWW','DDDW','DWDD','DDWD','WDDD','DDDD')

plot_me$pep <- factor(plot_me$pep, levels = seqs) %>% fct_rev()
plot_me$` ` <- plot_me$value

p1 <- plot_me %>% 
ggplot(., aes(position, pep, fill = ` `))+
  geom_tile(color = 'grey32', show.legend = T)+
  scale_fill_gradient2(midpoint = 33.5, high = '#00BFC4', low = '#F8766D', mid = 'white')+
  scale_x_continuous(breaks = 1:9)+
  theme_pubr(legend = 'top')+ylab('Motif')+xlab('Starting Position')

##save plot here
#ggsave('tetrapep-wd12.tiff',p1, height = 4, width = 3, dpi = 800, units = 'in')
#saveRDS(p1, 'tetrapep_wd12.rds')

#live percent calculation
x <- hold %>% select(starts_with('live'), -contains('percent'))
y <- hold %>% select(starts_with('total'))

sum(x %>% unlist()) #11979
sum(y %>% unlist()) #35712

11979/35712 #33.5 midpoint


###bonus plot on tetrapeptides
##add all positions together for over live%
hold <- hold %>% mutate(
  live_all = live_1 + live_2 + live_3 +live_4 +
    live_5 + live_6 + live_7 + live_8 + live_9,
  total_all = total_1 + total_2 + total_3 +total_4 +
    total_5 + total_6 + total_7 + total_8 + total_9
)

tetra_df <- hold %>% select(pep, live_all, total_all)
tetra_df$lp <- tetra_df$live_all/tetra_df$total_all * 100
tetra_df$set <- 'wd12'
