##compare 3 libraries
##data generated and saved in tetrapep x pos in hahn script
##   and tetrapep x position script

#packages
library(tidyverse)
library(ggpubr)

#get data from artificial tAD data
korn <- read_rds('Datasets/Yeast_TF_tiled.rds')

#only grab the tiled TF subset of this library
ro <- grep('tiles',korn$set)
korn <- korn[ro,]

#set up tetrapeptide search
aro <- '[WYF]'
acid <- '[DE]'
opts <- c(aro,acid)

combos <- expand.grid(opts,opts,opts,opts,stringsAsFactors = F)

#make regex for searching for these tetrapeptides
combos$pep <- paste(combos$Var1,combos$Var2,combos$Var3,combos$Var4,sep='')

#count tetrapeptides in each sequence
for(i in 1:nrow(combos)){
  korn[,combos$pep[i]] <- str_count(korn$aa_seq,combos$pep[i])
}

#split into live and die set from artificial tAD dataset.
live <- korn %>% filter(binary_stop=='live')
die <- korn %>% filter(binary_stop=='die')

l_sums2 <- apply(live[,10:25],2,sum)
d_sums2 <- apply(die[,10:25],2,sum)

#Create table based on counts and live percent per tetrapeptide
korn_tetrapep <- tibble(pep = names(l_sums2), live = l_sums2, die = d_sums2) %>% 
  mutate(
    total = live + die,
    lp = live/total * 100
  )

#Clean up the regex features
korn_tetrapep$pep <- gsub('\\[','',korn_tetrapep$pep)
korn_tetrapep$pep <- gsub('\\]','',korn_tetrapep$pep)
korn_tetrapep$pep <- gsub('[FEY]','',korn_tetrapep$pep)

#add set for combining all three libraries
korn_tetrapep$set <- 'korn'

#simplify table
korn_tetrapep <- korn_tetrapep %>% dplyr::select(pep,lp,set)

#tetrapap_3libraries is the saved tetrapeptide data from 
# gcn4 library and gal4 library, they are calculated in 
# other tetrapeptide scripts
tetra_df <- read_rds('Figure 2/tetrapep_3libraries.rds')
tetra_df <- rbind(tetra_df,korn_tetrapep)

seqs <- c('DDWW','DWWD','WWDW','DWWW','WDDW',
  'WDWW','DWDW','WWWD','WDWD','WWDD',
  'WWWW','DDDW','DWDD','DDWD','WDDD','DDDD')

tetra_df$pep <- factor(tetra_df$pep,
                       levels = seqs) %>% fct_rev()

tetra_df$library <- tetra_df$set

ro <- which(tetra_df$library == 'korn')
tetra_df$library[ro] <- 'kornberg'

tetra_df$pep2 <- gsub('W','WYF',tetra_df$pep)
tetra_df$pep2 <- gsub('D','DE',tetra_df$pep2)

tetra_df$pep2 <- gsub('ED','E-D',tetra_df$pep2)
tetra_df$pep2 <- gsub('FW','F-W',tetra_df$pep2)
tetra_df$pep2 <- gsub('EW','E-W',tetra_df$pep2)
tetra_df$pep2 <- gsub('FD','F-D',tetra_df$pep2)

p1 <- ggplot(tetra_df, aes(lp, pep2, color = library))+
  geom_point()+theme_pubr()+scale_color_manual(values = c('red','green','blue'))+
  geom_vline(xintercept = 16.8, color = 'red', linetype = 'dashed')+
  geom_vline(xintercept = 17, color = 'green', linetype = 'dashed')+
  geom_vline(xintercept = 33.5, color = 'blue', linetype = 'dashed')+
  ylab('Motif')+xlab('Functional tAD %')

#ggsave('tetrapep-3libraries.tiff',p1, height = 4, width = 5, dpi = 800, units = 'in')

#CHANGING WAY TO VIEW THIS GRAPH
tetra_df$avg <- 0
x <- expand.grid(c('W','D'),c('W','D'),
                 c('W','D'),c('W','D'),
                 stringsAsFactors = F)
x$motif <- paste(x$Var1,x$Var2,x$Var3,x$Var4,sep='')

for(i in 1:nrow(x)){
  ro <- which(tetra_df$pep == x$motif[i])
  tetra_df$avg[ro] <- sum(tetra_df$lp[ro])/3
}

hold <- tetra_df$pep[order(tetra_df$avg)]
these_ones <- hold[c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46)]

tetra_df$pep <- factor(tetra_df$pep,
                       levels = these_ones)

#fix library names
ro <- which(tetra_df$library=='kornberg')
tetra_df$library[ro] <- 'Artificial tAD'

ro <- which(tetra_df$library=='hahn')
tetra_df$library[ro] <- 'Gcn4 tAD'

ro <- which(tetra_df$library=='Gal4 tAD (wd12)')
tetra_df$library[ro] <- 'Gal4 tAD'

tetra_df$` ` <- tetra_df$library

p1 <- ggplot(tetra_df, aes(lp, pep, color = ` `))+
  geom_smooth(aes(group = ` `), se = F, show.legend = F, method = 'lm', 
              linetype = 'dashed', alpha = 1, size = 1.2)+
  geom_point(aes(shape = ` `), size = 2, alpha = 0.9)+theme_pubr()+
  scale_color_manual(values = c('#ff707c','#707cff','#56b820'))+
  ylab('Motif')+xlab('Functional tAD %')



#ggsave('Figure 2/tetrapep-3libraries.tiff',p1, height = 4, width = 4.5, dpi = 800, units = 'in')

