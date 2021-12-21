##looking for tetrapeptides in hahn and kornberg
library(tidyverse)

hahn <- read_rds('Datasets/Gcn4_random_library.rds')

######motif x position----
aro <- '[WYF]'
acid <- '[DE]'
x <- expand.grid(c(acid,aro),c(acid,aro),c(acid,aro),c(acid,aro))
x$motif <- paste(x$Var1,x$Var2,x$Var3,x$Var4,sep='')

#open position columns
for(i in 1:27){
  for(j in 1:nrow(x)){
    col <- paste('pos_',i,'_',x$motif[j],sep='')
    hahn[,col] <- 0
  }
}

for(i in 1:27){
  seqs <- str_sub(hahn$aa_seq,i,i+3)
  for(j in 1:nrow(x)){
    ro <- grep(x$motif[j], seqs)
    col <- paste('pos_',i,'_',x$motif[j],sep='')
    hahn[ro,col] <- 1
  }
}

#split into live and die
live <- hahn %>% filter(binary_stop == 'live')
die <- hahn %>% filter(binary_stop == 'die')

l_sums <- apply(live[,5:436],2,sum)
d_sums <- apply(die[,5:436],2,sum)

df <- l_sums/(d_sums + l_sums) * 100
x <- tibble(name = names(df), lp = df)

x$pos <- str_sub(x$name,1,6)
x$pos <- gsub('(_)$','',x$pos)
x$pos <- gsub('pos_','',x$pos)
x$pos <- as.numeric(x$pos)

x$motif <- str_sub(x$name,7,-1)
x$motif <- gsub('_','',x$motif)

x$yax <- gsub('\\[','',x$motif)
x$yax <- gsub('\\]','',x$yax)
x$yax <- gsub('E','',x$yax)
x$yax <- gsub('[YF]','',x$yax)

x$yax <- factor(x$yax,
                 levels = c('DDWW','DWWD','WWDW','DWWW','WDDW',
                            'WDWW','DWDW','WWWD','WDWD','WWDD',
                            'WWWW','DDDW','DWDD','DDWD','WDDD','DDDD')) %>% fct_rev()

x$` ` <- x$lp

p1 <- ggplot(x, aes(pos, yax, fill = ` `))+
  geom_tile(color = 'grey32', show.legend = T)+
  scale_fill_gradient2(midpoint = 16.8, high = '#00BFC4', low = '#F8766D', mid = 'white')+
  scale_x_continuous(breaks = seq(1,27,2))+
  theme_pubr(legend = 'top')+ylab('Motif')+xlab('Starting Position')

#ggsave('tetrapep-x-pos-hahn.tiff',p1, height = 4, width = 2.4, dpi = 800, units = 'in')
saveRDS(p1, 'tetrapep-hahn.rds')

#sum(l_sums)/(sum(l_sums)+sum(d_sums)) midpoint should be 16.8

#make bonus plot
name_hold <- names(l_sums)
name_hold <- gsub('pos_','',name_hold)
name_hold <- gsub('[0123456789]','',name_hold)
name_hold <- gsub('_','',name_hold)
name_hold <- gsub('\\[','',name_hold)
name_hold <- gsub('\\]','',name_hold)
name_hold <- gsub('[EFY]','',name_hold)

names(l_sums) <- name_hold

df <- tibble(pep = unique(name_hold), count = 0)
for(i in 1:nrow(df)){
  df$count[i] <- sum(l_sums[which(names(l_sums) == df$pep[i])])
}

#same for d_sums
name_hold <- names(d_sums)
name_hold <- gsub('pos_','',name_hold)
name_hold <- gsub('[0123456789]','',name_hold)
name_hold <- gsub('_','',name_hold)
name_hold <- gsub('\\[','',name_hold)
name_hold <- gsub('\\]','',name_hold)
name_hold <- gsub('[EFY]','',name_hold)

names(d_sums) <- name_hold

df$die_count <- 0
for(i in 1:nrow(df)){
  df$die_count[i] <- sum(d_sums[which(names(d_sums) == df$pep[i])])
}

df <- df %>% mutate(
  total = die_count + count,
  lp = count/total * 100,
  set = 'hahn'
)

df <- df %>% select(pep, lp, set)
tetra_df <- tetra_df %>% select(pep, lp, set)

tetra_df <- rbind(tetra_df, df)

