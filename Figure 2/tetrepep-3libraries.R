##tetrapep plot
##data generated in tetrapep x pos in hahn script
##   and tetrapep x position script

##could kornberg be added to this plot
korn <- read_rds('Datasets/Yeast_TF_tiled.rds')
ro <- grep('tiles',korn$set)
korn <- korn[ro,]

aro <- '[WYF]'
acid <- '[DE]'
opts <- c(aro,acid)

combos <- expand.grid(opts,opts,opts,opts,stringsAsFactors = F)

combos$pep <- paste(combos$Var1,combos$Var2,combos$Var3,combos$Var4,sep='')

combos$live_count <- 0
combos$die_count <- 0

for(i in 1:nrow(combos)){
  korn[,combos$pep[i]] <- str_count(korn$aa_seq,combos$pep[i])
}

live <- korn %>% filter(binary_stop=='live')
die <- korn %>% filter(binary_stop=='die')

l_sums2 <- apply(live[,10:25],2,sum)
d_sums2 <- apply(die[,10:25],2,sum)

korn_tetrapep <- tibble(pep = names(l_sums2), live = l_sums2, die = d_sums2) %>% 
  mutate(
    total = live + die,
    lp = live/total * 100
  )

korn_tetrapep$pep <- gsub('\\[','',korn_tetrapep$pep)
korn_tetrapep$pep <- gsub('\\]','',korn_tetrapep$pep)
korn_tetrapep$pep <- gsub('[FEY]','',korn_tetrapep$pep)

korn_tetrapep$set <- 'korn'

korn_tetrapep <- korn_tetrapep %>% select(pep,lp,set)

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

##whats live percent for kornberg
sum(korn_tetrapep$live)
sum(korn_tetrapep$total)
184/1085 #17%
