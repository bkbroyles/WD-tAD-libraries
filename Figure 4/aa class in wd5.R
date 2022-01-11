##Individual amino acids help
wd5 <- read_rds('Datasets/semi_random_SS.rds') %>% 
  filter(pattern == 'RdmWD5_20')

x <- str_split(wd5$aa_seq,'')

wd5$module <- ''

for(i in 1:length(x)){
  seq <- x[[i]][c(2,4,6,8,10,12,14,16,18,20)]
  mod <- paste(seq, collapse = '')
  wd5$module[i] <- mod
}

groups <- list(aliphatic = c('A','V','I','L','M'),
               aromatic = c('W','Y','F'),
               polar = c('S','T','N','Q'),
               C = 'C',G = 'G',P = 'P',
               basic = c('R','H','K'),
               acidic = c('D','E'))

for(i in 1:length(groups)){
  reggie <- paste(groups[[i]],collapse = '') %>% paste('[',.,']',sep='')
  wd5[,names(groups)[[i]]] <- str_count(wd5$module,reggie)
}

plotme <- tibble(aa = rep(names(groups),each = 11),
                 count = rep(0:10, 8), live_n = 0, die_n = 0)

for(i in 1:nrow(plotme)){
  ro <- which(wd5[,plotme$aa[i]]==plotme$count[i])
  lib <- wd5[ro,]
  y <- table(lib$binary_stop)
  plotme$live_n[i] <- y['live']
  plotme$die_n[i] <- y['die']
}

plotme$total <- plotme$live_n + plotme$die_n
plotme$lp <- plotme$live_n/plotme$total * 100

library(ggrepel)

##add label to end of lines
plotme$label <- NA
#just do it manually
plotme$label[9] <- 'aliphatic'
plotme$label[18] <- 'aromatic'
plotme$label[31] <- 'polar'
plotme$label[37] <- 'C'
plotme$label[49] <- 'G'
plotme$label[61] <- 'P'
plotme$label[72] <- 'basic'
plotme$label[82] <- 'acidic'

p2 <- plotme %>% filter(total > 9) %>% 
ggplot(., aes(count, lp, color = aa))+
  geom_hline(yintercept = 28.9, linetype = 'dashed')+
  geom_line(size = 1.5)+
  theme_pubr()+
  scale_x_continuous(breaks = 0:9)+
  geom_label_repel(aes(label = label), show.legend = F, min.segment.length = 0,
                  segment.color = 'grey', hjust = 0, direction = 'y',nudge_x = 1,
                  box.padding = 0.3, xlim = c(4,8))+
  ylab('Functional tAD %')+
  scale_color_manual(values = c("#F8766D", "#F0E442", "#E69F00", 
                                "#0072B2", "#56B4E9", "#999999", 
                                "#CC79A7", "#009E73"))

p2

#ggsave('wd5_aa_counts.tiff',p2, height = 4, width = 4, dpi = 800, units = 'in')
  