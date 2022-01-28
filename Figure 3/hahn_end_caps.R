##in hahn library can we find end caps
##given terminal D what are the best motifs before that

#packages
library(tidyverse)
library(ggpubr)
library(ggrepel)

hahn <- read_rds('Datasets/Gcn4_random_library.rds')

###just end with 'WWWW' or 'WWWD' or 'WWDW' or 'WWDD"

ends <- c('','W','WW','WWW','WWWW','WWWWW','WWWWWW') %>% gsub('W','[AVILMWYF]',.)
caps <- c('D','DD','DDD','W','WW','WWW')
caps <- gsub('D','[DE]',caps)
caps <- gsub('W','[WYF]',caps)


x <- expand.grid(ends, caps, stringsAsFactors = F)
x$motif <- paste(x$Var1,x$Var2,sep = '')

x$len <- rep(nchar(ends)/5,length(caps))

x$live_n <- 0
x$die_n <- 0

for(i in 1:nrow(x)){
  reggie <- paste('(',x$motif[i],')$',sep = '')
  ro <- grep(reggie, hahn$aa_seq)
  y <- table(hahn$binary_stop[ro])
  x$live_n[i] <- y['live']
  x$die_n[i] <- y['die']
}

x$total <- x$live_n + x$die_n
x$lp <- x$live_n/x$total * 100

x$end_cap <- x$Var2


#for aromatic stretches
x$len <- str_count(x$Var1,'\\[')
x$x <- paste(x$len,' ','[AVILMWYF]', sep = '')
ro <- which(x$x == '0 [AVILMWYF]')
x$x[ro] <- 'any'

x$x <- factor(x$x,
              levels = c('any','1 [AVILMWYF]','2 [AVILMWYF]','3 [AVILMWYF]',
                         '4 [AVILMWYF]','5 [AVILMWYF]','6 [AVILMWYF]'))

x$end_lab <- NA

ro <- which(x$x == '4 [AVILMWYF]')
x$end_lab[ro] <- c('[DE]','[DE][DE]','[DE][DE][DE]',
                   '[WYF]','[WYF][WYF]','[WYF][WYF][WYF]')

p2 <- x %>% filter(len < 5) %>% 
  ggplot(., aes(x, lp, color = end_cap))+
  geom_hline(yintercept = 3.62,linetype = 'dashed')+
  geom_text_repel(aes(label = end_lab),
                  nudge_x = 1, na.rm = T, show.legend = F, segment.color = 'grey',
                  min.segment.length = 0,box.padding = 0.1, hjust = 0,direction = 'y',
                  xlim = c(5,7))+
  geom_point(size = 2, show.legend = F)+
  geom_line(aes(group = end_cap), show.legend = F)+
  theme_pubr()+ylab('Functional tAD %')+xlab('Preceding amino acids')+
  ylim(0,15)+scale_color_manual(values = c('#ff6161','#e82323','#ff0000',
      '#ffb657','#eb9321','#ff9000'))+coord_cartesian(clip = 'off')+
  theme(plot.margin = unit(c(5.5,80,5.5,5.5), 'pt'),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust=1))#+
  #ggtitle(label = 'Gcn4 end caps')

p2

#ggsave('hahn_endcaps.tiff',p2, height = 4, width = 4, dpi = 800, units = 'in')

#saveRDS(p2,'Figure 3/fig3_hahn_endcaps_w-title.rds')
