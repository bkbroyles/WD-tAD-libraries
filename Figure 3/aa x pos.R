## aa x position graph for fig2

##packages
library(tidyverse)
library(ggpubr)

#load in wd12 data----
wd12 <- read_rds('Datasets/design_library.rds') %>% mutate(
  g_count = str_count(aa_seq, 'G')
) %>% filter(set == 'combinatorial', g_count == 8)

#need 24 aa x pos columns----
x <- expand.grid(c('W','D'), 1:12)
x$col <- paste(x$Var1,x$Var2,sep = '')

#initialize these columns
for(i in 1:nrow(x)){
  wd12[,x$col[i]] <- 0
}

#fill these 24 columns with 1 or 0
for(i in 1:nrow(wd12)){
  for(j in 9:20){
    let <- str_sub(wd12$aa_seq[i],j,j)
    col <- paste(let,j-8,sep = '')
    wd12[i,col] <- 1
  }
}

##Just live% for pos x aa----
hold <- wd12 %>% select('binary_stop',starts_with('W'),starts_with('D'))

###Plug in ML here real quick-----
#run_5ml(hold, my_lambda = 10^seq(-9,-3,length = 50)) %>% mean()


live <- hold %>% filter(binary_stop== 'live')
die <- hold %>% filter(binary_stop == 'die')

x <- map(live[,-1],sum) %>% unlist()
y <- map(die[,-1],sum) %>% unlist()
z <- x + y

yee <- x/z * 100

plotme <- tibble(aa = c(rep('W',12), rep('D',12)),
                 pos = c(1:12,1:12), check = names(yee),
                 live_percent = yee)

plotme$diff <- plotme$live_percent - 33.5

plotme$` ` <- plotme$aa
                 


p2 <- ggplot(plotme, aes(pos, live_percent, color = ` `))+
  geom_hline(yintercept = 33.5, linetype = 'dashed')+
  #geom_smooth(se=F, size = 4, alpha = 0.4, show.legend = F)+
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), size = 4,
              se = F, show.legend = F)+
  geom_point(shape = 21, color = 'black', aes(fill = ` `), size = 2)+
  theme_pubr()+scale_x_continuous(breaks = 1:12)+
  scale_color_manual(values = c('red','orange'))+
  scale_fill_manual(values = c('red','orange'))+
  xlab('Position')+ylab('Functional tAD %')

p2

#ggsave('aa-x-pos.tiff',p2, height = 4, width = 3, dpi = 800, units = 'in')
#saveRDS(p2, 'Figure 3/fig3_aaxpos.rds')

