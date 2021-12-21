##figure 3 scripts
##start with spot1d data

#load packages
library(tidyverse)
library(ggpubr)

x <- read_rds('Datasets/semi_random_SS.rds')

x <- x %>% filter(pattern == 'RdmWD5_20')

x$helical <- str_count(x$ss3,'H')/20 * 100



y <- table(x$helical, x$binary_stop)
plotme <- tibble(helical = row.names(y), die = y[,'die'], live = y[,'live'])
plotme$total <- plotme$live + plotme$die

plotme$lp <- plotme$live/plotme$total * 100

plotme$helical <- factor(plotme$helical,
                            levels = seq(0,85,5))

plotme$`Functional \n tAD %` <- plotme$lp

ggplot(plotme, aes(helical, lp, fill = `Functional \n tAD %`))+
  geom_hline(yintercept = 28.9, linetype = 'dashed')+
  geom_col(color = 'black')+theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Predicted Helical %')+
  ylab('Functional tAD %')+
  scale_fill_gradient2(midpoint = 28.9)

