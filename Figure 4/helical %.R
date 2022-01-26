##figure 3 scripts
##start with spot1d data

#load packages
library(tidyverse)
library(ggpubr)

x <- read_rds('Datasets/semi_random_library2.rds')

x <- x %>% filter(pattern.x == 'RdmWD5_20')

x$helical <- str_count(x$ss3,'H')/20 * 100



y <- table(x$helical, x$binary_stop.x)
plotme <- tibble(helical = row.names(y), die = y[,'die'], live = y[,'live'])
plotme$total <- plotme$live + plotme$die

plotme$lp <- plotme$live/plotme$total * 100

plotme$helical <- factor(plotme$helical,
                            levels = seq(0,85,5))

plotme$`Functional \n tAD %` <- plotme$lp

p1 <- ggplot(plotme, aes(helical, lp, fill = `Functional \n tAD %`))+
  geom_hline(yintercept = 19.6, linetype = 'dashed')+
  geom_col(color = 'black')+theme_pubr()+
  xlab('Predicted Helical %')+
  ylab('Functional tAD %')+
  scale_fill_gradient2(midpoint = 19.6, high = "#F8766D", low = "#00BFC4")+
  scale_x_discrete(breaks = seq(0,80,10))

ggsave('Figure 4/helicity_percent.tiff',p1, height = 4, width = 3, dpi = 800, units = 'in')
