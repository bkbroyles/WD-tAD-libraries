##align figures together
library(ggpubr)
library(tidyverse)

plot1 <- read_rds('Figure 2/fig2_mixing_wd12.rds')
plot2 <- read_rds('Figure 2/fig2_mixing_hahn.rds')

p1 <- ggarrange(plot1,plot2,ncol = 1)

#ggsave('mixing-both.tiff',p1, height = 4, width = 3, dpi = 800, units = 'in')
  