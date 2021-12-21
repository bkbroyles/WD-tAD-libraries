##align figures together
plot1 <- read_rds('mixing_wd12.rds')
plot2 <- read_rds('mixing_hahn.rds')

library(ggpubr)
library(tidyverse)

p1 <- ggarrange(plot1,plot2,ncol = 1)

ggsave('mixing-both.tiff',p1, height = 4, width = 3, dpi = 800, units = 'in')
  