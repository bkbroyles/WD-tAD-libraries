#align tetrapeptides
plot1 <- read_rds('Figure 2/fig2_tetrapep_wd12.rds')
plot2 <- read_rds('Figure 2/fig2_tetrapep-hahn.rds')

p1 <- ggarrange(plot1,plot2,ncol = 2, widths = c(1.2,2))


#ggsave('tetrapep-both.tiff',p1, height = 4, width = 7, dpi = 800, units = 'in')
