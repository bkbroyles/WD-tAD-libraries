##1 w and 1 d error bar plots

#libraries
library(tidyverse)
library(MASS)
library(ggpubr)

#load data----
hold <- read_rds('Datasets/design_library.rds')
ci_95 <- read_rds('Datasets/CI_95_dlib.rds')

hold <- hold %>% mutate(G = str_count(aa_seq,'G'),
                        D = str_count(aa_seq,'D'),
                        W = str_count(aa_seq,'W'))

hold <- hold %>% filter(G == 18, D == 1, W == 1)

hold$aa_seq

#delete later
hold$rep <- c(1,1,1,rep(2,9)) %>% as.factor()
hold$index <- c(1,2,3:12)

hold$low <- NA
hold$high <- NA

for(i in 1:nrow(hold)){
  ro <- which(ci_95$id_number == hold$id_number[i])
  hold$low[i] <- ci_95$low[ro]
  hold$high[i] <- ci_95$high[ro]
}

plot_1w1d <- ggplot(hold, aes(index, slope, color = rep))+
  geom_hline(yintercept = 0, size = 0.8)+
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, show.legend = F)+
  geom_point(size = 2, show.legend = F)+
  ylim(-1,1.8)+
  theme_pubr()+
  scale_x_continuous(breaks = 1:12)+ylab(' ')+
  scale_color_manual(values = c("#F8766D","#00BFC4"))+xlab("Index")

#ggsave('1w1d_CI.tiff', plot_1w1d, height = 4, width = 3, units = 'in', dpi = 400)
