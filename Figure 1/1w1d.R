##need to add WD repeats and DW repeats to fig 1
library(tidyverse)
library(ggpubr)

lib<- read_rds('Datasets/design_library.rds') %>% 
  mutate(
    g_count = str_count(aa_seq, 'G')
  )

lib <- lib %>% mutate(
  w = str_count(aa_seq, 'W'),
  d = str_count(aa_seq, 'D'),
  tots = g_count + w + d
)

hold <- lib %>% filter(w == 1, d == 1, tots == 20)
hold$index <- c(1,1,2:11)
hold$` ` <- c('DW','DW','DW', rep('WD',9))

plot_1w1d <- ggplot(hold, aes(index, slope, color = ` `))+
  geom_hline(yintercept = 0, size = 0.8)+
  geom_point(size = 2, show.legend = F)+
  ylim(-1,1.8)+
  theme_pubr()+
  scale_x_continuous(breaks = 1:12)+ylab(' ')+
  scale_color_manual(values = c("#F8766D","#00BFC4"))+xlab("Index")

#ggsave('1w1d.tiff', plot_1w1d, height = 4, width = 3, units = 'in', dpi = 800)


#saveRDS(plot_1w1d, 'Figure 1/fig1_1w1d.rds')

###ADDING BIOREP PERFORMANCE
 
ids <- hold$id_number
hold2 <- hold

counts <- read_rds('Process Design Library/Intermediate Datasets/design_library_raw_counts.rds')$A_screen

hold <- counts %>% filter(id_number %in% ids)

###back to scripts

df <- tibble(id_number = rep(hold$id_number, each = 5), aa_seq = rep(hold$aa_seq, each = 5),
             time0=0, time1=0, time2=0, time3=0, time4=0, rep = rep(1:5, nrow(hold)))

#big loop
for(i in 1:nrow(hold)){
  df$time0[(i*5)-4] <- hold$A0_bc1[i]
  df$time0[(i*5)-3] <- hold$A0_bc2[i]
  df$time0[(i*5)-2] <- hold$A0_bc3[i]
  df$time0[(i*5)-1] <- hold$A0_bc4[i]
  df$time0[(i*5)] <- hold$A0_bc5[i]
  
  df$time1[(i*5)-4] <- hold$A1_bc1[i]
  df$time1[(i*5)-3] <- hold$A1_bc2[i]
  df$time1[(i*5)-2] <- hold$A1_bc3[i]
  df$time1[(i*5)-1] <- hold$A1_bc4[i]
  df$time1[(i*5)] <- hold$A1_bc5[i]
  
  df$time2[(i*5)-4] <- hold$A2_bc1[i]
  df$time2[(i*5)-3] <- hold$A2_bc2[i]
  df$time2[(i*5)-2] <- hold$A2_bc3[i]
  df$time2[(i*5)-1] <- hold$A2_bc4[i]
  df$time2[i*5] <- hold$A2_bc5[i]
  
  df$time3[(i*5)-4] <- hold$A3_bc1[i]
  df$time3[(i*5)-3] <- hold$A3_bc2[i]
  df$time3[(i*5)-2] <- hold$A3_bc3[i]
  df$time3[(i*5)-1] <- hold$A3_bc4[i]
  df$time3[(i*5)] <- hold$A3_bc5[i]
  
  df$time4[(i*5)-4] <- hold$A4_bc1[i]
  df$time4[(i*5)-3] <- hold$A4_bc2[i]
  df$time4[(i*5)-2] <- hold$A4_bc3[i]
  df$time4[(i*5)-1] <- hold$A4_bc4[i]
  df$time4[(i*5)] <- hold$A4_bc5[i]
}

#need to find log2 transformation of these sequences
df <- df %>% mutate(
  time1 = log2(time1/time0),
  time2 = log2(time2/time0),
  time3 = log2(time3/time0),
  time4 = log2(time4/time0),
)


df$slope = NA
df$converged = NA

library(MASS)

for(i in 1:nrow(df)){
  df3 <- tibble(expr = df[i,c('time1','time2','time3','time4')] %>% unlist(), time = 1:4)
  res <- rlm(expr ~ 0 + time, data = df3)
  df$slope[i] <- res$coefficients
  df$converged[i] <- res$converged
}

df$slope <- df$slope + 1.05

df$label <- paste(df$id_number,'-',df$aa_seq,sep='')

###save columns to merge with agggregate slope info
df2 <- df %>% dplyr::select(id_number, aa_seq, rep, slope)

hold2 <- hold2 %>% dplyr::select(id_number, aa_seq, slope) %>% mutate(
  rep = 6
)

df2 <- rbind(df2, hold2)

df2$color <- ifelse(df2$rep==6,1,0) %>% as.factor()

df2 %>% filter(rep != 6) %>% 
ggplot(., aes(id_number, slope, color = color))+
  geom_hline(yintercept = 0, size = 0.8)+
  geom_boxplot(aes(group = id_number))+
  ylim(-2,2.5)+
  theme_pubr()
