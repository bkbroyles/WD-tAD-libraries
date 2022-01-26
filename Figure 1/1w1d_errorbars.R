##1 w and 1 d error bar plots

#libraries
library(tidyverse)
library(MASS)
library(ggpubr)

#load data----
hold <- read_rds('Process Design Library/Intermediate Datasets/A_screen_counts.rds')

hold <- hold %>% mutate(G = str_count(aa_seq,'G'),
                                D = str_count(aa_seq,'D'),
                                W = str_count(aa_seq,'W'))

#get day0sums -> day4sums----
day0sum <- hold %>% dplyr::select(starts_with('A0')) %>% unlist() %>% sum()
day1sum <- hold %>% dplyr::select(starts_with('A1')) %>% unlist() %>% sum()
day2sum <- hold %>% dplyr::select(starts_with('A2')) %>% unlist() %>% sum()
day3sum <- hold %>% dplyr::select(starts_with('A3')) %>% unlist() %>% sum()
day4sum <- hold %>% dplyr::select(starts_with('A4')) %>% unlist() %>% sum()

#only grab 1w and 1d sequences----
#this piece could be changed to investigate other sequence sets
hold <- hold %>% filter(G == 18, D == 1, W == 1)

##each sequence_id should have 5 rows, 1 for each biorep----
hold2 <- tibble(rows = 1:120, id_number = '', aa_seq = '',
                day0=0, day1=0, day2=0, day3=0, day4=0)

hold2$id_number <- rep(hold$id_number, each = 5)
hold2$aa_seq <- rep(hold$aa_seq, each = 5)

hold2$rep <- rep(1:5, 24)

hold2 <- hold2[,-1]

#add hold data to hold2----
for(i in 1:nrow(hold)){
  hold2[((i*5)-4),3:7] <- hold[i,c('A0_bc1','A1_bc1','A2_bc1','A3_bc1','A4_bc1')]
  hold2[((i*5)-3),3:7] <- hold[i,c('A0_bc2','A1_bc2','A2_bc2','A3_bc2','A4_bc2')]
  hold2[((i*5)-2),3:7] <- hold[i,c('A0_bc3','A1_bc3','A2_bc3','A3_bc3','A4_bc3')]
  hold2[((i*5)-1),3:7] <- hold[i,c('A0_bc4','A1_bc4','A2_bc4','A3_bc4','A4_bc4')]
  hold2[(i*5),3:7] <- hold[i,c('A0_bc5','A1_bc5','A2_bc5','A3_bc5','A4_bc5')]
}

#do these bioreps meet >=5 criteria----
hold2$include <- ifelse(hold2$day0 >= 5, 'yes', 'no')

#need to apply where stop----
hold2$where_stop <- 'day4'
for(i in 1:nrow(hold2)){
  lib <- hold2[i,] %>% dplyr::select(starts_with('day')) %>% unlist()
  x <- which(lib < 3) %>% names()
  if('day3' %in% x & 'day4' %in% x){
    hold2$where_stop[i] <- 'day3'
  }
  if('day2' %in% x & 'day3' %in% x){
    hold2$where_stop[i] <- 'day2'
  }
  if('day1' %in% x & 'day2' %in% x){
    hold2$where_stop[i] <- 'day1'
  }
}

#transform counts into tpm----
hold2$day0_tpm <- hold2$day0/day0sum * 1000000
hold2$day1_tpm <- hold2$day1/day1sum * 1000000
hold2$day2_tpm <- hold2$day2/day2sum * 1000000
hold2$day3_tpm <- hold2$day3/day3sum * 1000000
hold2$day4_tpm <- hold2$day4/day4sum * 1000000

##clean this up to look easier
hold3 <- hold2 %>% dplyr::select(id_number, aa_seq, contains(c('day0_tpm',
                                                              'day1_tpm',
                                                              'day2_tpm',
                                                              'day3_tpm',
                                                              'day4_tpm')),
                                 include, where_stop)


#log2fc transformation----
hold3$day1 <- log2(hold3$day1_tpm/hold3$day0_tpm)
hold3$day2 <- log2(hold3$day2_tpm/hold3$day0_tpm)
hold3$day3 <- log2(hold3$day3_tpm/hold3$day0_tpm)
hold3$day4 <- log2(hold3$day4_tpm/hold3$day0_tpm)

#find full slopes, then use where_stop to find adjusted slopes

#Regression to find slope----
regress1 <- function(rw) {
  df <- data.frame(expr = rw,
                   time = 1:4)
  res <- suppressWarnings(rlm(rw ~ 0 + time, data = df))
  return(res)
}

models <- apply(hold3[,c("day1","day2","day3","day4")],1,regress1)

hold3$slope <- sapply(models, function(i) i$coefficients)
hold3$converged <- sapply(models, function(i) i$converged)
hold3$mse <- 0
for(i in 1:nrow(hold3)){
  na_length <- which(is.na(hold3[i,c("day1","day2","day3","day4")])) %>% length()
  sse <- sum((models[[i]]$wresid)^2)
  hold3$mse[i] <- sse/(4 - na_length) 
}

#slopes with where stop filter----
for(i in 1:nrow(hold3)){
  x <- hold3$where_stop[i]
  if(x == 'day1'){
    hold3[i,c('day2','day3','day4')]<-NA
  }
  if(x == 'day2'){
    hold3[i,c('day3','day4')]<-NA
  }
  if(x == 'day3'){
    hold3[i,'day4']<-NA
  }
}

models <- apply(hold3[,c("day1","day2","day3","day4")],1,regress1)

hold3$slope2 <- sapply(models, function(i) i$coefficients)
hold3$converged2 <- sapply(models, function(i) i$converged)
hold3$mse2 <- 0
for(i in 1:nrow(hold3)){
  na_length <- which(is.na(hold3[i,c("day1","day2","day3","day4")])) %>% length()
  sse <- sum((models[[i]]$wresid)^2)
  hold3$mse2[i] <- sse/(4 - na_length) 
}

##add stop codon cut off to slopes
hold3$slope <- hold3$slope + 1.05
hold3$slope2 <- hold3$slope2 + 1.05

#grab data from design library to filter hold3---- 
dlib <- read_rds('Datasets/design_library.rds')

ro <- which(dlib$id_number %in% hold3$id_number)
dlib <- dlib[ro,]

ids <- dlib$id_number

#remove sequences that did not make final cut
hold4 <- hold3 %>% filter(id_number %in% ids)

hold4 <- hold4 %>% dplyr::select(id_number, aa_seq, slope, slope2, include)

#add dlib data to hold4
hold4$rep <- rep(1:5, 12)

x <- tibble(id_number = dlib$id_number, aa_seq = dlib$aa_seq,
            slope = dlib$slope, slope2 = dlib$slope, include = 'yes',
            rep = 6)

hold4 <- rbind(hold4, x)

#change rep to signify biorep vs final
hold4$rep <- ifelse(hold4$rep == 6, 1, 0) %>% as.factor()


#ggplot(hold4, aes(id_number, slope2, color = rep))+geom_point()
  
#add index to make it more like 1w 1d plot
hold4$index <- 0

hold4$index[grep('GGGGGGGGGGGGGGGGGGDW', hold4$aa_seq)] <- 1
hold4$index[grep('GGGGGGGGGDWGGGGGGGGG', hold4$aa_seq)] <- 2
hold4$index[grep('GGGGGGGGGGGGGGGGGGWD', hold4$aa_seq)] <- 3
hold4$index[grep('GGGGGGGGGGGGGGGGGWGD', hold4$aa_seq)] <- 4
hold4$index[grep('GGGGGGGGGGGGGGGGWGGD', hold4$aa_seq)] <- 5
hold4$index[grep('GGGGGGGGGGGGGGGWGGGD', hold4$aa_seq)] <- 6
hold4$index[grep('GGGGGGGGGGGGGGWGGGGD', hold4$aa_seq)] <- 7
hold4$index[grep('GGGGGGGGGGGGGWGGGGGD', hold4$aa_seq)] <- 8
#7 g sequence does not make it to final screen
hold4$index[grep('GGGGGGGGGGGGWGGGGGGD', hold4$aa_seq)] <- 9
hold4$index[grep('GGGGGGGGGGWGGGGGGGGD', hold4$aa_seq)] <- 10
hold4$index[grep('GGGGGGGGGWGGGGGGGGGD', hold4$aa_seq)] <- 11

#select max and min for each biorep per index for error bars
#using slope without early cutoff, it reflects final results better
hold <- hold4 %>% filter(rep == 1)
hold4 <- hold4 %>% filter(rep == 0)

hold$max <- 0
hold$min <- 0

for(i in 1:nrow(hold)){
  slopes <- hold4 %>% filter(index == hold$index[i]) %>% 
    dplyr::select(slope) %>% unlist()
  hold$max[i] <- max(slopes)
  hold$min[i] <- min(slopes)
}

hold$rep[1:3] <- 0


plot_1w1d <- ggplot(hold, aes(index, slope, color = rep))+
  geom_hline(yintercept = 0, size = 0.8)+
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.2, show.legend = F)+
  geom_point(size = 2, show.legend = F)+
  ylim(-1,1.8)+
  theme_pubr()+
  scale_x_continuous(breaks = 1:12)+ylab(' ')+
  scale_color_manual(values = c("#F8766D","#00BFC4"))+xlab("Index")

#ggsave('1w1d.tiff', plot_1w1d, height = 4, width = 3, units = 'in', dpi = 800)
