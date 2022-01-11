##process semi-random-library
library(Biobase)
library(tidyverse)
library(MASS)
hold <- read_rds("Process Semi-Random Library/processed_data_oneMisMatch.rds")

exp_hold <- exprs(hold) %>% as_tibble()
fdata_hold <- fData(hold) %>% as_tibble()

#combine exprs and fdata
semi_random <- cbind(exp_hold, fdata_hold) %>% as_tibble()

#fix T1 and T2 to A0 and H0
colnames(semi_random)[9:10] <- c('A0','H0')

#add pseudocount of 0.5 to library
semi_random[,1:10] <- semi_random[,1:10] + 0.5

#filter for counts above 5 - split into A and H screen
A_screen <- semi_random %>% filter(A0 > 5) %>% dplyr::select(-starts_with('H'))
H_screen <- semi_random %>% filter(H0 > 5) %>% dplyr::select(-matches('A[01234]'))

##add early stop
A_screen$where_stop <- 'A4'
for(i in 1:nrow(A_screen)){
  lib <- A_screen[i,] %>% dplyr::select(starts_with(c('A0','A1','A2','A3','A4'))) %>% unlist()
  hold <- which(lib < 3) %>% names()
  if('A3' %in% hold & 'A4' %in% hold){
    A_screen$where_stop[i] <- 'A3'
  }
  if('A2' %in% hold & 'A3' %in% hold){
    A_screen$where_stop[i] <- 'A2'
  }
  if('A1' %in% hold & 'A2' %in% hold){
    A_screen$where_stop[i] <- 'A1'
  }
}

#normalize counts to total for the day
A_sums <- map(A_screen[,1:5], sum)
H_sums <- map(H_screen[,1:5], sum)

for(i in 1:5){
  A_screen[,i] <- A_screen[,i]/A_sums[[i]] * 1000000
  H_screen[,i] <- H_screen[,i]/H_sums[[i]] * 1000000
}

#set counts to NA if they occur after where_stop
for(i in 1:nrow(A_screen)){
  if(A_screen$where_stop[i] == 'A1'){
    A_screen[i,c("A2","A3","A4")] <- NA
  }
  if(A_screen$where_stop[i] == 'A2'){
    A_screen[i,c("A3","A4")] <- NA
  }
  if(A_screen$where_stop[i] == 'A3'){
    A_screen[i,c("A4")] <- NA
  }
}


#compute log2 fc and run rlm----
regress1 <- function(rw) {
  df <- data.frame(expr = rw,
                   time = 1:4)
  res <- suppressWarnings(rlm(rw ~ 0 + time, data = df))
  return(res)
}

A_screen$A1_log <- log2(A_screen$A1/A_screen$A0)
A_screen$A2_log <- log2(A_screen$A2/A_screen$A0)
A_screen$A3_log <- log2(A_screen$A3/A_screen$A0)
A_screen$A4_log <- log2(A_screen$A4/A_screen$A0)
A_screen$A0_log <- log2(A_screen$A0/A_screen$A0)

#initialize columns
A_screen$A.coef <- 0
A_screen$A.conv <- 0
A_screen$A.rse <- 0
A_screen$A.rsquare <- 0

#do it in batches
ros <- 1:50000
ros <- 50001:100000
ros <- 100001:150000
ros <- 150001:nrow(A_screen)

aa.regress.R <- apply(A_screen[ros,c('A1_log', 'A2_log', 'A3_log', 'A4_log')] , 1, 
                      regress1)

A_screen$A.coef[ros] <- sapply(aa.regress.R, function(i) i$coefficients)
A_screen$A.conv[ros] <- sapply(aa.regress.R, function(i) i$converged)
A_screen$A.rse[ros] <- sapply(aa.regress.R, function(i) i$s)
A_screen$A.rsquare[ros] <- sapply(aa.regress.R, function(i) {
  sse<-sum((i$residuals)^2)
  ssr<- sum((i$fitted.values-mean(i$model[,1]))^2)
  sst<-sse+ssr
  r2<-1-(sse/sst)
  return(r2)
})

rm(aa.regress.R)


#save full library
#saveRDS(A_screen, 'Full_Semi_random_A_screen.rds')

A_screen <- read_rds('Full_Semi_random_A_screen.rds')

#remove sequences that didn't converge, and remove low r.square values
table(A_screen$A.conv) #484 did not converge
table(A_screen$A.rsquare < 0.3) #32524 below 

A_screen <- A_screen %>% filter(A.conv == 1)
A_screen <- A_screen %>% filter(A.rsquare >= 0.3)

##add binary stop
a_stops <- A_screen %>% filter(pattern == 'stop.first')

a_stops <- a_stops[order(a_stops$A.coef, decreasing = T),]

a_cut <- a_stops$A.coef[1:5] %>% mean()

A_screen <- A_screen %>% dplyr::select(aaseq, aaseq.full, pattern, pattern.temp, pattern.aa, A.coef, 
                                       A.rsquare)
colnames(A_screen)[6] <- 'slope'

A_screen$slope - a_cut

ggplot(A_screen, aes(slope))+geom_histogram()

A_screen$binary_stop <- ifelse(A_screen$slope > 0, 'live', 'die')

hold <- A_screen %>% filter(pattern.aa == 'RdmWD5_20')

