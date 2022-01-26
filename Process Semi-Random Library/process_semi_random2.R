##process_semi_random
library(Biobase)
library(tidyverse)
library(MASS)

Processed.Data.One.MisMatch <- read_rds("Process Semi-Random Library/processed_data_oneMisMatch.rds")
Counts.Data.One.MisMatch<-exprs(Processed.Data.One.MisMatch)

weedle<-fData(Processed.Data.One.MisMatch)
weedle<-cbind(weedle,Counts.Data.One.MisMatch)
weedle$A0<-weedle$T1
weedle$H0<-weedle$T2
weedle<-weedle[,-c(17,18)]
col_order<-c("Group.1","dnaseq","pattern.temp","aaseq.full","aaseq","pattern","pattern.aa","A0","A1","A2","A3","A4","H0","H1","H2","H3","H4")
weedle<-weedle[,col_order]

#This will create the new data sets for each of the screens
#We are requiring a starting count of at least 5 for these screens
Kakuna_A<-weedle[,c(1:12)]

#Next will need to add in the pseudo count for days 1:4
Kakuna_A$A0<-Kakuna_A$A0+0.5
Kakuna_A$A1<-Kakuna_A$A1+0.5
Kakuna_A$A2<-Kakuna_A$A2+0.5
Kakuna_A$A3<-Kakuna_A$A3+0.5
Kakuna_A$A4<-Kakuna_A$A4+0.5

#Need to set these variables before we remove them from the model 
Kakuna_A_timepoint_sums<-colSums(Kakuna_A[,8:12])

#remove sequences with low starting counts
Kakuna_A<-Kakuna_A%>%filter(A0>=5)

#Now we will get rid of low A counts the cutoff that is chosen is <3 counts 
x<-which(Kakuna_A$A1<3&Kakuna_A$A2<3)
Kakuna_A$A2[x]<-NA
Kakuna_A$A3[x]<-NA
Kakuna_A$A4[x]<-NA
y<-which(Kakuna_A$A2<3&Kakuna_A$A3<3)
Kakuna_A$A3[y]<-NA
Kakuna_A$A4[y]<-NA
z<-which(Kakuna_A$A3<3&Kakuna_A$A4<3)
Kakuna_A$A4[z]<-NA

#convert counts to transcripts/million
Kakuna_A$A0_per_million<-(Kakuna_A$A0/Kakuna_A_timepoint_sums[1])*1000000
Kakuna_A$A1_per_million<-(Kakuna_A$A1/Kakuna_A_timepoint_sums[2])*1000000
Kakuna_A$A2_per_million<-(Kakuna_A$A2/Kakuna_A_timepoint_sums[3])*1000000
Kakuna_A$A3_per_million<-(Kakuna_A$A3/Kakuna_A_timepoint_sums[4])*1000000
Kakuna_A$A4_per_million<-(Kakuna_A$A4/Kakuna_A_timepoint_sums[5])*1000000

#calculate log2FC from baseline
Kakuna_A$normalized_log_A0<-log2(Kakuna_A$A0_per_million/Kakuna_A$A0_per_million)
Kakuna_A$normalized_log_A1<-log2(Kakuna_A$A1_per_million/Kakuna_A$A0_per_million)
Kakuna_A$normalized_log_A2<-log2(Kakuna_A$A2_per_million/Kakuna_A$A0_per_million)
Kakuna_A$normalized_log_A3<-log2(Kakuna_A$A3_per_million/Kakuna_A$A0_per_million)
Kakuna_A$normalized_log_A4<-log2(Kakuna_A$A4_per_million/Kakuna_A$A0_per_million)

#grab log2fc data to run rlm
Kakuna_A_log<-Kakuna_A %>% dplyr::select(starts_with('normalized'))
x<-c("time0","time1","time2","time3","time4")
colnames(Kakuna_A_log)=x

#run rlm from MASS package and save error in regression model for later filtering
mycoefs<-c()
converged<-c()
t_value<-c()
r2<-c()
MSE<-c()
rse<-c()
for (i in (1:nrow(Kakuna_A_log))){
  myex<-data.frame(row.names = 0:4)
  myex$values<-Kakuna_A_log[i,] %>% unlist()
  myex$times<-0:4
  nalength<-length(which(is.na(myex$values)))
  mymodel<-suppressWarnings(rlm(values~times,data = myex))
  mycoefs[i]<-mymodel$coefficients[2]
  converged[i]<-mymodel$converged
  t_value[i]<-summary(mymodel)$coefficients[2,3]
  rse[i]<-mymodel$s
  sse<-sum((mymodel$wresid)^2)
  ssr<- sum((mymodel$fitted.values-mean(mymodel$model[,1]))^2)
  sst<-sse+ssr
  r2[i]<-1-(sse/sst)
  MSE[i]<-sse/(4-nalength)
}

Kakuna_A$slope<-mycoefs
Kakuna_A$converged<-converged
Kakuna_A$t_value<-t_value
Kakuna_A$MSE<-MSE
Kakuna_A$rse<-rse
Kakuna_A$log_r2<-r2

table((Kakuna_A$converged==TRUE))
table(Kakuna_A$t_value>=1.96)
table(Kakuna_A$rse<2.0)#This was around the same number we got with the r squared 


Kakuna_A_MSE_adj<-Kakuna_A%>%filter(MSE<=3.5)

Kakuna_A_stops<-Kakuna_A_MSE_adj%>%filter(Kakuna_A_MSE_adj$pattern=="stop.first")

Kakuna_A_stops_values<-sort(Kakuna_A_stops$slope,decreasing = T)
a_cut <- Kakuna_A_stops_values[1:5] %>% mean()

Kakuna_A_MSE_adj$slope <- Kakuna_A_MSE_adj$slope-a_cut

Kakuna_A_MSE_adj$binary_stop <- ifelse(Kakuna_A_MSE_adj$slope > 0, 'live', 'die')

hold <- Kakuna_A_MSE_adj %>% filter(pattern == 'RdmWD5_20')
table(hold$binary_stop)

colnames(Kakuna_A_MSE_adj)
hold <- Kakuna_A_MSE_adj %>% dplyr::select(aaseq, aaseq.full, pattern.temp,
                                           pattern.aa, pattern, slope, converged, binary_stop)

saveRDS(hold, 'semi_random_library2.rds')

#compare with other semi_random_library 
table(semi_random$pattern)
table(hold$pattern) #perfect only 2 extra for wd5

hold1 <- semi_random %>% filter(pattern == 'RdmWD5_20')
hold2 <- hold %>% filter(pattern == 'RdmWD5_20')
hold2$aa_seq <- hold2$aaseq

hold3 <- left_join(hold2, hold1, by = 'aa_seq')

x <- which(is.na(hold3$slope.y))
hold3 <- hold3[-x,]

hold3$helic <- str_count(hold3$ss3, 'H')

x <- table(hold3$helic, hold3$binary_stop.x)

y <- tibble(h = row.names(x) %>% as.numeric(), die = x[,1], live = x[,2])

y$total <- y$die + y$live

y$lp <- y$live/y$total

y$h <- y$h/20 * 100


ggplot(y, aes(h, lp))+geom_col()

x <- grep('WEDPWADPWPDGWPDLWADV',hold3$aa_seq)
hold3[x,]
