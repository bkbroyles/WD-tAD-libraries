##Process design library 7/12/21

#load packages
library(tidyverse)
library(Biobase)
library(xlsx)
library(MASS)

#load and merge NGS raw counts data with barcode/set data ----
load('Process Design Library/raw_data/design_lib_Data1.rData')
  
  #bioreps will have NGS counts, seq_data will have other sequence information
  bioreps <- exprs(esetCounts) %>% as_tibble()
  seq_data <- exprs(eset) %>% as_tibble()
  seq_data$id <- row.names(exprs(eset))

  #split up id column into more useful information
  seq_data$id_number <- ''
  seq_data$aa_seq <- ''
  seq_data$barcode <- ''

  for(i in 1:nrow(seq_data)){
    split <- strsplit(seq_data$id[i],'[|]')[[1]]
    seq_data$id_number[i] <- split[1]
    seq_data$aa_seq[i] <- split[2]
    seq_data$barcode[i] <- split[3]
  }

  #merge bioreps and seq_data
  bioreps$id_number <- seq_data$id_number
  bioreps$aa_seq <- seq_data$aa_seq
  bioreps$barcode <- seq_data$barcode

  #load barcode (barcode and set) data and merge with counts
  barcodes <- read_rds('Process Design Library/raw_data/barcodes.rds')
  bioreps <- left_join(bioreps, barcodes, by = 'barcode')

#split into A and H screen and fix column names ----
A_screen <- bioreps %>% dplyr::select(id_number, aa_seq, barcode, set,
                                      starts_with('P1'), 
                                      starts_with('R'))

#H_screen <- bioreps %>% dplyr::select(id_number, aa_seq, barcode, set,
#                                      starts_with('P2'),
#                                      starts_with('S'))

#reorder these columns
A_screen <- A_screen %>% dplyr::select(id_number, aa_seq, 
                                       barcode, set, contains('bc1'),contains('bc2'),
                                       contains('bc3'),contains('bc4'),contains('bc5'))

#H_screen <- H_screen %>% dplyr::select(id_number, aa_seq,
#                                       barcode, set, contains('bc1'),contains('bc2'),
#                                       contains('bc3'),contains('bc4'),contains('bc5'))

#rename R columns to be A, and S columns to be H
colnames(A_screen)[c(6:9,11:14,16:19,21:24,26:29)] <-  
  colnames(A_screen %>% dplyr::select(matches('R[1234]'))) %>%
  gsub('R', 'A', .)

colnames(A_screen)[c(5,10,15,20,25)] <-  
  colnames(A_screen %>% dplyr::select(contains('P'))) %>% 
  gsub('P1', 'A0', .)

#colnames(H_screen)[c(6:9,11:14,16:19,21:24,26:29)] <-  
#  colnames(H_screen %>% dplyr::select(matches('S[1234]'))) %>% 
#  gsub('S', 'H', .)

#colnames(H_screen)[c(5,10,15,20,25)] <- 
#  colnames(H_screen %>% dplyr::select(contains('P'))) %>% 
#  gsub('P2', 'H0', .)

####Processing Library Starts Here     #####
#add psuedocount of 0.5 to each time ----
A_screen$A0_bc1 <- A_screen$A0_bc1 + 0.5
A_screen$A0_bc2 <- A_screen$A0_bc2 + 0.5
A_screen$A0_bc3 <- A_screen$A0_bc3 + 0.5
A_screen$A0_bc4 <- A_screen$A0_bc4 + 0.5
A_screen$A0_bc5 <- A_screen$A0_bc5 + 0.5

A_screen$A1_bc1 <- A_screen$A1_bc1 + 0.5
A_screen$A1_bc2 <- A_screen$A1_bc2 + 0.5
A_screen$A1_bc3 <- A_screen$A1_bc3 + 0.5
A_screen$A1_bc4 <- A_screen$A1_bc4 + 0.5
A_screen$A1_bc5 <- A_screen$A1_bc5 + 0.5

A_screen$A2_bc1 <- A_screen$A2_bc1 + 0.5
A_screen$A2_bc2 <- A_screen$A2_bc2 + 0.5
A_screen$A2_bc3 <- A_screen$A2_bc3 + 0.5
A_screen$A2_bc4 <- A_screen$A2_bc4 + 0.5
A_screen$A2_bc5 <- A_screen$A2_bc5 + 0.5

A_screen$A3_bc1 <- A_screen$A3_bc1 + 0.5
A_screen$A3_bc2 <- A_screen$A3_bc2 + 0.5
A_screen$A3_bc3 <- A_screen$A3_bc3 + 0.5
A_screen$A3_bc4 <- A_screen$A3_bc4 + 0.5
A_screen$A3_bc5 <- A_screen$A3_bc5 + 0.5

A_screen$A4_bc1 <- A_screen$A4_bc1 + 0.5
A_screen$A4_bc2 <- A_screen$A4_bc2 + 0.5
A_screen$A4_bc3 <- A_screen$A4_bc3 + 0.5
A_screen$A4_bc4 <- A_screen$A4_bc4 + 0.5
A_screen$A4_bc5 <- A_screen$A4_bc5 + 0.5

#moving to H screen
#H_screen$H0_bc1 <- H_screen$H0_bc1 + 0.5
#H_screen$H0_bc2 <- H_screen$H0_bc2 + 0.5
#H_screen$H0_bc3 <- H_screen$H0_bc3 + 0.5
#H_screen$H0_bc4 <- H_screen$H0_bc4 + 0.5
#H_screen$H0_bc5 <- H_screen$H0_bc5 + 0.5

#H_screen$H1_bc1 <- H_screen$H1_bc1 + 0.5
#H_screen$H1_bc2 <- H_screen$H1_bc2 + 0.5
#H_screen$H1_bc3 <- H_screen$H1_bc3 + 0.5
#H_screen$H1_bc4 <- H_screen$H1_bc4 + 0.5
#H_screen$H1_bc5 <- H_screen$H1_bc5 + 0.5

#H_screen$H2_bc1 <- H_screen$H2_bc1 + 0.5
#H_screen$H2_bc2 <- H_screen$H2_bc2 + 0.5
#H_screen$H2_bc3 <- H_screen$H2_bc3 + 0.5
#H_screen$H2_bc4 <- H_screen$H2_bc4 + 0.5
#H_screen$H2_bc5 <- H_screen$H2_bc5 + 0.5

#H_screen$H3_bc1 <- H_screen$H3_bc1 + 0.5
#H_screen$H3_bc2 <- H_screen$H3_bc2 + 0.5
#H_screen$H3_bc3 <- H_screen$H3_bc3 + 0.5
#H_screen$H3_bc4 <- H_screen$H3_bc4 + 0.5
#H_screen$H3_bc5 <- H_screen$H3_bc5 + 0.5

#H_screen$H4_bc1 <- H_screen$H4_bc1 + 0.5
#H_screen$H4_bc2 <- H_screen$H4_bc2 + 0.5
#H_screen$H4_bc3 <- H_screen$H4_bc3 + 0.5
#H_screen$H4_bc4 <- H_screen$H4_bc4 + 0.5
#H_screen$H4_bc5 <- H_screen$H4_bc5 + 0.5

#Indicate bioreps with day0 counts < 5, filter out sequences with < 2 bioreps ----
A_screen <- A_screen %>% mutate(
  include1 = ifelse(A0_bc1 < 5, 'no', 'yes'),
  include2 = ifelse(A0_bc2 < 5, 'no', 'yes'),
  include3 = ifelse(A0_bc3 < 5, 'no', 'yes'),
  include4 = ifelse(A0_bc4 < 5, 'no', 'yes'),
  include5 = ifelse(A0_bc5 < 5, 'no', 'yes')
)

#H_screen <- H_screen %>% mutate(
#  include1 = ifelse(H0_bc1 < 5, 'no', 'yes'),
#  include2 = ifelse(H0_bc2 < 5, 'no', 'yes'),
#  include3 = ifelse(H0_bc3 < 5, 'no', 'yes'),
#  include4 = ifelse(H0_bc4 < 5, 'no', 'yes'),
#  include5 = ifelse(H0_bc5 < 5, 'no', 'yes')
#)

A_screen$multiple_biorep <- 'no'
A_screen$count <- 0
for(i in 1:nrow(A_screen)){
  x <- 0
  if(A_screen$include1[i] == 'yes')
    x <- x + 1
  if(A_screen$include2[i] == 'yes')
    x <- x + 1
  if(A_screen$include3[i] == 'yes')
    x <- x + 1
  if(A_screen$include4[i] == 'yes')
    x <- x + 1
  if(A_screen$include5[i] == 'yes')
    x <- x + 1
  A_screen$multiple_biorep[i] <- ifelse(x > 1, 'yes', 'no')
  A_screen$count[i] <- x
}


#H_screen$multiple_biorep <- 'no'
#H_screen$count <- 0
#for(i in 1:nrow(H_screen)){
#  x <- 0
#  if(H_screen$include1[i] == 'yes')
#    x <- x + 1
#  if(H_screen$include2[i] == 'yes')
#    x <- x + 1
#  if(H_screen$include3[i] == 'yes')
#    x <- x + 1
#  if(H_screen$include4[i] == 'yes')
#    x <- x + 1
#  if(H_screen$include5[i] == 'yes')
#    x <- x + 1
#  H_screen$multiple_biorep[i] <- ifelse(x > 1, 'yes', 'no')
#  H_screen$count[i] <- x
#}

#dlib <- list(A_screen, H_screen)
#names(dlib) <- c('A_screen', 'H_screen')
#saveRDS(dlib, 'design_library_raw_counts.rds')
#saveRDS(A_screen, 'A_screen_counts.rds')

#copy data frame for aggregate linear regression----
A_screen_agg <- A_screen
#H_screen_agg <- H_screen

#add up daily counts so I can tell where each sequences should stop
A_screen_agg$day0 <- 0
A_screen_agg$day1 <- 0
A_screen_agg$day2 <- 0
A_screen_agg$day3 <- 0
A_screen_agg$day4 <- 0

A_screen_agg$day0 <- apply(A_screen_agg %>% dplyr::select(starts_with('A0')), 1, sum)
A_screen_agg$day1 <- apply(A_screen_agg %>% dplyr::select(starts_with('A1')), 1, sum)
A_screen_agg$day2 <- apply(A_screen_agg %>% dplyr::select(starts_with('A2')), 1, sum)
A_screen_agg$day3 <- apply(A_screen_agg %>% dplyr::select(starts_with('A3')), 1, sum)
A_screen_agg$day4 <- apply(A_screen_agg %>% dplyr::select(starts_with('A4')), 1, sum)

#same for H_screen
#H_screen_agg$day0 <- 0
#H_screen_agg$day1 <- 0
#H_screen_agg$day2 <- 0
#H_screen_agg$day3 <- 0
#H_screen_agg$day4 <- 0

#H_screen_agg$day0 <- apply(H_screen_agg %>% dplyr::select(starts_with('H0')), 1, sum)
#H_screen_agg$day1 <- apply(H_screen_agg %>% dplyr::select(starts_with('H1')), 1, sum)
#H_screen_agg$day2 <- apply(H_screen_agg %>% dplyr::select(starts_with('H2')), 1, sum)
#H_screen_agg$day3 <- apply(H_screen_agg %>% dplyr::select(starts_with('H3')), 1, sum)
#H_screen_agg$day4 <- apply(H_screen_agg %>% dplyr::select(starts_with('H4')), 1, sum)

#INDICATE WHERE TO STOP (FALLS BELOW 3 COUNTS AND STAYS BELOW 3 COUNTS)
A_screen_agg$where_stop <- 'day4'
for(i in 1:nrow(A_screen_agg)){
 lib <- A_screen_agg[i,] %>% dplyr::select(starts_with('day')) %>% unlist()
 hold <- which(lib < 3) %>% names()
if('day3' %in% hold & 'day4' %in% hold){
   A_screen_agg$where_stop[i] <- 'day3'
 }
 if('day2' %in% hold & 'day3' %in% hold){
   A_screen_agg$where_stop[i] <- 'day2'
 }
 if('day1' %in% hold & 'day2' %in% hold){
   A_screen_agg$where_stop[i] <- 'day1'
 }
}

#H_screen_agg$where_stop <- 'day4'
#for(i in 1:nrow(H_screen_agg)){
#  lib <- H_screen_agg[i,] %>% dplyr::select(starts_with('day')) %>% unlist()
#  hold <- which(lib < 3) %>% names()
#  if('day3' %in% hold & 'day4' %in% hold){
#    H_screen_agg$where_stop[i] <- 'day3'
#  }
#  if('day2' %in% hold & 'day3' %in% hold){
#    H_screen_agg$where_stop[i] <- 'day2'
#  }
#  if('day1' %in% hold & 'day2' %in% hold){
#    H_screen_agg$where_stop[i] <- 'day1'
#  }
#}



#dlib <- list(A_screen_agg, H_screen_agg)
#names(dlib) <- c('A_screen', 'H_screen')
#saveRDS(dlib, 'design_library_agg_counts.rds')

#normalize by total counts of day regardless of biorep----
day0sum <- A_screen_agg %>% dplyr::select(starts_with('A0')) %>% unlist() %>% sum()
day1sum <- A_screen_agg %>% dplyr::select(starts_with('A1')) %>% unlist() %>% sum()
day2sum <- A_screen_agg %>% dplyr::select(starts_with('A2')) %>% unlist() %>% sum()
day3sum <- A_screen_agg %>% dplyr::select(starts_with('A3')) %>% unlist() %>% sum()
day4sum <- A_screen_agg %>% dplyr::select(starts_with('A4')) %>% unlist() %>% sum()
####MOVING SCRIPTS AFTER THIS POINT TO MAKE 1w and 1d error bar plot####

A_screen_agg$day0 <- A_screen_agg$day0/day0sum*1000000
A_screen_agg$day1 <- A_screen_agg$day1/day1sum*1000000
A_screen_agg$day2 <- A_screen_agg$day2/day2sum*1000000
A_screen_agg$day3 <- A_screen_agg$day3/day3sum*1000000
A_screen_agg$day4 <- A_screen_agg$day4/day4sum*1000000

#moving to H screen
#day0sum <- H_screen_agg %>% dplyr::select(starts_with('H0')) %>% unlist() %>% sum()
#day1sum <- H_screen_agg %>% dplyr::select(starts_with('H1')) %>% unlist() %>% sum()
#day2sum <- H_screen_agg %>% dplyr::select(starts_with('H2')) %>% unlist() %>% sum()
#day3sum <- H_screen_agg %>% dplyr::select(starts_with('H3')) %>% unlist() %>% sum()
#day4sum <- H_screen_agg %>% dplyr::select(starts_with('H4')) %>% unlist() %>% sum()

#H_screen_agg$day0 <- H_screen_agg$day0/day0sum*1000000
#H_screen_agg$day1 <- H_screen_agg$day1/day1sum*1000000
#H_screen_agg$day2 <- H_screen_agg$day2/day2sum*1000000
#H_screen_agg$day3 <- H_screen_agg$day3/day3sum*1000000
#H_screen_agg$day4 <- H_screen_agg$day4/day4sum*1000000


#remove sequence without multiple bioreps
A_screen_agg <- A_screen_agg %>% filter(count > 1)
#H_screen_agg <- H_screen_agg %>% filter(count > 1)

#set counts to NA if they occur after where_stop indicates
for(i in 1:nrow(A_screen_agg)){
  if(A_screen_agg$where_stop[i] == 'day1'){
    A_screen_agg[i,c("day2","day3","day4")] <- NA
  }
  if(A_screen_agg$where_stop[i] == 'day2'){
    A_screen_agg[i,c("day3","day4")] <- NA
  }
  if(A_screen_agg$where_stop[i] == 'day3'){
    A_screen_agg[i,c("day4")] <- NA
  }
}

#for(i in 1:nrow(H_screen_agg)){
#  if(H_screen_agg$where_stop[i] == 'day1'){
#    H_screen_agg[i,c("day2","day3","day4")] <- NA
#  }
#  if(H_screen_agg$where_stop[i] == 'day2'){
#    H_screen_agg[i,c("day3","day4")] <- NA
#  }
# if(H_screen_agg$where_stop[i] == 'day3'){
#    H_screen_agg[i,c("day4")] <- NA
#  }
#}

#dlib <- list(A_screen_agg, H_screen_agg)
#names(dlib) <- c('A_screen', 'H_screen')
#saveRDS(dlib, 'design_library_norm_counts.rds')

#compute log2 fc from day0 to each day
A_screen_agg$day1 <- log2(A_screen_agg$day1/A_screen_agg$day0)
A_screen_agg$day2 <- log2(A_screen_agg$day2/A_screen_agg$day0)
A_screen_agg$day3 <- log2(A_screen_agg$day3/A_screen_agg$day0)
A_screen_agg$day4 <- log2(A_screen_agg$day4/A_screen_agg$day0)
A_screen_agg$day0 <- log2(A_screen_agg$day0/A_screen_agg$day0)

#H_screen_agg$day1 <- log2(H_screen_agg$day1/H_screen_agg$day0)
#H_screen_agg$day2 <- log2(H_screen_agg$day2/H_screen_agg$day0)
#H_screen_agg$day3 <- log2(H_screen_agg$day3/H_screen_agg$day0)
#H_screen_agg$day4 <- log2(H_screen_agg$day4/H_screen_agg$day0)
#H_screen_agg$day0 <- log2(H_screen_agg$day0/H_screen_agg$day0)

#dlib <- list(A_screen_agg, H_screen_agg)
#names(dlib) <- c('A_screen', 'H_screen')
#saveRDS(dlib, 'design_library_log2fc_counts.rds')

#Robust linear regression on aggregate counts----
regress1 <- function(rw) {
  df <- data.frame(expr = rw,
                   time = 1:4)
  res <- suppressWarnings(rlm(rw ~ 0 + time, data = df))
  return(res)
}

models <- apply(A_screen_agg[,c("day1","day2","day3","day4")],1,regress1)

A_screen_agg$slope <- sapply(models, function(i) i$coefficients)
A_screen_agg$converged <- sapply(models, function(i) i$converged)
A_screen_agg$mse <- 0
for(i in 1:nrow(A_screen_agg)){
  na_length <- which(is.na(A_screen_agg[i,c("day1","day2","day3","day4")])) %>% length()
  sse <- sum((models[[i]]$wresid)^2)
  A_screen_agg$mse[i] <- sse/(4 - na_length) 
}
A_screen_agg$rmse <- sqrt(A_screen_agg$mse)
#models <- apply(H_screen_agg[,c("day1","day2","day3","day4")],1,regress1)

#H_screen_agg$slope <- sapply(models, function(i) i$coefficients)
#H_screen_agg$converged <- sapply(models, function(i) i$converged)
#H_screen_agg$mse <- 0
#for(i in 1:nrow(H_screen_agg)){
#  na_length <- which(is.na(H_screen_agg[i,c("day1","day2","day3","day4")])) %>% length()
#  sse <- sum((models[[i]]$wresid)^2)
#  H_screen_agg$mse[i] <- sse/(4 - na_length) 
#}

##adding live or die based on stop codon performance and adjust slope
A_screen <- A_screen_agg %>% dplyr::select(id_number, aa_seq, set,
                                           slope, converged, mse, rmse)

a_stops <- A_screen %>% filter(set == 'stop_codon') %>% 
  dplyr::select(slope) %>% unlist()

hist(a_stops)

a_stops <- sort(a_stops, decreasing = T)

a_cut <- mean(a_stops[1:5])

A_screen$binary_stop <- ifelse(A_screen$slope > a_cut, 'live','die')

A_screen$slope <- A_screen$slope - a_cut

#H_screen <- H_screen_agg %>% dplyr::select(id_number, aa_seq, set,
#                                           slope, converged, mse)

##save aggregate data and cleaned up data
#saveRDS(A_screen_agg, 'design_library_aggregate_with_log2_A.rds')
#saveRDS(H_screen_agg, 'design_library_aggregate_with_log2_H.rds')
#saveRDS(A_screen, 'design_library.rds')
#saveRDS(H_screen, 'design_library_H_screen.rds')


