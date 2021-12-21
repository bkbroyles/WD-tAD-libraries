##add binary stop
library(tidyverse)
a_screen <- read_rds('design_library.rds')
h_screen <- read_rds('Intermediate Datasets/design_library_H_screen.rds')

a_stops <- a_screen %>% filter(set == 'stop_codon') %>% dplyr::select(slope) %>% unlist()
h_stops <- h_screen %>% filter(set == 'stop_codon') %>% dplyr::select(slope) %>% unlist()

hist(a_stops)
hist(h_stops)

a_stops <- sort(a_stops, decreasing = T)
h_stops <- sort(h_stops, decreasing = T)

a_cut <- mean(a_stops[1:5])
h_cut <- mean(h_stops[1:5])

a_screen$binary_stop <- ifelse(a_screen$slope > a_cut, 'live','die')
h_screen$binary_stop <- ifelse(h_screen$slope > h_cut, 'live', 'die')

a_screen$slope <- a_screen$slope - a_cut
h_screen$slope <- h_screen$slope - h_cut

#save with binary_stop
saveRDS(a_screen, 'design_library.rds')
saveRDS(h_screen, 'Intermediate Datasets/design_library_H_screen.rds')


