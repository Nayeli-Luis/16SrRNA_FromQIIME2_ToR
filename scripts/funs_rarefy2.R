# Version 2 of funs_rarefy
# Author: Luis-Vargas M. Nayeli
# email: nayeli.luis@ciencias.unam.mx

library(tidyverse)
library(magrittr)

clean_data <- function(data){
  data %>% 
    pivot_longer(., cols = 2:101,  names_to = "depth", values_to = "asvs") %>% 
    separate(depth, c("depth", "iter"), sep = "_iter.") %>%
    separate(depth, c("depth", "depth_seqs")) %>%
    select(-c(depth)) -> data_long
  return(data_long)
}

asvs_mean <- function(data) {
  
  data %>% 
    select(asvs) %>% 
    summarize(mean_asv = mean(asvs)) %>% 
    pull(mean_asv)
} 