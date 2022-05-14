# **Parameters** of `clean_rarefy_data()`:
#   - `data` : A dataset with just one categorical column.
# - `catCol`:  1 or 101, the number represents de number of 
#     the categorical column, could be:  1 if your column are at the 
#     beginig o 101 if you column are at the end. 
# 
# **Output**: A list, each object of that list is a dataset. 
#     The number of datasets contained in the list depends on the 
#     number of samples or categories you have. 
#     Each dataset contains the 10 iterations from each read depth.


clean_rarefy_data <- function(data, catCol){
  
  library(tidyverse)
  library(magrittr)
  
  if(catCol == 1){
    data %<>% pivot_longer(., cols = 2:101, names_to = "depth", values_to = "sequences")
  }
  
  if(catCol == 101){
    data %<>% pivot_longer(., cols = 1:100, names_to = "depth", values_to = "sequences")
  }
  
  data %<>%
    separate(depth, c("depth", "iter"), sep = "_iter.") %>%
    separate(depth, c("depth", "Depth")) %>%
    select(-c(depth, iter)) %>%
    split(., .[,1])
}


rarefy <- function(df, var1, var2) {
  var1 <- enquo(var1); var2 <- enquo(var2);
  df %<>% group_by(!!var1) %>%
    dplyr::summarise(mean = mean(!!var2),
                     sd = sd(!!var2))
  
}

# **Parameters**: 
#   `L`: The list obtained in the last step.
# 
# **Output**: Another list.

rarefy_by <- function(L){
  for (i in seq_along(L)){
    L[[i]] <- rarefy(L[[i]], var1 = Depth, var2 = sequences) 
  }
  return(L)
} 

# **Parameters**:
#   `L`: The list obtained in the last step.

clean_rarefy_data2 <- function(L){
  obs <-  bind_rows(L, .id = "category") %>%
    as_tibble()
  x <- map_dbl(obs$Depth, as.numeric)
  obs %<>% 
    cbind(depth = x) %>%
    select(-Depth)
}



