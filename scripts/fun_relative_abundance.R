# Calculates the proportion of each taxa.
rel_ab <- function(dataTax, dataTotalSeq, iTax, iTotSeq){
  
  dataTax%<>%
    set_rownames(dataTax[,iTax]) %>%
    select(-iTax)
  
  dataTotalSeq %<>%
    set_rownames(dataTotalSeq[,iTotSeq]) %>%
    select(-iTotSeq)
  
  
  data.frame(mapply('/', dataTax, dataTotalSeq))
  
}

# Filter the taxa in function of a proportion
filter_abundances <- function(dataRelAb, filterValue){
  
  dataRelAb %>% 
    t() %>%
    as.data.frame() %>%
    filter_all(any_vars(. > filterValue))
}

# Add "Other"
other_row <- function(dataset, i){
  
  rownames(dataset) <- dataset[,i]
  dataset %<>% select(-i)
  
  total <- colSums(dataset[1:ncol(dataset)])
  other <-  1 - total
  
  dataset <- rbind(dataset, other)
  rownames(dataset)[nrow(dataset)] <- "Other"
  dataset <- rownames_to_column(dataset, var = "Phylum")
  dataset

}
