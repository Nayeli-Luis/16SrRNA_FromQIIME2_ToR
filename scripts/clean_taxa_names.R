# Clean taxa columns

clea_taxa_names <- function(data, colName, numberLevel){
  
  arg <- match.call()
  
  if(numberLevel == 1){
    
  data <- sub(pattern = "d__", replacement = "", data$Domain)
  }
  
  if(numberLevel == 2){
    
    data %<>% 
      separate(eval(arg$colName), c("Domain", "Phylum"), sep = ".p__") %>%
      mutate_all(na_if,"")
  }
  
  if(numberLevel == 3){
    
    data %<>% 
      separate(eval(arg$colName), c("Domain", "Phylum"), sep = "; p__") %>%
      separate(Phylm, c("Phylum", "Class"), sep = "; c__") %>%
      mutate_all(na_if,"")
  }
  
  if(numberLevel == 4){
    
    data %<>% 
      separate(eval(arg$colName), c("Domain", "Phylum"), sep = "; p__") %>%
      separate(Phylm, c("Phylum", "Class"), sep = "; c__") %>%
      separate(Class, c("Class", "Order"), sep = "; o__") %>%
      mutate_all(na_if,"")
  }
  
  if(numberLevel == 5){
    
    data %<>% 
      separate(eval(arg$colName), c("Domain", "Phylum"), sep = "; p__") %>%
      separate(Phylm, c("Phylum", "Class"), sep = "; c__") %>%
      separate(Class, c("Class", "Order"), sep = "; o__") %>%
      separate(Order, c("Order", "Family"), sep = "; f__") %>%
      mutate_all(na_if,"")
  }
  
  if(numberLevel == 6){
    
    data %<>% 
      separate(eval(arg$colName), c("Domain", "Phylum"), sep = "; p__") %>%
      separate(Phylm, c("Phylum", "Class"), sep = "; c__") %>%
      separate(Class, c("Class", "Order"), sep = "; o__") %>%
      separate(Order, c("Order", "Family"), sep = "; f__") %>%
      separate(Family, c("Family", "Genus"), sep = "; g__") %>%
      mutate_all(na_if,"")
  }
  
  
  if(numberLevel == 7){
    
    data %<>% 
      separate(eval(arg$colName), c("Domain", "Phylum"), sep = "; p__") %>%
      separate(Phylm, c("Phylum", "Class"), sep = "; c__") %>%
      separate(Class, c("Class", "Order"), sep = "; o__") %>%
      separate(Order, c("Order", "Family"), sep = "; f__") %>%
      separate(Family, c("Family", "Genus"), sep = "; g__") %>%
      separate(Genus, c("Genus", "Specie"), sep = "; s__") %>%
      mutate_all(na_if,"")
  }
  
  data$Domain <- sub (pattern = "d__", replacement = "", 
         
                                   data$Domain)
  
  data %>% na.omit()
  
}