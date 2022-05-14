# Análisis de especies indicadoras (IndVal)
# 10 julio 2021

library(tidyverse)
library(stats)
library(indicspecies)
library(gridExtra)
library(magrittr)

# Se utiliza a nivel de género, entonce se utiliza 'level-6.csv'

genus <- read.csv("level-6.csv", header = TRUE) %>%
  set_rownames(.$index) %>%
  select(-c(site, description, sample_type, index)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(., var = "Taxon")

# Separar taxa

genus %<>%
  separate(Taxon, c("Kingdom", "Phylum"), sep = ".p__") %>%
  separate(Phylum, c("Phylum", "Class"), sep = ".c__") %>%
  separate(Class, c("Class", "Order"), sep = ".o__") %>%
  separate(Order, c("Order", "Family"), sep = ".f__") %>%
  separate(Family, c("Family", "Genus"), sep = ".g__") %>%
  mutate_all(na_if,"") %>%
  select(-c(Kingdom:Family)) %>%
  na.omit() %>%
  filter(Genus != "uncultured") %>%# Eliminar los que digan uncultured
  set_rownames(.$Genus) %>%
  select(-Genus) %>%
  t() %>%
  as.data.frame()
  
# Queremos comparar por sitio no por muestra, entonces debemos elegir que 
# operación hacer. Creo que aquí sin conviene hacer un promedio. 

genus %<>%
  mutate(site = case_when(
    grepl("A", rownames(.)) ~ "A",
    grepl("B", rownames(.)) ~ "B",
    grepl("C", rownames(.)) ~ "C",
    grepl("D", rownames(.)) ~ "D"
  ))

# Preparar los datos para meterlos al indval
  
indval.matrix <- genus %>%
  select(-site) %>%
  as.matrix()
rownames(indval.matrix) <- NULL
sites <- as.vector(genus$site)

# Análisis

indval <- multipatt(indval.matrix, sites, duleg = FALSE, control = how(nperm=999))
summary(indval, indvalcomp = TRUE)

