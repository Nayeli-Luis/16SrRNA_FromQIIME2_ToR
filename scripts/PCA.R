#******************PCA fisicoquimicos******************************
library(tidyverse)
library(factoextra)

fq <- read.csv("envrionment_parameters.csv")
colnames(fq) <- c("season", "Tº", "CE", "OD", "pH", "TOC", "IC", "TN")


pca_fq <- prcomp(fq[,-c(1:2)], scale. = TRUE)
kc <- kmeans(fq[, 1:7], 2)

get_eig(pca_fq)
fviz_eig(pca_fq, addlabels = TRUE, hjust = -0.3) + 
  theme_minimal()

fviz_pca_ind(pca_fq)

fviz_pca_var(pca_fq, 
             col.var="contrib",
             gradient.cols = c("#00838F", "#1B5E20", "#CE6510"),
             repel = TRUE )

fviz_pca_biplot(pca_fq, 
                label = "var", 
                habillage = as.factor(kc$cluster))+
  labs(color = NULL) +
  ggtitle(" ") +
  guides(shape = FALSE) +
  scale_color_manual(labels = c("rainy", "dry"), values = c("#A23B72", "#0B6E4F")) +
  theme(
    text = element_text(size = 15), 
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

var <- get_pca_var(pca_fq)
var  
?fviz_pca_biplot
?get_pca_var


#Coordinates for the variables
head(var$coord)
#Correlations between variables and dimensions 
head(var$cor)
# Contributions of the wariables
head(var$contrib)

fviz_pca_var(pca_fq, 
             col.var="contrib",
             gradient.cols = c("#00838F", "#1B5E20", "#CE6510"),
             repel = TRUE ) 

#Supongo pa escoger los mas chilos
fviz_pca_var(pca_fam, col.var = "blue" , select.var = list(contrib = 5)) + 
  xlim(-0.15,0.15) + ylim(-0.1,0.1)


#Biplot para individuos y variables
biplot_fq <- fviz_pca_biplot(pca_fq, repel = FALSE,
                             col.var="contrib",
                             gradient.cols = c("#00838F", "#1B5E20", "#CE6510"), # Variables color
                             col.ind = "black" )# Individuals color

fviz_pca(pca_fq)

svg("biplot_fq.svg")
biplot_fq
dev.off()