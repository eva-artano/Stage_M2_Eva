#-------------------------------------------------------------------------------
# ENVIRONMENT
#-------------------------------------------------------------------------------

PATH_SAVE <- 'D:/Resultats_memoire/scRNAseq/WT/Saves'
source('D:/Resultats_memoire/Custom_Functions.R')

#-------------------------------------------------------------------------------
# K-MEANS CLUSTERING
#-------------------------------------------------------------------------------
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))


library(Seurat)
table <- OrderMatrix(data, 'Slingshot_1', slot = 'scale.data')

gene_list <- read.csv('D:/Resultats_memoire/scRNAseq/WT-KO/CompareExpr()/new_Candidats_compare_expr_TOX.txt', sep ="", header = FALSE)
gene_list <- gene_list$V1

## SETUP TABLES ----------------------------------------------------------------
scaled <- table
pseudotime <- scaled[,'Pseudotime']
color <- rep("grey", length(unique(candidats)))


# Filtering genes list to be in colnames and variables
candidats <- colnames(scaled)[colnames(scaled) %in% gene_list]

# Create tables 
scaled.candidats <- scaled[, c("Pseudotime", candidats)]
scaled_kmean <- scaled.candidats[,-1] #enlève la colonne Pseudotime
matrix_expr <- DrawExpr(scaled.candidats, feature.list = candidats, bin.number = 12,
                        scale = F, std = F, superposed = T, 
                        by.order = F, color = color, show.legend = FALSE)

## NUMBER OF CLUSTERS ----------------------------------------------------------

set.seed(975)

str(matrix_expr)
rownames(matrix_expr) <- matrix_expr$Gene

library(tidyr)
library(dplyr)

# On reformate : chaque ligne = un gène, colonnes = bins, valeurs = expression
expr_wide <- matrix_expr %>%
  select(Gene, Bin, Expr) %>%
  pivot_wider(names_from = Bin, values_from = Expr)

# Optionnel : mettre les noms de gènes en rownames
expr_mat <- as.data.frame(expr_wide)
rownames(expr_mat) <- expr_mat$Gene
expr_mat$Gene <- NULL

set.seed(975)
Tot.Wthnss <- c()
Clustering_list <- list()


for (i in 1:30) {
  km_res <- kmeans(expr_mat, centers = i)
  Clustering_list[[i]] <- km_res
  Tot.Wthnss <- c(Tot.Wthnss, km_res$tot.withinss)
}

ggplot(as.data.frame(Tot.Wthnss), aes(x= c(1:30), y = Tot.Wthnss)) +
  geom_line() + geom_point() + 
  scale_x_continuous(breaks = c(1:30)) +
  theme_classic() +
  xlab('Number of clusters k') + ylab('Total Within Sum of Square')

set.seed(975)
km_final <- kmeans(expr_mat, centers = 6)

gene_clusters <- data.frame(
  Gene = rownames(expr_mat),
  Cluster = as.factor(km_final$cluster)
)

# On split les gènes par cluster
gene_by_cluster <- split(gene_clusters$Gene, gene_clusters$Cluster)

# Création dynamique de variables : genes_cluster1, genes_cluster2, etc.
for (k in names(gene_by_cluster)) {
  assign(paste0("genes_cluster", k), gene_by_cluster[[k]])
}

#-------------------------------------------------------------------------------
# BOUCLE K-MEAN
#-------------------------------------------------------------------------------


set.seed(975)

for (k in 1:10) {
  # K-means clustering
  km_final <- kmeans(expr_mat, centers = k)
  
  # Assigner les gènes à leurs clusters
  gene_clusters <- data.frame(
    Gene = rownames(expr_mat),
    Cluster = as.factor(km_final$cluster)
  )
  
  # Split par cluster
  gene_by_cluster <- split(gene_clusters$Gene, gene_clusters$Cluster)
  
  # Pour chaque cluster, dessiner l'expression
  for (cluster_id in names(gene_by_cluster)) {
    genes <- gene_by_cluster[[cluster_id]]
    
    # Ajouter le gène "Tox"
    genes <- c(genes, "Tox")
    
    # Coloration : "Tox" sera en violet, le reste en gris
    color <- rep("grey", length(genes))
    genes <- sort(genes)
    tox_index <- which(genes == "Tox")
    if (length(tox_index) > 0) color[tox_index] <- "purple"
    
    # Titre du plot
    plot_title <- paste("k =", k, "- Cluster", cluster_id)
    
    # Dessin
    assign(
      paste0("Cluster_k", k, "_", cluster_id),
      DrawExpr(
        scaled,
        feature.list = genes,
        bin.number = 12,
        scale = FALSE,
        std = FALSE,
        superposed = TRUE,
        by.order = FALSE,
        color = color,
        lwd = 2,
        ylim = c(-1.5, 2.5),
        show.legend = FALSE,
        main = plot_title
      )
    )
  }
}
getwd()

Cluster1 <- unique(Cluster_k6_1$Gene)
write.table(Cluster1, file = "D:/Resultats_memoire/scRNAseq/kmean_TOX/Final_version/6/cluster1.txt" , sep = " ", quote = FALSE, row.names =  FALSE)

Cluster2 <- unique(Cluster_k6_2$Gene)
write.table(Cluster2, file = "D:/Resultats_memoire/scRNAseq/kmean_TOX/Final_version/6/cluster2.txt" , sep = " ", quote = FALSE, row.names =  FALSE)

Cluster3 <- unique(Cluster_k6_3$Gene)
write.table(Cluster3, file = "D:/Resultats_memoire/scRNAseq/kmean_TOX/Final_version/6/cluster3.txt" , sep = " ", quote = FALSE, row.names =  FALSE)

Cluster4 <- unique(Cluster_k6_4$Gene)
write.table(Cluster4, file = "D:/Resultats_memoire/scRNAseq/kmean_TOX/Final_version/6/cluster4.txt" , sep = " ", quote = FALSE, row.names =  FALSE)

Cluster5 <- unique(Cluster_k6_5$Gene)
write.table(Cluster5, file = "D:/Resultats_memoire/scRNAseq/kmean_TOX/Final_version/6/cluster5.txt" , sep = " ", quote = FALSE, row.names =  FALSE)

Cluster6 <- unique(Cluster_k6_6$Gene)
write.table(Cluster6, file = "D:/Resultats_memoire/scRNAseq/kmean_TOX/Final_version/6/cluster6.txt" , sep = " ", quote = FALSE, row.names =  FALSE)



#-------------------------------------------------------------------------------
# Figure propre
#-------------------------------------------------------------------------------



#
# Tentative de centroid pour le cluster1

x <- table_DC
bin.number <- 11
feature.list <- cluster3
by.order <- F
std <- F
scale <- F
superposed <- TRUE
color <- c(rep("grey",37), black)
lwd = 2
ylim = c(-1, 2)
xlab = "Pseudotime"
ylab = ""
main = "Expression DC"





suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# If x axis is setted to be cell ordered by rank and not by pseudotime value
if(by.order){
  x$Pseudotime <- 1:nrow(x)
}

# Subset table for more efficiency during following steps
x <- x[,colnames(x) %in% c('Pseudotime', feature.list)]

# Attributes Bin to each cell
x <- PseudotimeRepartition(x, bin.number)

# Initiate loop
x.dim         <- as.numeric(unique(x$Bin))
expr_table    <- data.frame()
absent_genes  <- c()

# Generate expr_table containing gene, bin, mean expr, std - - - - - - - - - -
for(feature in feature.list){
  if(feature %in% colnames(x)){
    for(b in x.dim){
      # Subset Bin
      values     <- x[x$Bin %in% b, feature]
      expr_table <- rbind(expr_table, c(feature, b, mean(values), sd(values)))
    }
    colnames(expr_table) <- c('Gene', 'Bin', 'Expr', 'Std')
    
    # Scale current gene mean expression if required - - - - - - - - - - - - -
    if(scale){
      scaled_expr <- as.numeric(expr_table$Expr[expr_table$Gene %in% feature])
      scaled_expr <- (scaled_expr-min(scaled_expr))/
        (max(scaled_expr)-min(scaled_expr))
      expr_table$Expr[expr_table$Gene %in% feature] <- scaled_expr
      # Modify parameters to correspond correctly to scaled plot
      ylim        <- c(0,1)
      std         <- F
    }
    
    # Draw individual curves if not superposed - - - - - - - - - - - - - - - -
    if(!superposed){
      p <- ggplot(expr_table[expr_table$Gene %in% feature,], 
                  aes(x = as.numeric(Bin), group = Gene)) +
        geom_line(linewidth = lwd, aes(y = as.numeric(Expr), color = Gene)) +
        scale_color_manual(values = color) +
        theme_classic() +
        ylim(ylim) 
      if(std){
        p <- p + geom_ribbon(
          aes(y = as.numeric(Expr), ymin = as.numeric(Expr) - as.numeric(Std), 
              ymax = as.numeric(Expr) + as.numeric(Std), fill = Gene), alpha = 0.1,) +
          scale_fill_manual(values = color)
      }
      p <- p + labs(x=xlab, y=ylab, title=ifelse(is.na(main), feature, main))
      print(p)
    }
  }else{
    absent_genes <- c(absent_genes, feature)
  }
}


expr_table <- expr_table %>%
  mutate(Expr = as.numeric(Expr)) 

########

centroids <- expr_table %>%
  group_by(Bin) %>%
  summarise(
    Gene = "ZZZCentroid",
    Expr = mean(Expr, na.rm = TRUE),
    Std = NA_real_
  )

expr_table <- expr_table %>%
  mutate(Std = as.numeric(Std))  # Conversion en numérique

# Fusion des données originales avec les centroïdes
df_final <- bind_rows(expr_table, centroids) %>%
  arrange(Bin, Gene)

# Affichage du résultat
print(df_final)

# Création du vecteur de couleurs
color <- rep("grey", length(unique(df_final$Gene)))  # Tous en gris
names(color) <- unique(df_final$Gene)  # Associer les noms des gènes

# Mettre "Centroid" en noir
color["ZZZCentroid"] <- "black"






df_final <- df_final %>%
  arrange(Gene != "ZZZCentroid")  # Met "Centroid" à la fin








# Draw merged plot if required - - - - - - - - - - - - - - - - - - - - - - - - 
if(superposed){
  p <- ggplot(df_final, 
              aes(x = as.numeric(Bin), group = Gene)) +
    geom_line(linewidth = lwd, aes(y = as.numeric(Expr), color = Gene)) +
    scale_color_manual(values = color) +
    theme_classic() +
    ylim(ylim) 
  if(std){
    p <- p + geom_ribbon(aes(y = as.numeric(Expr), 
                             ymin = as.numeric(Expr) - as.numeric(Std), 
                             ymax = as.numeric(Expr) + as.numeric(Std),
                             fill = Gene), alpha = 0.1,) +
      scale_fill_manual(values = color)
  }
  p <- p + labs(x=xlab, y=ylab, title=ifelse(is.na(main), '', main))
  print(p)
}
return(expr_table)







color1 <- rep("grey",37)
DrawExpr(table_ILC, bin.number = 11, cluster1, 
         by.order = FALSE, std = FALSE, scale = FALSE, superposed = TRUE, 
         color = color1, lwd = 2, ylim = c(-1.5, 2), 
         xlab = "Pseudotime", ylab = "", main = "Expression ILC")
























  
