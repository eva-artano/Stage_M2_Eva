#-------------------------------------------------------------------------------
# Environnement
#-------------------------------------------------------------------------------
LIST <- 'D:/Resultats_memoire/scRNAseq/WT/FT'
PATH <- 'D:/Resultats_memoire/scRNAseq/WT/tradeSeq/Data'
PATH_SAVE <- 'D:/Resultats_memoire/scRNAseq/WT/tradeSeq/Saves'

#-------------------------------------------------------------------------------
# IMPORT PACKAGES
#-------------------------------------------------------------------------------
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(clusterExperiment))

suppressPackageStartupMessages(library(slingshot))
`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------

#Importation liste de facteurs de transcription

TFlist <- read.csv(paste0(LIST, '/FT_VariableFeature.txt'),sep="", header=FALSE)

TFlist <- TFlist[,-2]
ft_names <- TFlist


# Counts
data <- readRDS(paste0(PATH, '/4_Filtered.rds'))
counts <- as.matrix(data@assays[["RNA"]]@layers[["counts"]])


cell_names <- Cells(data)
features_names <- Features(data)
rownames(counts) <- features_names
colnames(counts) <- cell_names
counts_FT <- counts[rownames(counts) %in% ft_names, ]


#Curves 
sce1 <- readRDS(paste0(PATH, '/Slingshot.rds'))
sce1_FT<- sce1[rownames(sce1) %in% ft_names, ]
crv <- as.SlingshotDataSet(sce1_FT)


sce <- readRDS(paste0(PATH, '/3_tradeSeq_ALL.rds'))
table(rowData(sce)$tradeSeq$converged)
gene_names_sce1 <- rownames(sce1)
rownames(sce) <- gene_names_sce1
sce_FT<- sce[rownames(sce) %in% ft_names, ]

#-------------------------------------------------------------------------------
# Discovering genes with different expression patterns
#-------------------------------------------------------------------------------

patternRes <- patternTest(sce_FT)
patternRes <- patternRes[order(patternRes$waldStat, decreasing = TRUE),]

patternRes$padjust <-p.adjust(patternRes$pvalue)
patternRes_filtered <- patternRes[patternRes$padjust < 0.05 & !is.na(patternRes$padjust), ]

getwd()
setwd('D:/Resultats_memoire/scRNAseq/WT/tradeSeq')

FT_list <- rownames(patternRes_filtered)
write.table(FT_list, "FT_list_patternRes_signi.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
diff <- setdiff(FT_list, ft_names)
print(diff)

patternRes_filtered$Gene <- row.names(patternRes_filtered)

Rank_tradeseq <- subset(patternRes_filtered, select = c(Gene, waldStat))
Rank_tradeseq$rank_tradeseq <- 1:nrow(Rank_tradeseq)

Rank_tradeseq <- Rank_tradeseq[,-2]
write.table(Rank_tradeseq, file = "D:/Resultats_memoire/scRNAseq/WT/tradeSeq/Rank_trade.txt" , sep = " ", quote = FALSE, row.names =  FALSE)



output_dir <- 'D:/Resultats_memoire/scRNAseq/WT/tradeSeq/TradeSeq_candidats/'

# Boucle pour tracer les deux types de graphiques pour chaque gène dans top_genes
for (gene in FT_list) {
  
  # Nom des fichiers pour les sauvegardes dans le dossier spécifié
  pseudotime_file <- paste0(output_dir, gene, "_p.png")

  # Tracer et sauvegarder le plot avec pseudotime
  png(pseudotime_file)
  print(plotSmoothers(sce_FT, counts_FT, gene = gene))
  dev.off()
  
}

#-------------------------------------------------------------------------------
# Distribution WaldTest en fonction des gènes
#-------------------------------------------------------------------------------

#All
library(ggplot2)

patternRes$rank <- rank(-patternRes$waldStat, ties.method = "first")

ggplot(patternRes, aes(x = rank, y = waldStat)) +
  geom_line(color = "#0072B2") +
  labs(
    title = "Distribution des valeurs de Wald test",
    x = "Rang des gènes (tri décroissant)",
    y = "Statistique de Wald"
  ) +
  theme_minimal()

#Only signi

patternRes_filtered$rank <- rank(-patternRes_filtered$waldStat, ties.method = "first")

ggplot(patternRes_filtered, aes(x = rank, y = waldStat)) +
  geom_line(color = "#0072B2") +
  labs(
    title = "Distribution des valeurs de Wald test",
    x = "Rang des gènes (tri décroissant)",
    y = "Statistique de Wald"
  ) +
  theme_minimal()

#-------------------------------------------------------------------------------
# TEST PLOT SMOOTHERS
#-------------------------------------------------------------------------------

curvesCols <- c("#8a26bd","#059025")

#plotSmoothers(sce_FT, counts_FT, gene = "Tox", curvesCols = curvesCols, border = FALSE, pointCol = pointCol)
plotSmoothers(sce_FT, counts_FT, gene = "Tox2", curvesCols = curvesCols, border = FALSE) +
  ggplot2::scale_color_manual(values= curvesCols)

plotSmoothers(sce_FT, counts_FT, gene="Tox2")


#-------------------------------------------------------------------------------
# DRAW ALL OUTLIERS
#-------------------------------------------------------------------------------

FT_list_signi <- rownames(patternRes_filtered)
FT_list <- rownames(patternRes)
  
FT_list_not_signi <- setdiff(FT_list, FT_list_signi)



output_dir <- 'D:/Resultats_memoire/scRNAseq/WT/tradeSeq/TradeSeq_candidats/tradeseq_not_signi/'

# Boucle pour tracer les deux types de graphiques pour chaque gène dans top_genes
for (gene in FT_list_not_signi) {
  
  # Nom des fichiers pour les sauvegardes dans le dossier spécifié
  pseudotime_file <- paste0(output_dir, gene, "_not_signi.png")
  
  # Tracer et sauvegarder le plot avec pseudotime
  png(pseudotime_file)
  curvesCols <- c("#FF9933", "#ffbbab")
  print(plotSmoothers(sce_FT, counts_FT, gene = gene), curvesCols = curvesCols, border = FALSE) +
    ggplot2::scale_color_manual(values=curvesCols)
  dev.off()
  
}




























