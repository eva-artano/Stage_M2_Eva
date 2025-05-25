#-------------------------------------------------------------------------------
# ENVIRONNEMENT
#-------------------------------------------------------------------------------

path <- 'D:/Resultats_memoire/scRNAseq/WT/Saves'
LIST <- 'D:/Resultats_memoire/scRNAseq/WT/FT'

#-------------------------------------------------------------------------------
# IMPORT PACKAGES
#-------------------------------------------------------------------------------
source("D:/Resultats_memoire/Custom_Functions.R")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(simspec))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(SingleCellExperiment))

#-------------------------------------------------------------------------------
# IMPORT DATA
#-------------------------------------------------------------------------------

data <- readRDS(paste0(path, '/5_Pseudotime.rds'))


test <- ProcessData(data, labels = c("DC", "ILC") , trajectory = c("Slingshot_2", "Slingshot_1"), type = 1 )

TFlist <- read.csv(paste0(LIST, '/FT_VariableFeature.txt'),sep="", header=FALSE)
TFlist <- TFlist[,-2]

valid_TF <- TFlist[TFlist %in% colnames(test)]

gene_significatifs <- c()
gene <- "Tox2"
result <- CompareExpr(x=test, bin.number=12, feature = gene , by.order = F,
                      scale = T, min.cells=5, stat ='t', color=c("#059025","#8a26bd"),
                      lwd=2, ylim=NULL,xlab='Pseudotime',ylab='Expression',
                      main = gene , show = T)


#-------------------------------------------------------------------------------
# Une seule pvalue

for (gene in valid_TF){
  result <- CompareExpr(x=test, bin.number=13, feature = gene, by.order = F,
                        scale = T, min.cells=5, stat ='t', color=rep(ColorBlind,10),
                        lwd=2, ylim=NULL,xlab='Pseudotime',ylab='Expression',
                        main = gene, show = F)
  
  result_df <- result
  
  if(any(result_df$DC_vs_ILC < 0.05, na.rm = TRUE)){
    gene_significatifs <- c(gene_significatifs, gene)
  }
  
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Au moins 2 pvalue consecutive


for (gene in valid_TF){
  result <- CompareExpr(x=test, bin.number=13, feature = gene, by.order = F,
                        scale = T, min.cells=5, stat ='t', color=rep(ColorBlind,10),
                        lwd=2, ylim=NULL, xlab='Pseudotime', ylab='Expression',
                        main = gene, show = F)
  
  result_df <- result
  
  pvals <- result_df$DC_vs_ILC
  if (any(pvals[-1] < 0.05 & pvals[-length(pvals)] < 0.05, na.rm = TRUE)){
    gene_significatifs <- c(gene_significatifs, gene)
  }
}



#-------------------------------------------------------------------------------
setwd('D:/Resultats_memoire/scRNAseq/WT/CompareExpr')
print(gene_significatifs)
diff <- setdiff(valid_TF, gene_significatifs)
print(diff)
write.table(gene_significatifs, "CompareTrajectory_signi.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(diff, "CompareTrajectory_not_signi.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

