#-------------------------------------------------------------------------------
# ENVIRONNEMENT
#-------------------------------------------------------------------------------

path <- 'D:/Resultats_memoire/scRNAseq/WT-KO/Saves'
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


test <- ProcessData(data, celltypes = c("overLIP", "TOX_KO") , trajectory = "Slingshot_1", type = 2 )

TOX_candidats <- read.csv(file="D:/Resultats_memoire/scRNAseq/WT/LEAP/ILC/candidate_ILC_TOX.txt", sep ="", header = FALSE)
TOX_candidats <- TOX_candidats$V1

gene <- "Rora"

result <- CompareExpr(x=test, bin.number=12, feature = gene, by.order = F,
                      scale = F, min.cells=5, stat ='t', color=c("#a1c8eb","#ff0000"),
                      lwd=2, ylim=NULL,xlab='Pseudotime',ylab='Expression',
                      main = gene, show = T)
print(result)

#-------------------------------------------------------------------------------
# FOR ALL GENE 
#-------------------------------------------------------------------------------
gene_significatifs <- c()

for (gene in TOX_candidats){
  result <- CompareExpr(x=test, bin.number=13, feature = gene, by.order = F,
                        scale = T, min.cells=5, stat ='t', color=rep(ColorBlind,10),
                        lwd=2, ylim=NULL, xlab='Pseudotime', ylab='Expression',
                        main = gene, show = F)
  
  result_df <- result
  
  pvals <- result_df$overLIP_vs_TOX
  if (any(pvals[-1] < 0.05 & pvals[-length(pvals)] < 0.05, na.rm = TRUE)){
    gene_significatifs <- c(gene_significatifs, gene)
  }
}


print(gene_significatifs)
write.table(gene_significatifs, "D:/Resultats_memoire/scRNAseq/WT-KO/COmpareExpr()/Candidats_compare_expr_TOX.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


#-------------------------------------------------------------------------------
# New version of Compareexpr
#-------------------------------------------------------------------------------
TOX_candidats <- read.csv(file="D:/Resultats_memoire/scRNAseq/WT/LEAP/ILC/candidate_ILC_TOX.txt", sep ="", header = FALSE)
TOX_candidats <- TOX_candidats$V1

gene_significatifs <- c()

# Une seule pvalue

for (gene in TOX_candidats){
  result <- CompareExpr_single_pvalue(x=test, bin.number=13, feature = gene, by.order = F,
                                      scale = T, min.cells=5, stat ='t', color=rep(ColorBlind,10),
                                      lwd=2, ylim=NULL,xlab='Pseudotime',ylab='Expression',
                                      main = gene, show = F)
  
  result_df <- result
  
  if(result <  0.05){
    gene_significatifs <- c(gene_significatifs, gene)
  }
  
}
print(gene_significatifs)
write.table(gene_significatifs, "D:/Resultats_memoire/scRNAseq/WT-KO/COmpareExpr()/new_Candidats_compare_expr_TOX.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

old <- read.csv(file="D:/Resultats_memoire/scRNAseq/WT-KO/COmpareExpr()/Candidats_compare_expr_TOX.txt", header = FALSE, sep ="")
old <- old$V1

new <- read.csv(file="D:/Resultats_memoire/scRNAseq/WT-KO/COmpareExpr()/new_Candidats_compare_expr_TOX.txt", header = FALSE, sep ="")
new <- new$V1

commun <- intersect(new,old)
print(commun)
