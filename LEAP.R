#-------------------------------------------------------------------------------
# LEAP
#-------------------------------------------------------------------------------

PATH_SAVE <- 'D:/Resultats_memoire/scRNAseq/WT/Saves'
path_save1 <- 'D:/Resultats_memoire/scRNAseq/WT/LEAP/ILC/'
path_save2 <- 'D:/Resultats_memoire/scRNAseq/WT/LEAP/DC/'
source('D:/Resultats_memoire/Custom_Functions.R')

suppressPackageStartupMessages(library(Seurat))
# DRAWEPR DC -------------------------------------------------------------------

data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))

DC  <- subset(data, Slingshot_2 %!in% NA)
ILC <- subset(data, Slingshot_1 %!in% NA)

table_ILC <- OrderMatrix(ILC,"Slingshot_1", min.cells=50, slot='data')
table_DC  <- OrderMatrix(DC,"Slingshot_2", min.cells=50, slot='data')


feature.list1 <- c("Nfil3","Id2")
color1 <-  ColorBlind[c(10,4)]

plot_ILC <- DrawExpr(table_ILC, bin.number = 13, feature.list1, 
                     by.order = FALSE, std = FALSE, scale = FALSE, superposed = TRUE, 
                     color = color1, lwd = 2, ylim = c(-1.5, 4), 
                     xlab = "Pseudotime", ylab = "", main = "Expression ILC")

plot_DC <- DrawExpr(table_DC, bin.number = 13, feature.list1, 
                    by.order = FALSE, std = FALSE, scale = FALSE, superposed = TRUE, 
                    color = color1, lwd = 2, ylim = c(-1.5, 4), 
                    xlab = "Pseudotime", ylab = "", main = "Expression DC")

# ILC---------------------------------------------------------------------------

MAC <- read.table(paste0(path_save1,"MAC_results_ILC.txt"), header = TRUE, sep = "")
index <- read.table(paste0(path_save1,"gene_index_ILC.txt"), header = TRUE, sep = "")

LEAP_ILC <- AnnotateLeap(MAC=MAC, index=index, write = FALSE, dir = getwd(), 
                         filename = 'results_indexed')
perm_ILC <- read.csv(paste0(path_save1,"perm_ILC.csv"), header = TRUE, sep=",")

LEAP_ILC = LEAP_ILC[abs(LEAP_ILC$Correlation) > 0.1700,]
LEAP_ILC = LEAP_ILC[LEAP_ILC$Lag != 0,]

#Nfil3
LEAP_ILC_NFIL3 = LEAP_ILC[LEAP_ILC$gene_row == "Nfil3",]

candidate_NFIL3_ILC <- LEAP_ILC_NFIL3$gene_col
write.table(candidate_NFIL3_ILC, file='D:/Resultats_memoire/scRNAseq/WT/LEAP/ILC/candidate_ILC_NFIL3.txt',quote=FALSE, row.names=FALSE,col.names=FALSE)

#Tox

LEAP_ILC_TOX = LEAP_ILC[LEAP_ILC$gene_row == "Tox",]

candidate_TOX_ILC <- LEAP_ILC_TOX$gene_col
write.table(candidate_TOX_ILC, file='D:/Resultats_memoire/scRNAseq/WT/LEAP/ILC/candidate_ILC_TOX.txt',quote=FALSE, row.names=FALSE,col.names=FALSE)

# DC----------------------------------------------------------------------------


MAC <- read.table(paste0(path_save2,"MAC_results_DC.txt"), header = TRUE, sep = "")
index <- read.table(paste0(path_save2,"gene_index_DC.txt"), header = TRUE, sep = "")

LEAP_DC <- AnnotateLeap(MAC=MAC, index=index, write = FALSE, dir = getwd(), 
                         filename = 'results_indexed')
perm_DC <- read.csv(paste0(path_save2,"perm_DC.csv"), header = TRUE, sep=",")

LEAP_DC = LEAP_DC[abs(LEAP_DC$Correlation) > 0.1600,]
LEAP_DC = LEAP_DC[LEAP_DC$Lag != 0,]



LEAP_DC_NFIL3 = LEAP_DC[LEAP_DC$gene_row == "Nfil3",]

candidate_DC_NFIL3<- LEAP_DC_NFIL3$gene_col
write.table(candidate_DC_NFIL3, file='D:/Resultats_memoire/scRNAseq/WT/LEAP/DC/candidate_DC_NFIL3.txt',quote=FALSE, row.names=FALSE,col.names=FALSE)

#LEAP_DC_ID2 = LEAP_DC[LEAP_DC$gene_row == "Id2",]
#LEAP_DC_TCF7 = LEAP_DC[LEAP_DC$gene_row == "Tcf7",]
#LEAP_DC_IRF8 = LEAP_DC[LEAP_DC$gene_row == "Irf8",]
#LEAP_DC_ZEB2 = LEAP_DC[LEAP_DC$gene_row == "Zeb2",]

#-------------------------------------------------------------------------------

#On croise les listes de gegne candidats ILC/DC
candidate_DC_NFIL3 = 'D:/Eva/Stage_M2/scRNAseq/Test/LEAP/DC/candidate_DC_NFIL3.txt'
candidate_ILC_NFIL3 = 'D:/Eva/Stage_M2/scRNAseq/Test/LEAP/ILC/candidate_ILC_NFIL3.txt'


candidate_DC_NFIL3<- read.csv(candidate_DC_NFIL3,sep="", header=FALSE)
candidate_ILC_NFIL3 <- read.csv(candidate_ILC_NFIL3, sep="", header=FALSE)
col1_DC_NFIL3 <- candidate_DC_NFIL3[,1]
col1_ILC_NFIL3 <- candidate_ILC_NFIL3[,1]
gene_commun <- intersect(col1_DC_NFIL3, col1_ILC_NFIL3 )
write.table(gene_commun, "D:/Eva/Stage_M2/scRNAseq/Test/LEAP/Gene_commun_NFIL3_DC_ILC.txt",quote=FALSE, row.names = FALSE, col.names=FALSE)

gene_diff_NFIL3_DC <- setdiff(col1_DC_NFIL3,col1_ILC_NFIL3)
write.table(gene_diff_NFIL3_DC, "Gene_diff_NFIL3_DC.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

gene_diff_NFIL3_ILC <- setdiff(col1DOWNTOX, col1DOWNALP)
gene_diff_TOX_DOWN <- intersect(gene_DOWN_diff,col1DOWNTOX)
write.table(gene_diff_TOX_DOWN, "Gene_DOWN_diff_TOX_KO.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)





























