#-------------------------------------------------------------------------------
# Preparation de l'environnement
#-------------------------------------------------------------------------------

path <- 'D:/Resultats_memoire/scRNAseq'
path_save <- 'D:/Resultats_memoire/scRNAseq/WT'


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

CellTypeOrder <- c('ALP', 'overLIP', 'TULIP', 'ILCpro')

dir.create(file.path(path_save, '/Figures'))
dir.create(file.path(path_save, '/Saves'))

PATH <- paste0(path, '/Raw/ILC/scRNA-seq_Raw')
PATH2 <- paste0(path,'/Raw/cDC1/raw_cDC1')

PATH_SAVE <- paste0(path_save, '/Saves')
PATH_FIG <- paste0(path_save, '/Figures')

#-------------------------------------------------------------------------------
# IMPORT DATA AND TURNING THEM INTO SEURAT OBJECT
#-------------------------------------------------------------------------------

#Experiment 1 ------------------------------------------------------------------

# All Lymphoid Progenitor
ALP <- Read10X(data.dir = paste0(PATH,'/1_ABSC_CLP')) %>%
  CreateSeuratObject(project = 'ALP', min.cells = 1, min.features = 0 )

#sEILP, cEILP and ILC progenitors
ILCpro <- Read10X(data.dir = paste0(PATH,'/3_ABSC_EILP_WT')) %>%
  CreateSeuratObject(project = 'ILCpro', min.cells = 1, min.features = 0)

# a4b7 intermediates
TULIP <- Read10X(data.dir = paste0(PATH,'/2_ABSC_TULIP')) %>%
  CreateSeuratObject(project = 'TULIP', min.cells = 1, min.features = 0)

#Experiment 2 ------------------------------------------------------------------

# Il7rLT+ XT
overLIP <-  Read10X(data.dir = paste0(PATH,'/2019_scRNAseq_data/WT')) %>%
  CreateSeuratObject(project = 'overLIP', min.cells = 1, min.features = 0)

# Spread datasets in groups of experiment
ALP@meta.data[['Experiment']]     <- 'Experiment 1'
ILCpro@meta.data[['Experiment']]  <- 'Experiment 1'
TULIP@meta.data[['Experiment']]   <- 'Experiment 1'
overLIP@meta.data[['Experiment']] <- 'Experiment 2'

# Merge all the experiments to perform the quality check 
data  <- merge(TULIP, c(ILCpro, ALP, overLIP), 
               add.cell.ids = c('TULIP', 'ILCpro', 'ALP', 'overLIP'), 
               project = 'scRNA_Integration')
data@meta.data[['Project']] <- 'scRNA_ILC'

# Organize cell type in the development order for the plots
data@active.ident <- factor(data@active.ident, levels = CellTypeOrder)

#Experiment 3 ------------------------------------------------------------------

#Run 1
CDP_1 <- Read10X(data.dir = paste0(PATH2,'/Run1')) %>%
  CreateSeuratObject(project = 'CDP', min.cells = 1, min.features = 0 )

#Run 2 
CDP_2 <- Read10X(data.dir = paste0(PATH2,'/Run2')) %>%
  CreateSeuratObject(project = 'CDP', min.cells = 1, min.features = 0)

CDP_1@meta.data[['Experiment']]   <- 'Experiment 3'
CDP_2@meta.data[['Experiment']] <- 'Experiment 3'

# Merge all the experiments to perform the quality check 
CDP  <- merge(CDP_1, CDP_2, 
              add.cell.ids = c('CDP',2), 
              project = 'scRNA_DC')

#On garde uniquement les cDC1
cell_names <- read.table(paste0(PATH_SAVE, "/cDC1_to_keep.txt"), header = FALSE)
cell_names <- cell_names$V1


# Créer un sous-ensemble de l'objet Seurat en fonction des cellules à garder
cDC1 <- CDP[, cell_names]

print(cell_names)

# Merge all
combined_data <- merge(data, cDC1 , project = 'scRNA_Combined')


#-------------------------------------------------------------------------------
# TEST
#-------------------------------------------------------------------------------

# based on : https://github.com/satijalab/seurat/issues/7191#issuecomment-1636319705

#extraire les noms des gènes et les noms des cellules dans des vecteurs

data1 <- JoinLayers(combined_data)
var <-LayerData(data1, assay="RNA", layer = "counts")
data1@assays[["RNA"]]@layers[["counts"]]@Dimnames <- dimnames(var)

#-------------------------------------------------------------------------------
# QUALITY CHECK 
#-------------------------------------------------------------------------------

# Calculate the percent of mitochondrial genes expression
data1$Percent.mt <- PercentageFeatureSet(data1, pattern = '^mt-')

plot0   <- VlnPlot(data1, features = c('nFeature_RNA','Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1,3,4)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot0, PATH_FIG)

data1@assays[["RNA"]]@cells
cell_names <- Cells(data1)
print(cell_names)
features_names <- Features(data1)
print(features_names)

# Elimination of cells that over-expressed mitochondrial genes (dead cells)
data.QC <- subset(data1, subset = ((Experiment == 'Experiment 1' & 
                                      nFeature_RNA > 200 & 
                                      nFeature_RNA < 5000 & Percent.mt < 2) |
                                     (Experiment != 'Experiment 1' & 
                                        nFeature_RNA > 1200 & 
                                        nFeature_RNA < 5000 & Percent.mt < 4) |
                                     (Experiment == 'Experiment 3' & 
                                        nFeature_RNA > 1000 & 
                                        nFeature_RNA < 4000 & Percent.mt < 6))) 

plot1   <- VlnPlot(data.QC, features = c('nFeature_RNA', 'Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1,3,4)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot1, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after quality check
saveRDS(data.QC, paste0(PATH_SAVE, '/1_AfterQC.rds'))
data.QC <- readRDS(paste0(PATH_SAVE, '/1_AfterQC.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#-------------------------------------------------------------------------------
# NORMALIZATION AND SCALING
#-------------------------------------------------------------------------------

data.QC <- NormalizeData(data.QC, normalization.method = 'LogNormalize')

# Cell Cycle Genes identification and scoring
s.genes       <- capitalize(tolower(cc.genes$s.genes))
g2m.genes     <- capitalize(tolower(cc.genes$g2m.genes))

# CellCycleScoring() dont work with Seurat V5 
# based on : https://github.com/satijalab/seurat/issues/7191#issuecomment-1636319705

#data.combined[["joined"]] <- JoinLayers(data.combined[["RNA"]])
#DefaultAssay(data.combined) <- "joined"

data.QC <- CellCycleScoring(data.QC, s.features = s.genes, 
                            g2m.features = g2m.genes, set.ident = TRUE)


# Calculate difference between S and G2M is described as more relevant
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Alternate Workflow
data.QC$CC.Difference <- data.QC$S.Score - data.QC$G2M.Score
data.QC <- ScaleData(data.QC, 
                     vars.to.regress = c('CC.Difference','Percent.mt'), 
                     verbose = TRUE, features = rownames(data.QC), 
                     
                     do.scale = TRUE)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after scaling step
saveRDS(data.QC, paste0(PATH_SAVE, '/2_Scaled.rds'))
data.combined <- readRDS(paste0(PATH_SAVE, '/2_Scaled.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#===============================================================================
## CLUSTERING SIMILARITY SPECTRUM INTEGRATION ----------------------------------
#===============================================================================

# Before Integration - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

data.combined <- FindVariableFeatures(data.combined, nfeatures = 3000)
data.combined <- RunPCA(data.combined, npcs = 50)

plot2         <- ElbowPlot(data.combined, 30) +
  labs(title = '', y = 'Standard deviation', x = 'Component number') + 
  theme(plot.title = element_text(hjust = 0.5))
writePlot(plot2, PATH_FIG)

data.combined <- RunUMAP(data.combined, dims = 1:16, metric = 'euclidean', 
                         n.neighbors = 50)

plot3         <- DimPlot(data.combined, group.by = 'Experiment', 
                         reduction = 'umap', cols = ColorBlind[c(1,3,10)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot3, PATH_FIG)

FeaturePlot(data.combined,feature="Cd74", reduction="umap")

FeaturePlot(data.combined,feature="Id2", reduction="umap")


plot4         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(7,10,3,1,9)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot4, PATH_FIG)

# After Integration Z-Transform  - - - - - - - - - - - - - - - - - - - - - - - -  

data          <- cluster_sim_spectrum(data.combined, label_tag='Experiment', 
                                      cluster_resolution = 0.6, 
                                      corr_method = 'spearman', lambda = 50, 
                                      reduction.name = 'cssz', 
                                      reduction.key = 'CSSZ_')
data          <- RunUMAP(data, reduction = 'cssz', 
                         dims = 1:ncol(Embeddings(data, 'cssz')), 
                         n.neighbors = 50, reduction.name='umap_cssz', 
                         reduction.key='UMAPCSSZ_')

plot5         <- DimPlot(data, group.by = 'Experiment', reduction = 'umap_cssz', 
                         cols = ColorBlind[c(1,3,10)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') 
writePlot(plot5, PATH_FIG)

plot6         <- DimPlot(data, group.by = 'orig.ident', reduction = 'umap_cssz', 
                         cols = ColorBlind[c(7,10,3,1,9)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') 
writePlot(plot6, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after integration
saveRDS(data, paste0(PATH_SAVE, '/3_Integrated.rds'))
data <- readRDS(paste0(PATH_SAVE, '/3_Integrated.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#-------------------------------------------------------------------------------
# CLUSTERING AND OUTLIER ELIMINATION 
#-------------------------------------------------------------------------------

data          <- FindNeighbors(data, reduction = 'umap_cssz', 
                               dims = 1:ncol(Embeddings(data, 'umap_cssz')))
data          <- FindClusters(data, resolution = 0.9)

plot7         <- DimPlot(data, reduction = 'umap_cssz', 
                         group.by = 'seurat_clusters', label = TRUE, 
                         repel = TRUE, cols = c(rep(ColorBlind,10)), 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2') 
writePlot(plot7, PATH_FIG)

#Keep cellnames cluster 3 and 26

data_3 <- subset(data, seurat_clusters %in% 3)
cellnames_3 <- colnames(data_3)

data_26 <- subset(data, seurat_clusters %in% 26)
cellnames_26 <- colnames(data_26)

cell_names_cluster_7 <- read.table("D:/Resultats_memoire/scRNAseq/Raw/cDC1/Saves/cells_cluster_7.txt", header = FALSE)
cell_names_cluster_7 <- cell_names_cluster_7$V1

commun_3 <- intersect (cell_names_cluster_7,cellnames_3)
commun_26 <- intersect (cell_names_cluster_7,cellnames_26)


plot8 <- FeaturePlot(data, feature = "Id2", reduction = 'umap_cssz', order = T)
writePlot(plot8, PATH_FIG)

plot9<-FeaturePlot(data, feature = "Zeb2", reduction = 'umap_cssz', order = T)
writePlot(plot9, PATH_FIG)


FeaturePlot(data, feature = "Nfil3", reduction = 'umap_cssz')
FeaturePlot(data, feature = "Tox", reduction = 'umap_cssz')

# Identification of cell types
clusterMarker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
top30         <- clusterMarker %>% group_by(cluster) %>% top_n(30, avg_log2FC) 
write.csv(top30, paste0(PATH_SAVE, '/WT_Top30markers.csv'), row.names = FALSE)

# For each cluster
for(cluster_id in unique(top30$cluster)) {
  cluster_data <- top30[top30$cluster == cluster_id, "gene" ]
  file_name <- paste0(PATH_SAVE, '/WT_Top30markers_cluster_', cluster_id, '.txt')
  write.table(cluster_data, file_name, row.names = FALSE, col.names = FALSE, quote = FALSE, sep="")}

# Lecture des noms de cellules du cluster 7
cell_names_cluster_7 <- read.table("D:/Resultats_memoire/scRNAseq/Raw/cDC1/Saves/cells_cluster_7.txt", header = FALSE)
cell_names_cluster_7 <- cell_names_cluster_7$V1

# Création de la colonne "ColorTag" dans meta.data avec des noms explicites
data$ColorTag <- ifelse(rownames(data@meta.data) %in% cell_names_cluster_7, "Cluster 7", "Autres")

# Choix des couleurs pour l'affichage
color_mapping <- c("Cluster 7" = "red", "Autres" = "gray")

# UMAP avec couleurs personnalisées et légende claire
DimPlot(data, group.by = "ColorTag", reduction = "umap_cssz", cols = color_mapping)


## DATA FILTERING --------------------------------------------------------------

# Remove mature cells (22,32), added B cells (27) and contaminants (28,29,30,31,33,36)
Inliers   <- subset(data, seurat_clusters %!in% c(3,13,22,23,24,25,26,27))
data  <- Inliers


Outliers_CDP <- subset(data, seurat_clusters %in% c(3,26))
cells_name_outliers_CDP <- colnames(Outliers_CDP)
write.table(cells_name_outliers_CDP, paste0(PATH_SAVE,'/outliers_CDP.txt'), row.names= FALSE, col.names = FALSE, quote = FALSE, sep="")


data  <- RunUMAP(data, reduction = 'cssz', 
                 dims = 1:ncol(Embeddings(data, 'cssz')), n.neighbors = 50, 
                 reduction.name='umap_cssz2', reduction.key='UMAPCSSZ2_')
data  <- FindNeighbors(data, reduction = 'umap_cssz2', 
                       dims = 1:ncol(Embeddings(data, 'umap_cssz2')))
data  <- FindClusters(data, resolution = 0.08)

plot11 <- DimPlot(data, reduction = 'umap_cssz2', 
                  group.by = 'seurat_clusters', 
                  label = TRUE, repel = TRUE, 
                  cols = c(rep(ColorBlind,10), "purple"), 
                  pt.size = 1) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot11, PATH_FIG)

plot12         <- DimPlot(data, group.by = 'orig.ident', reduction = 'umap_cssz2', 
                          cols = ColorBlind[c(12,10,3,1,9)], pt.size = 0.5, 
                          order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') 
writePlot(plot12, PATH_FIG)

FeaturePlot(data,feature="Icos", reduction="umap_cssz", order = T, pt.size = 0.5)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
saveRDS(data, paste0(PATH_SAVE, '/4_Filtered.rds'))
data <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


plot_density(data, features = "Flt3", pal='magma', size=2)
plot_density(data, features = "Tcf7", pal='magma', size=2)
plot_density(data, features = "Cd74", pal='magma', size=2)

cell_names_cluster_7 <- read.table("D:/Resultats_memoire/scRNAseq/Raw/cDC1/Saves/cells_cluster_7.txt", header = FALSE)
cell_names_cluster_7 <- cell_names_cluster_7$V1

#===============================================================================
## SLINGSHOT -------------------------------------------------------------------
#===============================================================================
#Based on : https://bustools.github.io/BUS_notebooks_R/slingshot.html

DimPlot(data, reduction = 'umap_cssz2', group.by = 'seurat_clusters', 
        label = TRUE, repel = TRUE, 
        cols = c(rep(ColorBlind,10)), pt.size = 0.9) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')


## PSEUDOTIME RECONSTRUCTION ---------------------------------------------------
#https://github.com/satijalab/seurat/issues/8248

sce                                     <- as.SingleCellExperiment(data)     
colData(sce)$Seurat_clusters            <- as.character(data@active.ident)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', 
                 reducedDim = toupper('umap_cssz2'), start.clus = 4, end.clus = c(0,3))


## ADDING PSEUDOTIMES TO DATA --------------------------------------------------
pseudotime                              <- slingPseudotime(sce)
for(i in 1:ncol(pseudotime)){
  data@meta.data[[paste0('Slingshot_', i)]] <- pseudotime[,i]
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Saving SCE file after processing
saveRDS(sce, paste0(PATH_SAVE, '/Slingshot.rds'))
sce <- readRDS(paste0(PATH_SAVE, '/Slingshot.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

## PLOTTING PSEUDOTIME TRAJECTORIES --------------------------------------------
red='umap_cssz2'
# Specifying gradient colors for plots
gradient         <- hcl.colors(100)

# Plot with curves
curves <- slingCurves(sce, as.df = TRUE) %>%
  dplyr::rename("x" = paste0(str_replace_all(
    toupper(red), "[[:punct:]]", ""), '_1'), 
    "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

#################


# First Trajectory
plot13 <- FeaturePlot(data, 'Slingshot_1', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) + # Remplacement de size par linewidth
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) 
writePlot(plot13, PATH_FIG)

# Second Trajectory
plot14 <- FeaturePlot(data, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient)
writePlot(plot14, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
saveRDS(data, paste0(PATH_SAVE, '/5_Pseudotime.rds'))
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
red='umap_cssz2'
# UMAP avec cellules colorées par 'orig.ident'
plot_umap_orig <- DimPlot(data, group.by = 'orig.ident', reduction = 'umap_cssz2', 
                          cols = ColorBlind[c(7,10,3,1,9)], pt.size = 1, 
                          order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') 


#extraire chaque trajectoire séparément
curves_df <- slingCurves(sce, as.df = TRUE)

# Séparer les deux trajectoires
curve1 <- curves_df %>% filter(Lineage == 1) %>% arrange(Order)
curve2 <- curves_df %>% filter(Lineage == 2) %>% arrange(Order)

curve1 <- curve1 %>%
  dplyr::rename("x" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_1'),
                "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

curve2 <- curve2 %>%
  dplyr::rename("x" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_1'),
                "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

plot_final <- plot_umap_orig +
  geom_path(data = curve1, aes(x = x, y = y), color = "#8a26bd", linewidth = 1.5) +
  geom_path(data = curve2, aes(x = x, y = y), color = "#059025", linewidth = 1.5)

# Sauvegarde
writePlot(plot_final, PATH_FIG)

#-------------------------------------------------------------------------------
# DRAW EXPRESSION 
#-------------------------------------------------------------------------------

data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))

ILC  <- subset(data, Slingshot_1 %!in% NA)
DC <- subset(data, Slingshot_2 %!in% NA)

table_ILC <- OrderMatrix(ILC,"Slingshot_1", min.cells=50, slot='scale.data')
table_DC  <- OrderMatrix(DC,"Slingshot_2", min.cells=50, slot='scale.data')

feature.list1 <- c("Nfil3","Tox","Tcf7","Gata3","Id2","Zbtb16")
feature.list2 <- c("Irf8","H2-Aa","Zeb2","Cd74","Ly6d","Batf3","Id2")

color1 <-  ColorBlind[c(10,9,1,8,2,3)]
color2 <- ColorBlind[c(11,12,1,9,4,8,1,13,5,10)]

plot_ILC <- DrawExpr(table_ILC, bin.number = 11, feature.list2, 
                     by.order = FALSE, std = FALSE, scale = FALSE, superposed = TRUE, 
                     color = color2, lwd = 2, ylim = c(-1.5, 2), 
                     xlab = "Pseudotime", ylab = "", main = "Expression ILC")

plot_DC <- DrawExpr(table_DC, bin.number = 11, feature.list2, 
                    by.order = FALSE, std = FALSE, scale = FALSE, superposed = TRUE, 
                    color = color2, lwd = 2, ylim = c(-1.5, 2), 
                    xlab = "Pseudotime", ylab = "", main = "Expression DC")

#-------------------------------------------------------------------------------
# Find variable feature
#-------------------------------------------------------------------------------
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))
data <- FindVariableFeatures(data, nfeatures = 3000)
top10 <- head(VariableFeatures(data), 10)


new <- VariableFeatures(data)
write.table(new, "new_variable_feature.txt",quote=FALSE, row.names = FALSE, col.names=FALSE)

#-------------------------------------------------------------------------------
# Find variable FT 
#-------------------------------------------------------------------------------
PATH <- '/Users/E192693Z/Desktop/Eva/Stage_M2/scRNAseq/Test/FT/'


VariableFeature <- read.csv(paste0(PATH, 'new_variable_feature.txt'),sep="", header=FALSE)
VariableFeature <- VariableFeature[,1]

FT_list <- read.csv(paste0(PATH, 'masterTFlist.Txt'), sep="", header = FALSE)
FT_list <- FT_list[,1]

FT <- intersect(FT_list, VariableFeature)
write.table(FT, "FT_VariableFeature.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#-------------------------------------------------------------------------------
# HEATMAP ? 
#-------------------------------------------------------------------------------
suppressPackageStartupMessages(library(pheatmap))

#Gene
feature <- read.csv(paste0(path_save, '/FT/CompareTrajectory_signi_two.txt'),sep="", header=FALSE)
feature <- feature[,1]

#Preparation des data
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))

DC  <- subset(data, Slingshot_2 %!in% NA)
ILC <- subset(data, Slingshot_1 %!in% NA)

table_ILC <- OrderMatrix(ILC,"Slingshot_1", min.cells=50, slot='scale.data')
table_DC  <- OrderMatrix(DC,"Slingshot_2", min.cells=50, slot='scale.data')

#

genes_to_plot <- intersect(feature, colnames(table_ILC))
data_to_plot <- table_ILC[, genes_to_plot, drop = FALSE]

pheatmap(t(data_to_plot), 
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         show_rownames = TRUE,  
         show_colnames = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(10))


#

genes_to_plot <- intersect(feature, colnames(table_DC))
data_to_plot <- table_DC[, genes_to_plot, drop = FALSE]

# Créer une annotation pour le pseudotime
annotation_col <- data.frame(Pseudotime = table_DC$Pseudotime)
rownames(annotation_col) <- rownames(table_DC)  

# Générer une palette de couleurs pour le pseudotime
pseudotime_colors <- hcl.colors(100)

# Associer les valeurs de pseudotime à une couleur
annotation_colors <- list(Pseudotime = pseudotime_colors)


pheatmap(t(data_to_plot),  
         cluster_rows = TRUE,  
         cluster_cols = FALSE,  
         show_rownames = TRUE,  
         show_colnames = FALSE, 
         annotation_col = annotation_col,  
         annotation_colors = annotation_colors,
         scale = "none",
         color = colorRampPalette(c("darkblue", "yellow"))(100))  


#-------------------------------------------------------------------------------
#Tentative de mettre DC et ILC sur la même heatmap

table_ILC <- OrderMatrix(ILC,"Slingshot_1", min.cells=50, slot='scale.data')
table_DC  <- OrderMatrix(DC,"Slingshot_2", min.cells=50, slot='scale.data')


cells_ILC <- rownames(table_ILC)
cells_DC <- rownames(table_DC)

common_cells <- intersect(cells_ILC, cells_DC)

table_DC_filtered <- table_DC[!rownames(table_DC) %in% common_cells, ]

table_DC_new <- rev(table_DC_filtered)

#valeur neg
table_DC_filtered$Pseudotime <- -table_DC_filtered$Pseudotime

merged_table <- rbind(table_ILC, table_DC_filtered)

merged_table <- merged_table[order(merged_table$Pseudotime, decreasing = FALSE), ]


genes_to_plot <- intersect(feature, colnames(merged_table))
data_to_plot <- merged_table[, genes_to_plot, drop = FALSE]

# Créer une annotation pour le pseudotime
annotation_col <- data.frame(Pseudotime = merged_table$Pseudotime)
rownames(annotation_col) <- rownames(merged_table)  

# Générer une palette de couleurs pour le pseudotime
pseudotime_colors <- hcl.colors(100)

# Associer les valeurs de pseudotime à une couleur
annotation_colors <- list(Pseudotime = pseudotime_colors)


pheatmap(t(data_to_plot),  
         cluster_rows = TRUE,  
         cluster_cols = FALSE,  
         show_rownames = TRUE,  
         show_colnames = FALSE, 
         annotation_col = annotation_col,  
         annotation_colors = annotation_colors,
         scale = "none",
         gaps_col = c(1185,2636),
         color = colorRampPalette(c("darkblue", "yellow"))(100))  

#############

table_ILC <- OrderMatrix(ILC,"Slingshot_1", min.cells=50, slot='scale.data')
table_DC  <- OrderMatrix(DC,"Slingshot_2", min.cells=50, slot='scale.data')


cells_ILC <- rownames(table_ILC)
cells_DC <- rownames(table_DC)

common_cells <- intersect(cells_ILC, cells_DC)

table_DC_filtered <- table_DC[!rownames(table_DC) %in% common_cells, ]
table_DC_reverse <- table_DC_filtered[order(table_DC_filtered$Pseudotime, decreasing = TRUE), ]

merged_table <- rbind(table_DC_reverse,table_ILC)


genes_to_plot <- intersect(feature, colnames(merged_table))
data_to_plot <- merged_table[, genes_to_plot, drop = FALSE]

# Créer une annotation pour le pseudotime
annotation_col <- data.frame(Pseudotime = merged_table$Pseudotime)
rownames(annotation_col) <- rownames(merged_table)  

# Générer une palette de couleurs pour le pseudotime
pseudotime_colors <- hcl.colors(100)

# Associer les valeurs de pseudotime à une couleur
annotation_colors <- list(Pseudotime = pseudotime_colors)

png("/Users/E192693Z/Desktop/Eva/Stage_M2/scRNAseq/Test/heatmap_Rplot2.png", width = 2113, height = 2000)
dev.list()

plot <- NULL


plot <- pheatmap(t(data_to_plot),  
         cluster_rows = TRUE,  
         cluster_cols = FALSE,  
         show_rownames = TRUE,  
         show_colnames = FALSE, 
         annotation_col = annotation_col,  
         annotation_colors = annotation_colors,
         scale = "none",
         gaps_col = c(1185,2636),
         cutree_rows=4,
         color = colorRampPalette(c("darkblue", "yellow"))(100))

dev.off()


#-------------------------------------------------------------------------------
# Autre test 

suppressPackageStartupMessages(library(pheatmap))

#Gene
feature <- read.csv(paste0(path_save, '/Figures/gene_cluster_heatmap.txt'),sep="", header=FALSE)
feature <- feature[,1]

#Preparation des data
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))

DC  <- subset(data, Slingshot_2 %!in% NA)
ILC <- subset(data, Slingshot_1 %!in% NA)

cell_id_DC <- data.frame(Ident = DC@meta.data[["orig.ident"]], row.names = rownames(DC@meta.data))
cell_id_ILC <- data.frame(Ident= ILC@meta.data[["orig.ident"]], row.names = rownames(ILC@meta.data))

table_ILC <- OrderMatrix(ILC,"Slingshot_1", min.cells=50, slot='scale.data')
table_DC  <- OrderMatrix(DC,"Slingshot_2", min.cells=50, slot='scale.data')

#Ajout de Ident aux DC
table_DC <- merge(table_DC, cell_id_DC, by = "row.names", all.x = TRUE)
rownames(table_DC) <- table_DC$Row.names
table_DC$Row.names <- NULL

table_DC_reverse <- table_DC[order(table_DC$Pseudotime, decreasing = TRUE), ]

#Ajout de Ident aux ILC 
table_ILC <- merge(table_ILC, cell_id_ILC, by = "row.names", all.x = TRUE)
rownames(table_ILC) <- table_ILC$Row.names
table_ILC$Row.names <- NULL

table_ILC <- table_ILC[order(table_ILC$Pseudotime), ]

merged_table <- rbind(table_DC_reverse,table_ILC)

genes_to_plot <- intersect(feature, colnames(merged_table))
data_to_plot <- merged_table[, genes_to_plot, drop = FALSE]

# Créer une annotation pour le pseudotime
annotation_col <- data.frame(Celltype = merged_table$Ident, Pseudotime = merged_table$Pseudotime)
rownames(annotation_col) <- rownames(merged_table)  

levels(as.factor(annotation_col$Celltype))


# Créer une annotation pour le type cellulaire 
ALP <- "#F0e442"
TULIP <- "#016ba0"
overLIP <- "#a1c8eb"
ILCpro <- "#FF9933"
CDP <- "#FF7070"

# Générer une palette de couleurs pour le pseudotime
pseudotime_colors <- hcl.colors(100)

# Associer les valeurs de pseudotime à une couleur
annotation_colors <- list(
  Pseudotime = pseudotime_colors,
  Celltype = c("ALP" = "#F0e442",
               "TULIP" = "#016ba0",
               "overLIP" = "#a1c8eb",
               "ILCpro" = "#FF9933",
               "CDP" = "#FF7070")
)



png("/Users/E192693Z/Desktop/Eva/Stage_M2/scRNAseq/Test/heatmap_finalV2.png", width = 2113, height = 1000)
dev.list()

plot <- NULL
plot <- pheatmap(t(data_to_plot),  
                 cluster_rows = FALSE,  
                 cluster_cols = FALSE,  
                 show_rownames = TRUE,  
                 show_colnames = FALSE, 
                 annotation_col = annotation_col,  
                 annotation_colors = annotation_colors,
                 scale = "none",
                 gaps_col = c(1185,2636,4087),
                 gaps_row = c(37,70),
                 color = colorRampPalette(c("darkblue","yellow"))(100))

dev.off()

#-------------------------------------------------------------------------------
# LEAP 
#-------------------------------------------------------------------------------
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))

ILC <- LeapMatrix(data,"Slingshot_1",min.cells=50, slot="data",write=T,dir=PATH_SAVE,suffix="ILC")
DC <- LeapMatrix(data,"Slingshot_2",min.cells=50, slot="data",write=T,dir=PATH_SAVE,suffix="DC")
getwd()

#-------------------------------------------------------------------------------
# CLUSTER 7
#-------------------------------------------------------------------------------

library(ggplot2)


library(ggplot2)
library(dplyr)

# Données
data <- data.frame(
  Cluster = c("3", "26"),
  Present = c(240, 34),
  Absent = c(39, 5)
)

# Calculs
data <- data %>%
  mutate(
    Total = Present + Absent,
    Pourcentage_Present = round((Present / Total) * 100)
  )

# Graphique
ggplot(data, aes(x = Cluster, y = Total)) +
  # Barre totale blanche
  geom_bar(stat = "identity", fill = "white", color = "black") +
  # Barre bleue représentant les absents
  geom_bar(aes(y = Absent), stat = "identity", fill = "darkblue") +
  # Nombre de présents affiché au milieu de la zone blanche
  geom_text(aes(y = Absent + (Present / 2), label = Present), color = "black", size = 4) +
  # Pourcentage de présents affiché juste sous le nombre
  geom_text(aes(y = Absent + (Present / 2) - 20, label = paste0(Pourcentage_Present, "%")), color = "black", size = 3.5) +
  # Total en haut de la barre
  geom_text(aes(y = Total + 15, label = paste0("n = ", Total)), size = 3.5) +
  theme_minimal(base_size = 14) +
  ylab("Nombre de gènes") +
  xlab("Cluster") +
  theme(
    panel.grid = element_blank(),     # Enlève le quadrillage
    legend.position = "none"
  )

#
library(ggplot2)
library(dplyr)
library(tidyr)

# Données
data <- data.frame(
  Cluster = c("3", "26"),
  Present = c(240, 34),
  Absent = c(39, 5)
)

# Format long pour ggplot
data_long <- data %>%
  pivot_longer(cols = c(Present, Absent), names_to = "Statut", values_to = "Nombre")

# Ordre pour empilement : Present (blanc) en dessous, Absent (bleu) au-dessus
data_long$Statut <- factor(data_long$Statut, levels = c("Absent", "Present"))

# Graphique
ggplot(data_long, aes(x = Cluster, y = Nombre, fill = Statut)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(
    values = c("Present" = "white", "Absent" = "darkblue"),
    labels = c("Absentes", "Présentes"),
    name = "Statut"
  ) +
  scale_y_continuous(
    breaks = seq(0, max(data$Present + data$Absent), by = 20),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 14) +
  xlab("Cluster") +
  ylab("Nombre de gènes") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "right"
  )








