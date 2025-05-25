#-------------------------------------------------------------------------------
# Preparation de l'environnement
#-------------------------------------------------------------------------------

path <- 'D:/Resultats_memoire/scRNAseq/Raw/cDC1/raw_cDC1'
path_save <- 'D:/Resultats_memoire/scRNAseq//Raw/cDC1'


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



# Create required directories
dir.create(file.path(path_save, '/Figures'))
dir.create(file.path(path_save, '/Saves'))

dir.create(file.path(paste0(path_save, '/Figures'), 'scRNA-seq'))
dir.create(file.path(paste0(path_save, '/Figures/scRNA-seq'), 'WT'))
dir.create(file.path(paste0(path_save, '/Figures/scRNA-seq/WT'), 'Markers'))
dir.create(file.path(paste0(path_save, '/Saves'), 'scRNA-seq'))
dir.create(file.path(paste0(path_save, '/Saves/scRNA-seq'), 'WT'))

# Set up path for figures and saves
PATH_FIG  <- paste0(path_save, '/Figures/scRNA-seq/WT')
PATH_SAVE <- paste0(path_save, '/Saves/scRNA-seq/WT')

#-------------------------------------------------------------------------------
# IMPORT DATA AND TURNING THEM INTO SEURAT OBJECT
#-------------------------------------------------------------------------------

#Run 1 -------------------------------------------------------------------------

cdc_1 <- Read10X(data.dir = paste0(path,'/Run1')) %>%
  CreateSeuratObject(project = 'cdc', min.cells = 1, min.features = 0 )

#Run 2 -------------------------------------------------------------------------

cdc_2 <- Read10X(data.dir = paste0(path,'/Run2')) %>%
  CreateSeuratObject(project = 'cdc', min.cells = 1, min.features = 0)


# Spread datasets in groups of experiment
cdc_1@meta.data[['Experiment']]     <- 'Experiment 1'
cdc_2@meta.data[['Experiment']] <- 'Experiment 2'

# Merge all the experiments to perform the quality check 
data  <- merge(cdc_1, cdc_2, 
               add.cell.ids = c('CDP',2), 
               project = 'scRNA_DC')
#-------------------------------------------------------------------------------
# TEST
#-------------------------------------------------------------------------------

# based on : https://github.com/satijalab/seurat/issues/7191#issuecomment-1636319705

#extraire les noms des gènes et les noms des cellules dans des vecteurs

data1 <- JoinLayers(data)
var <-LayerData(data1, assay="RNA", layer = "counts")
data1@assays[["RNA"]]@layers[["counts"]]@Dimnames <- dimnames(var)

#remplace matrice count par un df avec le nom des gènes et des cellules

#-------------------------------------------------------------------------------
# QUALITY CHECK 
#-------------------------------------------------------------------------------

# Calculate the percent of mitochondrial genes expression
data1$Percent.mt <- PercentageFeatureSet(data1, pattern = '^mt-')

plot0   <- VlnPlot(data1, features = c('nFeature_RNA','Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1, 3)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot0, PATH_FIG)

data1@assays[["RNA"]]@cells
cell_names <- Cells(data1)
print(cell_names)
features_names <- Features(data1)
print(features_names)

# Elimination of cells that over-expressed mitochondrial genes (dead cells)
data.QC <- subset(data1, subset = ((Experiment == 'Experiment 1' & 
                                      nFeature_RNA > 800 & 
                                      nFeature_RNA < 4500 & Percent.mt < 6) |
                                     (Experiment != 'Experiment 1' & 
                                        nFeature_RNA > 800 & 
                                        nFeature_RNA < 4500 & Percent.mt < 6))) 

plot1   <- VlnPlot(data.QC, features = c('nFeature_RNA', 'Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1, 3)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot1, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after quality check
saveRDS(data.QC, paste0(PATH_SAVE, '/1_AfterQC.rds'))
data.QC <- readRDS(paste0(PATH_SAVE, '/1_AfterQC.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#-------------------------------------------------------------------------------
# NORMALIZATION AND SCALING # STEP ON BIRD CLUSTER
#-------------------------------------------------------------------------------

data.combined <- NormalizeData(data.QC, normalization.method = 'LogNormalize')

# Cell Cycle Genes identification and scoring
s.genes       <- capitalize(tolower(cc.genes$s.genes))
g2m.genes     <- capitalize(tolower(cc.genes$g2m.genes))

# CellCycleScoring() dont work with Seurat V5 
# based on : https://github.com/satijalab/seurat/issues/7191#issuecomment-1636319705

#data.combined[["joined"]] <- JoinLayers(data.combined[["RNA"]])
#DefaultAssay(data.combined) <- "joined"

data.combined <- CellCycleScoring(data.combined, s.features = s.genes, 
                                  g2m.features = g2m.genes, set.ident = TRUE)


# Calculate difference between S and G2M is described as more relevant
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Alternate Workflow
data.combined$CC.Difference <- data.combined$S.Score - data.combined$G2M.Score
data.combined <- ScaleData(data.combined, 
                           vars.to.regress = c('CC.Difference','Percent.mt'), 
                           verbose = TRUE, features = rownames(data.combined), 
                           do.scale = TRUE)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after scaling step
saveRDS(data.combined, paste0(PATH_SAVE, '/2_Scaled.rds'))
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

data.combined <- RunUMAP(data.combined, dims = 1:18, metric = 'euclidean', 
                         n.neighbors = 50)

plot3         <- DimPlot(data.combined, group.by = 'Experiment', 
                         reduction = 'umap', cols = ColorBlind[c(1,10)], 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot3, PATH_FIG)

plot4         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(10,3,1,9)], 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot4, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after integration
saveRDS(data.combined, paste0(PATH_SAVE, '/3_Not_Integrated.rds'))
data <- readRDS(paste0(PATH_SAVE, '/3_Not_Integrated.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#-------------------------------------------------------------------------------
# CLUSTERING AND OUTLIER ELIMINATION 
#-------------------------------------------------------------------------------

data          <- FindNeighbors(data, reduction = 'umap', 
                               dims = 1:ncol(Embeddings(data, 'umap')))
data          <- FindClusters(data)

plot7         <- DimPlot(data, reduction = 'umap', 
                         group.by = 'seurat_clusters', label = TRUE, 
                         repel = TRUE, cols = c(rep(ColorBlind,10)), 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2') 
writePlot(plot7, PATH_FIG)


# Identification of cell types
clusterMarker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
top30         <- clusterMarker %>% group_by(cluster) %>% top_n(30, avg_log2FC) 
write.csv(top30, paste0(PATH_SAVE, '/WT_Top30markers.csv'), row.names = FALSE)

# For each cluster

for (cluster_id in unique(top30$cluster)) {
  cluster_data <- top30[top30$cluster == cluster_id, "gene"]
  
  # Construire le chemin proprement
  file_name <- paste0('WT_Top30markers_cluster_', cluster_id, '.txt')
  file_path <- file.path(PATH_SAVE, file_name)
  
  # Écrire le fichier
  write.table(cluster_data, file_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

#pre-cDC1

FeaturePlot(data, "Irf8", reduction = "umap")
FeaturePlot(data, "Batf3", reduction = "umap")
FeaturePlot(data, "Id2", reduction = "umap")
FeaturePlot(data, "Nfil3", reduction = "umap")
FeaturePlot(data, "Bcl6", reduction = "umap")

#pre-cDC2

FeaturePlot(data, "Irf4", reduction = "umap")
FeaturePlot(data, "Notch2", reduction = "umap")
FeaturePlot(data, "Klf4", reduction = "umap")

#pDC
FeaturePlot(data, "Tcf4", reduction = "umap")
FeaturePlot(data, "Zeb2", reduction = "umap")

## DATA FILTERING --------------------------------------------------------------
# Suppression des clusters aberrants
Inliers <- subset(data, seurat_clusters %!in% c(23,24,25))
data.combined <- Inliers

data.combined <- FindVariableFeatures(data.combined, nfeatures = 3000)
data.combined <- RunPCA(data.combined, npcs = 50)

plottest2         <- ElbowPlot(data.combined, 30) +
  labs(title = '', y = 'Standard deviation', x = 'Component number') + 
  theme(plot.title = element_text(hjust = 0.5))
writePlot(plottest2, PATH_FIG)

data.combined <- RunUMAP(data.combined, dims = 1:18, metric = 'euclidean', 
                         n.neighbors = 50)

plottest3         <- DimPlot(data.combined, group.by = 'Experiment', 
                         reduction = 'umap', cols = ColorBlind[c(1,10)], 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plottest3, PATH_FIG)



plottest4         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(10,3,1,9)], 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plottest4, PATH_FIG)


plottest4 <- DimPlot(data.combined,group.by = 'seurat_clusters',
  reduction = 'umap',cols = c(rep(ColorBlind, 10)),pt.size = 0.5,
  label = TRUE, label.size = 5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plottest4, PATH_FIG)


#pre-cDC1
FeaturePlot(data.combined, "Batf3", reduction = "umap")
FeaturePlot(data.combined, "Id2", reduction = "umap")
FeaturePlot(data.combined, "Zeb2", reduction = "umap")

VlnPlot(data.combined, features = "Id2",pt.size = 0)
VlnPlot(data.combined, features = "Zeb2", pt.size =0)
VlnPlot(data.combined, features = "Batf3", pt.size=0)

################################################################################

DimPlot(data.combined, reduction = 'umap', group.by = 'seurat_clusters', 
        label = TRUE, repel = TRUE, 
        cols = c(rep(ColorBlind,10)), pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')

################################################################################

#On garde que les cellules qui expriment Id2

data <- subset(data.combined, seurat_clusters %in% c(7,15,14))

cell_names_keep <- colnames(data)

write.table(cell_names_keep, file = "/Users/E192693Z/Desktop/cDC1_to_keep.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#On garde le nom des cellules du cluster 7 pour regarder si elles sont dans les CDP qu'on enlève les outlieres

data <- subset(data.combined, seurat_clusters %in% 7)
cellnames_cluster_7 <- colnames(data)

write.table(cellnames_cluster_7, file= "D:/Resultats_memoire/scRNAseq/Raw/cDC1/Saves/cells_cluster_7.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

















