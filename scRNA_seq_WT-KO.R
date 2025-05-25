#-------------------------------------------------------------------------------
# Preparation de l'environnement
#-------------------------------------------------------------------------------

path <- 'D:/Resultats_memoire/scRNAseq'
path_save <- 'D:/Resultats_memoire/scRNAseq/WT-KO'

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

CellTypeOrder <- c('ALP', 'TULIP', 'ILCpro', 'overLIP', 'TOX-KO', 'CDP')

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

# TOX-KO
TOX_KO <- Read10X(data.dir = paste0(PATH, '/2019_scRNAseq_data/TOX_KO')) %>%
  CreateSeuratObject(project = 'TOX_KO', min.cells = 1, min.features = 0)

# Spread datasets in groups of experiment
ALP@meta.data[['Experiment']]     <- 'Experiment 1'
ILCpro@meta.data[['Experiment']]  <- 'Experiment 1'
TULIP@meta.data[['Experiment']]   <- 'Experiment 1'
overLIP@meta.data[['Experiment']] <- 'Experiment 2'
TOX_KO@meta.data[['Experiment']]  <- 'Experiment 2'


# Merge all the experiments to perform the quality check 
data  <- merge(TULIP, c(ILCpro, ALP, overLIP, TOX_KO), 
               add.cell.ids = c('TULIP', 'ILCpro', 'ALP', 'overLIP', 'TOX_KO'), 
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
                                     (Experiment == 'Experiment 2' & 
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
                         reduction = 'umap', cols = ColorBlind[c(6,1,3)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot3, PATH_FIG)

FeaturePlot(data.combined,feature="Cd74", reduction="umap")

FeaturePlot(data.combined,feature="Id2", reduction="umap")


plot4         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(8,10,1,9,3,7)], 
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
                         cols = ColorBlind[c(6,1,3)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') 
writePlot(plot5, PATH_FIG)

plot6         <- DimPlot(data, group.by = 'orig.ident', reduction = 'umap_cssz', 
                         cols = ColorBlind[c(8,10,1,9,3,7)], pt.size = 0.5, 
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


plot8 <- FeaturePlot(data, feature = "Id2", reduction = 'umap_cssz')
writePlot(plot8, PATH_FIG)

plot9<-FeaturePlot(data, feature = "Zeb2", reduction = 'umap_cssz')
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

## DATA FILTERING --------------------------------------------------------------

# Remove mature cells (23,24), added B cells (20) and contaminants (25,28,29,5)
Inliers   <- subset(data, seurat_clusters %!in% c(23,24,20,25,28,29,5))
data  <- Inliers

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
                          cols = ColorBlind[c(8,10,1,9,3,7)], pt.size = 0.5, 
                          order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') 
writePlot(plot12, PATH_FIG)

FeaturePlot(data,feature="Id2", reduction="umap_cssz2")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
saveRDS(data, paste0(PATH_SAVE, '/4_Filtered.rds'))
data <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

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
                 reducedDim = toupper('umap_cssz2'), start.clus = 2, end.clus = c(3,6))


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

