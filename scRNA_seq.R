#-------------------------------------------------------------------------------
# Preparation de l'environnement
#-------------------------------------------------------------------------------

path <- 'D:/Resultats_memoire/scRNAseq/Raw/ILC/scRNA-seq_Raw'
path_save <- 'D:/Resultats_memoire/scRNAseq/Raw/ILC'

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
suppressPackageStartupMessages(library(Nebulosa))


CellTypeOrder <- c('ALP', 'overLIP', 'TULIP', 'ILCpro')

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

#Experiment 1 ------------------------------------------------------------------

# All Lymphoid Progenitor
ALP <- Read10X(data.dir = paste0(path,'/1_ABSC_CLP')) %>%
  CreateSeuratObject(project = 'ALP', min.cells = 1, min.features = 0 )

#sEILP, cEILP and ILC progenitors
ILCpro <- Read10X(data.dir = paste0(path,'/3_ABSC_EILP_WT')) %>%
  CreateSeuratObject(project = 'ILCpro', min.cells = 1, min.features = 0)

# a4b7 intermediates
TULIP <- Read10X(data.dir = paste0(path,'/2_ABSC_TULIP')) %>%
  CreateSeuratObject(project = 'TULIP', min.cells = 1, min.features = 0)


#Experiment 2 ------------------------------------------------------------------

# Il7rLT+ XT
overLIP <-  Read10X(data.dir = paste0(path,'/2019_scRNAseq_data/WT')) %>%
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
data@meta.data[['Project']] <- 'scRNA_WT_Integration'

# Organize cell type in the development order for the plots
data@active.ident <- factor(data@active.ident, levels = CellTypeOrder)

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
                                     nFeature_RNA > 200 & 
                                     nFeature_RNA < 5000 & Percent.mt < 2) |
                                    (Experiment != 'Experiment 1' & 
                                       nFeature_RNA > 1200 & 
                                       nFeature_RNA < 5000 & Percent.mt < 4))) 

plot1   <- VlnPlot(data.QC, features = c('nFeature_RNA', 'Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1, 3)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot1, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after quality check
saveRDS(data.QC, paste0(PATH_SAVE, '/1_AfterQC.rds'))
#data.QC <- readRDS(paste0(PATH, '/1_AfterQC.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#-------------------------------------------------------------------------------
# NORMALIZATION AND SCALING
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
#data.combined <- readRDS(paste0(PATH_SAVE, '/2_Scaled.rds'))
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
                         reduction = 'umap', cols = ColorBlind[c(1,10)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot3, PATH_FIG)

plot4         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(10,3,1,9)], 
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
                         cols = ColorBlind[c(1,10)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') + scale_x_reverse()
writePlot(plot5, PATH_FIG)

plot6         <- DimPlot(data, group.by = 'orig.ident', reduction = 'umap_cssz', 
                         cols = ColorBlind[c(10,3,1,9)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') + scale_x_reverse()
writePlot(plot6, PATH_FIG)

# Check cell cycle distribution
FeaturePlot(data, 'S.Score', reduction = 'umap_cssz') + scale_x_reverse()
FeaturePlot(data, 'G2M.Score', reduction = 'umap_cssz') + scale_x_reverse()

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
data          <- FindClusters(data)

plot7         <- DimPlot(data, reduction = 'umap_cssz', 
                         group.by = 'seurat_clusters', label = TRUE, 
                         repel = TRUE, cols = c(rep(ColorBlind,10)), 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2') + scale_x_reverse()
writePlot(plot7, PATH_FIG)


# Identification of cell types
clusterMarker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
 top30         <- clusterMarker %>% group_by(cluster) %>% top_n(30, avg_log2FC) 
write.csv(top30, paste0(PATH_SAVE, '/WT_Top30markers.csv'), row.names = FALSE)

# For each cluster
for(cluster_id in unique(top30$cluster)) {
  cluster_data <- top30[top30$cluster == cluster_id, "gene" ]
  file_name <- paste0(path, '/WT_Top30markers_cluster_', cluster_id, '.txt')
  write.table(cluster_data, paste0(PATH_SAVE, file_name), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="")}

FeaturePlot(data, "Icos", reduction = "umap_cssz")
FeaturePlot(data, "Bcl11b", reduction = "umap_cssz")
FeaturePlot(data, "Lgals3", reduction = "umap_cssz")
FeaturePlot(data, "Lyz2", reduction = "umap_cssz")

## KEY GENES REPRESENTATION ----------------------------------------------------

markers_genes <- c(
  'Nfil3',
  # B Cells
  'Pax5', 'Ebf1', 
  # Mast Cells
  'Cpa3', 'Gzmb', 'Kit', 
  # Mono/Macrophages
  'Lgals3', 'Lyz2', 
  # ILC2
  'Il12ra', 'Icos', 'Stab2', 'Bcl11b', 'Tox2'
)

for(g in markers_genes){
  if(g %in% rownames(data)){
    density <- plot_density(data, g, pal = 'magma', 
                            size = 2) + scale_x_reverse()
    density
    
    feature <- FeaturePlot(data, g, order = T, pt.size = 1, 
                           reduction = 'umap_cssz') + scale_x_reverse()
    feature
  }
}


## CONTAMINANT PERCENTAGE ------------------------------------------------------

# Overall outliers
Outliers  <- subset(data, seurat_clusters %in% c(15,16,17,21,22,24))
Inliers   <- subset(data, seurat_clusters %!in% c(15,16,17,21,22,24))

prop_out  <- table(factor(Outliers$orig.ident, levels = CellTypeOrder))
prop_in   <- table(factor(Inliers$orig.ident, levels = CellTypeOrder))
prop_all  <- (prop_out/(prop_in+prop_out))*100

plot8     <- barplot(prop_all[c(2,3)], col = ColorBlind[c(3,1)], 
                     border = 'white', ylab = 'Outlier percentage', 
                     main = 'Percentage of outliers' , ylim = c(0, 60))
y         <- prop_all[c(2,3)]
text(plot8, y+2, labels = paste(format(y, scientific=FALSE, digits = 2 ), '%'))
plot8     <- recordPlot()
writePlot(plot8, PATH_FIG)


# Outliers by cellype
prop_table           <- as.data.frame(matrix(0,4,4))
colnames(prop_table) <- c('ILC2', 'B_cells', 'Mast_Baso', 'Macro_Mono_Granulo')
rownames(prop_table) <- names(table(Outliers$orig.ident))

for(i in 1:length(Outliers$seurat_clusters)){
  if(Outliers$seurat_clusters[i] %in% c(15,21)){
    prop_table[Outliers$orig.ident[i],'ILC2'] <- prop_table[
      Outliers$orig.ident[i],'ILC2']+1
  }else if(Outliers$seurat_clusters[i] %in% 16){
    prop_table[Outliers$orig.ident[i],'B_cells'] <- prop_table[
      Outliers$orig.ident[i],'B_cells']+1
  }else if(Outliers$seurat_clusters[i] %in% 17){
    prop_table[Outliers$orig.ident[i],'Mast_Baso'] <- prop_table[
      Outliers$orig.ident[i],'Mast_Baso']+1
  }else if(Outliers$seurat_clusters[i] %in% c(22,24)){
    prop_table[Outliers$orig.ident[i],'Macro_Mono_Granulo'] <- prop_table[
      Outliers$orig.ident[i],'Macro_Mono_Granulo']+1
  }
}

prop_table <- rbind(prop_table, Total = colSums(prop_table))
perc_table <- prop_table[1:4,]/prop_table[rep(5,4),]*100
bar_table  <- data.frame(Percentage = c(perc_table[,1],perc_table[,2], 
                                        perc_table[,3], perc_table[,4]),
                         CellType = rep(colnames(perc_table), each = 4),
                         Sample = rep(rownames(perc_table)))

plot9 <- ggplot(bar_table[bar_table$Sample %in% c('TULIP', 'overLIP'),], 
                aes(fill=Sample, y=Percentage, x=CellType)) + 
  geom_bar(position='stack', stat='identity') + 
  ggtitle('Proportions of celltypes in outleirs') +
  scale_fill_manual(values= ColorBlind[c(3,1)]) + 
  theme_classic()
writePlot(plot9, PATH_FIG)


## DATA FILTERING --------------------------------------------------------------

# Remove mature cells (15,21), added B cells (16) and contaminants (17,22,24)
data  <- Inliers

data  <- RunUMAP(data, reduction = 'cssz', 
                 dims = 1:ncol(Embeddings(data, 'cssz')), n.neighbors = 50, 
                 reduction.name='umap_cssz2', reduction.key='UMAPCSSZ2_')
data  <- FindNeighbors(data, reduction = 'umap_cssz2', 
                       dims = 1:ncol(Embeddings(data, 'umap_cssz2')))
data  <- FindClusters(data, resolution = 0.1)


# Create inverted embedding to simplify representation
embed <- as.matrix(data.frame(
  -data@reductions[['umap_cssz2']]@cell.embeddings[,1],
  data@reductions[['umap_cssz2']]@cell.embeddings[,2]))

colnames(embed)           <- c('REVUMAPCSSZ2_1', 'REVUMAPCSSZ2_2')
data[['rev_umap_cssz2']]  <- CreateDimReducObject(embeddings = embed, 
                                                  key = 'REVUMAPCSSZ2_', 
                                                  assay = DefaultAssay(data))
data@reductions[['rev_umap_cssz2']]@global <- T


plot10 <- DimPlot(data, group.by = 'orig.ident', reduction = 'rev_umap_cssz2', 
                  cols = ColorBlind[c(10,3,1,9)], pt.size = 1.5, 
                  order = rev(CellTypeOrder)) + 
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot10, PATH_FIG)

plot11 <- DimPlot(data, reduction = 'rev_umap_cssz2', 
                  group.by = 'seurat_clusters', 
                  label = TRUE, repel = TRUE, 
                  cols = c(ColorBlind[c(3,9,2,7,8,10)], "purple"), 
                  pt.size = 1) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot11, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
saveRDS(data, paste0(PATH_SAVE, '/4_Filtered.rds'))
data <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

FeaturePlot(data, feature = "Flt3", reduction = 'rev_umap_cssz2', order = T, pt.size = 0.5)
FeaturePlot(data, feature = "Tcf7", reduction = 'rev_umap_cssz2', order = T, pt.size = 0.5)
FeaturePlot(data, feature = "Cd74", reduction = 'rev_umap_cssz2', order = T, pt.size = 0.5)


plot_density(data, features = "Flt3", pal='magma', size=2)
plot_density(data, features = "Tcf7", pal='magma', size=2)
plot_density(data, features = "Cd74", pal='magma', size=2)

# Set reductio nto use and plot axis limits
red    <- 'rev_umap_cssz2'
xlimit <- c(min(data@reductions[[red]]@cell.embeddings[,1]),
            max(data@reductions[[red]]@cell.embeddings[,1]))
ylimit <- c(min(data@reductions[[red]]@cell.embeddings[,2]),
            max(data@reductions[[red]]@cell.embeddings[,2]))      

#===============================================================================
## SLINGSHOT -------------------------------------------------------------------
#===============================================================================
#Based on : https://bustools.github.io/BUS_notebooks_R/slingshot.html

DimPlot(data, reduction = red, group.by = 'seurat_clusters', 
        label = TRUE, repel = TRUE, 
        cols = c(ColorBlind[c(3,9,2,7,8,10)], "purple"), pt.size = 0.9) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')


## PSEUDOTIME RECONSTRUCTION ---------------------------------------------------
#https://github.com/satijalab/seurat/issues/8248
library(Seurat)

sce                                     <- as.SingleCellExperiment(data)     
colData(sce)$Seurat_clusters            <- as.character(data@active.ident)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', 
                 reducedDim = toupper(red), start.clus = 2, end.clus = c(1,5))


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

# Specifying gradient colors for plots
gradient         <- hcl.colors(100)

# Plot with curves
curves <- slingCurves(sce, as.df = TRUE) %>%
  dplyr::rename("x" = paste0(str_replace_all(
    toupper(red), "[[:punct:]]", ""), '_1'), 
    "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

# First Trajectory
plot12 <- FeaturePlot(data, 'Slingshot_1', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            size = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot12, PATH_FIG)

# Second Trajectory
plot13 <- FeaturePlot(data, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            size = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot13, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
saveRDS(data, paste0(PATH_SAVE, '/5_Pseudotime.rds'))
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

ILC  <- subset(data, Slingshot_1 %!in% NA)
DC <- subset(data, Slingshot_2 %!in% NA)


#-------------------------------------------------------------------------------
# TradeSeq
#-------------------------------------------------------------------------------
# based on : https://statomics.github.io/tradeSeq/articles/tradeSeq.html

palette(brewer.pal(8, "Dark2"))

# Counts
data <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
counts <- as.matrix(data@assays[["RNA"]]@layers[["counts"]])

# Establish variable gene list
#suppressPackageStartupMessages(library(Hmisc))

#cc_genes <- capitalize(tolower(as.vector(unlist(cc.genes))))
#var.features <- rownames(data)
#var.features <- var.features[var.features %!in% cc_genes] 

# Subsetting counts matrix
#counts <- counts[rownames(counts) %in% var.features,]

# Curves
sce <- readRDS(paste0(PATH_SAVE, '/Slingshot.rds'))
crv <- as.SlingshotDataSet(sce)

# Cell types 
celltype <-data@meta.data[["orig.ident"]]

#-------------------------------------------------------------------------------
# Fit negative binomial model
#-------------------------------------------------------------------------------

#A faire sur le cluster bird 
set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
                   nGenes = 200, verbose = T)

#-------------------------------------------------------------------------------

set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)


sessionInfo()

library(devtools)

#-------------------------------------------------------------------------------
# NUMBER OF CELLS SPECIFIC FROM EACH TRAJECTORY
#-------------------------------------------------------------------------------

DC  <- subset(data, Slingshot_2 %!in% NA)
ILC <- subset(data, Slingshot_1 %!in% NA)

DC_cells <- Cells(DC)
ILC_cells <- Cells(ILC)

commun <- intersect (DC_cells, ILC_cells)


