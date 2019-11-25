############ This script performs Integration Analysis E11.5 vs E13.5 MNP cells using Seurat [v3.0.0]###############
############ Script carries out basic filtering, global scaling based normalization and scaling using Seurat Functions
############ Script regresses out Cellcycle difference between G2m and S phase using ScaleData Function
############ Script obtians integration anchors and then uses the same to integrate the datasets
############ Please Refer to Seurat [v3.0.0] manual for details on parameters and functions



#########Load R packages #########
library(Biobase)
library(Seurat)
library(cowplot)
library(plotly)
##################################


################################### Set Working Directory ###################################
setwd(WorkingDir)
#############################################################################################


################################### Read Counts Matrices and MetaFiles #####################
ctrl.data <- as.matrix(read.table(file = 'E11.5 Counts Matrix', sep = "\t", header=T, row.names=1, check.names=F))      #### Please provide E11.5MNP counts
stim.data <- as.matrix(read.table(file = "E13.5 Counts Matrix", sep = "\t", header=T, row.names=1, check.names=F))      #### Please provide E13.5MNP counts


meta <- as.matrix(read.table(file = "E11.5_MNP MetaData", sep = "\t", header=T, row.names=1, check.names=F))            #### Please provide E11.5MNP metafile
col_d <- data.frame(meta)
metactrl <- new("AnnotatedDataFrame", data = col_d)
rownames(metactrl) <- metactrl$Cells

meta <- as.matrix(read.table(file = "E13.5_MNP_MetaData.txt", sep = "\t", header=T, row.names=1, check.names=F))        #### Please provide E13.5MNP metafile
col_d <- data.frame(meta)
metastim <- new("AnnotatedDataFrame", data = col_d)
rownames(metastim) <- metastim$Cells
#############################################################################################


################# Setuo Seurat Object, Normalize, Find Variable Features and Integration #############
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "E11.5MNP")

CellsMeta = ctrl@meta.data
CellsMeta["Cluster"] <- metactrl$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
ctrl <- AddMetaData(ctrl, CellsMetaTrim)

CellsMeta = ctrl@meta.data
CellsMeta["Stage"] <- metactrl$Stage
CellsMetaTrim <- subset(CellsMeta, select = c("Stage"))
ctrl <- AddMetaData(ctrl, CellsMetaTrim)

ctrl$stim <- "E11.5MNP"
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)


stim <- CreateSeuratObject(counts = stim.data, project = "E13.5MNP")

CellsMeta = stim@meta.data
CellsMeta["Cluster"] <- metastim$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
stim <- AddMetaData(stim, CellsMetaTrim)

CellsMeta = stim@meta.data
CellsMeta["Stage"] <- metastim$Stage
CellsMetaTrim <- subset(CellsMeta, select = c("Stage"))
stim <- AddMetaData(stim, CellsMetaTrim)

stim$stim <- "E13.5MNP"
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

mnp.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
mnp.combined <- IntegrateData(anchorset = mnp.anchors, dims = 1:20)

DefaultAssay(mnp.combined) <- "integrated"
###########################################################################################################


################# Cell-Cycle Regression #############
cc.genes <- readLines(con = CellCycleGenesFile)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
mnp.combined <- CellCycleScoring(object = mnp.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mnp.combined@meta.data$CC.Difference <- mnp.combined@meta.data$S.Score - mnp.combined@meta.data$G2M.Score
#####################################################


############################ Perform PCA, UMAP/tSNE, Clustering and Visualization ###########################
mnp.combined_aftercc <- ScaleData(mnp.combined, vars.to.regress='CC.Difference', verbose = TRUE)

mnp.combined_aftercc <- RunPCA(mnp.combined_aftercc, npcs = 30, verbose = TRUE) # tSNE and Clustering
mnp.combined_aftercc <- RunUMAP(mnp.combined_aftercc, reduction = "pca", dims = 1:20)
mnp.combined_aftercc <- RunTSNE(mnp.combined_aftercc, reduction = "pca", dims = 1:20)
mnp.combined_aftercc <- FindNeighbors(mnp.combined_aftercc, reduction = "pca", dims = 1:20)
mnp.combined_aftercc <- FindClusters(mnp.combined_aftercc)

# Visualization using UMAP
p1 <- DimPlot(mnp.combined_aftercc, reduction = "umap", group.by = "stim")
p2 <- DimPlot(mnp.combined_aftercc, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
ggsave('UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval.png', width=15, height=15)

DimPlot(mnp.combined_aftercc, reduction = "umap", split.by = "stim")
ggsave('UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_Comparison.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_Comparison.html')

DimPlot(mnp.combined_aftercc, reduction = "umap", group.by="Cluster")
ggsave('UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithOriginalClusters.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithOriginalClusters.html')

DimPlot(mnp.combined_aftercc, reduction = "umap", group.by="Phase")
ggsave('UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_CellCyclePhase.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_CellCyclePhase.html')

DimPlot(mnp.combined_aftercc, reduction = "umap", group.by="Stage")
ggsave('UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithStage.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'UMAP_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithStage.html')

# Visualization using tSNE
p1 <- DimPlot(mnp.combined_aftercc, reduction = "tsne", group.by = "stim")
p2 <- DimPlot(mnp.combined_aftercc, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)
ggsave('tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval.png', width=15, height=15)

DimPlot(mnp.combined_aftercc, reduction = "tsne", split.by = "stim")
ggsave('tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_Comparison.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_Comparison.html')

DimPlot(mnp.combined_aftercc, reduction = "tsne", group.by="Cluster")
ggsave('tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithOriginalClusters.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithOriginalClusters.html')

DimPlot(mnp.combined_aftercc, reduction = "tsne", group.by="Phase")
ggsave('tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_CellCyclePhase.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_CellCyclePhase.html')

DimPlot(mnp.combined_aftercc, reduction = "tsne", group.by="Stage")
ggsave('tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithStage.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'tSNE_IntegrationPlot_E11.5_vs_E13.5_MNP_afterCCremoval_WithStage.html')

Idents(mnp.combined_aftercc) <- mnp.combined_aftercc@meta.data$seurat_clusters

### Identify conserved cell type markers
DefaultAssay(mnp.combined_aftercc) <- "RNA"
#######################################################################################################################


################################## Marker Analysis and Saving R object ####################################################################
Markers <- FindAllMarkers(mnp.combined_aftercc)
write.table(Markers,"Markers_All_IntegrationClusters_afterCCremoval.txt",sep="\t",quote=F)

save(mnp.combined_aftercc,file = 'E11.5_vs_E13.5_MNP_IntegrationAnalysis_afterCCremoval.Robj')
###########################################################################################################################################