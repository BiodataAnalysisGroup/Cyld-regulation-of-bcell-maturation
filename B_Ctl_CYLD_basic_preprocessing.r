#This is the script to analyze CYLDflx/flx scRNA-seq data from B cells from bone marrow

library(Seurat)
library(tidyseurat)
library(SingleR)
library("bluster")
library(Nebulosa)
setwd("./CYLD_project")


#1. Load Seurat object and perform basic QC analysis
counts = read.delim("Matrix_files/BM_Ctl_CYLD_project.gz", h=T, row.names=1) # also in 10.5281/zenodo.10074123
B_Ctl <- CreateSeuratObject(counts=counts, project="BM_Ctl", min.cells=5, min.features=200)
sample_name = "BM_Ctl"

#Getting the mitochondrial percentage in our dataset
B_Ctl[["percent.mt"]] <- PercentageFeatureSet(B_Ctl, pattern = "^mt-")
VlnPlot(B_Ctl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
head(B_Ctl@meta.data, 5)
plot1 <- FeatureScatter(B_Ctl, feature1="nCount_RNA", feature2="percent.mt", cols=NULL)
plot2 <- FeatureScatter(B_Ctl, feature1="nCount_RNA", feature2="nFeature_RNA")
plot1 + plot2
summary(B_Ctl$percent.mt)
summary(B_Ctl$nFeature_RNA)

#Subsetting to remove cells with unwanted characteristics (possible doublets or empty droplets, dead cells)
B_Ctl <- subset(B_Ctl, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
dim(B_Ctl@assays$RNA@data)

cc.genes <- readLines(con="mmusculus_cell_cycle_genes.txt") # 10.5281/zenodo.10074123
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
B_Ctl <- CellCycleScoring(B_Ctl, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
head(B_Ctl[[]])

#Use SCTransform as an improved approach for Normalization, Scaling, detection of Variable Features and removal of confounding sources of variation (mitoch., cell cycle)
B_Ctl <- B_Ctl %>% 
SCTransform(vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), vst.flavor = "v2")

#Pushing the PCs to a high number based on SCTransform specifications
B_Ctl <- RunPCA(B_Ctl, npcs = 30, features=VariableFeatures(B_Ctl))
DimPlot(B_Ctl, reduction = "pca")
B_Ctl <- RunUMAP(B_Ctl, dims=1:30, reduction = "pca", verbose = FALSE)
ElbowPlot(B_Ctl)
nPCs = 30 
res = 0.4 
B_Ctl <- B_Ctl %>% 
FindNeighbors(dims=1:nPCs, reduction = "pca") %>% 
FindClusters(resolution=res, algorithm = 1)
DimPlot(B_Ctl, label=T)



#2. Perform cell annotation with SingleR
blueprint <- celldex::ImmGenData()
cell_type_df <-
  GetAssayData(B_Ctl, slot = 'counts', assay = "SCT") %>%
  log1p() %>%
  Matrix::Matrix(sparse = TRUE) %>%
  SingleR::SingleR(
    ref = blueprint,
    labels = blueprint$label.fine, ##We want the exhaustive list, thus $label.fine instead of $label.main
    method = "single"
  ) %>%
  as.data.frame() %>%
  as_tibble(rownames = "cell") %>%
  dplyr::select(cell, labels)
# Join UMAP and cell type info
B_Ctl <-
  B_Ctl %>%
  left_join(cell_type_df, by = "cell")
# Reorder columns
B_Ctl %>%
  select(cell, labels, everything())
B_Ctl %>%
  ggplot(aes(UMAP_1, UMAP_2, UMAP_3, color = labels)) +
  geom_point()
B_Ctl %>%
  ggplot(aes(UMAP_1, UMAP_2, UMAP_3, color = seurat_clusters)) +
  geom_point()
sce <- as.SingleCellExperiment(DietSeurat(B_Ctl))
sce
Imm.ref <- celldex::ImmGenData()
Imm.fine <- SingleR(test = sce,assay.type.test = 1,ref = Imm.ref,labels = Imm.ref$label.fine)
plotScoreHeatmap(Imm.fine)
# http://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html
#Comparing to unsupervised clustering

pairwiseRand(sce$seurat_clusters, Imm.fine$labels, mode="index")
tab <- table(cluster=sce$seurat_clusters, label=Imm.fine$labels) 
pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing.



#3. Subset only to B cells
## Remove Non-B cells
nonBclones <- c("Neutrophils (GN)", "Neutrophils (GN.ARTH)", "Neutrophils (GN.URAC)", 
                "DC (DC.PDC.8+)", "Macrophages (MFIO5.II+480INT)", 
                "Macrophages (MFIO5.II+480LO)", 
                "T cells (T.8EFF.OT1.48HR.LISOVA)", 
                "Stem cells (GMP)", "Stem cells (MLP)", "Stem cells (SC.MEP)", 
                "Macrophages (MF.II+480LO)", "Macrophages (MF.RP)", "Stem cells (SC.CDP)", "DC (DC.PDC.8-)", "Monocytes (MO.6C-II+)")

b_srt <- B_Ctl %>% filter(!(labels %in% nonBclones))
# Basic preprocessing, again...
b_srt <- CellCycleScoring(b_srt, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
b_srt <- SCTransform(b_srt, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), vst.flavor = "v2", verbose = FALSE)
b_srt <- RunPCA(b_srt, features=VariableFeatures(b_srt))
b_srt <- RunUMAP(b_srt, dims=1:30, reduction = "pca", verbose = FALSE, n.epochs = 200, min.dist = 0.1)
nPCs = 30
res = 0.5
b_srt <- FindNeighbors(b_srt, dims=1:nPCs)
b_srt <- FindClusters(b_srt, resolution=res)
b_srt <- RunUMAP(b_srt, dims=1:30, reduction = "pca", verbose = FALSE, n.epochs = 200) #, n.epochs = 200, min.dist = 0.1, n.neighbors = 20)
DimPlot(b_srt, group.by = "seurat_clusters", label = TRUE)
DimPlot(b_srt, group.by = "labels", label = TRUE)

#Change final idents
DimPlot(b_srt, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(b_srt, reduction = "umap", label = TRUE, group.by = "labels")
levels(b_srt)
# new.cluster.ids <- c("Mature B", "Pre B", "Immature B", "Pro B") ## Align 4 Seurat clusters with these B cell groupings
# names(new.cluster.ids) <- levels(b_srt)
# b_srt <- RenameIdents(b_srt, new.cluster.ids)
b_srt$Idents <- Idents(b_srt)
p1 <- DimPlot(b_srt, group.by= "seurat_clusters") + ggtitle("Seurat clusters")
p2 <- DimPlot(b_srt, group.by= "labels") + ggtitle("SingleR granular annotation")
p3 <- Dimplot(b_srt, group.by = "Idents", label=T) + ggtitle("Non-granular annotation")
p4 <- FeaturePlot(b_srt, features = "Cyld")
saveRDS(b_srt, "B_control.rds")

#4. Sanity check for B cell identities (https://doi.org/10.1016/j.it.2022.01.003)
#Pre pro-B cell
Nebulosa::plot_density(b_srt, features = c("Irf8", "Tcf4", "Bst2"), joint = TRUE, pal = "viridis") + ggtitle('Pre pro-B cell')
# Early-Pro B cell
Nebulosa::plot_density(b_srt, features = c("Ebf1", "Bok", "Ifitm2", "Ifitm3"), joint = TRUE, pal = "viridis") + ggtitle('Early-Pro B cell')
#Late-Pro B cell
Nebulosa::plot_density(b_srt, features = c("Vpreb1", "Igll1"), joint = TRUE, pal = "viridis") + ggtitle('Late-Pro B cell')
#Pre B cell
Nebulosa::plot_density(b_srt, features = c("Vpreb1", "Igll1", "Nrgn", "Pax5"), joint = TRUE, pal = "viridis") + ggtitle('Pre B cell')
#Pre B cell [2]--Large pre-B cells
Nebulosa::plot_density(b_srt, features = c("Bub1", "Dlgap5", "Neil3", "Ncapg", "Med10"), joint = TRUE, pal = "viridis") + ggtitle('Large pre-B cells')
#Pre B cell [3]--Small pre-B cells
Nebulosa::plot_density(b_srt, features = c("Ckap2", "Cep55", "Cdca5", "Ska1", "Pax5"), joint = TRUE, pal = "viridis") + ggtitle('Small pre-B cells')
#Mature B cells
Nebulosa::plot_density(b_srt, features = c("Ms4a1", "Bank1", "Ltb", "Ctsh", "Ly6a", "Fcer2a", "Bcl2", "Fcrl1", "H2-DMa"), joint = TRUE, pal = "viridis") + ggtitle('Mature B cells')


#5. Plotting for Figure S5
plot <- ((p1 + p2) / (p3 + p4)) + patchwork::plot_annotation(tag_levels = "A")
ggsave("figureS5.png", plot = plot, device = "png", width = 10, height = 14)
