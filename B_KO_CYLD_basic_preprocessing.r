#This is the script to analyze Mbe1Cre-CYLDflx/flx scRNA-seq data from B cells from bone marrow


library(Seurat)
library(tidyseurat)
library(SingleR)
library("bluster")
library(Nebulosa)
setwd("./CYLD_project")


#1. Load Seurat object and perform basic QC analysis and downsample
counts = read.delim("BM_KO_CYLD.gz", h=T, row.names=1)
B_cyld <- CreateSeuratObject(counts=counts, project="BM_ΚΟ", min.cells=5, min.features=200) 
sample_name = "BM_ΚΟ"
B_cyld

#Getting the mitochondrial percentage in our dataset because dead cells can half an increased mitochondrial percentage
B_cyld[["percent.mt"]] <- PercentageFeatureSet(B_cyld, pattern = "^mt-")
VlnPlot(B_cyld, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
head(B_cyld@meta.data, 5)
plot1 <- FeatureScatter(B_cyld, feature1="nCount_RNA", feature2="percent.mt", cols=NULL)
plot2 <- FeatureScatter(B_cyld, feature1="nCount_RNA", feature2="nFeature_RNA")
plot1 + plot2
summary(B_cyld$percent.mt)
summary(B_cyld$nFeature_RNA)

#Subsetting to remove cells with unwanted characteristics (possible doublets or empty droplets, dead cells)
B_cyld <- subset(B_cyld, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
dim(B_cyld@assays$RNA@data)

cc.genes <- readLines(con="mmusculus_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
B_cyld <- CellCycleScoring(B_cyld, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
head(B_cyld[[]])
B_cyld

#Basic-basic pre-processing
B_cyld <- NormalizeData(B_cyld)
B_cyld <- FindVariableFeatures(B_cyld)
B_cyld

#Downsample to 700 cells
B_atoms <- LeverageScoreSampling(object = B_cyld, num.cells = 700)

#Use SCTransform as an improved approach for Normalization, Scaling, detection of Variable Features and removal of confounding sources of variation (mitoch., cell cycle)
B_atoms <- CellCycleScoring(B_atoms, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
head(B_atoms[[]])
B_atoms <- B_atoms %>% 
SCTransform(vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), vst.flavor = "v2")

#Pushing the PCs to a high number based on SCTransform specifications
B_atoms <- RunPCA(B_atoms, npcs = 30, features=VariableFeatures(B_atoms))
DimPlot(B_atoms, reduction = "pca")

B_atoms <- RunUMAP(B_atoms, dims=1:30, reduction = "pca", verbose = FALSE)
ElbowPlot(B_atoms)
nPCs = 30 
# res = 1, this is the initial value of this script
res = 0.5
B_atoms <- B_atoms %>% 
FindNeighbors(dims=1:nPCs, reduction = "pca") %>% 
FindClusters(resolution=res, algorithm = 1)
DimPlot(B_atoms, label=T)


#2. Perform cell annotation with SingleR
blueprint <- celldex::ImmGenData()
cell_type_df <-
  GetAssayData(B_atoms, slot = 'counts', assay = "SCT") %>%
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
B_atoms <-
  B_atoms %>%
  left_join(cell_type_df, by = "cell")
# Reorder columns
B_atoms %>%
  tidyseurat::select(cell, labels, everything())
# B_atoms %>%
#   ggplot(aes(UMAP_1, UMAP_2, UMAP_3, color = first.labels)) +
#   geom_point()
# B_atoms %>%
#   ggplot(aes(UMAP_1, UMAP_2, UMAP_3, color = seurat_clusters)) +
#   geom_point()
sce <- as.SingleCellExperiment(DietSeurat(B_atoms2))
sce
Imm.ref <- celldex::ImmGenData()
Imm.fine <- SingleR(test = sce,assay.type.test = 1,ref = Imm.ref,labels = Imm.ref$label.fine)
plotScoreHeatmap(Imm.fine)
# http://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html
#heatmap of the per-cell and label scores
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
                "Macrophages (MF.II+480LO)", "Macrophages (MF.RP)", "Stem cells (SC.CDP)", 
                "DC (DC.PDC.8-)", "Monocytes (MO.6C-II+)", "DC (DC.8-)", "Stem cells (SC.MDP)")

b_srt <- B_atoms %>% filter(!(first.labels %in% nonBclones))
DimPlot(b_srt, label=T, group.by="first.labels")
# Basic preprocessing, again...
b_srt <- CellCycleScoring(b_srt, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
b_srt <- SCTransform(b_srt, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), vst.flavor = "v2", verbose = FALSE)
b_srt <- RunPCA(b_srt, features=VariableFeatures(b_srt))
b_srt <- RunUMAP(b_srt, dims=1:30, reduction = "pca", verbose = FALSE, n.epochs = 300 , min.dist = 0.1) #, n.neighbors = 20)
nPCs = 30
#res=0.8 ##that was the initial value of this parameter for this script
res = 0.5
b_srt <- FindNeighbors(b_srt, dims=1:nPCs)
b_srt <- FindClusters(b_srt, resolution=res)
b_srt@active.ident
DimPlot(b_srt, group.by = "seurat_clusters", label = TRUE)
DimPlot(b_srt, group.by = "first.labels", label = TRUE)
DimPlot(b_srt)
b_srt@active.ident

#Change idents-manual clustering
DimPlot(b_srt, label = T, reduction = "umap", group.by = "labels")
plot = DimPlot(b_so, label = T, reduction = "umap")
plot
b_srt = CellSelector(plot = plot, object = b_srt) #manual clustering
levels(b_srt)
DimPlot(b_srt, label = T, reduction = "umap")
new.cluster.ids <- c("Mature B", "Pre B", "Immature B", "Pro B", "Pro/Pre B")
names(new.cluster.ids) <- levels(b_srt)
b_so <- RenameIdents(b_so, new.cluster.ids)
DimPlot(b_so, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(b_so, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "labels")
b_srt$Idents <- Idents(b_srt)
p1 <- DimPlot(b_srt, group.by= "seurat_clusters") + ggtitle("Seurat clusters")
p2 <- DimPlot(b_srt, group.by= "labels") + ggtitle("SingleR granular annotation")
p3 <- Dimplot(b_srt, group.by = "Idents", label=T) + ggtitle("Non-granular annotation")
p4 <- FeaturePlot(b_srt, features = "Cyld")
saveRDS(b_srt, "B_CYLDKO.rds")

#4. Sanity check for cell identities (https://doi.org/10.1016/j.it.2022.01.003)
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


#5. Plotting for Figure S6
plot <- ((p1 + p2) / (p3 + p4)) + patchwork::plot_annotation(tag_levels = "A")
ggsave("figureS5.png", plot = plot, device = "png", width = 10, height = 14)
