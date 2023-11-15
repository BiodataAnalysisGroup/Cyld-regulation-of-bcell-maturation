#This is the script for pre- and post- scFEA analysis through Seurat | Figures S9, S10, S11
# https://doi.org/10.1101/gr.271205.120

library(Seurat)
library(tidyseurat)
library(patchwork)
library(ggplot2)
setwd("./CYLD_project")

#1. Load Seurat objects for Control and Cyldko
#Control (Cyldflx/flx)
b_ctl <- readRDS("Seurat_objects/B_control.rds")
DefaultAssay(b_ctl) <- "SCT"

Data <- b_ctl@assays$SCT@data
write.csv(Data_ctl, "Data.csv")
Idents_Ctl <- Idents(b_ctl)
save(Idents_Ctl, file = "Data/scFEA_fluxomics/case/mouse_cell_ident.RData")

#Cyldko (Mbe1-Cre Cyldflx/flx)
b_cyldko <- readRDS("Seurat_objects/B_CYLDKO.rds")
DefaultAssay(b_cyldko) <- "SCT"

Data_CYLD <- b_cyldko@assays$SCT@data
write.csv(Data_cyldko, "Data_CYLD.csv")
Idents_Cyldko <- Idents(b_cyldko)
save(Idents_Cyldko, file = "Data/scFEA_fluxomics/case/mouse_cell_ident_CYLD.RData")


#--------------------------------------------------------------------------------------
#2. Go to Python to launch scFEA and run the fluxomic analysis (scFEA_CYLD_commands.sh)
#--------------------------------------------------------------------------------------


#3. Post scFEA analysis in Seurat
# add flux as a new assay


#Control (Cyldflx/flx)
flux <- read.csv("Data/scFEA_fluxomics/control/mouse_control_flux.csv", header = T, row.names = 1)
View(mouse_control_flux)
flux <- data.matrix(flux)
flux <- t(flux)
b_ctl <- readRDS("Seurat_objects/B_control.rds")
DefaultAssay(b_srt) <- "SCT"
Bctl[["FLUX"]] <- CreateAssayObject(counts = flux)
DefaultAssay(Bctl) <- 'FLUX'
obj <- Bctl
obj <- FindVariableFeatures(obj,  selection.method = "vst", nfeatures = 2000, verbose = T, assay = "FLUX")
obj <- ScaleData(obj, features = rownames(obj), assay = 'FLUX', verbose = TRUE)
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 10, 
              reduction.name = 'pca.flux', verbose = T, assay = "FLUX")
DimPlot(obj, reduction = "pca.flux")
obj <- FindNeighbors(obj, dims = 1:10, verbose = F, reduction = "pca.flux")
obj <- FindClusters(obj, resolution = 0.2, verbose = F)
obj <- RunUMAP(obj, dims = 1:10, assay = 'FLUX', reduction.name = "umap.flux", verbose = T, reduction = "pca.flux")
saveRDS(obj, "Seurat_objects/BCtl_scFEA.rds")
plot1 <- DimPlot(obj, reduction = "umap.flux", label = T) + ggtitle('UMAP of Flux with scFEA clusters')
plot1
plot2 <- DimPlot(obj, reduction = "umap.flux", group.by = "Idents") + ggtitle('UMAP of Flux with cell annotation')
plot2

plot1 | plot2
plot3 <- DimPlot(BCtl_scFEA, reduction = "umap", group.by = "Idents", label = T)

Idents(obj) <- obj$Idents
obj.markers <- FindAllMarkers(obj, only.pos = F, logfc.threshold = 0.1, 
                              features = VariableFeatures(object = obj), assay = 'FLUX', 
                              slot = "scale.data", test.use = "wilcox")
obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_diff) -> top10
plot4 <- DoHeatmap(obj, features = top10$gene) + NoLegend()

#Inspecting Glycolysis
Supermodules <- read_csv("Data/scFEA_fluxomics/mouse_metabolic_map.csv")
Supermodules$Supermodule_id
Glycolysis <- with(Supermodules, Modules[Supermodule_id == 1])
plot5 <- DoHeatmap(obj, features = Glycolysis, assay = 'FLUX', slot = 'scale.data', 
                group.by = 'Idents', size = 4) + NoLegend()


#Cyldko (Mbe1-Cre Cyldflx/flx)
flux <- read.csv("Data/scFEA_fluxomics/case/mouse_case_flux.csv", header = T, row.names = 1)
flux <- data.matrix(flux)
flux <- t(flux)
b_cyldko <- readRDS("Seurat_objects/B_CYLDKO.rds")
DefaultAssay(b_cyldko) <- "SCT"
b_cyldko[["FLUX"]] <- CreateAssayObject(counts = flux)
DefaultAssay(b_cyldko) <- 'FLUX'
obj <- Bctl
obj <- FindVariableFeatures(obj,  selection.method = "vst", nfeatures = 2000, verbose = T, assay = "FLUX")
obj <- ScaleData(obj, features = rownames(obj), assay = 'FLUX', verbose = TRUE)
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 10, 
              reduction.name = 'pca.flux', verbose = T, assay = "FLUX")
DimPlot(obj, reduction = "pca.flux")
obj <- FindNeighbors(obj, dims = 1:10, verbose = F, reduction = "pca.flux")
obj <- FindClusters(obj, resolution = 0.2, verbose = F)
obj <- RunUMAP(obj, dims = 1:10, assay = 'FLUX', reduction.name = "umap.flux", verbose = T, reduction = "pca.flux")
saveRDS(obj, "Seurat_objects/Bcyldko_scFEA.rds")
plot6 <- DimPlot(obj, reduction = "umap.flux", label = T) + ggtitle('UMAP of Flux with scFEA clusters')
plot6
plot7 <- DimPlot(obj, reduction = "umap.flux", group.by = "Idents") + ggtitle('UMAP of Flux with cell annotation')
plot7

plot6 | plot7
plot8 <- DimPlot(BCtl_scFEA, reduction = "umap", group.by = "Idents", label = T)

Idents(obj) <- obj$Idents
obj.markers <- FindAllMarkers(obj, only.pos = F, logfc.threshold = 0.1, 
                              features = VariableFeatures(object = obj), assay = 'FLUX', 
                              slot = "scale.data", test.use = "wilcox")
obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_diff) -> top10
plot9 <- DoHeatmap(obj, features = top10$gene) + NoLegend()

#Inspecting Glycolysis
# Supermodules <- read_csv("Data/scFEA_fluxomics/mouse_metabolic_map.csv")
# Supermodules$Supermodule_id
# Glycolysis <- with(Supermodules, Modules[Supermodule_id == 1])
plot10 <- DoHeatmap(obj, features = Glycolysis, assay = 'FLUX', slot = 'scale.data', 
                   group.by = 'Idents', size = 4) + NoLegend()



 
#Figure S9
plotS9 <- (plot1 | plot2) / (plot6 | plot7) + patchwork::plot_annotation(tag_levels = "A")
ggsave("figureS9.png", plot = plotS9, device = "png", width = 10, height = 14)


#Figure S10
plotS10 <-  (plot4 / plot9) + patchwork::plot_annotation(tag_levels = "A")


#Figure S11
plotS11 <- (plot5 + plot10) + patchwork::plot_annotation(tag_levels = "A")
