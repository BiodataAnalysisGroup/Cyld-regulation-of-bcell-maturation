## Scenario 2 scMINER

library(scMINER) # version 0.1.0
library(Seurat)
library(tidyseurat)
library(openxlsx)
library(tidyverse)
library(SingleCellExperiment)

rm(list = ls())
gc()

seurat_bcells <- readRDS("cyldIntegrated_STrv2.rds")

feature.data<-data.frame(rownames(seurat_bcells@assays$integrated@data))
colnames(feature.data) <- "geneSymbol"
rownames(feature.data) <- feature.data$geneSymbol

meta.data <- seurat_bcells@meta.data
meta.data$CellNames <- rownames(meta.data)

meta.data$orig.ident <- str_remove_all(meta.data$orig.ident, "BM_")

eset <- CreateSparseEset(data=seurat_bcells@assays$integrated@data,
                       meta.data = meta.data,
                       feature.data = feature.data,
                       add.meta = T)

# Create SJARACNe Objects

generateSJARACNeInput(
  input_eset = eset, funcType = "TF", # Considered all genes not only TFs, "TF" is necessary to run the algorithm 
  ref = "mm",  #mouse
  wd.src = "SJAR/SJAR_TF",  #Output directory
  group_name = "orig.ident") # KO and Ctl

## Run SJARACNe with bash

## Let's use B cell as an example
# For Ctl
# sjaracne local -e Ctl_3000_3000_629.txt \
  # -g all_genes.txt \
  # -o output_all_genes_Ctl \
  # -n 100 \
  # -pc 0.01

# For KO
# sjaracne local -e ΚΟ_2999_2999_695.txt \
# -g all_genes.txt \
# -o output_all_genes_KO \
# -n 100 \
# -pc 0.01

acs.lupin <- GetActivityFromSJARACNe(
  SJARACNe_output_path ="SJAR/SJAR_TF", # The folder TF is by default the folder of output
  SJARACNe_input_eset = eset,
  activity.method="unweighted", # we highly recommend using 'unweighted' as activity calculation method
  activity.norm = FALSE, 
  group_name = "orig.ident", # which group was used to partition expression profiles
  save_network_file = FALSE, 
  functype = "tf",# whether or not save network for each group
  save_path = NULL) 

saveRDS(acs.lupin, file = "acs.lupin.rds")

DAG_result <- get.DA(input_eset = acs.lupin, group_name = "orig.ident")

saveRDS(DAG_result, file = "DAG_result_tf.rds")
write.xlsx(DAG_result, "DAG.xlsx")
