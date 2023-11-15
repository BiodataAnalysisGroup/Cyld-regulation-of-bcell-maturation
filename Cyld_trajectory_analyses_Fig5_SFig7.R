#Trajectory analysis for Figure 5, Supplementary Figure 7
library(Seurat)
library(SCORPIUS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
setwd("./CYLD_project")

#1.SCORPIUS to plot BMCtl Seurat object
b_srt <- readRDS("Seurat_objects/B_control.rds")
DefaultAssay(b_srt) <- "SCT"

#Run SCORPIUS on UMAP co-ordinates
space <- Embeddings(b_srt, reduction = "umap") 
expression <- t(as.matrix(b_srt@assays$SCT@scale.data))
idents <- b_srt$Idents

#Calculate Trajectory Inference
draw_trajectory_plot(space, progression_group = idents, 
                     progression_group_palette = c("Pro B" = "magenta", 
                                                   "Pre B" = "royalblue", 
                                                   "Immature B" = "red", 
                                                   "Mature B" = "green"), 
                     contour = FALSE)
traj <- infer_trajectory(space)
p1 <- draw_trajectory_plot(
  space, 
  progression_group = idents, 
  progression_group_palette = c("Pro B" = "magenta", "Pre B" = "royalblue", 
                                "Immature B" = "red", "Mature B" = "green"),
  path = traj$path,
  contour = FALSE
)
p1

# Dif. expressed genes
gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
gene_sel <- gimp[1:100,]
expr_sel <- expression[,gene_sel$gene]
draw_trajectory_heatmap(expr_sel, traj$time, idents, progression_group_palette = 
                          c("Pro B" = "magenta", "Pre B" = "royalblue", "Immature B" = "red", 
                            "Mature B" = "green"))

# Make modules of trajectory-related genes
modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
p2 <- draw_trajectory_heatmap(expr_sel, rever$time, seuratclusters, modules, show_labels_row = T,
                              fontsize_row=6, progression_group_palette =
                                c("Pro B" = "magenta", "Pre B" = "royalblue", "Immature B" = "red",
                                  "Mature B" = "green"))
p2

# # SCORPIUS is trajectory-agnostic so it is possible that we need to reverse the direction of the trajectory
# # Reverse trajectory if Pro B cells are not at Pseudotime 0
# rever <- reverse_trajectory(traj)
# p1 <- draw_trajectory_plot(
#   space,
#   progression_group = idents,
#   progression_group_palette = c("Pro B" = "magenta", "Pre B" = "royalblue",
#                                 "Immature B" = "red", "Mature B" = "green"),
#   path = traj$path,
#   contour = FALSE
# ) & ggtitle("Trajectory")
# p1
# gimp <- gene_importances(expression, rever$time, num_permutations = 0, num_threads = 8)
# gene_sel <- gimp[1:100,]
# expr_sel <- expression[,gene_sel$gene]
# draw_trajectory_heatmap(expr_sel, rever$time, idents, progression_group_palette =
#                           c("Pro B" = "magenta", "Pre B" = "royalblue", "Immature B" = "red",
#                             "Mature B" = "green"))
# 
# 
# 
# #Make modules of trajectory-related genes
# modules <- extract_modules(scale_quantile(expr_sel), rever$time, verbose = T)
# p2 <- draw_trajectory_heatmap(expr_sel, rever$time, idents, modules, show_labels_row = T,
#                         fontsize_row=6, progression_group_palette =
#                           c("Pro B" = "magenta", "Pre B" = "royalblue", "Immature B" = "red",
#                             "Mature B" = "green"))
# p2


#Plot pseudotime on Seurat object "b_srt"
b_srt <- AddMetaData(object=b_srt, 
                     metadata = rever$time,  #or rever$time if you have reversed above
                     col.name="pseudotime")
p3 <- draw_trajectory_plot(
  space,
  progression_group = rever$time, 
  path = traj$path,
  contour = FALSE
) & scale_color_viridis_c() & ggtitle("Pseudotime") & guides(color = guide_legend(title = "Pseudotime"))
p3


# Focus on the top-20 most predictive genes of the trajectory
#a
# For control case (CYLDflx/flx)
SCORPIUS_modules_B_Ctl <- read_excel("Data/SCORPIUS_trajectories/SCORPIUS_modules_B_Ctl.xlsx")
SCORPIUS_modules_B_Ctl <- SCORPIUS_modules_B_Ctl[order(SCORPIUS_modules_B_Ctl$orig_index),]
features <- SCORPIUS_modules_B_Ctl$feature[1:20]

expr_sel_test <- expression[,features]
draw_trajectory_heatmap(expr_sel, traj$time, idents)
gene_sel_S <-features

#b.
# expr_sel_S <- t(expr_sel)
expr_sel_S <- t(expr_sel_test)
expr_sel_S <- as.matrix(expr_sel_S)

#c.
ident_Seurat <- as.character(b_srt$Idents)
traj_time <- as.numeric(traj$time)
# traj_time <- as.numeric(rever$time) # depending on the trajectory direction
matrix <- rbind(expr_sel_S,traj_time,ident_Seurat)

#d.
expr_S_matrix <- as.data.frame(t(matrix))

#e.
S_df <- pivot_longer(expr_S_matrix, gene_sel_S, names_to = 'feature', values_to = 'expr')
S_df$expr <- as.numeric(as.character(S_df$expr )) #Some of the columns changed to characters after `pivot_longer` for some reason.
S_df$traj_time <- as.numeric(as.character(S_df$traj_time)) 

#f.
p <- ggplot(S_df, mapping = aes(x=traj_time, y=expr, color=ident_Seurat)) + geom_jitter(size=1) + theme_classic() + xlab('traj_time') + ylab('Expression') + theme(plot.title = element_text(size=16, hjust =  0.5, face = 'bold'), strip.text = element_text(size=12, face = 'bold'),strip.background = element_rect(size = 0)) + guides(color = guide_legend(override.aes = list(linetype = 'blank'))) + scale_y_log10() + facet_wrap(~feature,scales = "free_y") 
p_col <- p + scale_color_manual(values=c("red", "green", "royalblue", "magenta"))
pline1 <- p_col + geom_smooth(aes(color = expr), method = 'gam', se=F, color = 'black') 
pline1






#2.SCORPIUS to plot BM CyldKO (Mbe1 Cre Cyldflx/flx) Seurat object
b_srt <- readRDS("Seurat_objects/B_CYLDKO.rds")
DefaultAssay(b_srt) <- "SCT"
#Run SCORPIUS on UMAP co-ordinates
space <- Embeddings(b_srt, reduction = "umap") 
expression <- t(as.matrix(b_srt@assays$SCT@scale.data))
idents <- b_srt@active.ident

#Calculate Trajectory Inference
draw_trajectory_plot(space, progression_group = idents, 
                     progression_group_palette = c("Pro B" = "magenta",
                                                   "Pro/Pre B" = "orange",
                                                   "Pre B" = "royalblue", 
                                                   "Immature B" = "red", 
                                                   "Mature B" = "green"), 
                     contour = FALSE)
traj <- infer_trajectory(space)
p4 <- draw_trajectory_plot(
  space, 
  progression_group = idents, 
  progression_group_palette = c("Pro B" = "magenta", "Pro/Pre B" = "orange", 
                                "Pre B" = "royalblue", 
                                "Immature B" = "red", "Mature B" = "green"),
  path = traj$path,
  contour = FALSE
) + ggtitle("Trajectory")
p4

# dif. expressed genes
gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
gene_sel <- gimp[1:100,]
expr_sel <- expression[,gene_sel$gene]
draw_trajectory_heatmap(expr_sel, traj$time, idents, progression_group_palette = 
                          c("Pro B" = "magenta", "Pro/Pre B" = "orange", "Pre B" = "royalblue", "Immature B" = "red", 
                            "Mature B" = "green"))

#Make modules of trajectory-related genes
modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
p5 <- draw_trajectory_heatmap(expr_sel, traj$time, idents, modules, show_labels_row = T, 
                              fontsize_row=6, progression_group_palette = 
                                c("Pro B" = "magenta", "Pro/Pre B" = "orange", "Pre B" = "royalblue", "Immature B" = "red", 
                                  "Mature B" = "green"))
p5


# # reverse trajectory if Pro B cells are not at Pseudotime 0
# rever <- reverse_trajectory(traj)
# p4 <- draw_trajectory_plot(
#   space,
#   progression_group = idents, 
#   progression_group_palette = c("Pro B" = "magenta", "Pre B" = "royalblue", 
#                                 "Immature B" = "red", "Mature B" = "green"),
#   path = traj$path,
#   contour = FALSE
# ) & ggtitle("Trajectory")
# p4
# gimp <- gene_importances(expression, rever$time, num_permutations = 0, num_threads = 8)
# gene_sel <- gimp[1:100,]
# expr_sel <- expression[,gene_sel$gene]
# draw_trajectory_heatmap(expr_sel, rever$time, seuratclusters, progression_group_palette = 
#                           c("Pro B" = "magenta", "Pre B" = "royalblue", "Immature B" = "red", 
#                             "Mature B" = "green"))
# 
# #Make modules of trajectory-related genes
# modules <- extract_modules(scale_quantile(expr_sel), rever$time, verbose = FALSE)
# p5 <- draw_trajectory_heatmap(expr_sel, rever$time, seuratclusters, modules, show_labels_row = T, 
#                         fontsize_row=6, progression_group_palette = 
#                           c("Pro B" = "magenta", "Pre B" = "royalblue", "Immature B" = "red", 
#                             "Mature B" = "green"))
# p5

#Plot pseudotime on Seurat object "b_srt"
b_srt <- AddMetaData(object=b_srt, 
                     metadata = traj$time,  #or rever$time if you have reversed above
                     col.name="pseudotime")
p6 <- draw_trajectory_plot(
  space,
  progression_group = traj$time, 
  path = traj$path,
  contour = FALSE
) & scale_color_viridis_c() & ggtitle("Pseudotime") & guides(color = guide_legend(title = "Pseudotime"))
p6

# Focus on the top-20 most predictive genes of the trajectory
#a
# For cyldko case (Mbe1CreCyldflx/flx)
SCORPIUS_modules_B_cyldko <- read_excel("Data/SCORPIUS_trajectories/SCORPIUS_modules_B_cyldko.xlsx")
SCORPIUS_modules_B_cyldko <- SCORPIUS_modules_B_cyldko[order(SCORPIUS_modules_B_cyldko$orig_index),]
features <- SCORPIUS_modules_B_cyldko$feature[1:20]

expr_sel_test <- expression[,features]
draw_trajectory_heatmap(expr_sel, traj$time, idents)
gene_sel_S <-features

#b.
# expr_sel_S <- t(expr_sel)
expr_sel_S <- t(expr_sel_test)
expr_sel_S <- as.matrix(expr_sel_S)

#c.
ident_Seurat <- as.character(b_srt$Idents)
traj_time <- as.numeric(traj$time)
# traj_time <- as.numeric(rever$time) # depending on the trajectory direction
matrix <- rbind(expr_sel_S,traj_time,ident_Seurat)

#d.
expr_S_matrix <- as.data.frame(t(matrix))

#e.
S_df <- pivot_longer(expr_S_matrix, gene_sel_S, names_to = 'feature', values_to = 'expr')
S_df$expr <- as.numeric(as.character(S_df$expr )) #Some of the columns changed to characters after `pivot_longer` for some reason.
S_df$traj_time <- as.numeric(as.character(S_df$traj_time)) 

#f.
p <- ggplot(S_df, mapping = aes(x=traj_time, y=expr, color=ident_Seurat)) + geom_jitter(size=1) + theme_classic() + xlab('traj_time') + ylab('Expression') + theme(plot.title = element_text(size=16, hjust =  0.5, face = 'bold'), strip.text = element_text(size=12, face = 'bold'),strip.background = element_rect(size = 0)) + guides(color = guide_legend(override.aes = list(linetype = 'blank'))) + scale_y_log10() + facet_wrap(~feature,scales = "free_y") 
p_col <- p + scale_color_manual(values=c("red", "green", "royalblue", "magenta", "orange"))
pline2 <- p_col + geom_smooth(aes(color = expr), method = 'gam', se=F, color = 'black') 
pline2





# 3. Plotting

# For Figure 5: p1, p2, p3, p4, p5, p6

# For Supl. Figure 7: pline1, pline2