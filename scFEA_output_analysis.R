#'
#'This R script produces the histographs for each cell type per module.
#'
#'Input: scFEA flux predictions file, file with all barcodes and cell types.
#'
#'Output: histographs for each cell type per module
#'

rm(list = ls())
gc()

library(data.table)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(reshape)
library(ggridges)
library(ggplot2)


# load predicted flux---------------------------------------------------------
path = "case/mouse_case_flux.csv"
data_c <- read.csv(path)
data_c0 <- as.matrix(data_c[,-1])
#first column of data_c becomes row names of data_c0
rownames(data_c0) <- as.character(data_c[,1])
#transpose
data_c0 <- t(data_c0)
ppp_all <-c()

#load factors with barcodes and cell types
load('case/mouse_cell_ident_CYLD.RData')

#keep the elements of cell_id from column names of data_c0
yyy <- Idents_CYLD[colnames(data_c0)]
for(ii in 1:nrow(data_c0)){
    xxx <- data_c0[ii,]
    final_df <- cbind(paste('X', 1:length(xxx), sep=''), xxx, yyy)
    final_df <- as.data.frame(final_df)
    final_df[,2] <- as.numeric(final_df[,2])
    colnames(final_df) <- c('var', 'flux', 'cellType')
    pp <- sd(final_df$flux)/abs(mean(final_df$flux))
    ppp_all <- c(ppp_all, pp)
}
#take those modules with npn zero flux (here set as 1e-10, re-adjust?)
tg_ids <- which(ppp_all > 1e-10)

# load mouse module info
load('mouse_module_info.RData')

#create list of modules number to be plotted e.g. 1 for M_1
plot_modules = c(71, 1, 2, 3, 4, 6, 78, 33, 15, 5, 7,8,9,10,11,12,13,18,92,91,48,49,51)
plot_path_save = sub("_flux.csv.*", "_", path)  

for (ss in plot_modules){
#insert module number of interest to get the row name position
# jj is position/row of Module ss
jj = which(tg_ids == ss) 
#check if Module ss exists in tg_ids
if(identical(jj, integer(0))){ 
    print(paste0("Module M_", ss," did not pass threshold."))
}else{
if(length(tg_ids) > 0){
        #jj = 20 # check module 2, runs on all modules
        xxx <- data_c0[paste0("M_",tg_ids[jj]), ]
        final_df <- data.frame(var = paste('X', 1:length(xxx), sep=''),
                               flux = xxx,
                               cellType = yyy)
        #mouse
        title <- mouse_module_info[paste0("M_",tg_ids[jj]), 'M_name']

        
        aa <- ggplot(final_df, aes(x = flux, y = cellType, fill = cellType)) +
            geom_density_ridges() +
            theme_ridges() +
            theme(legend.position = 'none') +
            ggtitle(title) +
            theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(
            filename = paste0( plot_path_save ,"module_", ss,".pdf"),
            width = 13, height = 10, units = "in"
        )
        }
    }  
}


