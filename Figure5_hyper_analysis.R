#Using the hypeR package, performing pathway enrichment to capture biological insights for gene modules
#from SCORPIUS trajectory analysis

library(hypeR)
library(readxl)
setwd("./CYLD_project")


#1. Data preparation for control case
SCORPIUS_modules_B_Ctl <- read_excel("Data/SCORPIUS_trajectories/SCORPIUS_modules_B_Ctl.xlsx")
M1 <- with(SCORPIUS_modules_B_Ctl, feature[module == "1"])
M2 <- with(SCORPIUS_modules_B_Ctl, feature[module == "2"])
M3 <- with(SCORPIUS_modules_B_Ctl, feature[module == "3"])
M4 <- with(SCORPIUS_modules_B_Ctl, feature[module == "4"])
M5 <- with(SCORPIUS_modules_B_Ctl, feature[module == "5"])



#2. Actually perform hypeR analysis
hyperSCENIC <- list(M1=as.character(M1), M2=as.character(M2),
                    M3=as.character(M3), M4=as.character(M4),
                    M5=as.character(M5))
genesets <- msigdb_gsets("Mus musculus", clean=TRUE, category = "C2", subcategory = "CP:REACTOME")
print(genesets)
mhyp <- hypeR(hyperSCENIC, genesets, test="hypergeometric", background=30000)
p1 <- hyp_dots(mhyp, merge=TRUE, fdr=0.05, title="SCORPIUS modules BMCtl_Reactome") #include in Figure 5
p2 <- hyp_dots(mhyp, merge=TRUE, pval=0.05, title="SCORPIUS modules BMCtl_Reactome")
ggsave("SCORPIUS_BCtl_Reactome_fdr.png", plot = p1, device = "png", width = 10, height = 14)
ggsave("SCORPIUS_BCtl_Reactome_pvalue.png", plot = p2, device = "png", width = 10, height = 14)

genesets <- msigdb_gsets("Mus musculus", clean=TRUE, category = "C5", subcategory = "GO:BP")
print(genesets)
mhyp <- hypeR(hyperSCENIC, genesets, test="hypergeometric")
p3 <- hyp_dots(mhyp, merge=TRUE, fdr=0.05, title="SCORPIUS modules BMCtl_GO:BP") #include in Figure 5
p4 <- hyp_dots(mhyp, merge=TRUE, pval=0.05, title="SCORPIUS modules BMCtl_GO:BP")
ggsave("SCORPIUS_BCtl_GOBP_fdr.png", plot = p3, device = "png", width = 10, height = 14)
ggsave("SCORPIUS_BCtl_GOBP_pvalue.png", plot = p4, device = "png", width = 10, height = 14)



#3. Data preparation for cyldko case 
SCORPIUS_modules_B_cyldko <- read_excel("Data/SCORPIUS_trajectories/SCORPIUS_modules_B_cyldko.xlsx")
M1 <- with(SCORPIUS_modules_B_cyldko, feature[module == "1"])
M2 <- with(SCORPIUS_modules_B_cyldko, feature[module == "2"])
M3 <- with(SCORPIUS_modules_B_cyldko, feature[module == "3"])



#4. Actually perform hypeR analysis
hyperSCENIC <- list(M1=as.character(M1), M2=as.character(M2),
                    M3=as.character(M3))
genesets <- msigdb_gsets("Mus musculus", clean=TRUE, category = "C2", subcategory = "CP:REACTOME")
print(genesets)
mhyp <- hypeR(hyperSCENIC, genesets, test="hypergeometric")
p5 <- hyp_dots(mhyp, merge=TRUE, fdr=0.05, title="SCORPIUS modules BMCYLDKO_Reactome") #include in Figure 5
p6 <- hyp_dots(mhyp, merge=TRUE, pval=0.05, title="SCORPIUS modules BMCYLDKO_Reactome")
ggsave("SCORPIUS_Bcyldko_Reactome_fdr.png", plot = p5, device = "png", width = 10, height = 14)
ggsave("SCORPIUS_Bcyldko_Reactome_pvalue.png", plot = p6, device = "png", width = 10, height = 14)

genesets <- msigdb_gsets("Mus musculus", clean=TRUE, category = "C5", subcategory = "GO:BP")
print(genesets)
mhyp <- hypeR(hyperSCENIC, genesets, test="hypergeometric")
p7 <- hyp_dots(mhyp, merge=TRUE, fdr=0.05, title="SCORPIUS modules BMCYLDKO_GO:BP") #include in Figure 5
p8 <- hyp_dots(mhyp, merge=TRUE, pval=0.05, title="SCORPIUS modules BMCYLDKO_Reactome")
ggsave("SCORPIUS_Bcyldko_GOBP_fdr.png", plot = p7, device = "png", width = 10, height = 14)
ggsave("SCORPIUS_Bcyldko_GOBP_pvalue.png", plot = p8, device = "png", width = 10, height = 14)
