ProB <- with(SCENIC_BMCtl_TFs, feature[module == "ProB"])
PreB <- with(SCENIC_BMCtl_TFs, feature[module == "PreB"])
Immature_B <- with(SCENIC_BMCtl_TFs, feature[module == "ImmatureB"])
Mature_B<- with(SCENIC_BMCtl_TFs, feature[module == "MatureB"])


hyperSCENIC <- list(ProB=as.character(ProB), PreB=as.character(PreB),
                    Immature_B=as.character(Immature_B), Mature_B=as.character(Mature_B))

genesets <- msigdb_gsets("Mus musculus", clean=TRUE, category = "H")
                         #subcategory = "CP:REACTOME")
print(genesets)
mhyp <- hypeR(hyperSCENIC, genesets, test="hypergeometric")
hyp_dots(mhyp, merge=TRUE, pval=0.05, title="SCENIC modules BMCtl", val = "pval", sizes = T)
hyp_dots(mhyp, merge=TRUE, fdr=0.05, title="SCENIC modules BMCtl", val = "fdr", sizes = T)



#BMCyldKO
ProB <- with(SCENIC_BMCyldKO_TFs, feature[module == "ProB"])
PreB <- with(SCENIC_BMCyldKO_TFs, feature[module == "PreB"])
ProPreB <- with(SCENIC_BMCyldKO_TFs, feature[module == "ProPreB"])
Immature_B <- with(SCENIC_BMCyldKO_TFs, feature[module == "ImmatureB"])
Mature_B<- with(SCENIC_BMCyldKO_TFs, feature[module == "MatureB"])


hyperSCENIC <- list(ProB=as.character(ProB), PreB=as.character(PreB), ProPreB=as.character(ProPreB),
                    Immature_B=as.character(Immature_B), Mature_B=as.character(Mature_B))

genesets <- msigdb_gsets("Mus musculus", clean=TRUE, category = "H")
#subcategory = "CP:REACTOME")
print(genesets)
mhyp <- hypeR(hyperSCENIC, genesets, test="hypergeometric", background=30000)
hyp_dots(mhyp, merge=TRUE, pval=0.05, title="SCENIC modules BMCyldKO", val = "pval", sizes = T)
hyp_dots(mhyp, merge=TRUE, fdr=0.05, title="SCENIC modules BMCyldKO", val = "fdr", sizes = T)

?hypeR()
?hyp_dots()
?hypeR::hyp_emap()
hypeR::hyp_emap(hyp_obj = mhyp, similarity_metric = c("jaccard_similarity", "overlap_similarity"),
                similarity_cutoff = 0.2, val = "fdr")
hypeR::hyp_hmap(hyp_obj = mhyp)
