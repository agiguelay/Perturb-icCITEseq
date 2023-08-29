##########################################################
#                                                        #
#        Genome-wide CRISPR screen in human T cells      #
#              reveals regulators of FOXP3               #
#                                                        #
##########################################################
#                                                        #
#                  Module definition                     #
#                                                        #
##########################################################

rm(list = ls())

# Libraries --------------------------------------------------------------------

library(umap)
library(circlize)

# Load data --------------------------------------------------------------------

beta = read.table("Beta_matrix_LR4.txt", check.names = F)
beta = beta[-c(1:7),]
ncol(beta)*nrow(beta)
sum(beta != 0)

# Gene filtering ---------------------------------------------------------------
beta = beta[rowSums(abs(beta)) != 0 , colSums(abs(beta)) != 0]
beta = beta[, apply(abs(beta), 2, function(x) sum(x >= 0.08)) > 3] 

# Target modules ---------------------------------------------------------------
target.mod = kmeans(beta, centers = 20, nstart = 500)
target.mod = data.frame(Target = names(target.mod$cluster), Module = target.mod$cluster)
target.mod$Gene = sapply(target.mod$Target, function(x) strsplit(x, "_")[[1]][2])
target.mod = target.mod[order(target.mod$Module),]

target.umap = umap(beta, n_neighbors = 5, min_dist = 0.01)
df.umap = data.frame(UMAP1 = target.umap$layout[,1], UMAP2 = target.umap$layout[,2])
target.mod = cbind(target.mod, df.umap[target.mod$Target,])

write.csv2(target.mod, "TargetModules_kmeans20.csv")

# Gene modules -----------------------------------------------------------------
gene.mod = kmeans(t(beta), centers = 20, nstart = 500, iter.max = 20)
gene.mod = data.frame(Gene = names(gene.mod$cluster), Module = gene.mod$cluster)
gene.mod = gene.mod[order(gene.mod$Module),]

gene.umap = umap(t(beta), n_neighbors = 20, min_dist = 0.1)
df.umap = data.frame(UMAP1 = gene.umap$layout[,1], UMAP2 = gene.umap$layout[,2])
gene.mod = cbind(gene.mod, df.umap[gene.mod$Gene,])

write.csv2(gene.mod, "GeneModules_kmeans20.csv")

