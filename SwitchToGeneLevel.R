##########################################################
#                                                        #
#        Genome-wide CRISPR screen in human T cells      #
#              reveals regulators of FOXP3               #
#                                                        #
##########################################################
#                                                        #
#               Switch to the gene level                 #
#                                                        #
##########################################################

rm(list = ls())

# Libraries --------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(cape)

# Compute the residual matrix --------------------------------------------------
coeff = fread(paste0("Beta_matrix_LR2.txt"), sep='\t', data.table = F)
rownames(coeff) = coeff$V1
coeff = coeff[grep("gRNA", rownames(coeff)),]
coeff = coeff[,-1]
coeff = coeff[rowSums(coeff)!=0 | grepl("_NTC-", rownames(coeff)),]
coeff = coeff[rowSums(coeff)!=0,]

corr.mtx = cor(t(coeff))

# Null distribution of inter-target correlations -------------------------------
pairs = data.frame(gRNA1 = as.vector(sapply(rownames(corr.mtx), function(x) rep(x, nrow(corr.mtx)))),
                   gRNA2 = rep(rownames(corr.mtx), nrow(corr.mtx)),
                   Correlation = as.vector(corr.mtx))

pairs$Target1 = target[pairs$gRNA1]
pairs$Target2 = target[pairs$gRNA2]

pairs.null = pairs[pairs$Target1 != pairs$Target2,]
pairs.null = pairs.null[as.character(pairs.null$Target1) > as.character(pairs.null$Target2),] # Because the mtx is symetrical

pairs.target = pairs[pairs$Target1 == pairs$Target2,]
pairs.target = pairs.target[pairs.target$gRNA1 > pairs.target$gRNA2,]
calc_emp_p <- function(obs_dist, null_dist){
  p_fun <- ecdf(null_dist)
  emp_p <- 1-unlist(lapply(obs_dist, p_fun))
  return(emp_p)
}

pairs.target$emp.pval = calc_emp_p(obs_dist = pairs.target$Correlation, null_dist = pairs.null$Correlation)
pairs.target$adj.pval = p.adjust(pairs.target$emp.pval, method = "BH", n = length(pairs.target$emp.pval))

write.table(pairs.target, paste0("Pvalues_gRNAs.txt"), sep = "\t")
sel = pairs.target[pairs.target$adj.pval <= 0.5,]

