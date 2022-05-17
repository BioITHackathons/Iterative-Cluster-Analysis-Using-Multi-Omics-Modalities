# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---
install.packages("pvclust")
install.packages("tidyverse")

library(parallel)
library (pvclust)
library(tidyverse)

setwd ("/Users/deslattesmaysa2/projects/Iterative-Cluster-Analysis-Using-Multi-Omics-Modalities/data")
getwd()

tumor_gene_expression_file<- "tumor-gene-expression-rsem-tpm-collapsed.tsv"

tumor_gene_expression <- read_table(tumor_gene_expression_file)
tumor_gene_expression[1:4,1:4]
dim(tumor_gene_expression)

tge <- column_to_rownames(tumor_gene_expression,var="gene_id")
tge[1:4,1:4]
dim(tge)
typeof(tge)

gcm <- t(tge)
gcm[1:4,1:4]
dim(gcm)
typeof(gcm)
is.matrix(gcm)

ttge <- t(tge)
typeof(ttge)
is.matrix(ttge)


# +
# makeCluster is not allowing more than 10 cores to be made
# -

cl <- makeCluster(detectCores()-2)
cl
# +
# run the clustering by gene with the gcm matrix
# -
tinygcm<-gcm[1:1000,1:1000]
tinygcm.pv <- parPvclust(cl, tinygcm, method.hclust="ward.D2",
                        method.dist="minkowski", use.cor="pairwise.complete.obs",
                        nboot=10, r=seq(.5,1.4,by=.1))


# plot a 6x3 inches cluster dendogram and draw the au=95%
pdf("/sbgenomics/output-files/cluster_by_gene_plot.pdf", width=6, height=3)
plot(gcm.pv)
pvrect(gcm.pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)
dev.off()

pick <- pvpick(gcm.pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)

cluster_assignments <- do.call('rbind',lapply(1:length(pick$clusters),function(i){
  cbind(rep(i,length(pick$clusters[[i]])),pick$clusters[[i]])
}))
colnames(cluster_assignments) <- c("cluster","gene")

write.csv(cluster_assignments, "/sbgenomics/output-files/cluster_by_gene_assignments.csv",row.names=FALSE)


ttge.pv <- parPvclust(cl, ttge, method.hclust="ward.D2",
           method.dist="minkowski", use.cor="pairwise.complete.obs",
           nboot=1000, r=seq(.5,1.4,by=.1))


# 6x3 inches
pdf("/sbgenomics/output-files/cluster_by_tissue_plot.pdf", width=6, height=3)
plot(ttge.pv)
pvrect(ttge.pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)
dev.off()

pick <- pvpick(ttge.pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)
cluster_assignments <- do.call('rbind',lapply(1:length(pick$clusters),function(i){
    cbind(rep(i,length(pick$clusters[[i]])),pick$clusters[[i]])
}))
colnames(cluster_assignments) <- c("cluster","biospecimen")

write.csv(cluster_assignments, "/sbgenomics/output-files/cluster_by_tissue_assignments.csv",row.names=FALSE)



