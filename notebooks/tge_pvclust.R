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

library(parallel)
library (pvclust)
library(tidyverse)
tumor_gene_expression_file<- "/sbgenomics/project-files/tumor-gene-expression-rsem-tpm-collapsed.tsv"

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

numCores <-detectCores()
numCores

# +
# makeCluster is not allowing more than 10 cores to be made
# -

cl <- makeCluster(10)

# +
#"parPvclust" has been integrated into pvclust (with "parallel" option).
# It is available for back compatibility but will be unavailable in the future.”
# -

ttge_t.pv <- parPvclust(cl, ttge, method.hclust="ward.D2",
           method.dist="minkowski", use.cor="pairwise.complete.obs",
           nboot=10, r=seq(.5,1.4,by=.1))

plot(ttge_t.pv)

pvrect(ttge_t.pv)

pick <- pvpick(ttge.pv, alpha=0.80, pv="au", type="geq", max.only=TRUE)

cluster_assignments <- do.call('rbind',lapply(1:length(pick$clusters),function(i){
    cbind(rep(i,length(pick$clusters[[i]])),pick$clusters[[i]])
}))
colnames(cluster_assignments) <- c("cluster","biospecimen")

cluster_assignments



write.csv(cluster_assignments, "cluster_assignments.csv",row.names=FALSE)


