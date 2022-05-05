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

tumor_gene_expression_file<- "/sbgenomics/project-files/tumor-gene-expression-rsem-tpm-collapsed.tsv"

tumor_gene_expression <- read.table(tumor_gene_expression_file, header=TRUE, row.names=1, sep="\t")

head (tumor_gene_expression)


dim(tumor_gene_expression)

patient_cluster_matrix <- t(tumor_gene_expression)

patient_cluster_matrix[1:4,1:4]

reduced_pcm <- patient_cluster_matrix[1:100, 1:100]

head(reduced_pcm)

reduced_pcm_t <- t(reduced_pcm)

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

pcm_t.pv <- parPvclust(cl, reduced_pcm_t, method.hclust="average",
           method.dist="correlation", use.cor="pairwise.complete.obs",
           nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE,
           init.rand=NULL, iseed=NULL, quiet=FALSE)

p <-plot(reduced_pcm_t.pv)

x <- seplot(reduced_pcm_t.pv, identify=TRUE)

pvrect(p)

reduced_pcm_t.pv$labels

pick <- pvpick(reduced_pcm_t.pv, alpha=0.60, pv="au", type="geq", max.only=TRUE)

cluster_assignments <- do.call('rbind',lapply(1:length(pick$clusters),function(i){
    cbind(rep(i,length(pick$clusters[[i]])),pick$clusters[[i]])
}))
colnames(cluster_assignments) <- c("cluster","biospecimen")

cluster_assignments



write.csv(cluster_assignments, "cluster_assignments.csv",row.names=FALSE)


