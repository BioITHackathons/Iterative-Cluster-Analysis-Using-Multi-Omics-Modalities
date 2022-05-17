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

current_dir     <- getwd()
current_dir
clone_data_dir  <- "Iterative-Cluster-Analysis-Using-Multi-Omics-Modalities/data"
clone_data_dir
tge_file        <- "tumor-gene-expression-rsem-tpm-collapsed.tsv"
tge_file
tge_with_dir    <- paste(paste(current_dir, clone_data_dir,sep="/"),tge_file, sep="/")
tge_with_dir

tumor_gene_expression <- read_table(tumor_gene_expression_file)
tumor_gene_expression[1:4,1:4]
dim(tumor_gene_expression)

tge <- column_to_rownames(tumor_gene_expression,var="gene_id")
tge[1:4,1:4]
dim(tge)
typeof(tge)

gcm <- t(tge)
dim(gcm)
typeof(gcm)
is.matrix(gcm)

ttge <- t(tge)
tge <- t(ttge)
typeof(tge)
is.matrix(tge)

gcm[1:4,1:4]
tge[1:4,1:4]

# +
# makeCluster is not allowing more than 10 cores to be made
# -

cl <- makeCluster(detectCores()-2)
cl
tinytge<- tge[1:1000,1:1000]

tinytge.pv <- parPvclust(cl, tinytge, method.hclust="ward.D2",
           method.dist="minkowski", use.cor="pairwise.complete.obs",
           nboot=10, r=seq(.5,1.4,by=.1))
plot(tinytge.pv)
pvrect(tinytge.pv)

pick <- pvpick(tinytge.pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)


cluster_assignments <- do.call('rbind',lapply(1:length(pick$clusters),function(i){
  cbind(rep(i,length(pick$clusters[[i]])),pick$clusters[[i]])}))
colnames(cluster_assignments) <- c("cluster","gene")
#

output_dir      <- "/sbgenomics/output-files"
output_dir
outfilename <- "tiny_tissue_by_gene_experiment_cluster.csv"
output_dir_file <- paste(output_dir, outfilename, sep="/")
output_dir_file

write.csv(cluster_assignments, output_dir_file,row.names=FALSE)

tgenboot10.pv <- parPvclust(cl, tge, method.hclust="ward.D2",
                         method.dist="minkowski", use.cor="pairwise.complete.obs",
                         nboot=10, r=seq(.5,1.4,by=.1))

pick <- pvpick(tgenboot10.pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)


cluster_assignments <- do.call('rbind',lapply(1:length(pick$clusters),function(i){
  cbind(rep(i,length(pick$clusters[[i]])),pick$clusters[[i]])}))
colnames(cluster_assignments) <- c("cluster","gene")
#

output_dir      <- "/sbgenomics/output-files"
output_dir
outfilename <- "tge_nboot10_tissue_by_gene_experiment_cluster.csv"
output_dir_file <- paste(output_dir, outfilename, sep="/")
output_dir_file

write.csv(cluster_assignments, output_dir_file,row.names=FALSE)install.packages("pvclust")
