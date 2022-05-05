# Iterative cluster analysis using multi-omics modalities and interpretation with the data translator

This repository contains the source code and results of the iterative clustering of multiomics data and interpretation with the Biomedical Data Translator project. This project began at the 2022 Bio-IT FAIR Data Hackathon. We use the [Gabriella Miller Kids First Data Resource Center](https://kidsfirstdrc.org/) data supported by the NIH Common Fund--this resource contains data from over 11,000 samples, including DNA and RNA as well as clinical information. 

In this project, we initially focused on clustering of gene expression profiles from RNA-Seq data collected from pediatric tumor samples. We then create a simple and interpretable predictive model to determine the gene expression signatures that differentiate the clusters from one another. To gain additional translational insight into the clusters we sought to annotate the important genes from each cluster with data from the [NCATS Biomedical Data Translator](https://ncats.nih.gov/translator) ([github org](https://github.com/NCATSTranslator)). Our analyses were executed on the [Cavatica](https://www.cavatica.org/) cloud-based data analysis and sharing  platform. 

## Summary

Our workflow consisted of the following core steps:

* Wrangle the NICHD Kids First data
    * See the following notebooks
* Perform unsupervised clustering of gene expression gene expression using [pvclust](https://cran.r-project.org/web/packages/pvclust/pvclust.pdf), which is hierarchichal clustering approach that implements a bootstrapping method for assessing statistical significance of clusters
* Develop a classification model using [xgboost](https://xgboost.readthedocs.io/en/stable/index.html) to predict the cluster assignments from the gene expression data. The xgboost model provides feature importance metrics that we use the 
* Annotate results by querying the NCATS Biomedical Data Translator

## Future directions

### Data types
In the future we hope to integrate additional omics data modalities available through the Kids First Data Resrouce-such as somatic mutation calls from tumor sequencing, HPO phenotypes and patient clinical characteristics-and additional disease states like the [INCLUDE](https://www.nih.gov/include-project) project focused on Down Syndrome.

### Methodology
In the future we would like to explore the use of feature selection methods, such as recursive feature elimination, to reduce the number of genes required to make cluster predictions. We would then iterate on the clustering process to see if pvclust performs better on the reduced feature set. 

## Technical details

### Platforms
* Cloud platform: [Cavatica](https://www.cavatica.org/)
* Biomedical Data Translator: [github org](https://github.com/NCATSTranslator)

### Dependencies

* Python
    * [scikit learn](https://scikit-learn.org/stable/)
    * [xgboost](https://xgboost.readthedocs.io/en/stable/index.html)
* R
    * [pvclust](https://cran.r-project.org/web/packages/pvclust/pvclust.pdf)