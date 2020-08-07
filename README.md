# hif2gene: Dense, high-resolution mapping of cells and tissues from pathology images for the interpretable prediction of molecular phenotypes in cancer

## Overview
> `hif2gene` is the primary code repository for reproducing analyses in Diao, Chui, and Wang et al. 2020: "Dense, high-resolution mapping of cells and tissues from pathology images for the interpretable prediction of molecular phenotypes in cancer"

## Installation
- Clone this repo to your local machine using https://github.com/jamesdiao/hif2gene.git

## Contents
1. `/scripts` contains all Python and R code used to produce content in `/figures` from `/data`
2. `/data` contains all raw and cached data objects needed to reproduce the analysis from scratch
3. `/figures` contains all primary figures as individual vectorized `.pdf` 

## Version and package requirements 
- R version: 3.6.2
    - `caret`
    - `cluster`
    - `cowplot`
    - `data.table`
    - `devtools`
    - `dplyr`
    - `ggplot2`
    - `ggpubr`
    - `mclust`
    - `openxlsx`
    - `readxl`
    - `stringr`
- Python version: 3.7.4
    - `collections`
    - `copy`
    - `group_lasso`
    - `ipython`
    - `lifelines`
    - `math`
    - `matplotlib`
    - `numpy`
    - `pandas`
    - `pickle`
    - `plotly`
    - `random`
    - `scipy`
    - `seaborn`
    - `sklearn`
    - `statsmodels`
    - `sys`
    - `umap`
    - `warnings`

## Abstract 

While computational methods have made substantial progress in improving the accuracy and throughput of pathology workflows for diagnosis, prognosis, and the prediction of genomic features, a lack of interpretability remains a significant barrier to clinical integration. In this study, we present a novel approach for prediction of clinically-relevant molecular phenotypes from histopathology whole-slide images (WSIs) using human-interpretable image features (HIFs). Our method leverages >1.6 million annotations from board-certified pathologists across >5,700 WSIs to train deep learning models for high-resolution tissue classification and cell detection across entire WSIs in five cancer types. Combining cell- and tissue-type models enables computation of 607 HIFs that comprehensively capture specific and biologically-relevant characteristics of multiple tumors. We demonstrate that these correlate with well-known markers of the tumor microenvironment (TME) and can predict diverse molecular signatures, including immune checkpoint protein expression and homologous recombination deficiency (HRD). Our HIF-based approach therefore provides a novel, quantitative, and interpretable window into the composition and spatial architecture of the TME.

