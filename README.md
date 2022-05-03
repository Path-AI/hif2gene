# hif2gene: Dense, high-resolution mapping of cells and tissues from cancer pathology images for the interpretable prediction of molecular phenotypes

## Overview
> `hif2gene` is the primary code repository for reproducing analyses in Diao, Chui, and Wang et al. 2020: "Dense, high-resolution mapping of cells and tissues from pathology images for the interpretable prediction of molecular phenotypes in cancer"

You can read the full publication in Nature: [Human-interpretable image features derived from densely mapped cancer pathology slides predict diverse molecular phenotypes](https://www.nature.com/articles/s41467-021-21896-9)

## Installation
- Clone this repo to your local machine using https://github.com/Path-AI/hif2gene.git
- Estimated install time: <1 minute
- Estimated run time: varies between scripts (<1 minute to ~1 hour)

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

Computational methods have made substantial progress in improving the accuracy and throughput of pathology workflows for diagnostic, prognostic, and genomic prediction. Still, lack of interpretability remains a significant barrier to clinical integration. We present a novel approach for predicting clinically-relevant molecular phenotypes from whole-slide histopathology images using human-interpretable image features (HIFs). Our method leverages >1.6 million annotations from board-certified pathologists across >5,700 samples to train deep learning models for cell and tissue classification that can exhaustively map whole-slide images at two and four micron-resolution. Cell- and tissue-type model outputs are combined into 607 HIFs that quantify specific and biologically-relevant characteristics across five cancer types. We demonstrate that these HIFs correlate with well-known markers of the tumor microenvironment and can predict diverse molecular signatures (AUROC 0.601-0.864), including expression of four immune checkpoint proteins and homologous recombination deficiency, with performance comparable to ‘black-box’ methods. Our HIF-based approach provides a comprehensive, quantitative, and interpretable window into the composition and spatial architecture of the tumor microenvironment.

## License

Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
