# Differential expression of single-cell RNA-seq data using Tweedie models

This repository contains R codes to reproduce the real data analyses described in the [Tweedie single-cell paper](https://www.biorxiv.org/content/10.1101/2021.03.28.437378v3). Check out the [associated software (Tweedieverse)](https://github.com/himelmallick/Tweedieverse) for a collection of hands-on tutorials on how to use Tweedieverse with various omics data types. Check out the related [BenchmarkSingleCell](https://github.com/himelmallick/BenchmarkSingleCell) repository that contains the codes to carry out the comprehensive simulation study described in the paper.

# Directory Guide

The structure of the repository is as follows:

* `Codes/`: Preprocessing and analysis scripts to carry out the real data analyses described in the [Tweedie single-cell paper](https://www.biorxiv.org/content/10.1101/2021.03.28.437378v3).
* `Data/`: Gene expression data matrices, with the exception of PBMC and Klein datasets that exceed GitHub’s file size limit but can be downloaded from [here](https://www.dropbox.com/s/63gnlw45jf7cje8/pbmc3k_final.rds?dl=1) and [here](https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/14634416/klein_indrops_control_GSM1599501.h5ad). 
* `Library/`: Helper functions.
* `Output/`: Output files generated by executing the analysis scripts in `Codes/`.

**Each script must specify a working directory which is generally the parent directory containing the above folders.**

# Analysis Guide

The source for most analyses can be found in the following scripts. These scripts call individuals data files from the `Data/` folder (the PBMC and Klein datasets must be downloaded in the `Data/` folder before executing the scripts) and are generally written to the `Output/` folder.

* 1: `Tweedie_index_negative_control`: Script for Tweedie index parameter estimation on negative control datasets as described in **Section 1**.
* 2: `Tweedie_index_positive_control`: Script for Tweedie index parameter estimation on positive control datasets as described in **Section 1**.
* 3: `analysis_DE.R`: Script for DE analyses using six methods (DESeq2, edgeR, CPLM, ZICP, scREHurdle, and MAST) on four datasets (Brain, Kidney, PBMC, Petropoulos) as described in **Section 3.3**.
* 4: `analysis_DEmock.R`: Script for mock DE analyses using six methods (DESeq2, edgeR, CPLM, ZICP, scREHurdle, and MAST) on four datasets (Kidney, PBMC, Klein, Svensson) as described in **Section 3.4**.

# Citation

To cite **`Tweedieverse`** in publications, please use:

Mallick H et al. (2021). [Differential Expression of Single-cell RNA-seq Data using Tweedie Models](https://www.biorxiv.org/content/10.1101/2021.03.28.437378v1). bioRxiv, <https://doi.org/10.1101/2021.03.28.437378>.

To cite the **`Tweedieverse`** software, please use:

Mallick H et al. (2021). [Tweedieverse - A Unified Statistical Framework for Differential Analysis of Multi-omics Data](https://github.com/himelmallick/Tweedieverse). R package, <https://github.com/himelmallick/Tweedieverse>.
