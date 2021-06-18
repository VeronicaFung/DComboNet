# DComboNet

## Introduction

` DComboNet` is an R package for personalized anti-cancer drug combination prediction based on multi-level data integration. There are two main prediction models contained in the package. The level-one model is for generalized anti-cancer drug combination effectiveness prediction and level-two model is for cancer sample specific drug combination prediction. The two-level model based on a network-based method which integrates five subnetwork including drug-drug, drug-gene, drug-pathway, gene-gene and pathway-pathway association networks. Random walk with restart(RWR) algorithm is used to capture global proximity between drugs and the result is based on the rank returned via RWR algorithm. ` DComboNet` also provide clues for the potential mechanisms of drug combinations by extracting the top ranked genes/drugs between predicted drug combinations.

This tutorial provides the instruction of the main usage of this package fitting different scenario. This tutorial will lead to know the basic usage of the prediction, the description of prepared data how to extend the network construction to fit your own dataset. This tutorial will not present detail description of all functions contain in the packagem but you can easily learn those in help document with the R package. 

## Package installation

` DComboNet` has been upload in Github and can be install in R as follow:

```{r, eval=FALSE}
  install.packages("devtools")
  devtools::install_github("VeronicaFung/DComboNet")

  library(DComboNet)
```

Due to the large size of drug-induced gene expression profiles for Level 2 model, we upload only cell line data for repeating the results in our manuscripts, for more cell lines, we prepare an extra data folder in [links]<links>. 

After install the package, you should download this compressed data set, unzip it and put in a reachable path. This path should be used as input for parameter ` load_dir` for functions in ` DComboNet`. For more details about using level-two model with provided cancer sample specific data or your own data table, you can check the tutorial below.

## Tutorial  

After installing package, you can use sample codes (describe below) for  testing.
The detail usages of ` DComboNet`, can be found here: [Instruction](https://veronicafung.github.io/DComboNet/DComboNet-vignette.html) or DComboNet-vignette.pdf if the website preview failed.

## Test codes

Under  folder ` ./test/`, we provide testing codes, which include:
  1. `test_L1_LOOCV.R`: Codes used for Leave-One-Combination-Out Cross Validation (LOOCOCV) for level-one model. 

  2. `test_L2_LOOCV.R`: Codes used for Leave-One-Combination-Out Cross Validation (LOOCOCV) for level-two models. In `test_L2_LOOCV.R`, Level-one model is also provided on the task of predicting drug combination in the context of given cancer-cell line. Here, we provide hepatocellular carcinoma cell line HepG2 (which we show in the main result of the manuscript) and Breast cancer cell line MCF7.
    
  3. `test_OCILY3.R`: Codes to conduct prediction task on OCI-LY3 cell line (from DREAM challenge 7). This code can be used to reproduce the independent validation result in the manuscript.
   
  4. folder ` ./test/mlr_codes/`: machine learning based models we construct for model comparison with DComboNet level-one model. Note that to run these model, you should install Rpackage ` mlr`, ` mlrMBO` and machine learning model packages which listed in supplementary methods in our manuscript.
    
  5. ` p1_features_test.R` and ` plot.R`: main codes used to generate the figures, may need to modify the path to fetch necessary data files.

*Note: for runing the test scripts, model performance evaluation R packages (e.g. ` ROCR`) are also require, please review the codes before running test.*

## Data availability
All data files used in level-one and level-two models (including data for HEPG2 specific model and OCI-LY3 specific model) are incorporated in the packages via this link: https://github.com/VeronicaFung/DComboNet/. 

Data under `/data/` folder includes:
  1. drug and drug pair lists:
      druglistlib.RData: 
      drugpair_lib.RData: 

  2. data files used to calculate pharmacological similarities:
      atc_net.RData: 
      drug_atc_lib.RData: 
      fingerprint_lib.RData: 
      drug_sideeffect.RData: 

  3. data files used to construct drug-gene, drug-pathway, gene-gene and pathway-pathway interaction networks:
      cancer_gene_L1.RData: 
      druggene_lib.RData: 
      drugtarget_lib.RData: 
      inBio_PPI_net.RData: 
      pathway_pathway_interaction.RData: 
  4. extra data files used in level-two data:
      cellline_name.RData: names of cell lines which is integrated in the package.

Data under `/inst/extdata/` folder includes:
   `/drug_net/`: pharmacological similarties and integrated pharmacological score for Level one model, level-two model on HepG2, OCI-LY3 and MCF7 cell lines.
   `/LINCS_data/`: Differentially expressed genes of HepG2 (under the folder `GSE92742`), OCI-LY3 (under the folder `OCI-LY3`) and MCF7  (under the folder `GSE70138`) cell lines.

`OCI-LY3.zip` provided other necessary data for testing model on OCI-LY3 cancer cell line.
`results_figures_manuscripts.zip` contains results and figures files used in manuscript.

## Environment

` DComboNet` is developed in R version 3.5.1. For running the test code in Instruction, you can check the "sessionInfo" chapter for more detailed session information.

## Citation

We have been submit our manuscript to bioaxiv, if you use `DComboNet` in your publication(s), please cite:
Feng F, Zhang Z, Ding G, et al. Personalized anti-cancer drug combination prediction by an Integrated Multi-level Network[J]. bioRxiv, 2020.


