# Bi-EB
![Bi_EB_example](https://user-images.githubusercontent.com/53017373/126046429-469fb8d6-1504-42d8-8dc9-fd3451db268c.png)
#### Bicluster for (a) breast cancer luminal subtype (b) breast cancer basal-like subtype. Red color shows higher probability and green shows lower probability of belonging to bicluster. 
#### Clinical evidence: these important drug targets in breast cancer, ESR1, PGR, HER2, EGFR and AR have a high similarity in mRNA and protein variation in both tumors and cell lines [1][2]. GATA3 and RP56KB1 are two promising drug targets for breast cancer [1][2].
[1] Jiang GL, Zhang SJ, Yazdanparast A, Li M, Vikram Pawar A, Liu YL, Inavolu SM, Cheng LJ. Comprehensive comparison of molecular portraits between cell lines and tumors in breast cancer. BMC genomics, 2016, 17(7), 281-301. <p>
[2] Yazdanparast A, Li L, Radovich M, Cheng LJ. Signal translational efficiency between mRNA expression and antibody-based protein expression for breast cancer and its subtypes from cell lines to tissue. International Journal of Computational Biology and Drug Design , 2018, 11 (1-2), 67-89.


#### Example data
![Bi_EB_example_data](https://github.com/lijcheng12/Bi-EB/blob/main/Example%20data%20for%20Bi-EB.xlsx) <p>
![Bi_EB_synthetic_data](https://github.com/lijcheng12/Bi-EB/blob/main/synthetic_data.xlsx) <p>
##### Supplementary data is for systematic patterns of the (gene expression/Protein amount) ratio absed on the Cancer Genomics Atlas (TCGA) and the Cancer Cell Line Encyclopedia (CCLE) breast cancer data 
  
#### Example code
![Bi_EB_example_code](https://github.com/lijcheng12/Bi-EB/blob/main/Bi-EB_Example.R)

## Introduction
The novelty Bi-EB algorithm can search the coherent and flexible co-regulation patterns across multi-omics data both in tumors and cancer cells. Transparent probabilistic interpretation and ratio strategy for omics data is first time proposed to detect the co-regulation patterns of drug targets and identify their associated molecular functions. 

## Author summary
The genome molecular features shared between cell lines and tumors give us insight into discovering potential drug targets for cancer patients. Our previous studies demonstrate that these important drug targets in breast cancer, ESR1, PGR, HER2, EGFR, and AR have a high similarity in mRNA and protein variation in both tumors and cell lines [1-2]. Based on previous studies we made specific hypothesis that there exist translational gene sets that are characterized by highly correlated molecular profiles among RNA, and proteins. There are translational gene sets that are shared between tumor tissues and cancer cell lines. These gene sets show similar pattern in a subgroup of cell line and tissue samples. In this study, we aim to integrate cell line and tissue RNA and protein profiles to characterize drug-able target expression alterations across both RNA and protein data by using bi-clustering method. Here we developed a biclustering method based on empirical bayesian (Bi-EB), to detect the local pattern of integrated omics data both in cancer cells and tumors. We adopt a data driven statistics strategy by using Expected-Maximum (EM) algorithm to extract the foreground bicluster pattern from its background noise data in an iterative search. Our novel Bi-EB statistical model has better chance to detect co-current patterns of gene and protein expression variation than the existing biclustering algorithms and seek the drug targets’ co-regulated modules.

## API Link
### https://github.com/lijcheng12/Bi-EB/README.

## Features
-Muliple omics data integration <p>
-Doing biclusters among multiple omics data and specises.<p>
-Easy to use. We provide an example how to use Bi-EB, including example data and code.<p>
-systematic patterns of the (gene expression/Protein amount) ratio is found 

## Contact with

lijun.cheng@osumc.edu
