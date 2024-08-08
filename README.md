# BIGf90  
## (B)reeding (I)nsight (G)enomics f90
  
BIGf90 is a wrapper package for BLUPf90 family of programs, please use the following citations:

#### BIGf90 Reference:
  Chinchilla-Vargas J, Taniguti C, Sandercock A, Breeding Insight Team (2024). BIGf90: Breeding Insight Genomics R front face to blupf90 modules. R package version 0.3.2, https://github.com/josuechinchilla/BIGf90
  
#### BLUPF90 Reference:
  Misztal, I., S. Tsuruta, D.A.L. Lourenco, I. Aguilar, A. Legarra, and Z. Vitezica. 2014. Manual for BLUPF90 family of programs: http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all2.pdf")
  
### Overview
This package has functions to:
* Process parameter files with Renumf90.
* Calculate EBVs through BLUP, GBLUP,ssGBLUP with Blupf90+
* Calculate adjusted phenotypes using Predictf90
* Perform cross-validation analyses of predictions produced with Blupf90+.
* Estimate variance components and run other analyses that use Gibbsf90+ and Postgibbsf90.

Please note that as of version 0.3.0 functions are written to work on both unix and windows enviroments.
When running on windows always run RStudio as administrator to avoid issues with BIGf90 functions.

### To install package:  
install.packages("devtools") #If not already installed   
library(devtools)  
devtools::install_github("josuechinchilla/BIGf90")  
library("BIGf90")  
