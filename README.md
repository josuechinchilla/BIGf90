# BIGf90  
## (B)reeding (I)nsight (G)enomics f90

<!-- badges: start -->
[![Development Status](https://img.shields.io/badge/status-active%20development-yellow)](https://github.com/Breeding-Insight/BIGf90)
[![R-CMD-check](https://github.com/Breeding-Insight/BIGf90/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Breeding-Insight/BIGf90/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%3E%3D%204.3-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![GitHub issues](https://img.shields.io/github/issues/Breeding-Insight/BIGf90)](https://github.com/Breeding-Insight/BIGf90/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/Breeding-Insight/BIGf90)](https://github.com/Breeding-Insight/BIGf90/pulls)
[![GitHub Release](https://img.shields.io/github/v/release/Breeding-Insight/BIGf90?include_prereleases)](https://github.com/Breeding-Insight/BIGf90/releases/latest)
<!-- badges: end -->
  
BIGf90 is a wrapper package for BLUPf90 family of programs, please use the following citations:

#### BIGf90 Reference:
  Chinchilla-Vargas J, Taniguti C, Sandercock A, Breeding Insight Team (2024). BIGf90: Breeding Insight Genomics R front face to blupf90 modules. R package version 0.4.0, https://github.com/Breeding-Insight/BIGf90
  
#### BLUPF90 Reference:
  Misztal, I., S. Tsuruta, D.A.L. Lourenco, I. Aguilar, A. Legarra, and Z. Vitezica. 2014. Manual for BLUPF90 family of programs: http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all2.pdf")

*The BLUPF90 programs are free for research, but their use should be acknowledged in publications. For commercial use,
please contact Ignacy Misztal (ignacy@uga.edu) or Daniela Lourenco (danilino@uga.edu). Provide your name, company, and purpose of use in the email.*   
  
### Overview
This package has functions to:
* Write RENUMF90 parameter (.par) files with write_par().
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
devtools::install_github("Breeding-Insight/BIGf90")  
library("BIGf90")  
