# BIGf90  
#### (B)reeding (I)nsight (G)enomics f90
  
An R front face to run blupf90 modules.
#### Overview
This package has functions to:
* Process parameter files with Renumf90.
* Calculate EBVs through BLUP, GBLUP,ssGBLUP with Blupf90+
* Calculate adjusted phenotypes using Predictf90
* Perform cross-validation analyses of predictions produced with Blupf90+.
* Estimate variance components and run other analyses that use Gibbsf90+ and Postgibbsf90.

Please note that as of version 0.3.0 functions are written to work on both unix and windows enviroments.
When running on windows always run RStudio as administrator to avoid issues with BIGf90 functions.

#### To install package:  
install.packages("devtools") #If not already installed   
library(devtools)  
devtools::install_github("josuechinchilla/BIGf90")  
library("BIGf90")  
