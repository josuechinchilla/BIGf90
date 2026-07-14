<!-- badges: start -->
[![Development Status](https://img.shields.io/badge/status-active%20development-yellow)](https://github.com/Breeding-Insight/BIGf90)
[![R-CMD-check](https://github.com/Breeding-Insight/BIGf90/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Breeding-Insight/BIGf90/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%3E%3D%204.3-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![GitHub issues](https://img.shields.io/github/issues/Breeding-Insight/BIGf90)](https://github.com/Breeding-Insight/BIGf90/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/Breeding-Insight/BIGf90)](https://github.com/Breeding-Insight/BIGf90/pulls)
[![GitHub Release](https://img.shields.io/github/v/release/Breeding-Insight/BIGf90?include_prereleases)](https://github.com/Breeding-Insight/BIGf90/releases/latest)
<!-- badges: end -->

<!-- Add a hex logo here, e.g.:
<div align="center">
  <img width="250" height="250" alt="BIGf90 logo" src="LOGO_URL" />
</div>
-->

# BIGf90
### (B)reeding (I)nsight (G)enomics f90 — an R front end to the BLUPf90 programs

BIGf90 is an R package developed by [Breeding Insight](https://breedinginsight.org/) that wraps the [BLUPf90 family of programs](http://nce.ads.uga.edu/wiki/doku.php) for genetic evaluation. It provides R functions to write RENUMF90 parameter files, run the BLUPf90 modules, estimate variance components, compute (genomic) breeding values, and run K-fold cross-validation — for both pedigree and genomic (GBLUP/ssGBLUP) models.

It provides functions to:

* Write RENUMF90 parameter (`.par`) files with `write_par()`.
* Build a BLUPf90 genotype (`.geno`) file from a PLINK `.ped` with `write_geno()`.
* Process parameter files with RENUMF90 (`run_renum()`).
* Calculate EBVs through BLUP, GBLUP, and ssGBLUP with `blupf90+` (`run_blup()`, `clean_ebvs()`).
* Calculate adjusted phenotypes with `predictf90` (`run_predict()`).
* Estimate variance components with `gibbsf90+` / `postgibbsf90` (`run_gibbs()`, `run_postgibbs()`).
* Run K-fold cross-validation of predictions with `bf90_cv()`.

### Installation

Install the development version from GitHub with `remotes`:

```R
install.packages("remotes")
remotes::install_github("Breeding-Insight/BIGf90", dependencies = TRUE)
library(BIGf90)
```

You also need the BLUPf90 executables (`renumf90`, `blupf90+`, `predictf90`, `gibbsf90+`, `postgibbsf90`) in a folder you point the functions to — download them from the [BLUPF90 wiki](http://nce.ads.uga.edu/wiki/doku.php).

##### Note: BIGf90 is currently in development. Please report any bugs or issues on the GitHub Issues page. On Windows, run RStudio as administrator to avoid BLUPf90 file-permission issues.

### Getting started

See the tutorial vignette (`vignettes/BIGf90.Rmd`) for a full worked pipeline: parameter file → variance components → breeding values → cross-validation, with and without genotypes.

### Funding

BIGf90 development is supported by Breeding Insight, a USDA-funded initiative based at the University of Florida - IFAS.

## Citation

If you use BIGf90 in your research, please cite both the package and the BLUPf90 programs:

**BIGf90:** Chinchilla-Vargas J, Taniguti C, Sandercock A, Breeding Insight Team (2024). BIGf90: Breeding Insight Genomics R front face to blupf90 modules. R package version 0.5.0, https://github.com/Breeding-Insight/BIGf90

**BLUPF90:** Misztal, I., S. Tsuruta, D.A.L. Lourenco, I. Aguilar, A. Legarra, and Z. Vitezica. 2014. Manual for BLUPF90 family of programs. http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all2.pdf

*The BLUPF90 programs are free for research, but their use should be acknowledged in publications. For commercial use, please contact Ignacy Misztal (ignacy@uga.edu) or Daniela Lourenco (danilino@uga.edu). Provide your name, company, and purpose of use in the email.*
