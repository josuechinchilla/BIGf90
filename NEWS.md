# BIGf90 0.5.0

* New function write_geno() to build a BLUPf90 genotype (.geno) file from a PLINK .ped
* Reads the --recode12 layout (two allele codes per locus) and collapses each pair to a 0/1/2 dosage (count of count_allele), writing a missing allele as missing_code (default 5); set alleles_per_locus = 1 for input that is already one dosage per locus
* Writes the renumf90-ready contiguous format (ID left-justified, dosages column-aligned), with an optional .map locus-count check and an optional BLUPf90 marker-map output
* Added a testthat suite covering write_geno()


# BIGf90 0.4.0

* New function write_par() to build the raw RENUMF90 parameter (.par) file from R instead of writing it by hand
* Effects are given as a named list keyed by a label; each effect gives its data-file column (col, or pos as an alias), type (cross/cov) and class (numer/alpha); an optional comment overrides the EFFECT #... label
* random marks the random effect by its column (or label); the random effect is written last automatically so RANDOM follows it, with the pedigree FILE / FILE_POS block after
* Supports genomic models (snp_file), permanent-environment / repeatability terms (optional = "pe") and (CO)VARIANCES starting values
* Multi-trait models: give a trait vector with matrix RESIDUAL_VARIANCE / (CO)VARIANCES, and effect columns repeat per trait
* Reproduces the sealice test .par files (pedigree and genomic)
* Added a testthat suite covering write_par() and the create_folds() helper


# BIGf90 0.3.1

* Users now can have executable files, input files, and output files in different locations
* Warnings were added in case of run_gibbs and run_postgibbs functions that requires inputs to be together with renum results
* Code changed to always use global path even if it is referring to working directory
* Verbose argument added to give option to print or not information on screen while running
* Default values added
* Output directories need to exist. The function will stop if they don´t
* Automatic installation checks included
* NEWS file included
* GitHub example
