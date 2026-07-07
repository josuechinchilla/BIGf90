# BIGf90 0.4.0

* New function write_par() to create the raw RENUMF90 parameter (.par) file from R instead of writing it by hand
* Effects are given as a named list: the name is the effect label (written as EFFECT #name), with col for the data-file column, type (cross/cov) and class (numer/alpha); pos is accepted as an alias for col
* Optional per-effect comment overrides the #... label
* random_effect marks the animal effect and adds the pedigree FILE / FILE_POS block; optional snp_file adds genotypes for genomic models
* Reproduces the sealice test .par files (pedigree and genomic)

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
