# BIGf90 1.0.0

* First stable release, consolidating the write_par()/write_geno() builders, the individual module wrappers (run_pregs, run_postgs, run_predf, ...), the GWAS and cross-validation workflows (run_gwas, run_cva), heritability (calc_h2) and the GWAS plots
* New function manhattan_plot() to draw a Manhattan plot straight from the postGSf90 outputs (snp_sol or windows_variance) produced by run_postgs()/run_gwas(): plots the proportion of genetic variance explained per SNP or per non-overlapping window (the ssGBLUP/wssGBLUP standard) or -log10(p) when p-values are available; alternates colours by chromosome, with point highlighting and PNG/PDF output
* manhattan_plot() threshold lines: a numeric cutoff (e.g. 0.02) for the variance plot, or, for the p-value plot, any of "bonferroni" (-log10(alpha/m)), "fdr" (Benjamini-Hochberg at level alpha) and "meff" (Bonferroni on a supplied effective number of tests) drawn as labelled lines. Multiple-testing lines are refused on the variance plot (variance explained is not a test statistic)
* manhattan_plot() and qq_plot() gain gc_correct: single-parameter genomic-control correction that divides every chi-square by the genome-wide inflation factor lambda before recomputing p (Devlin and Roeder 1999); lambda is reported either way
* New function qq_plot() to draw a QQ plot of the exact-GWAS p-values (a^2/var(a) ~ chi-square, Aguilar et al. 2019), overall or one panel per chromosome (per_chromosome = TRUE), reporting the genomic inflation factor lambda (median chi-square / 0.4549) for each
* Both GWAS plotters accept a chromosomes = argument to subset chromosomes (e.g. 1:28 to drop a rank-coded / unmapped chromosome)
* New function mcmc_diagnostics() for post-Gibbs MCMC convergence diagnostics via the coda package (Suggests): reads gibbsf90+ gibbs_samples (default, robust to a failed postgibbsf90) or postgibbsf90 postgibbs_samples, builds a coda mcmc object of the (co)variance components plus derived heritability, and reports autocorrelation, effective sample size, and the Geweke, Raftery-Lewis and Heidelberger-Welch diagnostics, writing a multi-page PDF (trace/density, autocorrelation, Geweke, cross-correlation, normal-QQ) and a text summary. Based on the CODA post-Gibbs workflow of Vallejo et al.
* mcmc_diagnostics() now also returns a decision, not just plots: a per-parameter table (posterior mean/SD/median, 95% HPD interval, ESS, Geweke z, Raftery-Lewis dependence factor, and a converged/borderline/failed status), an overall converged flag (with a warning when any parameter fails), suggested run_gibbs() re-run settings (gibbs_iter/gibbs_burn/gibbs_keep from the Raftery-Lewis diagnostic), and the heritability posterior with a 95% HPD interval
* Internal tidy-up: the package hooks (globalVariables, startup message) and the shared GWAS/reader helpers now live in a single R/utils.R (previous zzz.R and gwas_utils.R removed); no user-visible change
* run_gwas() gains a snp_pvalue argument (default FALSE): set TRUE to add OPTION snp_p_value to the blupf90+ and postGSf90 steps, so snp_sol carries the SNP-solution variance (column var_a_hat) needed by qq_plot() and by manhattan_plot(statistic = "pvalue")
* run_gwas() now reads snp_sol with header detection: newer postGSf90 writes a header row, which was being read as data - this corrupted the re-estimated SNP weights on the 2nd iteration (a stray "weight" label made blupf90+ read 0 weights) and returned one spurious row. wssGBLUP iterations past the first are now correct
* New function write_geno() to build a BLUPf90 genotype (.geno) file from a PLINK .ped
* Reads the --recode12 layout (two allele codes per locus) and collapses each pair to a 0/1/2 dosage (count of count_allele), writing a missing allele as missing_code (default 5); set alleles_per_locus = 1 for input that is already one dosage per locus
* Writes each line as the ID left-justified to a fixed width then the contiguous dosage string, so every row's genotypes start at the same column (BLUPf90 fixed format); optional .map locus-count check and optional BLUPf90 marker-map output
* Added a testthat suite covering write_geno()
* New function run_pregs() wrapping preGSf90: genotype QC and construction/inversion of the genomic relationship matrix (G, A22, GimA22i) for GBLUP/ssGBLUP
* New function run_postgs() wrapping postGSf90: back-solves SNP effects (ssGWAS), window variances, SNP weights and Manhattan plots after a genomic run_blup()
* run_pregs() and run_postgs() return the paths of the outputs they produced, so higher-level workflow functions can chain the wrappers
* New function run_predf() wrapping predf90: predicts direct genomic values (DGV) for new/young genotyped animals from the SNP effects (snp_pred) produced by run_postgs(); command-line driven, with optional reliabilities via acc = TRUE
* New workflow function run_gwas(): a weighted single-step GBLUP (wssGBLUP) GWAS pipeline that chains run_renum -> run_pregs -> iterated run_blup + run_postgs, re-estimating the SNP weights from the SNP effects each iteration (iterations = 1 is ssGBLUP; 2 reproduces the second-iteration weights of Vallejo et al. 2024). Returns the parsed SNP-solution table and the window-variance path; supports 1-Mb windows (windows_mbp) or SNP-count windows (windows_snp)
* New function calc_h2() to estimate narrow-sense heritability from a Gibbs run: by default from postgibbsf90's posterior means (postmean), or with from_postmean = FALSE straight from the gibbsf90+ samples (gibbs_samples, returning the posterior mean h2 + SD) - the samples path avoids postgibbsf90's diagnostic step, which can crash on some chains. Provides the h2 for run_cva
* run_renum() now also symlinks the SNP map_file into the renumbered directory (previously only the genotype file and its XrefID), so downstream genomic steps find the map
* run_blup() gains a par_file argument (default "renf90.par") so workflows can point it at a custom parameter file
* run_cva() now, for genomic runs, draws CV folds only from genotyped animals (testing on genotyped+phenotyped animals) while all phenotyped animals train and G is built on all genotyped animals - matching the reference complete-data CVA; genotyped IDs are read from the renumf90 <geno>_XrefID file
* run_cva() gains a genotyped_only argument (default FALSE): set TRUE to restrict the whole analysis to the genotype+phenotype intersection - phenotyped-only records are dropped from training and genotyped-only animals are removed from the genotype file so G uses the intersection. Ignored for pedigree-only runs Genotyped IDs are read from the renumf90 <geno>_XrefID file; the argument is ignored for pedigree-only runs
* run_cva() now reports, for genomic runs, how many phenotype records are genotyped (used for testing) versus phenotyped-only and genotyped-only (not used) - a warning when phenotyped animals lack a genotype, otherwise an informational message


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
* Output directories need to exist. The function will stop if they don't
* Automatic installation checks included
* NEWS file included
* GitHub example
