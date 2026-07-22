#' Run K-fold cross-validation analysis (CVA)
#'
#' This function runs a K-fold cross-validation analysis (CVA) using blupf90 modules.
#'
#' This function sets up and runs a K-fold cross-validation analysis (CVA) using blupf90+ and predictf90.
#' The function run_renumf90 needs to be used beforehand to process a .par file created by the user.
#' Using the phenotype corrected for all fixed effects (y*), the function reports predictive ability as the correlation between y* and the ebvs, accuracy as that correlation divided by the square-root of the narrow-sense heritability, and bias as the regression of y* on the ebvs.
#'
#' @param missing_value_code code used in the .par file after OPTION MISSING to indicate missing phenotype, if this option is no use, this value must be 0.
#' @param random_effect_col Column where random effect is located, found under RANDOM_GROUP in the renf90.par file.
#' @param h2 estimate of narrow-sense heritabilty.This value is use to calculate accuracy of ebvs
#' @param num_runs Number of independent cross-validation runs to be performed.
#' @param num_folds Number of folds to be generated within each independent run.
#' @param genotyped_only logical (default FALSE); only relevant for genomic runs (a SNP_file in the parameter file). When FALSE (the default, matching the reference complete-data CVA), all phenotyped animals are used for TRAINING and G is built on all genotyped animals, but only genotyped animals are TESTED (placed in the CV folds). When TRUE, the analysis is restricted to animals that have BOTH a genotype and a phenotype: phenotyped-only records are dropped from training and genotyped-only animals are removed from the genotype file, so G is built on the intersection. Ignored for pedigree-only runs.
#' @param output_table_name Name of the final tab-separated out-up file. This field should be in quotes "".
#' @param path_2_execs path to a folder that holds all blupf90 executables that ill be used (blupf90+,predictf90). This field should be in quotes "".
#' @param input_files_dir directory containing files renf90.par renf90.fields renf90.inb renf90.tables renf90.dat
#' @param output_files_dir path to directory to store the output files
#' @param seed set seed for the stochastic process
#' @param verbose logical defining if information will be printed on the console
#'
#' @return a tab-separated file that includes accuracy and bias estimates of ebvs.
#' @import dplyr
#' @examples
#' ## Example for a CVA with 5 independent runs dividing the data in 10 folds.
#'
#'
#' # run_cva(path_2_execs = "/Users/johndoe/Desktop/bf90_execs/",
#' #      missing_value_code = -999,
#' #      random_effect_col= 3,
#' #      h2 = 0.5,
#' #      num_runs = 5,
#' #      num_folds = 10,
#' #      output_table_name = "example_run",
#' #      input_files_dir = ".",
#' #      output_files_dir = ".")
#'
#'
#' @export
run_cva <- function(missing_value_code = NULL,
                    random_effect_col = NULL,
                    h2 = NULL,
                    num_runs = NULL,
                    num_folds = NULL,
                    path_2_execs = ".",
                    input_files_dir = ".",
                    output_files_dir = ".",
                    output_table_name = NULL,
                    genotyped_only = FALSE,
                    seed = 101919,
                    verbose = TRUE) {
  
  input_files_dir <- normalizePath(input_files_dir)
  
  # Checks
  output_files_dir <- normalizePath(output_files_dir)
  if(file.exists(output_files_dir)){
    check_files <- list.files(output_files_dir)
    if(length(check_files) > 0) warning(paste("Directory", output_files_dir, "is not empty. Some files may be replaced."))
  } else {
    stop(paste("Directory", output_files_dir, "does not exist. Create it before running the function."))
  }
  
  if(is.null(missing_value_code)) stop("Specify missing_value_code")
  if(is.null(random_effect_col)) stop("Specify random_effect_col")
  if(is.null(h2)) stop("Specify h2")
  if(is.null(num_runs)) stop("Specify num_runs")
  if(is.null(num_folds)) stop("Specify num_folds")
  
  if (!file.exists(file.path(input_files_dir,"renf90.par"))) stop("File 'renf90.par' not found at: ", input_files_dir)
  if (!file.exists(file.path(input_files_dir,"renf90.fields"))) stop("File 'renf90.fields' not found at: ", input_files_dir)
  if (!file.exists(file.path(input_files_dir,"renf90.tables"))) stop("File 'renf90.tables' not found at: ", input_files_dir)
  if (!file.exists(file.path(input_files_dir,"renf90.dat"))) stop("File 'renf90.dat' not found at: ", input_files_dir)
  
  # Check files in the parameter files
  parfile <- readLines(file.path(input_files_dir,"renf90.par"))
  
  snpfile <- grep("SNP_file", parfile)
  if(length(snpfile) != 0) {
    snp_file_name <- sapply(strsplit(parfile[snpfile], " "), function(x) x[length(x)])
    if(!file.exists(file.path(input_files_dir, snp_file_name))) stop(paste("File", file.path(input_files_dir, snp_file_name),
                                                                           "especified in parameter file",file.path(input_files_dir,"renf90.par"),
                                                                           "line", snpfile, "does not exist."))
    snp_file_name <- normalizePath(file.path(input_files_dir, snp_file_name))
  } else snp_file_name <- NULL

  # cleaned-SNP and saved G-inverse filenames (built once on the full data, reused per fold)
  if(!is.null(snp_file_name)){
    snp_base  <- base::basename(snp_file_name)      # e.g. sealice_match.geno
    clean_snp <- base::paste0(snp_base, "_clean")   # produced by saveCleanSNPs
  }
  
  pedfile <- grep(" FILE", parfile) + 1
  if(length(pedfile) != 0) {
    renf90_ped_name <- gsub(" ", "", parfile[pedfile])
    if(!file.exists(file.path(input_files_dir, renf90_ped_name))) stop(paste("File", file.path(input_files_dir, renf90_ped_name),
                                                                             "especified in parameter file",file.path(input_files_dir,"renf90.par"),
                                                                             "line", pedfile, "does not exist."))
    renf90_ped_name <- normalizePath(file.path(input_files_dir, renf90_ped_name))
  } else renf90_ped_name <- NULL
  
  if(is.null(output_table_name)) stop("Define output table name.")
  
  # set seed for reproducibility
  base::set.seed(seed)
  
  path_2_execs <- normalizePath(path_2_execs)
  
  # Inform parameters and directories set
  if(verbose){
    cat("Parameters set:\n",
        "  missing_value_code = ", missing_value_code,"\n",
        "  random_effect_col = ", random_effect_col,"\n",
        "  h2 = ",h2 ,"\n",
        "  num_runs = ", num_runs,"\n",
        "  num_folds =", num_folds,"\n",
        "Directories:\n",
        "  input files:", input_files_dir, "\n",
        "  output files:", output_files_dir, "\n",
        "  executable files:", path_2_execs)
  }
  
  # define working directory
  wd_path <- base::getwd()
  on.exit(setwd(wd_path), add = TRUE)   # restore the working directory even if the function errors
  
  #Assign .exes or not based on OS
  if (.Platform$OS.type == "unix") {
    predict = "predictf90"
    blup = "blupf90+"
  } else if (.Platform$OS.type == "windows") {
    predict = "predictf90.exe"
    blup = "blupf90+.exe"
  }
  
  # Run BF90 programs for the whole dataset
  setwd(input_files_dir)

  # --- Genotyped-animal handling (genomic runs) -------------------------------
  # Default: only genotyped animals are TESTED (put in folds); all phenotyped
  # animals stay in TRAINING. With genotyped_only = TRUE only genotyped+phenotyped
  # animals are kept: phenotyped-only records and genotyped-only individuals (in the
  # genotype file) are both removed, so training, testing and G use the intersection.
  swap_datafile <- function(par_lines, new_dat) {          # point DATAFILE at a filtered data file
    i <- base::grep("^DATAFILE", par_lines)[1]
    if(!base::is.na(i) && i < base::length(par_lines)) par_lines[i + 1] <- new_dat
    par_lines
  }
  geno_ids       <- NULL
  data_file_name <- "renf90.dat"
  cv_snp         <- if(!is.null(snp_file_name)) snp_base else NULL   # genotype file to use (default: the original)
  if(!is.null(snp_file_name)){
    xref <- base::paste0(snp_base, "_XrefID")                        # renumf90: renumbered <-> original genotyped IDs
    if(base::file.exists(xref)){
      xref_tab <- utils::read.table(xref, header = FALSE, stringsAsFactors = FALSE)   # V1 renumbered, V2 original
      geno_ids <- base::as.character(xref_tab[[1]])
      dat_tab  <- utils::read.table("renf90.dat", sep = " ", header = FALSE)               # sep=" " matches bf90_phenos (renf90.dat has a leading space); same row order as readLines
      dat_id   <- base::as.character(dat_tab[, -1, drop = FALSE][[random_effect_col + 1]]) # animal id per phenotype record
      is_geno  <- dat_id %in% geno_ids
      n_pheno  <- base::length(dat_id)
      n_used   <- base::sum(is_geno)                                 # genotyped + phenotyped (used for testing)
      n_phenoonly <- n_pheno - n_used                               # phenotyped, no genotype
      n_genoonly  <- base::sum(!(geno_ids %in% base::unique(dat_id)))  # genotyped, no phenotype

      # Report how many records are used vs not used
      info <- base::paste0("run_cva (genomic): of ", n_pheno, " phenotype records, ", n_used,
                           if(genotyped_only) " are genotyped and used for training and testing" else " are genotyped and used for testing",
                           "; ", n_phenoonly, " phenotyped records have no genotype (",
                           if(genotyped_only) "dropped" else "kept for training, not tested",
                           "); ", n_genoonly, " genotyped animals have no phenotype (",
                           if(genotyped_only) "removed from the genotypes" else "used in G only", ").")
      if(n_phenoonly > 0 || (genotyped_only && n_genoonly > 0)) base::warning(info) else if(verbose) base::message(info)

      if(genotyped_only){
        pheno_renum <- base::unique(dat_id[is_geno])                 # renumbered ids of genotyped + phenotyped animals
        base::writeLines(base::readLines("renf90.dat")[is_geno], "renf90_geno.dat")   # keep only their phenotype records
        data_file_name <- "renf90_geno.dat"
        if(n_genoonly > 0){                                          # drop genotyped-but-not-phenotyped from the genotypes
          kx        <- base::as.character(xref_tab[[1]]) %in% pheno_renum
          keep_orig <- base::as.character(xref_tab[[2]][kx])
          gl        <- base::readLines(snp_file_name)
          gid       <- base::sub("[[:space:]].*$", "", base::trimws(gl))              # first field = original id
          base::writeLines(gl[gid %in% keep_orig], "cv_geno.geno")
          utils::write.table(xref_tab[kx, , drop = FALSE], "cv_geno.geno_XrefID",
                             row.names = FALSE, col.names = FALSE, quote = FALSE)
          cv_snp <- "cv_geno.geno"
        }
      }
      clean_snp <- base::paste0(cv_snp, "_clean")                    # reflect the genotype file actually used
    } else {
      base::warning("Genotype cross-reference '", xref, "' not found in ", input_files_dir,
                    "; testing on all phenotyped animals. Run run_renum() to create it.")
    }
  } else if(genotyped_only){
    base::message("run_cva: 'genotyped_only = TRUE' ignored - no SNP_file (pedigree-only run).")
  }

  # For genomic analyses, make sure the G-inverse (Gi) and cleaned SNPs are saved on
  # this full-data run so every fold can reuse them (built once here).
  blup_parfile <- "renf90.par"
  if(!is.null(snp_file_name)){
    blup_par <- base::readLines("renf90.par")
    if(data_file_name != "renf90.dat") blup_par <- swap_datafile(blup_par, data_file_name)
    if(!is.null(cv_snp) && cv_snp != snp_base) blup_par <- base::gsub(snp_base, cv_snp, blup_par, fixed = TRUE)  # build G on the intersection genotypes
    if(!any(base::grepl("saveGInverse",  blup_par))) blup_par <- c(blup_par, "OPTION saveGInverse")
    if(!any(base::grepl("saveCleanSNPs", blup_par))) blup_par <- c(blup_par, "OPTION saveCleanSNPs")
    base::writeLines(blup_par, "renf90_blup.par")
    blup_parfile <- "renf90_blup.par"
  }
  output <- execute_command(command = paste0(file.path(path_2_execs, blup)," ", blup_parfile), logfile = "run_blup.log")

  # predictf90 on a copy of renf90.par that adjusts the phenotype for all effects
  # except the random (animal) effect -> yhat = corrected phenotype (y*)
  predict_par <- base::readLines("renf90.par")
  if(data_file_name != "renf90.dat") predict_par <- swap_datafile(predict_par, data_file_name)
  predict_par <- c(predict_par, base::paste("OPTION include_effects", random_effect_col))
  base::writeLines(predict_par, "renf90_predict.par")
  command_predict <- paste0(file.path(path_2_execs, predict), " ", "renf90_predict.par")
  output <- execute_command(command = command_predict, logfile = "run_predict.log")
  
  if(is.null(snp_file_name)){
    files_res <- c("run_blup.log", "bvs.dat", "bvs2.dat", "yhat_residual", "solutions")
  } else {
    files_res <- c("run_blup.log", "bvs.dat", "bvs2.dat", "yhat_residual", "solutions",
                   "freqdata.count", "freqdata.count.after.clean", "Gen_call_rate", "Gen_conflicts",
                   "run_predict.log", "sum2pq",
                   "Gi", clean_snp, base::paste0(clean_snp, "_XrefID"))   # reused by every fold
  }

  for(i in 1:length(files_res)) if(file.exists(files_res[i])) file.rename(from = files_res[i], to = file.path(output_files_dir,files_res[i]))
  
  # Prepare files for each BLUP run
  renf90 <- base::readLines(paste0("renf90.par"))
  #ped_file <- dirname(renf90_ped_name)
  fields_file <- base::readLines(paste0("renf90.fields"))
  tables_file <- base::readLines(paste0("renf90.tables"))

  # Genomic: point each fold at the G-inverse + cleaned SNPs built on the full data,
  # so folds skip genotype QC and G construction (identical results, much faster).
  if(!is.null(snp_file_name)){
    renf90 <- base::gsub("saveGInverse",  "readGInverse",       renf90, fixed = TRUE)
    renf90 <- base::gsub("saveCleanSNPs", "no_quality_control", renf90, fixed = TRUE)
    if(!any(base::grepl("readGInverse",       renf90))) renf90 <- c(renf90, "OPTION readGInverse")
    if(!any(base::grepl("no_quality_control", renf90))) renf90 <- c(renf90, "OPTION no_quality_control")
    renf90 <- base::gsub(snp_base, clean_snp, renf90, fixed = TRUE)   # SNP_file -> cleaned set
    renf90 <- renf90[!base::grepl("map_file", renf90, ignore.case = TRUE)]   # folds run plain BLUP; the SNP map is not needed
    reuse_files <- base::file.path(output_files_dir, c("Gi", clean_snp, base::paste0(clean_snp, "_XrefID")))
  } else reuse_files <- NULL
  
  # Read and preprocess the phenotype data
  bf90_phenos <- utils::read.table(data_file_name, sep = " ", header = FALSE) %>%   # training set (all animals, or genotyped-only when genotyped_only = TRUE)
    dplyr::select(-1)
  
  # Pool of animals eligible to be TESTED (masked): genotyped animals for genomic
  # runs, all animals for pedigree-only runs. Non-genotyped animals are never
  # masked, so they always stay in the training set.
  fold_pool <- bf90_phenos
  if(!is.null(geno_ids))
    fold_pool <- fold_pool[base::as.character(fold_pool[[random_effect_col + 1]]) %in% geno_ids, , drop = FALSE]

  # Shuffle the eligible animals and create folds
  data_shuffled <- base::lapply(1:num_runs, function(x) fold_pool[base::sample(base::nrow(fold_pool)), ] %>%
                                  dplyr::select((random_effect_col + 1)))
  
  folds <- base::lapply(data_shuffled, function(x) create_folds(x, num_folds))
  mutated_data <- base::lapply(folds, function(f) mutate_folds(bf90_phenos, f, num_folds, missing_value_code, random_effect_col + 1))  # id col = random_effect_col+1
  
  setwd(output_files_dir)
  for (run in 1:num_runs) {
    for (fold in 1:num_folds) {
      dir_path <- base::sprintf("run%d/fold%d", run, fold)
      create_cv_datasets(run, 
                         fold, 
                         data_frame =  mutated_data[[run]][[fold]],
                         dir_path, 
                         renf90, 
                         renf90_ped_name,
                         input_files_dir,
                         reuse_files)
    }
  }
  
  # Run BLUPf90+ for each fold and collect EBVs
  ebvs_for_cv_runs <- base::list()
  for (run in 1:num_runs) {
    ebvs_for_cv_list <- base::list()
    for (fold in 1:num_folds) {
      base::setwd(base::file.path(output_files_dir, base::sprintf("run%d/fold%d", run, fold)))
      command <- base::paste0(file.path(path_2_execs, blup), base::sprintf(" renf90_run%d_fold%d.par", run, fold))
      logfile <- base::sprintf("blup_fold%d_run%d.log", fold, run)
      execute_command(command = command, logfile = logfile)
      
      data_file <- base::sprintf("renf90_run%d_fold%d.dat", run, fold)
      masked_ids <- utils::read.table(data_file) %>%
        dplyr::filter(V1 == missing_value_code) %>%
        dplyr::select(random_effect_col + 1) %>%   # id column
        base::unlist()
      
      ebvs_for_cv <- utils::read.table("solutions", header = FALSE, skip = 1) %>%
        dplyr::filter(V2 == random_effect_col & V3 %in% masked_ids) %>%
        dplyr::select(3, 4)
      
      ebvs_for_cv_list[[fold]] <- ebvs_for_cv
      utils::write.table(ebvs_for_cv, file = base::sprintf("ebvs_for_cv_run%d_fold%d.dat", run, fold), row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    }
    ebvs_for_cv_runs[[base::sprintf("ebvs_for_cv_run%d", run)]] <- base::do.call(base::rbind, ebvs_for_cv_list)
  }
  base::setwd(output_files_dir)
  
  # Calculate correlations and bias
  corrected_phenos <- utils::read.table(paste0(output_files_dir, "/yhat_residual"), header = FALSE) %>% dplyr::select(1, 2)
  
  ystar_correlations <- base::numeric(num_runs)
  bias_list <- base::numeric(num_runs)
  accuracy_list <- base::numeric(num_runs)
  
  for (i in 1:num_runs) {
    ystar <- dplyr::inner_join(ebvs_for_cv_runs[[base::sprintf("ebvs_for_cv_run%d", i)]], corrected_phenos, by = c("V3" = "V1"))
    ystar_correlations[i] <- base::round(stats::cor(ystar$V4, ystar$V2), 3)
    
    model <- stats::lm(V2 ~ V4, data = ystar)   # bias (LR b1): corrected phenotype (V2) regressed on EBV (V4)
    bias_list[i] <- base::round(stats::coefficients(model)["V4"], 3)
    
    accuracy_list[i] <- base::round(ystar_correlations[i] / base::sqrt(h2), 3)
  }
  
  # Collecting results for verbose output
  verbose_results <- data.frame(
    Metric = rep(c("predictive ability", "bias", "accuracy"), each = num_runs),
    Run = paste("Run", rep(1:num_runs, times = 3)),
    Value = c(ystar_correlations, bias_list, accuracy_list)
  )
  
  # Calculate averages
  ystar_accuracy <- base::round(base::mean(ystar_correlations), 3)
  y_corrected_accuracy <- base::round(ystar_accuracy / base::sqrt(h2), 3)
  average_bias <- base::round(base::mean(bias_list), 3)
  
  # Create data frame with only the desired averaged metrics
  summary_data <- base::data.frame(
    Metric = c("predictive ability", "bias", "accuracy"),
    Run = "Average",
    Value = c(ystar_accuracy, average_bias, y_corrected_accuracy)
  )
  
  # Open the output file for writing
  output_file <- file(output_table_name, "w")
  
  # Write the parameters to the file
  writeLines("***** Parameters used for CV Analysis *****", output_file)
  writeLines(paste("  missing_value_code =", missing_value_code), output_file)
  writeLines(paste("  random_effect_col =", random_effect_col), output_file)
  writeLines(paste("  h2 =", h2), output_file)
  writeLines(paste("  num_runs =", num_runs), output_file)
  writeLines(paste("  num_folds =", num_folds), output_file)
  
  # Write section title and detailed per-run results
  writeLines("\n***** Predictive Ability, Bias, and Accuracy Per Run *****", output_file)
  utils::write.table(verbose_results, file = output_file, row.names = FALSE, quote = FALSE, append = TRUE)
  
  # Write section title and summary results
  writeLines("\n***** Predictive Ability, Bias, and Accuracy Results (Averages) *****", output_file)
  utils::write.table(summary_data, file = output_file, row.names = FALSE, quote = FALSE, append = TRUE)
  
  # Close the output file
  close(output_file)
  
  # Print the results per run if verbose is TRUE
  if (verbose) {
    base::cat("\n*****Predictive Ability, Bias, and Accuracy Per Run*****\n")
    print(verbose_results, row.names = FALSE)
  }
  
  # Print out summarized results regardless of Verbose
  base::cat("\n*****Predictive Ability, Bias, and Accuracy Results*****\n")
  base::cat("predictive ability: ", ystar_accuracy, "\n", sep = "")
  base::cat("bias: ", average_bias, "\n", sep = "")
  base::cat("accuracy: ", y_corrected_accuracy, "\n", sep = "")
  
  # Set the working directory back to the original path
  setwd(wd_path)
}
