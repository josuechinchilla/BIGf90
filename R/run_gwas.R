#' Run a single-step GWAS (ssGBLUP / wssGBLUP) workflow
#'
#' Orchestrates a genome-wide association analysis with the BLUPf90 programs, using
#' weighted single-step GBLUP (wssGBLUP). It chains the individual wrappers:
#' \code{\link{run_renum}} renumbers the model, \code{\link{run_pregs}} does the genomic
#' QC and saves the cleaned SNPs, then the SNP effects are estimated iteratively --
#' \code{\link{run_blup}} solves the (weighted) genomic evaluation and
#' \code{\link{run_postgs}} back-solves the SNP effects and their weights, which are fed
#' back in for the next iteration. With \code{iterations = 1} all SNP weights are 1, i.e.
#' plain ssGBLUP; with more iterations the weights become SNP-specific (wssGBLUP). This
#' follows Wang et al. (2012).
#'
#' The model is defined by the parameter file you pass -- build it with
#' \code{\link{write_par}}, including a \code{SNP_file}, a \code{map_file} (needed for
#' 1-Mb windows and Manhattan plots) and the variance components (estimate them first,
#' e.g. with \code{\link{run_gibbs}}).
#'
#' @param path_2_execs folder with the BLUPf90 executables (renumf90, preGSf90, blupf90+, postGSf90). In quotes "".
#' @param raw_par_file the raw model parameter file (as fed to renumf90); its directory
#'   should also hold the data, pedigree, genotype and map files. In quotes "".
#' @param output_files_dir directory to run in and store all outputs (created if needed).
#' @param iterations number of wssGBLUP iterations. 1 = ssGBLUP (equal weights);
#'   2 (default) reproduces the second-iteration weights used by Vallejo et al. (2024).
#' @param windows_mbp optional window size in Mb for 1-Mb sliding windows
#'   (OPTION windows_variance_mbp); needs a map_file. NULL to skip.
#' @param windows_snp optional window size in number of adjacent SNPs
#'   (OPTION windows_variance); use when no map is available. NULL to skip.
#' @param which_weight optional postGSf90 OPTION which_weight value (the SNP-weight
#'   formula). NULL uses the program default.
#' @param snp_pvalue logical (default FALSE); if TRUE add OPTION snp_p_value to the
#'   blupf90+ and postGSf90 steps so exact-GWAS p-values are produced (snp_sol column
#'   var_a_hat), enabling p-value Manhattan and QQ plots. Slower and memory-heavy.
#' @param return_iteration which iteration's SNP solutions to return. NULL = last.
#' @param verbose logical; if TRUE prints progress and the program logs.
#'
#' @return (invisibly) a list with: \code{snp_sol} (the returned iteration's SNP-solution
#'   table as a named data frame), \code{windows_variance} (path to that iteration's
#'   window-variance file, if produced), \code{iterations}, \code{return_iteration} and
#'   \code{run_dir}. Each iteration also leaves snp_sol_iterN / windows_variance_iterN in
#'   the run directory.
#' @references Wang H et al. (2012) Genet Res 94:73-83. Vallejo RL et al. (2024)
#'   Aquaculture 586:740819. Misztal I et al. (2015) BLUPF90 family of programs.
#' @examples
#' \dontrun{
#' run_gwas(path_2_execs     = "/path/to/bf90_execs/",
#'          raw_par_file     = "gwas_model.par",
#'          output_files_dir = "gwas_run",
#'          iterations = 2, windows_mbp = 1)
#' }
#'
#' @export
run_gwas <- function(path_2_execs,
                     raw_par_file,
                     output_files_dir,
                     iterations = 2,
                     windows_mbp = NULL,
                     windows_snp = NULL,
                     which_weight = NULL,
                     snp_pvalue = FALSE,
                     return_iteration = NULL,
                     verbose = TRUE) {

  if(missing(raw_par_file)) stop("Provide the model parameter file in 'raw_par_file'.")
  raw_par_file <- normalizePath(raw_par_file)
  if(!file.exists(raw_par_file)) stop("Parameter file not found: ", raw_par_file)
  if(iterations < 1) stop("'iterations' must be >= 1.")
  input_dir <- dirname(raw_par_file)
  if(!dir.exists(output_files_dir)) dir.create(output_files_dir, recursive = TRUE)
  output_files_dir <- normalizePath(output_files_dir)

  say     <- function(...) if(verbose) base::message(...)
  add_opt <- function(lines, opt) if(any(grepl(opt, lines, fixed = TRUE))) lines else c(lines, paste("OPTION", opt))

  ## 1. renumf90 --------------------------------------------------------------
  say("run_gwas [renum]: renumbering ...")
  run_renum(path_2_execs, raw_par_file = raw_par_file, output_files_dir = output_files_dir, verbose = FALSE)
  inb <- file.path(input_dir, "renf90.inb")
  if(file.exists(inb)) file.rename(inb, file.path(output_files_dir, "renf90.inb"))

  par_path  <- file.path(output_files_dir, "renf90.par")
  par_lines <- readLines(par_path)

  # bring the SNP map into the run dir (renumf90 does not copy it)
  map_line <- grep("map_file", par_lines, ignore.case = TRUE, value = TRUE)
  map_name <- if(length(map_line)) basename(trimws(sub(".*map_file", "", map_line[1]))) else NULL
  if(!is.null(map_name) && !file.exists(file.path(output_files_dir, map_name)) && file.exists(file.path(input_dir, map_name)))
    file.copy(file.path(input_dir, map_name), file.path(output_files_dir, map_name))

  snp_line <- grep("SNP_file", par_lines, value = TRUE)
  if(length(snp_line) == 0) stop("No SNP_file in the parameter file - run_gwas needs a genomic model.")
  snp_base  <- basename(trimws(sub(".*SNP_file", "", snp_line[1])))
  clean_snp <- paste0(snp_base, "_clean")

  cur <- getwd(); on.exit(setwd(cur), add = TRUE)
  setwd(output_files_dir)

  ## 2. preGSf90: QC + cleaned SNP set ---------------------------------------
  say("run_gwas [preGS]: genotype QC ...")
  pre_par <- add_opt(par_lines, "saveCleanSNPs")
  pre_par <- pre_par[!grepl("OPTION[[:space:]]+missing", pre_par, ignore.case = TRUE)]  # preGSf90 rejects OPTION missing
  writeLines(pre_par, "renf90.par")
  run_pregs(path_2_execs, input_files_dir = output_files_dir, output_files_dir = output_files_dir, verbose = FALSE)
  if(!file.exists(clean_snp)) stop("preGSf90 did not produce the cleaned SNP file '", clean_snp,
                                   "'. Check run_pregs.log in ", output_files_dir)

  # number of clean SNPs = length of the genotype string on the first line
  first <- strsplit(trimws(readLines(clean_snp, n = 1)), "\\s+")[[1]]
  nsnp  <- nchar(first[length(first)])

  ## build the per-iteration parameter files (cleaned SNPs, weighted G) -------
  base_par <- gsub(snp_base, clean_snp, par_lines, fixed = TRUE)                    # point at cleaned SNPs
  if(!is.null(map_name)) base_par <- gsub(map_name, paste0(map_name, "_clean"), base_par, fixed = TRUE)  # and the matching cleaned map
  base_par <- gsub("saveCleanSNPs", "no_quality_control", base_par, fixed = TRUE)
  base_par <- add_opt(base_par, "no_quality_control")
  base_par <- c(base_par, "OPTION weightedG w")
  blup_par <- add_opt(base_par, "saveGInverse")                                    # blup: build+save weighted G-inverse
  if(snp_pvalue) blup_par <- add_opt(blup_par, "snp_p_value")                      # store MME-inverse elements for exact-GWAS p-values
  writeLines(blup_par, "blup.par")

  postgs_par <- add_opt(base_par, "readGInverse")                                  # postGS: reuse the G-inverse
  postgs_par <- postgs_par[!grepl("OPTION[[:space:]]+missing", postgs_par, ignore.case = TRUE)]  # postGSf90 rejects OPTION missing
  if(snp_pvalue) postgs_par <- add_opt(postgs_par, "snp_p_value")                  # compute p-values from the stored MME-inverse elements
  if(!is.null(windows_mbp))  postgs_par <- c(postgs_par, paste("OPTION windows_variance_mbp", windows_mbp))
  if(!is.null(windows_snp))  postgs_par <- c(postgs_par, paste("OPTION windows_variance", windows_snp))
  if(!is.null(which_weight)) postgs_par <- c(postgs_par, paste("OPTION which_weight", which_weight))
  writeLines(postgs_par, "postgs.par")

  ## 3. wssGBLUP iterations ---------------------------------------------------
  writeLines(as.character(rep(1, nsnp)), "w")                                      # iteration 1: equal weights (ssGBLUP)
  snp_sol_iters <- vector("list", iterations)
  for(i in seq_len(iterations)){
    say(sprintf("run_gwas [iter %d/%d]: blupf90+ + postGSf90 ...", i, iterations))
    if(verbose) run_blup(path_2_execs, par_file = "blup.par")
    else invisible(utils::capture.output(run_blup(path_2_execs, par_file = "blup.par")))
    run_postgs(path_2_execs, input_files_dir = output_files_dir, output_files_dir = output_files_dir,
               par_file = "postgs.par", verbose = FALSE)
    if(!file.exists("snp_sol")) stop("postGSf90 did not produce 'snp_sol' at iteration ", i,
                                     ". Check run_postgs.log in ", output_files_dir)
    ss <- .read_snp_sol("snp_sol")                                                # header-tolerant: keeps weights numeric
    file.copy("snp_sol", sprintf("snp_sol_iter%d", i), overwrite = TRUE)
    if(file.exists("windows_variance")) file.copy("windows_variance", sprintf("windows_variance_iter%d", i), overwrite = TRUE)
    snp_sol_iters[[i]] <- ss
    if(i < iterations && ncol(ss) >= 7) writeLines(as.character(ss[[7]]), "w")     # new weights = SNP weight column
  }

  ## return -------------------------------------------------------------------
  ri <- if(is.null(return_iteration)) iterations else return_iteration
  snp_sol <- snp_sol_iters[[ri]]                                                  # already named by .read_snp_sol()
  wv <- file.path(output_files_dir, sprintf("windows_variance_iter%d", ri))
  say("run_gwas: done. ", iterations, " iteration(s); returning iteration ", ri, " (",
      nrow(snp_sol), " SNPs).")
  invisible(list(snp_sol          = snp_sol,
                 windows_variance = if(file.exists(wv)) wv else NULL,
                 iterations       = iterations,
                 return_iteration = ri,
                 run_dir          = output_files_dir))
}
