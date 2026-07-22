#' Run preGSf90
#'
#' This function runs preGSf90 to process the genomic information for the BLUPf90 family.
#'
#' preGSf90 is the genomic pre-processor: it runs quality control on the genotypes
#' and constructs and inverts the genomic relationship matrix (G) and the pedigree
#' relationship matrix for genotyped animals (A22), producing \code{GimA22i}
#' (\eqn{G^{-1} - A22^{-1}}) for single-step GBLUP. It is driven by the renf90.par
#' parameter file, which must contain a \code{SNP_file} line (and, for map-aware QC
#' and downstream ssGWAS, a \code{map_file}); the QC and matrix-saving behaviour is
#' controlled by OPTION lines in the parameter file. Run \code{\link{run_renum}}
#' beforehand so the renf90 files and the genotype \code{_XrefID} are present.
#'
#' @param path_2_execs path to a folder that holds the preGSf90 executable. This field should be in quotes "".
#' @param input_files_dir directory with renf90.par, the genotype file and its _XrefID. Defaults to ".".
#' @param output_files_dir directory to store the results. Defaults to \code{input_files_dir}.
#' @param par_file name of the parameter file to run. Defaults to "renf90.par".
#' @param verbose logical; if TRUE prints the command and the log.
#'
#' @return (invisibly) a named character vector of the key output files that were
#'   produced (those depend on the OPTIONs used), e.g. \code{GimA22i},
#'   \code{freqdata.count}, \code{Gen_conflicts} and the cleaned SNP files. Names are
#'   the file names; values are the full paths. Intended for use by downstream steps.
#' @examples
#'
#' \donttest{
#'  # run_pregs(path_2_execs   = "/path/to/bf90_execs/",
#'  #           input_files_dir = "renumbered")
#' }
#'
#' @export
run_pregs <- function(path_2_execs,
                      input_files_dir = ".",
                      output_files_dir = input_files_dir,
                      par_file = "renf90.par",
                      verbose = TRUE) {

  # Checks
  if(file.exists(output_files_dir)){
    output_files_dir <- normalizePath(output_files_dir)
  } else {
    stop(paste("Directory", output_files_dir, "does not exist. Create it before running the function."))
  }
  path_2_execs <- normalizePath(path_2_execs)
  cur_dir <- getwd()
  on.exit(setwd(cur_dir), add = TRUE)   # restore working dir even on error

  # OS-specific executable name
  pregs <- if (.Platform$OS.type == "windows") "preGSf90.exe" else "preGSf90"
  pregs_exec <- file.path(path_2_execs, pregs)
  if (!file.exists(pregs_exec)) stop("Executable not found at: ", pregs_exec)
  if (!file.exists(file.path(input_files_dir, par_file)))
    stop("Parameter file not found: ", file.path(input_files_dir, par_file))

  input_files_dir <- normalizePath(input_files_dir)

  # preGSf90 reads the parameter-file name from standard input
  temp_input_file <- tempfile()
  writeLines(par_file, temp_input_file)
  command_pregs <- if (.Platform$OS.type == "windows")
    paste0("type ", temp_input_file, " | \"", pregs_exec, "\"") else
    paste0("cat ", temp_input_file, " | ", pregs_exec)
  if(verbose) cat("Running command:", command_pregs, "\n")

  setwd(input_files_dir)
  execute_command(command = command_pregs, logfile = "run_pregs.log")
  unlink(temp_input_file)

  # Collect the outputs that were produced (which exist depends on the OPTIONs used)
  fixed <- c("GimA22i", "Gi", "Ginv", "G", "A22", "A22i",
             "freqdata.count", "freqdata.count.after.clean",
             "Gen_call_rate", "Gen_conflicts", "Gen_conflicts_all", "sum2pq", "run_pregs.log")
  produced <- unique(c(fixed[file.exists(fixed)], list.files(pattern = "_clean")))  # <geno>_clean, _clean_XrefID

  # Move results to output_files_dir when it differs from the run directory
  if(output_files_dir != input_files_dir)
    for(f in produced) file.rename(from = f, to = file.path(output_files_dir, f))

  if(verbose){
    logp <- file.path(output_files_dir, "run_pregs.log")
    if(file.exists(logp)){ cat("Log file content:\n"); cat(readLines(logp), sep = "\n") }
  }

  paths <- file.path(output_files_dir, produced)
  names(paths) <- produced
  invisible(paths)
}
