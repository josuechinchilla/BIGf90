#' Run predf90
#'
#' This function runs predf90 to predict direct genomic values (DGV) for animals
#' from their genotypes and previously estimated SNP effects.
#'
#' predf90 computes the direct genomic value (\eqn{\hat{u} = Z\hat{a}}) for young or
#' newly genotyped animals, using the SNP effects in the \code{snp_pred} file produced
#' by \code{\link{run_postgs}} (postGSf90). Unlike the other application programs it is
#' driven by command-line options rather than a parameter file: the genotypes of the
#' animals to predict are supplied through \code{snp_file}. With \code{acc = TRUE} it
#' also computes reliabilities, which requires \code{OPTION snp_p_value} in the
#' preceding \code{\link{run_blup}} and \code{OPTION snp_var} in postGSf90 (so the
#' \code{snp_var} file is present).
#'
#' @param path_2_execs path to a folder that holds the predf90 executable. This field should be in quotes "".
#' @param snp_file genotype file of the animals to be predicted (passed to predf90 as
#'   \code{--snpfile}); same format as the SNP file used by preGSf90. This field should be in quotes "".
#' @param input_files_dir directory with the \code{snp_pred} file (from run_postgs) and
#'   the genotype file. Defaults to ".".
#' @param output_files_dir directory to store the results. Defaults to \code{input_files_dir}.
#' @param acc logical; if TRUE adds \code{--acc} to compute reliabilities (needs the
#'   \code{snp_var} file from run_postgs() and OPTION snp_p_value in run_blup()).
#' @param acc_type accuracy type used with \code{acc}: 1.0 for dairy cattle (reliability)
#'   or 0.5 for beef cattle (BIF accuracy). Defaults to 1.0.
#' @param verbose logical; if TRUE prints the command and the log.
#'
#' @return (invisibly) a named character vector of the files predf90 produced (detected
#'   as the files new to the run directory, plus run_predf.log). Names are the file
#'   names; values are the full paths.
#' @examples
#'
#' \donttest{
#'  # run_predf(path_2_execs    = "/path/to/bf90_execs/",
#'  #           snp_file        = "young_animals.geno",
#'  #           input_files_dir = "renumbered")
#' }
#'
#' @export
run_predf <- function(path_2_execs,
                      snp_file,
                      input_files_dir = ".",
                      output_files_dir = input_files_dir,
                      acc = FALSE,
                      acc_type = 1.0,
                      verbose = TRUE) {

  # Checks
  if(missing(snp_file) || is.null(snp_file))
    stop("Provide the genotype file of the animals to predict in 'snp_file'.")
  if(file.exists(output_files_dir)){
    output_files_dir <- normalizePath(output_files_dir)
  } else {
    stop(paste("Directory", output_files_dir, "does not exist. Create it before running the function."))
  }
  path_2_execs <- normalizePath(path_2_execs)
  cur_dir <- getwd()
  on.exit(setwd(cur_dir), add = TRUE)   # restore working dir even on error

  # OS-specific executable name
  predf <- if (.Platform$OS.type == "windows") "predf90.exe" else "predf90"
  predf_exec <- file.path(path_2_execs, predf)
  if (!file.exists(predf_exec)) stop("Executable not found at: ", predf_exec)
  if (!file.exists(file.path(input_files_dir, snp_file)))
    stop("Genotype file not found: ", file.path(input_files_dir, snp_file))
  if (!file.exists(file.path(input_files_dir, "snp_pred")))
    warning("No 'snp_pred' file in ", input_files_dir,
            "; predf90 needs the SNP effects from a prior run_postgs() in the same directory.")
  if (acc && !file.exists(file.path(input_files_dir, "snp_var")))
    warning("acc = TRUE but no 'snp_var' file in ", input_files_dir,
            "; it is produced by run_postgs() with OPTION snp_var (plus OPTION snp_p_value in run_blup()).")

  input_files_dir <- normalizePath(input_files_dir)

  # predf90 is driven by command-line options (no parameter file)
  command_predf <- paste0(predf_exec, " --snpfile ", snp_file)
  if(acc) command_predf <- paste0(command_predf, " --acc --acc_type ", acc_type)
  if(verbose) cat("Running command:", command_predf, "\n")

  setwd(input_files_dir)
  before <- list.files()
  execute_command(command = command_predf, logfile = "run_predf.log")
  produced <- unique(c(setdiff(list.files(), before), "run_predf.log"))   # files new to the run dir + the log
  produced <- produced[file.exists(produced)]

  # Move results to output_files_dir when it differs from the run directory
  if(output_files_dir != input_files_dir)
    for(f in produced) file.rename(from = f, to = file.path(output_files_dir, f))

  if(verbose){
    logp <- file.path(output_files_dir, "run_predf.log")
    if(file.exists(logp)){ cat("Log file content:\n"); cat(readLines(logp), sep = "\n") }
  }

  paths <- file.path(output_files_dir, produced)
  names(paths) <- produced
  invisible(paths)
}
