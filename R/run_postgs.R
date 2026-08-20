#' Run postGSf90
#'
#' This function runs postGSf90 to extract SNP effects after a genomic evaluation.
#'
#' postGSf90 is the genomic post-processor: within the single-step GBLUP framework it
#' back-solves the SNP effects from the GEBVs (ssGWAS), and can compute the variance
#' explained by SNP windows, SNP weights (for a weighted G) and, with the relevant
#' OPTIONs, p-values and Manhattan plots. It is run AFTER a genomic \code{\link{run_blup}}
#' in the same directory, so the \code{solutions} file, the cleaned genotypes and the
#' genomic setup from \code{\link{run_pregs}} are available. The parameter file should
#' carry a \code{map_file} (for SNP positions / Manhattan plots) and the desired
#' postGSf90 OPTIONs (e.g. \code{windows_variance}, \code{Manhattan_plot}).
#'
#' @param path_2_execs path to a folder that holds the postGSf90 executable. This field should be in quotes "".
#' @param input_files_dir directory with renf90.par, the solutions file and the genomic files. Defaults to ".".
#' @param output_files_dir directory to store the results. Defaults to \code{input_files_dir}.
#' @param par_file name of the parameter file to run (may be a postGSf90-specific par
#'   with map_file and the postGS OPTIONs). Defaults to "renf90.par".
#' @param verbose logical; if TRUE prints the command and the log.
#'
#' @return (invisibly) a named character vector of the key output files that were
#'   produced, e.g. \code{snp_sol}, \code{windows_variance}, \code{snp_pred} and any
#'   Manhattan-plot files. Names are the file names; values are the full paths.
#' @examples
#' \dontrun{
#' run_postgs(path_2_execs    = "/path/to/bf90_execs/",
#'            input_files_dir = "renumbered")
#' }
#'
#' @export
run_postgs <- function(path_2_execs,
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
  postgs <- if (.Platform$OS.type == "windows") "postGSf90.exe" else "postGSf90"
  postgs_exec <- file.path(path_2_execs, postgs)
  if (!file.exists(postgs_exec)) stop("Executable not found at: ", postgs_exec)
  if (!file.exists(file.path(input_files_dir, par_file)))
    stop("Parameter file not found: ", file.path(input_files_dir, par_file))
  if (!file.exists(file.path(input_files_dir, "solutions")))
    warning("No 'solutions' file in ", input_files_dir,
            "; postGSf90 needs the solutions from a prior genomic run_blup() in the same directory.")

  input_files_dir <- normalizePath(input_files_dir)

  # postGSf90 reads the parameter-file name from standard input
  temp_input_file <- tempfile()
  writeLines(par_file, temp_input_file)
  command_postgs <- if (.Platform$OS.type == "windows")
    paste0("type ", temp_input_file, " | \"", postgs_exec, "\"") else
    paste0("cat ", temp_input_file, " | ", postgs_exec)
  if(verbose) cat("Running command:", command_postgs, "\n")

  setwd(input_files_dir)
  execute_command(command = command_postgs, logfile = "run_postgs.log")
  unlink(temp_input_file)

  # Collect the outputs that were produced (which exist depends on the OPTIONs used)
  fixed <- c("snp_sol", "chrsnp", "chrsnpvar", "windows_variance", "windows_segment",
             "snp_pred", "sum2pq", "run_postgs.log")
  plots <- list.files(pattern = "^(Sft|Vft|Pft|manplot)|\\.gnuplot$")   # graphic/plot files
  produced <- unique(c(fixed[file.exists(fixed)], plots))

  # Move results to output_files_dir when it differs from the run directory
  if(output_files_dir != input_files_dir)
    for(f in produced) file.rename(from = f, to = file.path(output_files_dir, f))

  if(verbose){
    logp <- file.path(output_files_dir, "run_postgs.log")
    if(file.exists(logp)){ cat("Log file content:\n"); cat(readLines(logp), sep = "\n") }
  }

  paths <- file.path(output_files_dir, produced)
  names(paths) <- produced
  invisible(paths)
}
