#' Run renumf90
#'
#' This function runs renumf90.
#'
#' This function runs renumf90 to process the raw parameter(.par) file to be used with the blupf90 suite of programs.
#' The outputs will be the standard output files produced by renumf90. A log file called run_renum.log is also produced.
#'
#' @param path_2_execs path to a folder that holds the renumf90 executable. This field should be in quotes "".
#' @param raw_par_file name of the .par file that will be processed.  This field should be in quotes "".
#' @examples
#' ## Example
#'
#' # run_renum(path_2_execs = "/Users/johndoe/Desktop/bf90_execs/",
#'  raw_par_file = "weight_2022_no_cov_cv.par")
#'
#' @export
run_renum <- function(path_2_execs, raw_par_file) {

  # Function to run commands on the terminal and log output
  mac_terminal_command <- function(command, logfile) {
    base::system(paste(command, "2>&1 | tee -a", logfile))
  }

  raw_file <- paste(base::getwd(), raw_par_file, sep = "/")
  #run mac_terminal_command
  mac_terminal_command(command = paste0(path_2_execs, "renumf90 ", raw_par_file), logfile = "run_renum.log")

}
