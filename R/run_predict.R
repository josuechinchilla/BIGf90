#' Run predictf90
#'
#' This function runs predict90.
#'
#' This function runs predictf90 using a pre-processed parameter file called renf90.par to calculate adjusted phenotypes.
#' This function is written to use a renf90.par file and a solutions file. Therefore run_renum and run_blup should be ran before using this function.
#' The output files are be the standard output files produced by predictf90. A log file called run_predict.log is also produced.
#'
#' @param path_2_execs path to a folder that holds the renumf90 executable. This field should be in quotes "".
#' @examples
#' ## Example
#'
#' # run_predict(path_2_execs = "/Users/johndoe/Desktop/bf90_execs/")
#'
#' @export
run_predict <- function(path_2_execs) {

  # Function to run commands on the terminal and log output
  mac_terminal_command <- function(command, logfile) {
    base::system(paste(command, "2>&1 | tee -a", logfile))
  }

  #run mac_terminal_command
  mac_terminal_command(command = paste0(path_2_execs, "predictf90 renf90.par"), logfile = "run_predict.log")

}
