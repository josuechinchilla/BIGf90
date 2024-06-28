#' Run blupf90+
#'
#' This function runs blupf90+ using a pre-processed parameter file called renf90.par.
#'
#' This function runs blupf90+ using a parameter file named renf90.par.
#' Since this function depends only on the renf90.par parameter file, the only input needed from the user is a path where the blupf90+ executable is located. A log file called run_blup.log is also produced.
#'
#' @param path_2_execs path to a folder that holds the renumf90 executable. This field should be in quotes "".
#' @examples
#' ## Example
#'
#' # run_blup(path_2_execs = "/Users/johndoe/Desktop/bf90_execs/")
#'
#' @export
run_blup <- function(path_2_execs) {

  # Function to run commands on the terminal and log output
  mac_terminal_command <- function(command, logfile) {
    base::system(paste(command, "2>&1 | tee -a", logfile))
  }

  #run mac_terminal_command
  mac_terminal_command(command = paste0(path_2_execs, "blupf90+ renf90.par"), logfile = "run_blup.log")

}
