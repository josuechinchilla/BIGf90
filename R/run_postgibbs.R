#' Run postgibbsf90
#'
#' This function runs postgibbsf90 using the renf90.par file used to run gibbsf90+
#'
#' This function runs postgibbsf90 to provide diagnostic statistics,posterior means and standard deviations for analyses performed through gibbsf90+.
#' This function is written to run using the same renf90.par file used with run_gibbs and its output files are be the standard output files produced by postgibbsf90.
#' The user will have to enter the number of samples to burn and the number used to thin samples. A log file called run_postgibbs.log is also produced.
#'
#' @param path_2_execs path to a folder that holds the renumf90 executable. This field should be in quotes "".
#' @param postgibbs_burn number of samples to be discarded at the begining of the Gibbs sampler
#' @param postgibbs_keep the interval to save samples (thinning). Entering a 1 means all samples are kept.#'
#' @examples
#' ## Example
#'
#' # run_postgibbs( path_2_execs = "/Users/johndoe/Desktop/bf90_execs/",
#' # postgibbs_burn =1,
#' # postgibbs_keep = 100)
#'
#' @export
run_postgibbs <- function(path_2_execs, postgibbs_burn, postgibbs_keep) {
  # Define the mac_terminal_command function
  mac_terminal_command <- function(command, logfile) {
    base::system(paste(command, "2>&1 | tee -a", logfile))
  }

  # Construct the command with the given inputs
  command_postgibbs <- paste0("printf 'renf90.par \n", postgibbs_burn, "\n", postgibbs_keep, "\n 0 ' | ", path_2_execs, "postgibbsf90")

  # Run the command and log the output
  mac_terminal_command(command = command_postgibbs, logfile = "run_postgibbs.log")
}
