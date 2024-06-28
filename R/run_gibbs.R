#' Run gibbsf90+
#'
#' This function runs gibbsf90+.
#'
#' This function runs gibbsf90+ using input from the user and a renf90.par file.
#' This function is written to use a renf90.par paramererfile, therefore the function run_renumf90 needs to be used beforehand to process a .par file created by the user.
#' The user will have to enter the number of samples in the MCMC chain, the number of samples to burn and the number used to thin samples. A log file called run_gibbs.log is also produced.
#'
#' @param path_2_execs path to a folder that holds the renumf90 executable. This field should be in quotes "".
#' @param gibbs_iter number of samples in the Gibbs sampler.
#' @param gibbs_burn number of samples to be discarded at the begining of the Gibbs sampler
#' @param gibbs_keep the interval to save samples (thinning). Entering a 1 means all samples are kept.
#' @examples
#' ## Example
#'
#' # run_gibbs( path_2_execs = "/Users/johndoe/Desktop/bf90_execs/",
#' # gibbs_iter = 250000,
#' # gibbs_burn = 20000
#' # gibbs_keep = 1)
#'
#' @export
run_gibbs <- function(path_2_execs, gibbs_iter, gibbs_burn, gibbs_keep) {
  # Define the mac_terminal_command function
  mac_terminal_command <- function(command, logfile) {
    base::system(paste(command, "2>&1 | tee -a", logfile))
  }

  # Construct the command with the given inputs
  command_gibbs <- paste0("printf 'renf90.par \n", gibbs_iter, "\n", gibbs_burn, "\n", gibbs_keep, "' | ", path_2_execs, "gibbsf90+")

  # Run the command and log the output
  mac_terminal_command(command = command_gibbs, logfile = "run_gibbs.log")
}
