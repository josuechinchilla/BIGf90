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
  # Function to run commands on the terminal and log output
  execute_command <- function(command, logfile) {
    if (.Platform$OS.type == "unix") {
      output <- system(paste(command, "2>&1 | tee -a", logfile), intern = TRUE)
    } else if (.Platform$OS.type == "windows") {
      output <- system(paste("cmd /c", command, ">", logfile, "2>&1"), intern = TRUE)
    } else {
      stop("Unsupported OS type")
    }
    return(output)
  }

  #Assign .exe or not based on OS
  if (.Platform$OS.type == "unix") {
    gibbs = "gibbsf90+"
  } else if (.Platform$OS.type == "windows") {
    gibbs = "gibbsf90+.exe"
  }
  gibbs_exec <- paste0(path_2_execs, gibbs)

  # Check if executable exists
  if (!file.exists(gibbs_exec)) {
    stop("Executable not found at: ", gibbs_exec)
  }
  if (!file.exists("renf90.par")) {
    stop("Parameter file not found: renf90.par")
  }

  # Create a temporary file with the inputs
  temp_input_file <- tempfile()
  writeLines(c("renf90.par", gibbs_iter, gibbs_burn, gibbs_keep), temp_input_file)

  # Construct the command to use the input file
  if (.Platform$OS.type == "unix") {
    command_gibbs <- paste0("cat ", temp_input_file, " | ", gibbs_exec)
  } else if (.Platform$OS.type == "windows") {
    command_gibbs <- paste0("type ", temp_input_file, " | \"", gibbs_exec, "\"")
  }

  # Debugging: Print the constructed command
  cat("Running command:", command_gibbs, "\n")

  # Run the command and log the output
  output <- execute_command(command = command_gibbs, logfile = "run_gibbs.log")

  # Remove the temporary file
  unlink(temp_input_file)

  # Debugging: Print the output
  cat("Command output:", output, "\n")

  # Capture and print the log file content
  if (file.exists("run_gibbs.log")) {
    cat("Log file content:\n")
    cat(readLines("run_gibbs.log"), sep = "\n")
  } else {
    cat("Log file not created.\n")
  }
}
