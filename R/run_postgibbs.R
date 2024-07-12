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
    postgibbs = "postgibbsf90"
  } else if (.Platform$OS.type == "windows") {
    postgibbs = "postgibbsf90+.exe"
  }
  postgibbs_exec <- paste0(path_2_execs, postgibbs)

  # Check if executable exists

  if (!file.exists(postgibbs_exec)) {
    stop("Executable not found at: ", postgibbs_exec)
  }
  if (!file.exists("renf90.par")) {
    stop("Parameter file not found: renf90.par")
  }

  # Create a temporary file with the inputs
  temp_input_file <- tempfile()
  writeLines(c("renf90.par", postgibbs_burn, postgibbs_keep, "0"), temp_input_file)

  # Construct the command to use the input file
  if (.Platform$OS.type == "unix") {
    command_postgibbs <- paste0("cat ", temp_input_file, " | ", postgibbs_exec)
  } else if (.Platform$OS.type == "windows") {
    command_postgibbs <- paste0("type ", temp_input_file, " | \"", postgibbs_exec, "\"")
  }

  # Debugging: Print the constructed command
  cat("Running command:", command_postgibbs, "\n")

  # Run the command and log the output
  output <- execute_command(command = command_postgibbs, logfile = "run_postgibbs.log")

  # Remove the temporary file
  unlink(temp_input_file)

  # Debugging: Print the output
  cat("Command output:", output, "\n")

  # Capture and print the log file content
  if (file.exists("run_postgibbs.log")) {
    cat("Log file content:\n")
    cat(readLines("run_postgibbs.log"), sep = "\n")
  } else {
    cat("Log file not created.\n")
  }
}
