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
    predict = "predictf90"
  } else if (.Platform$OS.type == "windows") {
    predict = "predictf90.exe"
  }

  # Check if executable exists
  if (!file.exists(paste0(path_2_execs, predict))) {
    stop("Executable not found at: ", paste0(path_2_execs, predict))
  }
  if (!file.exists("renf90.par")) {
    stop("Parameter file not found: renf90.par")
  }

  # Construct the command
  command_predict <- paste0(path_2_execs, predict, " renf90.par")

  # Run the command and log the output
  output <- execute_command(command = command_predict, logfile = "run_predict.log")

  # Capture and print the log file content
  if (file.exists("run_predict.log")) {
    cat("Log file content:\n")
    cat(readLines("run_predict.log"), sep = "\n")
  } else {
    cat("Log file not created.\n")
  }
}
