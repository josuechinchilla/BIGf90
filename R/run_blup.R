#' Run blupf90+
#'
#' This function runs blupf90+ using a pre-processed parameter file called renf90.par.
#'
#' This function runs blupf90+ using a parameter file named renf90.par.
#' Since this function depends only on the renf90.par parameter file, the only input needed from the user is a path where the blupf90+ executable is located. A log file called run_blup.log is also produced.
#'
#' @param path_2_execs path to a folder that holds the renumf90 executable. This field should be in quotes "".
#' @param par_file name of the parameter file to run. Defaults to "renf90.par".
#'
#' @return No return value, called for side effects: runs blupf90+ using the parameter file, producing the solutions file and the run_blup.log log file.
#' @examples
#' ## Example
#'
#' # run_blup(path_2_execs = "/Users/johndoe/Desktop/bf90_execs/")
#'
#' @export
run_blup <- function(path_2_execs, par_file = "renf90.par") {
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
    blup = "blupf90+"
  } else if (.Platform$OS.type == "windows") {
    blup = "blupf90+.exe"
  }

  # Check if executable exists
  if (!file.exists(file.path(path_2_execs, blup))) {
    stop("Executable not found at: ", file.path(path_2_execs, blup))
  }
  if (!file.exists(par_file)) {
    stop("Parameter file not found: ", par_file)
  }

  # Run the command
  output <- execute_command(command = paste0(file.path(path_2_execs, blup), " ", par_file), logfile = "run_blup.log")

  # Capture and print the log file content
  if (file.exists("run_blup.log")) {
    cat("Log file content:\n")
    cat(readLines("run_blup.log"), sep = "\n")
  } else {
    cat("Log file not created.\n")
  }

}
