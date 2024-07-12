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
  if (!file.exists(paste0(path_2_execs, blup))) {
    stop("Executable not found at: ", paste0(path_2_execs, blup))
  }
  if (!file.exists("renf90.par")) {
    stop("Parameter file not found: renf90.par")
  }

  # Run the command
  output <- execute_command(command = paste0(path_2_execs, blup, " renf90.par"), logfile = "run_blup.log")

  # Capture and print the log file content
  if (file.exists("run_blup.log")) {
    cat("Log file content:\n")
    cat(readLines("run_blup.log"), sep = "\n")
  } else {
    cat("Log file not created.\n")
  }

}
