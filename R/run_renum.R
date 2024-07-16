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
#' 
#' \donttest{
#'  #run_renum(path_2_execs = "path/bf90_execs/", 
#'  #raw_par_file = "weight_2022_no_cov_cv.par")
#' }
#' 
#' 
#' @export
run_renum <- function(path_2_execs, raw_par_file) {

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
    renum = "renumf90"
  } else if (.Platform$OS.type == "windows") {
    renum = "renumf90.exe"
  }

  raw_file <- paste(base::getwd(), raw_par_file, sep = "/")
  # Construct the command
  command_renum <- paste0(path_2_execs, renum," ",raw_par_file)

  # Check if executable and parameter files exist
  if (!file.exists(paste0(path_2_execs, renum))) {
    stop("Executable not found at: ", paste0(path_2_execs, renum))
  }
  if (!file.exists(raw_file)) {
    stop("Parameter file not found at: ", raw_file)
  }

  # Run the command and log the output
  output <- execute_command(command = command_renum, logfile = "run_renum.log")

  # Capture and print the log file content
  if (file.exists("run_renum.log")) {
    cat("Log file content:\n")
    cat(readLines("run_renum.log"), sep = "\n")
  } else {
    cat("Log file not created.\n")
  }
}
