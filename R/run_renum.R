#' Run renumf90
#'
#' This function runs renumf90.
#'
#' This function runs renumf90 to process the raw parameter(.par) file to be used with the blupf90 suite of programs.
#' The outputs will be the standard output files produced by renumf90. A log file called run_renum.log is also produced.
#'
#' @param path_2_execs path to a folder that holds the renumf90 executable. This field should be in quotes "".
#' @param raw_par_file name of the .par file that will be processed.  This field should be in quotes "".
#' @param output_files_dir path to the folder to store the output files renadd03.ped renf90.dat renf90.fields renf90.inb renf90.par renf90.tables run_renum.log
#' @param verbose logical if TRUE prints log information
#' @examples
#'
#' \donttest{
#'  #run_renum(path_2_execs = "path/bf90_execs/",
#'  #input_files_dir = "weight_2022_no_cov_cv.par")
#' }
#'
#'
#' @export
run_renum <- function(path_2_execs = ".",
                      raw_par_file = NULL,
                      snp_file_name = NULL,
                      output_files_dir = NULL,
                      verbose = TRUE) {
  
  # Checks
  if(is.null(output_files_dir)) stop("Define a output directory in output_files_dir")
  
  if(file.exists(output_files_dir)){
    check_files <- list.files(output_files_dir)
    if(length(check_files) > 0) warning(paste("Directory", output_files_dir, "is not empty. Some files may be replaced."))
    output_files_dir <- normalizePath(output_files_dir)
  } else {
    stop(paste("Directory '", output_files_dir, "' does not exist. Create it before running the function."))
  }
  
  path_2_execs <- normalizePath(path_2_execs)
  
  if(is.null(raw_par_file)) stop("Define raw_par_file.")
  raw_file <- normalizePath(raw_par_file)
  input_files_dir <- dirname(raw_file)
  
  if (!file.exists(raw_file)) {
    stop("Parameter file not found at: ", raw_file)
  }
  
  #Assign .exe or not based on OS
  if (.Platform$OS.type == "unix") {
    renum = "renumf90"
  } else if (.Platform$OS.type == "windows") {
    renum = "renumf90.exe"
  }
  
  cur_dir <- getwd() # save working directory location
  
  # Construct the command
  command_renum <- paste0(file.path(path_2_execs, renum)," ", paste0("'",raw_file, "'"))
  
  # Check if executable and parameter files exist
  if (!file.exists(file.path(path_2_execs, renum))) {
    stop("Executable not found at: ", paste0("'",path_2_execs, renum,"'"))
  }
  
  # Run the command and log the output
  setwd(input_files_dir)
  output <- execute_command(command = command_renum, logfile = "run_renum.log")
  
  # Capture and print the log file content
  if(verbose){
    if (file.exists("run_renum.log")) {
      cat("Log file content:\n")
      cat(readLines("run_renum.log"), sep = "\n")
    } else {
      cat("Log file not created.\n")
    }
  }
  
  # Move generated files to working directory
  files_res <- c("renadd03.ped", "renf90.dat", "renf90.fields", "renf90.inb", "renf90.par" ,"renf90.tables", "run_renum.log")
  for(i in 1:length(files_res)) file.rename(from = files_res[i], to = file.path(output_files_dir,files_res[i]))
  
  # copy snp file
  if(!is.null(snp_file_name) & !file.exists(file.path(output_files_dir, basename(snp_file_name)))){
    snp_file_name <- normalizePath(snp_file_name)
    file.symlink(base::file.path(snp_file_name), base::file.path(output_files_dir, basename(snp_file_name)))
    file.symlink(paste0(base::file.path(snp_file_name),"_XrefID"), base::file.path(output_files_dir, paste0(basename(snp_file_name), "_XrefID")))
  }
  
  # Return to past working directory
  setwd(cur_dir)
}





