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
#'
#' @return No return value, called for side effects: runs renumf90, producing the renumbered renf90.* files, the renumbered pedigree (renadd0X.ped) and the genotype cross-reference (_XrefID) file.
#' @examples
#'
#' \donttest{
#'  #run_renum(path_2_execs = "path/bf90_execs/",
#'  #          raw_par_file = "my_analysis.par",
#'  #          output_files_dir = "results")
#' }
#'
#'
#' @export
run_renum <- function(path_2_execs = ".",
                      raw_par_file = NULL,
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
  
  #### Check par files ####
  parfile <- base::readLines(raw_file)
  
  ## Check DATAFILE
  datafile <- parfile[grep("DATAFILE", parfile) + 1]
  if(length(datafile) == 0) stop(paste("DATAFILE was not defined in the parameter file", raw_file))
  
  # Check if file exist
  if(!file.exists(file.path(input_files_dir, datafile))) stop(paste("File:", datafile, "defined in the file", raw_file,"line", grep("DATAFILE", parfile) + 1,"does not exist."))
  
  ## Check pedfile
  pedfile <- parfile[which(grepl("^FILE$", parfile) | grepl("^FILE ", parfile) | grepl("^FILE#", parfile)) + 1]
  
  if(length(pedfile) != 0){
    pedfile <- pedfile[1]   # first match only -> scalar (avoids length>1 in conditions)
    # Check if file exist
    if(!file.exists(file.path(input_files_dir, pedfile))) stop(paste("File:", pedfile, "defined in the file", raw_file,"line", which(grepl("^FILE$", parfile) | grepl("^FILE ", parfile) | grepl("^FILE#", parfile)) + 1,"does not exist."))
  } else pedfile <- "not defined"
  
  ## Check snpfile
  snpfile <- parfile[grep("SNP_FILE", parfile) + 1]
  
  if(length(snpfile) != 0){
    snpfile <- snpfile[1]   # first match only -> scalar (avoids length>1 in conditions)
    # Check if file exist
    if(!file.exists(file.path(input_files_dir, snpfile))) stop(paste("File:", snpfile, "defined in the file", raw_file,"line", grep("SNP_FILE", parfile) + 1, "does not exist."))
  } else snpfile <- "not defined"
  
  if(verbose){
    cat(paste("DATAFILE:", datafile, "\n"))
    cat(paste("FILE (pedigree file):", pedfile, "\n"))
    cat(paste("SNP_FILE:", snpfile, "\n"))
  }
  
  #####
  
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
  result <- if(file.exists("run_renum.log")) readLines("run_renum.log") else character()
  
  if(verbose){
    if (file.exists("run_renum.log")) {
      cat("Log file content:\n")
      cat(result, sep = "\n")
    } else {
      cat("Log file not created.\n")
    }
  }
  
  if(any(grepl("EFFECT: not found", result))) {
    files_res <- c("run_renum.log")
    file.rename(from = files_res, to = file.path(output_files_dir,files_res))
    stop(paste("Effect not found. Check ",file.path(output_files_dir,files_res)))
  } else {
    # Move generated files to working directory
    renadd <- base::list.files(pattern = "^renadd[0-9]+\\.ped$")   # renumf90 names ped by effect position
    if (file.exists("renf90.inb")) {   # renumf90 just wrote it in the current dir; files not moved yet
      files_res <- c(renadd, "renf90.dat", "renf90.fields", "renf90.inb", "renf90.par", "renf90.tables", "run_renum.log")
    } else {
      files_res <- c(renadd, "renf90.dat", "renf90.fields", "renf90.par" ,"renf90.tables", "run_renum.log")
    }
    
    for(i in 1:length(files_res)) file.rename(from = files_res[i], to = file.path(output_files_dir,files_res[i]))
  }
  
  # copy snp file
  if(snpfile != "not defined" && !file.exists(file.path(output_files_dir, basename(snpfile)))){
    snpfile <- normalizePath(snpfile)
    file.symlink(base::file.path(snpfile), base::file.path(output_files_dir, basename(snpfile)))
    file.symlink(paste0(base::file.path(snpfile),"_XrefID"), base::file.path(output_files_dir, paste0(basename(snpfile), "_XrefID")))
  }
  
  # copy map file (OPTION map_file), if referenced in the parameter file and present
  map_idx <- grep("map_file", parfile, ignore.case = TRUE)
  if(length(map_idx) != 0){
    map_name <- basename(trimws(sub(".*map_file", "", parfile[map_idx[1]])))
    map_src  <- file.path(input_files_dir, map_name)
    if(file.exists(map_src) && !file.exists(file.path(output_files_dir, map_name)))
      file.symlink(normalizePath(map_src), file.path(output_files_dir, map_name))
  }

  # Return to past working directory
  setwd(cur_dir)
}





