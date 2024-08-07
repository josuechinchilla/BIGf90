#' Function to run commands on the terminal and log output
#' 
#' @param command comment line used to run executable file
#' @param logfile logfile name
#' 
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

#' 
create_folds <- function(data, num_folds) {
  n <- base::nrow(data)
  fold_size <- n %/% num_folds
  folds <- base::vector("list", num_folds)
  for (i in 1:num_folds) {
    start_index <- (i - 1) * fold_size + 1
    end_index <- if (i == num_folds) n else i * fold_size
    folds[[i]] <- data[start_index:end_index, ]
  }
  return(folds)
}

# Function to mutate phenotypes to missing value for given folds
mutate_folds <- function(phenos, folds, num_folds, missing_value_code) {
  base::lapply(1:num_folds, function(i) {
    phenos %>%
      dplyr::mutate(V2 = base::ifelse(V5 %in% folds[[i]], missing_value_code, V2))
  })
}

# Create cross-validation datasets
create_cv_datasets <- function(run, fold, data_frame, dir_path, renf90, renf90_ped_name, input_files_dir, snp_file_name) {
  file_name <- base::sprintf("renf90_run%d_fold%d.dat", run, fold)
  utils::write.table(data_frame, file = file_name, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  if (!base::file.exists(dir_path)) {
    base::dir.create(dir_path, recursive = TRUE)
  }
  base::file.rename(file_name, base::file.path(dir_path, file_name))
  
  modified_content <- base::gsub("renf90.dat", base::sprintf("renf90_run%d_fold%d.dat", run, fold), renf90)
  base::writeLines(modified_content, base::file.path(dir_path, base::sprintf("renf90_run%d_fold%d.par", run, fold)))
  
  # Create symbolic links instead of copying files
  file.symlink(base::file.path(renf90_ped_name), base::file.path(dir_path, basename(renf90_ped_name)))
  file.symlink(base::file.path(input_files_dir, "renf90.fields"), base::file.path(dir_path, "renf90.fields"))
  file.symlink(base::file.path(input_files_dir, "renf90.inb"), base::file.path(dir_path, "renf90.inb"))
  file.symlink(base::file.path(input_files_dir, "renf90.tables"), base::file.path(dir_path, "renf90.tables"))
  
  if (!is.null(snp_file_name)) {
    Xref_file <- paste0(snp_file_name, "_XrefID")
    file.symlink(base::file.path(snp_file_name), base::file.path(dir_path, basename(snp_file_name)))
    file.symlink(base::file.path(Xref_file), base::file.path(dir_path, basename(Xref_file)))
  }
}
