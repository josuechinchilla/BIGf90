#' Function to run commands on the terminal and log output
#' 
#' @param command comment line used to run executable file
#' @param logfile logfile name
#' @noRd
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
#' @noRd
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
mutate_folds <- function(phenos, folds, num_folds, missing_value_code, id_col) {
  base::lapply(1:num_folds, function(i) {
    phenos %>%
      dplyr::mutate(V2 = base::ifelse(phenos[[id_col]] %in% base::unlist(folds[[i]]), missing_value_code, V2))  # id_col = position of id column; unlist -> compare IDs, not a data.frame
  })
}

# Create cross-validation datasets
create_cv_datasets <- function(run, fold, data_frame, dir_path, renf90, renf90_ped_name, input_files_dir, reuse_files) {
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
  if (base::file.exists(base::file.path(input_files_dir, "renf90.inb")))   # renf90.inb is optional
    file.symlink(base::file.path(input_files_dir, "renf90.inb"), base::file.path(dir_path, "renf90.inb"))
  file.symlink(base::file.path(input_files_dir, "renf90.tables"), base::file.path(dir_path, "renf90.tables"))
  
  # Genomic: link the pre-built G-inverse (Gi) + cleaned SNPs so the fold reuses them
  if (!is.null(reuse_files)) {
    for (f in reuse_files) file.symlink(f, base::file.path(dir_path, base::basename(f)))
  }
}
