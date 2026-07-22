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


## ---------------------------------------------------------------------------
## Package hooks
## ---------------------------------------------------------------------------

# Register the NSE column names used across the package so R CMD check does not
# flag them as undefined global variables.
utils::globalVariables(c("V1", "V2", "V3", "V4", "trait.effect"))

.onAttach <- function(libname, pkgname){
  msg <- paste0("\nBIGf90 is a wrapper package for BLUPf90 family of programs, please use the following citations:\n\n",
 "BIGf90 Reference: \n",
  "Chinchilla-Vargas J, Taniguti C, Sandercock A, Breeding Insight Team (2024). BIGf90:
  Breeding Insight Genomics R front face to blupf90 modules. R package version 1.0.0,
  https://github.com/Breeding-Insight/BIGf90 \n\n",
  "BLUPF90 Reference: \n",
  "Misztal, I., S. Tsuruta, D.A.L. Lourenco, I. Aguilar, A. Legarra, and Z. Vitezica. 2014.\n",
  "Manual for BLUPF90 family of programs: http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all2.pdf \n\n",
  "The BLUPF90 programs are free for research, but their use should be acknowledged in publications. For commercial use,
please contact Ignacy Misztal (ignacy@uga.edu) or Daniela Lourenco (danilino@uga.edu). Provide your name, company, and purpose of use in the email."
               )

  packageStartupMessage(msg)
}


## ---------------------------------------------------------------------------
## Internal helpers for the GWAS readers and plot functions
## ---------------------------------------------------------------------------

# Canonical postGSf90 snp_sol column names (see BLUPF90 manual, POSTGSF90 outputs).
.snp_sol_names <- c("trait", "effect", "snp", "chr", "pos", "snp_effect",
                    "weight", "variance_explained", "var_a_hat")

# Name (the first columns of) a snp_sol data frame with the canonical names.
.name_snp_sol <- function(d) {
  k <- min(ncol(d), length(.snp_sol_names))
  names(d)[seq_len(k)] <- .snp_sol_names[seq_len(k)]
  d
}

# Resolve x to a single file path: a directory -> file.path(x, fname); a file -> x.
.resolve_path <- function(x, fname) {
  if(!is.character(x) || length(x) != 1)
    stop("Provide a directory or a file path (a single character string).")
  if(dir.exists(x)) file.path(x, fname) else x
}

# Read a postGSf90 'snp_sol' file, tolerating the optional header line that newer
# postGSf90 versions write. Detection: if the first token of line 1 is not numeric
# it is treated as a header and skipped.
.read_snp_sol <- function(path) {
  if(!file.exists(path)) stop("snp_sol file not found: ", path)
  first_tok  <- strsplit(trimws(readLines(path, n = 1)), "\\s+")[[1]][1]
  has_header <- is.na(suppressWarnings(as.numeric(first_tok)))
  d <- utils::read.table(path, header = FALSE, skip = if(has_header) 1L else 0L,
                         stringsAsFactors = FALSE)
  .name_snp_sol(d)
}

# Read a postGSf90 'windows_variance' file and derive, per window, the chromosome
# and the midpoint position. Columns (manual): 6 = start "Chr_Position",
# 7 = end "Chr_Position", last = variance explained.
.read_windows <- function(path) {
  if(!file.exists(path)) stop("windows_variance file not found: ", path)
  w <- utils::read.table(path, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  sp <- strsplit(as.character(w[[6]]), "_")
  ep <- strsplit(as.character(w[[7]]), "_")
  head_num <- function(z) suppressWarnings(as.numeric(z[1]))            # chromosome
  tail_num <- function(z) suppressWarnings(as.numeric(z[length(z)]))    # position
  data.frame(chr = vapply(sp, head_num, numeric(1)),
             pos = (vapply(sp, tail_num, numeric(1)) + vapply(ep, tail_num, numeric(1))) / 2,
             variance_explained = suppressWarnings(as.numeric(w[[ncol(w)]])),
             stringsAsFactors = FALSE)
}

# Per-SNP chi-square(1) statistic for exact ssGWAS: a^2 / var(a) (Aguilar et al. 2019).
# Needs OPTION snp_p_value in the run (otherwise var_a_hat is zero).
.snp_chisq <- function(d) {
  if(is.null(d$var_a_hat) || all(d$var_a_hat == 0, na.rm = TRUE))
    stop("No p-values available (snp_sol column 'var_a_hat' is absent or all zero). ",
         "Re-run run_gwas(..., snp_pvalue = TRUE) to add OPTION snp_p_value.")
  d$snp_effect^2 / d$var_a_hat
}

# Genomic inflation factor (Devlin & Roeder 1999): median observed chi-square over
# the null median, qchisq(0.5, 1) = 0.4549.
.lambda_gc <- function(chisq)
  stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, df = 1)

# Benjamini-Hochberg FDR p-value threshold: the largest p(i) with p(i) <= (i/m) q.
# Returns NA if no test passes at level q.
.fdr_threshold <- function(p, q = 0.05) {
  p <- sort(p[is.finite(p)]); m <- length(p)
  if(m == 0) return(NA_real_)
  pass <- which(p <= (seq_len(m) / m) * q)
  if(length(pass)) p[max(pass)] else NA_real_
}

# Read a gibbsf90+ 'gibbs_samples' file into a matrix of (co)variance-component
# samples. Layout: line 1 = "<code> <n_params> ..."; the next n_params lines are
# descriptors "pos eff1 eff2 trt1 trt2"; then per stored sample a PAIR of lines
# (a "round n_params" line followed by the n_params values). Returns the value
# matrix (named by effect), the round numbers, and the raw descriptors.
.read_gibbs_matrix <- function(path) {
  if(!file.exists(path)) stop("Gibbs sample file not found: ", path)
  lines <- readLines(path)
  if(length(lines) < 5) stop("Gibbs sample file '", path, "' has no samples.")
  hdr <- suppressWarnings(as.integer(strsplit(trimws(lines[1]), "\\s+")[[1]]))
  np  <- hdr[2]
  if(is.na(np) || np < 1) stop("Could not read the parameter count from ", path, ".")
  desc <- lapply(lines[2:(1 + np)], function(x) suppressWarnings(as.integer(strsplit(trimws(x), "\\s+")[[1]])))
  body <- lines[-(1:(1 + np))]
  round_lines <- body[seq(1, length(body), by = 2)]
  val_lines   <- body[seq(2, length(body), by = 2)]
  n <- min(length(round_lines), length(val_lines))
  rounds <- suppressWarnings(as.integer(sub("\\s.*", "", trimws(round_lines[seq_len(n)]))))
  vals <- do.call(rbind, lapply(strsplit(trimws(val_lines[seq_len(n)]), "\\s+"),
                                function(x) suppressWarnings(as.numeric(x))))
  nm <- vapply(desc, function(d) {
    if(length(d) >= 3 && !is.na(d[2]) && d[2] > 0 && d[2] == d[3]) paste0("var_eff", d[2])
    else if(length(d) >= 2 && !is.na(d[2]) && d[2] == 0)          "residual_var"
    else if(length(d) >= 3 && d[2] != d[3])                       paste0("cov_", d[2], "_", d[3])
    else                                                          "param"
  }, character(1))
  colnames(vals) <- make.unique(nm[seq_len(ncol(vals))])
  ok <- stats::complete.cases(vals) & is.finite(rounds)
  list(values = vals[ok, , drop = FALSE], rounds = rounds[ok], desc = desc)
}

# Read a postgibbsf90 'postgibbs_samples' file (tabular: sample#, iteration#,
# n_params, then the n_params (co)variance components) into a value matrix + rounds.
.read_postgibbs_matrix <- function(path) {
  if(!file.exists(path)) stop("postgibbs_samples file not found: ", path)
  if(length(readLines(path, n = 1)) == 0)
    stop("postgibbs_samples is empty at ", path,
         " (postgibbsf90 may have failed; use source = \"gibbs\").")
  tab <- utils::read.table(path, header = FALSE)
  np  <- suppressWarnings(as.integer(tab[1, 3]))
  if(is.na(np) || np < 1 || ncol(tab) < 3 + np)
    stop("Unexpected postgibbs_samples layout in ", path, " (use source = \"gibbs\").")
  vc <- as.matrix(tab[, 4:(3 + np), drop = FALSE])
  colnames(vc) <- paste0("var_", seq_len(ncol(vc)))
  colnames(vc)[ncol(vc)] <- "residual_var"          # BLUPf90 writes the residual last
  list(values = vc, rounds = tab[[2]], desc = NULL)
}

# Build the named list of Manhattan reference-line heights from 'threshold'.
# Numeric values are drawn as-is (variance or p-value scale); the character options
# "bonferroni" / "fdr" / "meff" are computed on the -log10(p) scale and only valid
# for the per-SNP p-value plot. Returns a named list (name -> y height).
.manhattan_lines <- function(threshold, d, pv_mode, alpha, meff) {
  if(is.null(threshold)) return(list())
  if(isTRUE(threshold)) threshold <- "bonferroni"
  out <- list()
  if(is.numeric(threshold)){
    for(v in threshold) out[[sprintf("y = %.3g", v)]] <- v
    return(out)
  }
  if(is.character(threshold)){
    if(!pv_mode)
      stop("Named thresholds (\"bonferroni\"/\"fdr\"/\"meff\") apply only to statistic = ",
           "\"pvalue\", by = \"snp\". For a variance plot pass a numeric cutoff.")
    m <- nrow(d)
    for(t in tolower(threshold)){
      if(t == "bonferroni")      out[["Bonferroni"]]      <- -log10(alpha / m)
      else if(t == "meff"){
        if(is.null(meff)) stop("threshold = \"meff\" needs the effective number of tests via meff = .")
        out[["Meff Bonferroni"]] <- -log10(alpha / meff)
      } else if(t == "fdr"){
        pf <- .fdr_threshold(d$p, alpha)
        if(is.finite(pf)) out[["FDR (BH)"]] <- -log10(pf)
        else message("No SNP passes FDR at q = ", alpha, "; FDR line omitted.")
      } else stop("Unknown threshold: ", t)
    }
    return(out)
  }
  stop("'threshold' must be NULL, a number, or a character vector of ",
       "\"bonferroni\" / \"fdr\" / \"meff\".")
}
