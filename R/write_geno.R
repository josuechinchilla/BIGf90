#' Build a BLUPf90 genotype (.geno) file from a PLINK .ped
#'
#' Converts a PLINK \code{.ped} into the genotype file that the BLUPf90 programs
#' read, with each line is the animal ID, left-justified in a fixed-width field, then a
#' contiguous string of one 0/1/2 allele-dosage per SNP (missing written as
#' \code{missing_code}, 5 by default).
#'
#' By default the input is the \code{.ped} layout, with TWO allele codes
#' (1/2, and 0 for a missing allele) per locus; each pair is collapsed to a dosage
#' equal to the count of \code{count_allele} (so 1 1 -> 0, 1 2 and 2 1-> 1, 2 2 -> 2, and
#' any 0 -> missing). Set \code{alleles_per_locus = 1} if the genotypes are already
#' one 1/2 dosage per locus.
#'
#' The \code{.ped} must have NO header: a set of leading, non-genotype columns
#' (the six PLINK columns FID IID PAT MAT SEX PHENOTYPE, by default) followed by
#' the genotypes. The animal ID is taken from \code{id_col} (column 2, the PLINK
#' IID, by default) and must match the IDs in your data and pedigree files -- note
#' that PLINK output made from a VCF sometimes carries the usable ID in the FID
#' column (\code{id_col = 1}) instead, so check your file.
#'
#' @param ped path to the PLINK .ped file (no header, whitespace-separated).
#'   This field should be in quotes "".
#' @param file path/name of the .geno file to create. This field should be in quotes "".
#' @param map optional path to the matching PLINK .map file. When given it is used
#'   to check that the number of loci equals the number of markers (its number of
#'   rows) and to write a BLUPf90 marker map (see \code{map_out}).
#' @param map_out optional path to write a BLUPf90 marker map (snp chr pos), built
#'   from \code{map}. When \code{map} is provided and \code{map_out} is NULL, the
#'   map is written automatically to the same path as \code{file} with a .mapout
#'   extension. Set to NULL explicitly to skip only when \code{map} is also NULL.
#' @param id_col leading column holding the animal ID. Defaults to 2 (the PLINK IID).
#' @param n_lead_cols number of non-genotype leading columns before the genotypes.
#'   Defaults to 6 (the PLINK FID IID PAT MAT SEX PHENOTYPE columns).
#' @param alleles_per_locus 2 for the \code{--recode12} layout (two allele codes
#'   per locus, the default) or 1 for genotypes already coded as one 1/2 dosage.
#' @param count_allele allele code counted toward the dosage when
#'   \code{alleles_per_locus = 2}. Defaults to 2 (so 1 1 -> 0, 2 2 -> 2), matching
#'   the usual \code{--recode12} recoding.
#' @param missing value(s) that denote a missing genotype in addition to the
#'   built-ins (allele code 0 when \code{alleles_per_locus = 2}, and NA); these
#'   are written as \code{missing_code}. Defaults to NA.
#' @param missing_code code written for missing genotypes. Defaults to 5 (BLUPf90).
#' @param id_width fixed width for the left-justified ID field. Defaults to the
#'   longest ID + 1, so every row's dosage string starts at the same column.
#' @param overwrite logical; if FALSE (default) the function stops when \code{file}
#'   already exists.
#'
#' @return (invisibly) a character vector of the animal IDs written, in file order.
#' @examples
#' \dontrun{
#' # PLINK --recode12 ped (6 lead cols, two allele codes per locus, ID in col 2):
#' write_geno(ped  = "ped012.ped",
#'            file = "bf90_geno.txt",
#'            map  = "ped012.map")   # marker map written to bf90_geno.mapout
#' }
#'
#' @export
write_geno <- function(file = NULL,
                       ped = NULL,
                       map = NULL,
                       map_out = NULL,
                       id_col = 2,
                       n_lead_cols = 6,
                       alleles_per_locus = 2,
                       count_allele = 2,
                       missing = NA,
                       missing_code = 5,
                       id_width = NULL,
                       overwrite = FALSE) {
  if(is.null(ped))  stop("Define the input .ped file in 'ped'.")
  if(is.null(file)) stop("Define the output .geno file name in 'file'.")
  if(!base::file.exists(ped)) stop(paste("Input ped file", ped, "was not found."))
  if(!overwrite && base::file.exists(file))
    stop(paste("File", file, "already exists. Set overwrite = TRUE to replace it."))
  if(id_col > n_lead_cols)
    stop("'id_col' must be within the leading columns (id_col <= n_lead_cols).")
  if(!alleles_per_locus %in% c(1, 2))
    stop("'alleles_per_locus' must be 1 (dosage) or 2 (allele codes).")
  raw <- utils::read.table(ped, header = FALSE, colClasses = "character",
                           stringsAsFactors = FALSE, check.names = FALSE)
  if(ncol(raw) <= n_lead_cols)
    stop(paste0("The ped has only ", ncol(raw), " column(s); expected more than ",
                n_lead_cols, " (n_lead_cols) leading columns plus the genotypes."))
  ids <- raw[[id_col]]
  g   <- as.matrix(raw[, (n_lead_cols + 1):ncol(raw), drop = FALSE])
  miss_set <- as.character(missing)
  in_set <- function(m, s) matrix(m %in% s, nrow(m), ncol(m))
  if(alleles_per_locus == 2) {
    if(ncol(g) %% 2 != 0)
      stop(paste0("alleles_per_locus = 2 expects two allele columns per locus, but found ",
                  ncol(g), " genotype column(s) (odd). Check 'n_lead_cols' or use plink --recode12."))
    A <- g[, c(TRUE, FALSE), drop = FALSE]
    B <- g[, c(FALSE, TRUE), drop = FALSE]
    stray <- base::setdiff(base::unique(c(A, B)), c("0", "1", "2", miss_set))
    if(length(stray) > 0)
      stop(paste0("Expected PLINK --recode12 allele codes {1,2} (0 = missing); found: ",
                  paste(utils::head(stray, 5), collapse = ", "),
                  ". Use plink --recode12, or set alleles_per_locus = 1 for dosage input."))
    cnt  <- as.character(count_allele)
    dose <- (A == cnt) + (B == cnt)
    geno <- matrix(as.character(dose), nrow(dose), ncol(dose))
    miss <- in_set(A, c("0", miss_set)) | in_set(B, c("0", miss_set))
    geno[miss] <- as.character(missing_code)
  } else {
    geno <- g
    miss <- is.na(geno) | in_set(geno, miss_set)
    geno[miss] <- as.character(missing_code)
    stray <- base::setdiff(base::unique(base::as.vector(geno)), c("0", "1", "2", as.character(missing_code)))
    if(length(stray) > 0)
      stop(paste0("Found dosage code(s) not in {0,1,2,", missing_code, "}: ",
                  paste(utils::head(stray, 5), collapse = ", "),
                  ". For a two-allele PLINK .ped use alleles_per_locus = 2."))
  }
  n_snp <- ncol(geno)
  if(!is.null(map)) {
    if(!base::file.exists(map)) stop(paste("Map file", map, "was not found."))
    map_tab <- utils::read.table(map, header = FALSE, colClasses = "character",
                                 stringsAsFactors = FALSE, check.names = FALSE)
    if(nrow(map_tab) != n_snp)
      stop(paste0("Locus count mismatch: ", n_snp, " locus/loci in 'ped' but ",
                  nrow(map_tab), " marker(s) in 'map'. Check 'n_lead_cols'/'alleles_per_locus'."))
    if(is.null(map_out) || isTRUE(map_out)) map_out <- sub("\\.[^.]+$", ".mapout", file)
    map_lines <- paste(map_tab[[2]], map_tab[[1]], map_tab[[4]])
    base::writeLines(map_lines, con = map_out)
    message("write_geno: wrote marker map to '", map_out, "'.")
  }
  geno_str <- apply(geno, 1, paste0, collapse = "")
  if(is.null(id_width)) id_width <- max(nchar(ids)) + 1
  lines <- paste0(formatC(ids, flag = "-", width = id_width), geno_str)
  base::writeLines(lines, con = file)
  message("write_geno: wrote ", length(ids), " animals x ", n_snp,
          " SNPs to '", file, "'.")
  invisible(ids)
}