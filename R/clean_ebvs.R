#' Clean ebvs
#'
#' This function formats BLUP solutions.
#'
#' This function cleans and formats the raw solutions file produced by blupf90+ (or gibbsf90+) by keeping the
#' solutions of the random (animal) effect (the ebvs) and matching the renumf90-processed ID back to the
#' original ID for each individual. For a single-trait model it writes two columns, ID and EBV. For a
#' multi-trait model (\code{multivariate} > 1) it keeps the trait dimension and writes one EBV column per
#' trait (ID, EBV_1, EBV_2, ...).
#'
#' @param random_effect_col Column where random effect is located, found under RANDOM_GROUP in the renf90.par file.
#' @param solutions_output_name name for the output file.  This field should be in quotes "".
#' @param multivariate number of traits in the model. Default 1 gives the single-trait output (ID, EBV). A value
#'   greater than 1 keeps every trait and writes one EBV column per trait (ID, EBV_1, ...); it is also checked
#'   against the number of traits actually present in the solutions file.
#' @param trait_names optional character vector of names for the per-trait EBV columns (multi-trait only),
#'   in trait order. Defaults to EBV_1, EBV_2, ...
#' @return a tab-separated file with the original id and ebv(s): two columns for a single trait, or ID plus one
#'   column per trait for a multi-trait model.
#' @import dplyr
#' @examples
#'
#' \donttest{
#'   #clean_ebvs(3, "my_clean_ebvs")                    # single trait -> ID, EBV
#'   #clean_ebvs(3, "my_clean_ebvs", multivariate = 3)  # 3 traits    -> ID, EBV_1, EBV_2, EBV_3
#' }
#'
#' @export
clean_ebvs <- function(random_effect_col, solutions_output_name, multivariate = 1, trait_names = NULL) {

  # Check if necessary files exist
  if (!file.exists("solutions")) {
    stop("Solutions file not found: solutions")
  }
  if (!file.exists("renf90.inb")) {
    stop("ID file not found: renf90.inb")
  }

  # Read in the solutions file positionally (skip the header line): 1 trait, 2 effect, 3 level (renumbered id),
  # 4 solution. Reading positionally also tolerates the extra SD column that gibbsf90+ writes in final_solutions.
  sols <- tryCatch({
    s <- utils::read.table("solutions", header = FALSE, skip = 1)
    colnames(s)[1:4] <- c("trait", "effect", "level", "solution")
    s[s$effect == random_effect_col, c("trait", "level", "solution"), drop = FALSE]   # keep only the random effect (ebvs)
  }, error = function(e) {
    stop("Error reading solutions file: ", e$message)
  })
  if (nrow(sols) == 0) stop("No solutions found for effect ", random_effect_col,
                            ". Check random_effect_col against RANDOM_GROUP in renf90.par.")

  # Read in the .inb file: column 1 = original id, column 3 = renumbered id (matches 'level')
  ids <- tryCatch({
    m <- utils::read.table("renf90.inb", header = FALSE, row.names = NULL)[, c(1, 3)]
    colnames(m) <- c("ID", "level")
    m
  }, error = function(e) {
    stop("Error reading ID file: ", e$message)
  })

  # Join ids and solutions
  final_output <- tryCatch({
    if (multivariate == 1) {                                           # single trait -> ID + EBV
      dplyr::left_join(ids, sols, by = "level") %>%
        dplyr::select("ID", "solution") %>%
        dplyr::rename("EBV" = "solution")
    } else {                                                           # multi-trait -> ID + one EBV column per trait
      traits <- base::sort(base::unique(sols$trait))
      if (length(traits) != multivariate)
        warning("multivariate = ", multivariate, " but ", length(traits), " trait(s) found in the solutions file.")
      if (is.null(trait_names)) trait_names <- base::paste0("EBV_", traits)
      out <- base::unique(ids["ID"])
      for (k in base::seq_along(traits)) {
        e <- dplyr::left_join(ids, sols[sols$trait == traits[k], c("level", "solution")], by = "level")
        e <- stats::setNames(e[, c("ID", "solution")], c("ID", trait_names[k]))
        out <- dplyr::left_join(out, e, by = "ID")
      }
      out
    }
  }, error = function(e) {
    stop("Error joining IDs and solutions: ", e$message)
  })

  # Write output
  tryCatch({
    utils::write.table(final_output, solutions_output_name, row.names = FALSE, quote = FALSE)
    cat("Output written to:", solutions_output_name, "\n")
  }, error = function(e) {
    stop("Error writing output file: ", e$message)
  })

  invisible(final_output)
}
