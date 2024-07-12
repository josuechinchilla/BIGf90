#' Clean ebvs
#'
#' This function formats BLUP solutions.
#'
#' This function cleans and formats the raw solutions file produced by blupf90+ by removing the solutions to all effects other than the random effect (ebvs).
#' I also matches the processed ID assigned by renumf90 to the original ID for each individual tested producing a interpretable output file with two columns: ID and EBV.
#'
#' @param random_effect_col Column where random effect is located, found under RANDOM_GROUP in the renf90.par file.
#' @param solutions_output_name name for the output file.  This field should be in quotes "".
#' @return a tab-separated file that includes the original id and ebv for all individuals for which an EBV was produced.
#' @import dplyr
#' @examples
#' ## Example
#'
#'clean_ebvs(3,
#'"my_clean_ebvs")
#'
#' @export
clean_ebvs <- function(random_effect_col, solutions_output_name) {

  # Check if necessary files exist
  if (!file.exists("solutions")) {
    stop("Solutions file not found: solutions")
  }
  if (!file.exists("renf90.inb")) {
    stop("ID file not found: renf90.inb")
  }

  # Read in the solutions file
  sols <- tryCatch({
    utils::read.table("solutions", header = TRUE, row.names = NULL) %>%
      dplyr::select(2, 3, 4) %>% # keep cols for effect number, id and solution
      base::subset(trait.effect == random_effect_col) # keep only the solutions for the random effect (ebvs)
  }, error = function(e) {
    stop("Error reading solutions file: ", e$message)
  })

  # Read in the .inb file
  ids <- tryCatch({
    utils::read.table("renf90.inb", header = FALSE, row.names = NULL) %>%
      dplyr::select(1, 3) # keep only columns for original id and processed id
  }, error = function(e) {
    stop("Error reading ID file: ", e$message)
  })

  # Join ids and solutions
  final_output <- tryCatch({
    dplyr::left_join(ids, sols, by = c("V3" = "level")) %>%
      dplyr::select(1, 4) %>%
      dplyr::rename("ID" = "V1", "EBV" = "solution")
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
}
