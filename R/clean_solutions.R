#' Clean solutions
#'
#' This function cleans and formats BLUP solutions.
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
#'clean_solutions(3,
#'"my_clean_ebvs")
#'
#' @export
clean_solutions <- function(random_effect_col, solutions_output_name) {
  # read in the solutions file
  sols <- utils::read.table("solutions", header = TRUE, row.names = NULL) %>%
    dplyr::select(2, 3, 4) %>% # keep cols for effect number, id and solution
    base::subset(trait.effect == random_effect_col) # keep only the solutions for the random effect (ebvs)

  # read in the .inb file
  ids <- utils::read.table("renf90.inb", header = FALSE, row.names = NULL) %>%
    dplyr::select(1, 3) # keep only columns for original id and processed id

  # join ids and solutions
  final_output <- dplyr::left_join(ids, sols, by = c("V3" = "level")) %>%
    dplyr::select(1, 4) %>%
    dplyr::rename("ID" = "V1", "EBV" = "solution")

  # write output out
  utils::write.table(final_output, solutions_output_name, row.names = FALSE, quote = FALSE)
}
