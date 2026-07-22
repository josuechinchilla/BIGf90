#' Write a RENUMF90 parameter (.par) file
#'
#' This function writes the raw parameter (.par) file to be used with RENUMf90.
#'
#' This function assembles the raw parameter file that you feed to \code{\link{run_renum}}
#' , NOT the renf90.par file that RENUMF90 produces. 
#' It is built from the standard RENUMF90 blocks in order:
#' DATAFILE, TRAITS, FIELDS_PASSED TO OUTPUT, WEIGHT(S), RESIDUAL_VARIANCE, the EFFECT
#' blocks, the random-effect block, an optional SNP_FILE, an optional (CO)VARIANCES block
#' and any OPTION lines.
#'
#' Effects are supplied as a named list. Each list name is the effect's label, written as the
#' "#..." comment on its EFFECT line (a per-effect \code{comment} overrides it). Inside each
#' element give the data-file column with
#' \code{col} (\code{pos} is accepted as an alias), the \code{type} ("cross" or "cov") and, for
#' cross-classified effects, the \code{class} ("numer" or "alpha"). For multi-trait models \code{col}
#' may be a vector with one column per trait; a single value is repeated across traits (so two
#' traits give, e.g., \code{2 2 cross alpha}).
#'
#' Mark the random effect with \code{random}, giving its column number (e.g. \code{random = 1})
#' or its label. Because RENUMF90 requires RANDOM to directly follow its effect, the random effect
#' is always written last (it is moved to the end automatically, with a message, if you list it
#' earlier). RANDOM (with \code{random_type}) and an optional OPTIONAL line (e.g.
#' \code{optional = "pe"} for a permanent-environment / repeatability term) belong to that effect
#' block; they are followed by the separate pedigree block FILE / FILE_POS, then an optional
#' SNP_FILE and an optional (CO)VARIANCES block. The order of OPTION lines is irrelevant to RENUMF90.
#'
#' @param file path/name of the .par file to create. This field should be in quotes "".
#' @param datafile name of the data file written under DATAFILE. This field should be in quotes "".
#' @param traits column position(s) of the trait(s). For multi-trait give a vector, e.g. c(6, 7).
#' @param effects named list of effects, in model order. Each list name is the effect's label
#'   (written as "EFFECT #name"). Inside each
#'   element supply \code{col} (the effect's data-file column; \code{pos} is an alias, and a vector
#'   gives one column per trait), \code{type} ("cross"/"cov") and, for cross effects, \code{class}
#'   ("numer"/"alpha"). An optional \code{comment} overrides the "#..." label.
#' @param random which effect is the random effect, given as its column number (e.g.
#'   \code{random = 1}) or its label. Matching is by \code{col}, so the effect's "#..." label can
#'   be anything you like. This effect is written last (RENUMF90 requires RANDOM to follow it),
#'   with RANDOM immediately after it. Leave NULL for a fixed-effects-only model.
#' @param random_effect deprecated alias for \code{random} (name-based); use \code{random}.
#' @param pedigree_file pedigree file written under FILE. Required when a random effect is set.
#' @param pedigree_pos integer vector written under FILE_POS (animal, sire, dam, then two 0s).
#'   Defaults to c(1, 2, 3, 0, 0).
#' @param random_type keyword written under RANDOM. Defaults to "animal".
#' @param optional keyword(s) written under an OPTIONAL line after RANDOM (e.g. "pe" for a
#'   permanent-environment / repeatability term). NULL to omit.
#' @param residual_variance starting residual (co)variance written under RESIDUAL_VARIANCE. A
#'   scalar for one trait, or a matrix for multi-trait (one row per line).
#' @param covariances starting (co)variance for the random effect, written under (CO)VARIANCES.
#'   Scalar or matrix. NULL to omit.
#' @param snp_file genotype file written under SNP_FILE (inside the random-effect block). Leave
#'   NULL for a pedigree-only model.
#' @param missing_value value written as "OPTION missing <value>". Leave NULL to omit.
#' @param options character vector of extra OPTION lines given WITHOUT the leading "OPTION "
#'   (e.g. c("callrate 0.9", "cat 0 0")).
#' @param fields_passed column position(s) written under FIELDS_PASSED TO OUTPUT (NULL = blank).
#' @param weights column position written under WEIGHT(S) (NULL = blank).
#' @param overwrite logical; if FALSE (default) the function stops when \code{file} already exists.
#'
#' @return (invisibly) the character vector of lines written to \code{file}.
#' @examples
#'
#' \donttest{
#'  # two-trait repeatability model: animal + permanent environment, genomic
#'  # write_par(file     = "act.par",
#'  #           datafile = "act_act_auc.csv",
#'  #           traits   = c(6, 7),
#'  #           residual_variance = diag(2),
#'  #           effects  = list(lact_group         = list(col = 2, class = "alpha"),
#'  #                           breed_group        = list(col = 3, class = "alpha"),
#'  #                           DIM_group          = list(col = 4, class = "alpha"),
#'  #                           contemporary_group = list(col = 5, class = "alpha"),
#'  #                           animal             = list(col = 1, class = "alpha")),
#'  #           random        = 1,          # the animal effect is column 1 (label can be anything)
#'  #           optional      = "pe",
#'  #           pedigree_file = "ponderosa.ped",
#'  #           snp_file      = "ponderosa_bulls.geno",
#'  #           covariances   = matrix(c(1, 0.1, 0.1, 1), 2, 2),
#'  #           options       = "cat 0 0")
#' }
#'
#' @export
write_par <- function(file = NULL,
                      datafile = NULL,
                      traits = NULL,
                      effects = NULL,
                      random = NULL,
                      random_effect = NULL,
                      pedigree_file = NULL,
                      pedigree_pos = c(1, 2, 3, 0, 0),
                      random_type = "animal",
                      optional = NULL,
                      residual_variance = NULL,
                      covariances = NULL,
                      snp_file = NULL,
                      missing_value = NULL,
                      options = NULL,
                      fields_passed = NULL,
                      weights = NULL,
                      overwrite = FALSE) {

  # Checks
  if(is.null(file))              stop("Define the output .par file name in 'file'.")
  if(is.null(datafile))          stop("Define the data file name in 'datafile'.")
  if(is.null(traits))            stop("Define the trait column position(s) in 'traits'.")
  if(is.null(residual_variance)) stop("Define a starting value in 'residual_variance'.")
  if(is.null(effects) || !is.list(effects) || length(effects) == 0)
    stop("'effects' must be a non-empty named list (see ?write_par).")
  if(is.null(names(effects)) || any(names(effects) == ""))
    stop("Every element of 'effects' must be named; the names are used as the effect labels.")
  rnd <- if(!is.null(random)) random else random_effect      # random effect: its col number or its label
  if(!is.null(rnd) && is.null(pedigree_file))
    stop("'pedigree_file' is required when a random effect is set (see 'random').")
  if(!overwrite && file.exists(file))
    stop(paste("File", file, "already exists. Set overwrite = TRUE to replace it."))

  n_traits <- length(traits)

  ## Format a scalar/vector/matrix (co)variance as one line per matrix row
  var_lines <- function(x) {
    if(is.matrix(x)) apply(x, 1, function(r) paste(r, collapse = " ")) else paste(x, collapse = " ")
  }

  ## Format one EFFECT definition line: one column per trait, then cross/cov [class]
  effect_line <- function(e, cols, nm) {
    type <- if(is.null(e[["type"]])) "cross" else e[["type"]]
    if(!type %in% c("cross", "cov")) stop(paste0("Effect '", nm, "': 'type' must be 'cross' or 'cov'."))
    posdef <- paste(cols, collapse = " ")
    if(type == "cov") return(paste(posdef, "cov"))
    cls <- if(is.null(e[["class"]])) "numer" else e[["class"]]
    if(!cls %in% c("numer", "alpha")) stop(paste0("Effect '", nm, "': 'class' must be 'numer' or 'alpha'."))
    paste(posdef, "cross", cls)
  }

  # Assemble the file, block by block, in RENUMF90 order
  lines <- c("DATAFILE", datafile,
             "TRAITS",   paste(traits, collapse = " "),
             "FIELDS_PASSED TO OUTPUT", if(is.null(fields_passed)) "" else paste(fields_passed, collapse = " "),
             "WEIGHT(S)",              if(is.null(weights))       "" else paste(weights, collapse = " "),
             "RESIDUAL_VARIANCE",      var_lines(residual_variance))

  # Resolve each effect's data-file column(s) and find the random effect. 'col' (or 'pos' alias)
  # may be a vector, one column per trait; a single value is repeated across traits.
  cols   <- list()
  rnd_nm <- NULL
  for(nm in names(effects)) {
    e   <- base::as.list(effects[[nm]])
    col <- if(!is.null(e[["col"]])) e[["col"]] else e[["pos"]]                  # data-file column ('pos' = alias)
    if(is.null(col))
      stop(paste0("Effect '", nm, "': supply its data-file column with col = (e.g. col = 1)."))
    col <- suppressWarnings(as.integer(col))
    if(any(is.na(col)) || any(col < 1))
      stop(paste0("Effect '", nm, "': 'col' must be a positive integer (one per trait)."))
    if(length(col) == 1 && n_traits > 1) col <- rep(col, n_traits)             # same column across traits
    if(length(col) != n_traits)
      stop(paste0("Effect '", nm, "': give one column (recycled) or one per trait (", n_traits, ")."))
    cols[[nm]] <- col
    # match by col number (label-independent) or by name (alias); keep the first match
    if(!is.null(rnd) && is.null(rnd_nm) &&
       ((is.numeric(rnd) && rnd %in% col) || identical(as.character(rnd), nm))) rnd_nm <- nm
  }
  if(!is.null(rnd) && is.null(rnd_nm))
    stop(paste0("No effect matches 'random' = ", rnd,
                ". Give the col number (or label) of the random effect."))

  # RENUMF90 requires RANDOM to directly follow its effect, so the random effect is written last.
  ord <- names(effects)
  if(!is.null(rnd_nm)) {
    if(ord[length(ord)] != rnd_nm)
      message("write_par: moving random effect '", rnd_nm, "' to the end (RANDOM must follow it).")
    ord <- c(setdiff(ord, rnd_nm), rnd_nm)
  }

  # EFFECT blocks in order; RANDOM (+ OPTIONAL) belong to the random effect, then the separate
  # pedigree block FILE / FILE_POS, then SNP_FILE and (CO)VARIANCES.
  snp_done <- FALSE
  for(nm in ord) {
    e      <- base::as.list(effects[[nm]])
    cmt    <- if(is.null(e[["comment"]])) nm else e[["comment"]]                # "#..." label (defaults to the name)
    header <- if(is.na(cmt) || !nzchar(cmt)) "EFFECT" else paste0("EFFECT #", cmt)
    lines  <- c(lines, header, effect_line(e, cols[[nm]], nm))

    if(identical(nm, rnd_nm)) {
      lines <- c(lines, "RANDOM", random_type)                                 # RANDOM (+ OPTIONAL): part of the effect block
      if(!is.null(optional)) lines <- c(lines, "OPTIONAL", paste(optional, collapse = " "))
      lines <- c(lines, "FILE", pedigree_file, "FILE_POS", paste(pedigree_pos, collapse = " "))  # separate pedigree block
      if(!is.null(snp_file))    { lines <- c(lines, "SNP_FILE", snp_file); snp_done <- TRUE }
      if(!is.null(covariances))   lines <- c(lines, "(CO)VARIANCES", var_lines(covariances))
    }
  }

  # SNP_FILE for a genomic model with no random pedigree effect (rare)
  if(!is.null(snp_file) && !snp_done) lines <- c(lines, "SNP_FILE", snp_file)

  # OPTION lines (order is irrelevant to RENUMF90)
  if(!is.null(missing_value)) lines <- c(lines, paste("OPTION missing", missing_value))
  if(!is.null(options))       lines <- c(lines, paste("OPTION", options))

  base::writeLines(lines, con = file)
  invisible(lines)
}
