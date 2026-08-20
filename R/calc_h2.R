#' Estimate heritability from Gibbs variance components
#'
#' Computes narrow-sense heritability from a Gibbs variance component run. By default it
#' reads the posterior-mean variance components from postgibbsf90's \code{postmean} file
#' (\code{\link{run_postgibbs}}). If \code{from_postmean = FALSE} it instead reads the
#' raw variance-component samples written by gibbsf90+ (\code{gibbs_samples},
#' \code{\link{run_gibbs}}) and returns the posterior mean of the per-sample heritability
#' -- which also avoids postgibbsf90's convergence-diagnostic step, that can crash on some
#' chains and leave \code{postmean} empty. Heritability is \eqn{h^2 = \sigma^2_a /
#' \sigma^2_p}, with \eqn{\sigma^2_a} the additive genetic variance of the animal effect
#' and \eqn{\sigma^2_p} the sum of all variance components. This is for single-trait models.
#'
#' @param random_effect the effect number of the additive (animal) genetic effect (the
#'   same number used for random_effect_col elsewhere).
#' @param input_files_dir directory holding the postmean / gibbs_samples file. Defaults to ".".
#' @param from_postmean logical (default TRUE); if TRUE read postgibbsf90's posterior
#'   means (postmean), if FALSE read the gibbsf90+ samples (gibbs_samples) directly.
#' @param postmean_file name of the posterior-mean file. Defaults to "postmean".
#' @param samples_file name of the Gibbs sample file. Defaults to "gibbs_samples".
#' @param verbose logical; if TRUE prints the variance components and h2.
#'
#' @return (invisibly) a list with \code{h2}, \code{h2_sd} (posterior SD of h2, only when
#'   reading samples; NA otherwise), \code{sigma_a}, \code{sigma_e}, \code{sigma_p} and
#'   \code{n_samples}.
#' @examples
#' \dontrun{
#' # after run_gibbs(...) and run_postgibbs(...):
#' h <- calc_h2(random_effect = 3, input_files_dir = "renumbered")           # from postmean
#' h <- calc_h2(random_effect = 3, input_files_dir = "renumbered",
#'              from_postmean = FALSE)                                        # from the samples
#' }
#'
#' @export
calc_h2 <- function(random_effect,
                    input_files_dir = ".",
                    from_postmean = TRUE,
                    postmean_file = "postmean",
                    samples_file = "gibbs_samples",
                    verbose = TRUE) {

  if(missing(random_effect))
    stop("Provide the effect number of the animal (genetic) effect in 'random_effect'.")

  if(from_postmean){
    ## ---- posterior means from postgibbsf90 'postmean' -----------------------
    f <- file.path(input_files_dir, postmean_file)
    if(!file.exists(f))
      stop("Posterior-mean file not found: ", f,
           " (run run_postgibbs() first, or set from_postmean = FALSE to read the Gibbs samples).")
    lines <- readLines(f)
    vars <- list(); cur <- NA_character_
    for(ln in lines){
      if(grepl("G matrix for effect", ln)){
        cur <- paste0("effect_", trimws(sub(".*=", "", ln))); vars[[cur]] <- numeric(0)
      } else if(grepl("R matrix", ln, ignore.case = TRUE)){
        cur <- "residual"; vars[[cur]] <- numeric(0)
      } else if(!is.na(cur)){
        v <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]])); v <- v[!is.na(v)]
        if(length(v)) vars[[cur]] <- c(vars[[cur]], v)
      }
    }
    if(length(vars) == 0)
      stop("Could not parse any variance components from ", f,
           " (it may be empty if postgibbsf90 failed). Set from_postmean = FALSE to read the Gibbs samples.")
    gen_key <- paste0("effect_", random_effect)
    if(is.null(vars[[gen_key]]))
      stop("No 'G matrix for effect = ", random_effect, "' in ", f, ". Found: ",
           paste(names(vars), collapse = ", "), ".")
    first   <- vapply(vars, function(x) x[1], numeric(1))
    sigma_a <- unname(first[gen_key])
    sigma_e <- if("residual" %in% names(first)) unname(first["residual"]) else NA_real_
    sigma_p <- sum(first)
    out <- list(h2 = sigma_a / sigma_p, h2_sd = NA_real_, sigma_a = sigma_a,
                sigma_e = sigma_e, sigma_p = sigma_p, n_samples = NA_integer_)

  } else {
    ## ---- posterior-mean h2 from gibbsf90+ 'gibbs_samples' -------------------
    f <- file.path(input_files_dir, samples_file)
    if(!file.exists(f)) stop("Gibbs sample file not found: ", f, " (run run_gibbs() first).")
    lines <- readLines(f)
    if(length(lines) < 5) stop("Gibbs sample file '", f, "' has no samples.")
    hdr1 <- suppressWarnings(as.integer(strsplit(trimws(lines[1]), "\\s+")[[1]])); np <- hdr1[2]
    if(is.na(np) || np < 1) stop("Could not read the parameter count from ", f, ".")
    desc <- lapply(lines[2:(1 + np)], function(x) suppressWarnings(as.integer(strsplit(trimws(x), "\\s+")[[1]])))
    body      <- lines[-(1:(1 + np))]
    val_lines <- body[seq(2, length(body), by = 2)]
    vals <- do.call(rbind, lapply(strsplit(trimws(val_lines), "\\s+"),
                                  function(x) suppressWarnings(as.numeric(x))))
    vals <- vals[stats::complete.cases(vals), , drop = FALSE]
    if(nrow(vals) == 0 || ncol(vals) != np)
      stop("Could not parse the variance-component samples from ", f, ".")
    gen_col <- which(vapply(desc, function(d) length(d) >= 3 && d[2] == random_effect && d[3] == random_effect, logical(1)))
    res_col <- which(vapply(desc, function(d) length(d) >= 2 && d[2] == 0, logical(1)))
    if(length(gen_col) == 0)
      stop("No variance parameter for effect ", random_effect, " in ", f, ".")
    gen_col <- gen_col[1]
    sigma_p_s <- rowSums(vals); h2_s <- vals[, gen_col] / sigma_p_s
    out <- list(h2 = mean(h2_s), h2_sd = stats::sd(h2_s), sigma_a = mean(vals[, gen_col]),
                sigma_e = if(length(res_col)) mean(vals[, res_col[1]]) else NA_real_,
                sigma_p = mean(sigma_p_s), n_samples = nrow(vals))
  }

  if(verbose){
    src <- if(from_postmean) "postmean (posterior means)" else sprintf("%d Gibbs samples", out$n_samples)
    cat(sprintf("Variance components from %s:\n", src))
    cat(sprintf("  sigma_a (effect %s) = %g\n", random_effect, out$sigma_a))
    cat(sprintf("  sigma_e             = %g\n", out$sigma_e))
    cat(sprintf("  sigma_p             = %g\n", out$sigma_p))
    if(is.na(out$h2_sd)) cat(sprintf("h2 = %g\n", out$h2))
    else                 cat(sprintf("h2 = %g  (posterior SD %g)\n", out$h2, out$h2_sd))
  }
  invisible(out)
}
