#' MCMC convergence diagnostics for a Gibbs variance-component run
#'
#' Runs post-Gibbs diagnostics on the MCMC samples from \code{\link{run_gibbs}} /
#' \code{\link{run_postgibbs}} using the \pkg{coda} package, following the CODA
#' post-Gibbs workflow of Vallejo et al. Beyond the diagnostic plots it returns a
#' decision: a per-parameter convergence table, an overall pass/fail flag, a
#' ready-to-use set of re-run settings for \code{\link{run_gibbs}}, and the
#' heritability posterior with a 95\% HPD interval.
#'
#' It builds a \code{coda} \code{mcmc} object from the sampled (co)variance
#' components (and, when the genetic effect is identified via \code{random_effect},
#' the derived heritability). For each parameter it reports the posterior mean, SD
#' and median, the 95\% HPD interval, the effective sample size (ESS), the Geweke z
#' and the Raftery-Lewis dependence factor, and a \code{status} of "converged",
#' "borderline" or "failed". A \code{warning} is emitted if any parameter has not
#' converged. It also returns suggested \code{run_gibbs()} settings (burn-in, total
#' iterations, thinning) from the Raftery-Lewis diagnostic, and the heritability
#' posterior with a 95\% HPD interval.
#'
#' By default it reads gibbsf90+'s \code{gibbs_samples} (available even when
#' postgibbsf90's diagnostic step fails and leaves \code{postgibbs_samples} empty);
#' set \code{source = "postgibbs"} to read \code{postgibbs_samples} instead. A
#' multi-page PDF (trace/density, autocorrelation, Geweke and normal-QQ) and a text
#' summary (mcmc_diagnostics.pdf / .txt) are written to \code{output_dir}.
#'
#' @param input_files_dir directory holding the sample file. Defaults to ".".
#' @param random_effect the effect number of the additive (animal) genetic effect,
#'   used to label the genetic variance and derive heritability. NULL skips the
#'   labelling / heritability (diagnostics still run on the raw components).
#' @param source which sample file to read: "gibbs" (default, gibbsf90+
#'   \code{gibbs_samples}) or "postgibbs" (postgibbsf90 \code{postgibbs_samples}).
#' @param genetic_col for \code{source = "postgibbs"} only (no effect labels in that
#'   file): which variance column is the genetic one. Defaults to 1.
#' @param output_dir directory for the PDF and text outputs. NULL uses \code{input_files_dir}.
#' @param verbose logical; if TRUE prints the verdict, suggested settings and paths.
#'
#' @return (invisibly) a list with: \code{table} (the per-parameter data frame:
#'   mean, sd, median, hpd_low, hpd_high, ess, geweke_z, raftery_I, status),
#'   \code{converged} (overall logical), \code{suggested} (a list of run_gibbs
#'   gibbs_iter / gibbs_burn / gibbs_keep, or NULL), \code{heritability} (mean,
#'   median, 95\% HPD, MCSE, SD; or NULL), the raw \code{mcmc} object, the
#'   \code{effective_size} / \code{geweke} / \code{raftery} diagnostics, and the
#'   \code{stats_file} / \code{plot_file} paths.
#' @references Vallejo RL et al. (2024) Aquaculture 586:740819. Plummer M et al.
#'   (2006) CODA: convergence diagnosis and output analysis for MCMC. R News 6:7-11.
#'   Geweke J (1992). Raftery AE, Lewis SM (1992).
#' @examples
#'
#' \donttest{
#'  # after run_gibbs(...) [and optionally run_postgibbs(...)]:
#'  # d <- mcmc_diagnostics(input_files_dir = "renumbered", random_effect = 3)
#'  # d$converged; d$table; d$suggested        # verdict, table, re-run recipe
#' }
#'
#' @export
mcmc_diagnostics <- function(input_files_dir = ".",
                             random_effect = NULL,
                             source = c("gibbs", "postgibbs"),
                             genetic_col = 1,
                             output_dir = NULL,
                             verbose = TRUE) {

  if(!requireNamespace("coda", quietly = TRUE))
    stop("Package 'coda' is required for mcmc_diagnostics(). Install it with install.packages(\"coda\").")
  source <- match.arg(source)
  if(is.null(output_dir)) output_dir <- input_files_dir
  ess_min <- 100                                          # ESS threshold for the verdict

  ## ---- read samples --------------------------------------------------------
  if(source == "gibbs"){
    g <- .read_gibbs_matrix(file.path(input_files_dir, "gibbs_samples"))
    V <- g$values
    if(!is.null(random_effect)){
      gi <- which(vapply(g$desc, function(d)
        length(d) >= 3 && !is.na(d[2]) && d[2] > 0 && d[2] == d[3] && d[2] == random_effect, logical(1)))
      if(length(gi)) colnames(V)[gi[1]] <- "genetic_var"
    }
  } else {
    g <- .read_postgibbs_matrix(file.path(input_files_dir, "postgibbs_samples"))
    V <- g$values
    if(genetic_col >= 1 && genetic_col <= ncol(V)) colnames(V)[genetic_col] <- "genetic_var"
  }
  rounds <- g$rounds
  if(nrow(V) < 10) stop("Too few samples (", nrow(V), ") for diagnostics.")

  ## ---- derived heritability + drop constant columns ------------------------
  is_var <- grepl("var", colnames(V)) & !grepl("^cov", colnames(V))
  if("genetic_var" %in% colnames(V) && sum(is_var) >= 2)
    V <- cbind(V, heritability = V[, "genetic_var"] / rowSums(V[, is_var, drop = FALSE]))
  const <- apply(V, 2, function(z) !is.finite(stats::sd(z)) || stats::sd(z) == 0)
  if(any(const)){
    if(verbose) base::message("mcmc_diagnostics: dropping constant parameter(s): ",
                              paste(colnames(V)[const], collapse = ", "))
    V <- V[, !const, drop = FALSE]
  }
  if(ncol(V) == 0) stop("No varying parameters left to diagnose.")

  ## ---- coda mcmc object ----------------------------------------------------
  thin  <- if(length(rounds) > 1) stats::median(diff(rounds)) else 1
  if(!is.finite(thin) || thin < 1) thin <- 1
  start <- if(length(rounds)) rounds[1] else 1
  S  <- coda::mcmc(V, start = start, thin = thin)
  nm <- colnames(V); np <- length(nm); n <- nrow(V)

  ## ---- diagnostics ---------------------------------------------------------
  smry <- summary(S)
  ess  <- coda::effectiveSize(S)
  gew  <- coda::geweke.diag(S)
  raf  <- try(coda::raftery.diag(S), silent = TRUE)
  st <- smry$statistics; if(is.null(dim(st))) { st <- t(as.matrix(st)); rownames(st) <- nm }
  qt <- smry$quantiles;  if(is.null(dim(qt))) { qt <- t(as.matrix(qt)); rownames(qt) <- nm }
  hpd <- coda::HPDinterval(S, prob = 0.95)
  gz  <- gew$z

  ## Raftery-Lewis burn-in (M), total (N), dependence factor (I)
  rI <- stats::setNames(rep(NA_real_, np), nm); rM <- NA_real_; rN <- NA_real_
  raf_ok <- !inherits(raf, "try-error") && !is.null(raf$resmatrix)
  if(raf_ok){
    rm <- raf$resmatrix; if(is.null(dim(rm))) rm <- matrix(rm, nrow = 1, dimnames = list(nm, names(rm)))
    rI[rownames(rm)] <- rm[, 4]; rM <- max(rm[, 1], na.rm = TRUE); rN <- max(rm[, 2], na.rm = TRUE)
  }

  ## ---- per-parameter verdict + table ---------------------------------------
  status <- vapply(seq_len(np), function(i){
    z <- abs(gz[nm[i]]); I <- rI[nm[i]]; e <- ess[nm[i]]
    if((is.finite(z) && z > 3) || (is.finite(I) && I > 5) || (is.finite(e) && e < ess_min/2)) "failed"
    else if(is.finite(z) && z < 2 && (!is.finite(I) || I < 5) && is.finite(e) && e >= ess_min) "converged"
    else "borderline"
  }, character(1))
  tab <- data.frame(parameter = nm,
                    mean = st[nm, "Mean"], sd = st[nm, "SD"], median = qt[nm, "50%"],
                    hpd_low = hpd[nm, 1], hpd_high = hpd[nm, 2],
                    ess = round(ess[nm], 1),
                    geweke_z = round(gz[nm], 2), raftery_I = round(rI[nm], 2),
                    status = status, row.names = NULL, stringsAsFactors = FALSE)
  converged <- all(status == "converged")

  ## ---- suggested run_gibbs re-run (Raftery-Lewis) --------------------------
  suggested <- NULL
  if(raf_ok && is.finite(rN))
    suggested <- list(gibbs_iter = as.integer(ceiling(rN)),
                      gibbs_burn = as.integer(ceiling(rM)),
                      gibbs_keep = as.integer(max(ceiling(max(rI, na.rm = TRUE)), as.numeric(thin))))

  ## ---- heritability posterior with a 95% HPD interval ----------------------
  heritability <- NULL
  if("heritability" %in% nm)
    heritability <- list(mean = st["heritability", "Mean"], median = qt["heritability", "50%"],
                         hpd_low = hpd["heritability", 1], hpd_high = hpd["heritability", 2],
                         mcse = st["heritability", "Time-series SE"], sd = st["heritability", "SD"])

  if(!converged)
    warning("MCMC not converged for: ",
            paste(tab$parameter[tab$status != "converged"], collapse = ", "),
            ". See $table and $suggested.", call. = FALSE)

  ## ---- text summary --------------------------------------------------------
  sp  <- file.path(output_dir, "mcmc_diagnostics.txt")
  con <- file(sp, "w")
  writeLines(c("BIGf90 MCMC diagnostics",
               sprintf("samples: %d | thin: %d | overall: %s",
                       n, thin, if(converged) "CONVERGED" else "NOT CONVERGED"),
               "", "== Per-parameter table =="), con)
  utils::capture.output(print(tab, row.names = FALSE), file = con)
  if(!is.null(suggested))
    writeLines(c("", "== Suggested run_gibbs() settings ==",
                 sprintf("gibbs_iter = %d ; gibbs_burn = %d ; gibbs_keep = %d",
                         suggested$gibbs_iter, suggested$gibbs_burn, suggested$gibbs_keep)), con)
  if(!is.null(heritability))
    writeLines(c("", "== Heritability ==",
                 sprintf("mean = %.4f ; median = %.4f ; 95%% HPD = [%.4f, %.4f]",
                         heritability$mean, heritability$median,
                         heritability$hpd_low, heritability$hpd_high)), con)
  close(con)

  ## ---- diagnostic plots (single multi-page PDF) ----------------------------
  pp <- file.path(output_dir, "mcmc_diagnostics.pdf")
  grDevices::pdf(pp, width = 8, height = 6)
  op <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(op); grDevices::dev.off() }, add = TRUE)
  graphics::plot(S)
  coda::autocorr.plot(S, auto.layout = TRUE)
  try(coda::geweke.plot(S, auto.layout = TRUE), silent = TRUE)
  nc <- ceiling(sqrt(np)); graphics::par(mfrow = c(ceiling(np / nc), nc), mar = c(4, 4, 2, 1))
  for(j in seq_len(np)){
    stats::qqnorm(V[, j], main = nm[j], pch = 19, cex = 0.4); stats::qqline(V[, j], col = "red")
  }

  if(verbose){
    base::message(sprintf("mcmc_diagnostics: %s (%d/%d parameters). ESS(min) = %.0f.",
                          if(converged) "CONVERGED" else "NOT CONVERGED",
                          sum(status == "converged"), np, min(ess)))
    if(!is.null(suggested))
      base::message(sprintf("  suggested run_gibbs: gibbs_iter=%d, gibbs_burn=%d, gibbs_keep=%d",
                            suggested$gibbs_iter, suggested$gibbs_burn, suggested$gibbs_keep))
    base::message("  wrote ", sp, " and ", pp)
  }

  invisible(list(table = tab, converged = converged, suggested = suggested,
                 heritability = heritability, mcmc = S, effective_size = ess,
                 geweke = gew, raftery = raf, stats_file = sp, plot_file = pp))
}
