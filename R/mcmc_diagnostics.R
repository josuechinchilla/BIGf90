#' MCMC convergence diagnostics for a Gibbs variance-component run
#'
#' Runs post-Gibbs diagnostics on the MCMC samples from \code{\link{run_gibbs}} /
#' \code{\link{run_postgibbs}} using the \pkg{coda} package. It returns a
#' per-parameter convergence table, an overall pass/fail flag, a ready-to-use set
#' of re-run settings for \code{\link{run_gibbs}}, and the heritability posterior
#' with a 95\% HPD interval.
#'
#' It builds a \code{coda} mcmc object from the sampled (co)variance components
#' (and, when the genetic effect is identified via random_effect, the derived
#' heritability). For each parameter it reports the posterior mean, SD and median,
#' the 95\% HPD interval, the effective sample size (ESS), the Geweke z and the
#' Raftery-Lewis dependence factor, and a \code{status} of "converged",
#' "borderline" or "failed". A \code{warning} is emitted if any parameter has not
#' converged. It also returns suggested \code{run_gibbs()} settings (burn-in, total
#' iterations, thinning) from the Raftery-Lewis diagnostic, and the heritability
#' posterior with a 95\% HPD interval.
#'
#' By default it reads gibbsf90+'s \code{gibbs_samples} (available even when
#' postgibbsf90's diagnostic step fails and leaves \code{postgibbs_samples} empty);
#' set \code{source = "postgibbs"} to read \code{postgibbs_samples} instead. A
#' multi-page PDF (trace/density, autocorrelation, Geweke and normal-QQ) is written
#' to \code{output_dir}; the normal-QQ plots are also rendered to the active device
#' (RStudio Plots pane when running interactively).
#'
#' @param input_files_dir directory holding the sample file. Defaults to ".".
#' @param random_effect the effect number of the additive (animal) genetic effect,
#'   used to label the genetic variance and derive heritability. NULL skips the
#'   labelling / heritability (diagnostics still run on the raw components).
#' @param source which sample file to read: "gibbs" (default, gibbsf90+
#'   \code{gibbs_samples}) or "postgibbs" (postgibbsf90 \code{postgibbs_samples}).
#' @param genetic_col for \code{source = "postgibbs"} only (no effect labels in that
#'   file): which variance column is the genetic one. Defaults to 1.
#' @param output_dir directory for the PDF and text outputs. NULL uses
#'   \code{input_files_dir}.
#' @param separate_plots logical (default FALSE); if TRUE, also write each diagnostic set
#'   (trace/density, autocorrelation, Geweke, normal-QQ) to its own file, in addition to the
#'   combined \code{mcmc_diagnostics.pdf}, for publication or other single-panel use.
#' @param plot_format file type for the per-set files when \code{separate_plots = TRUE}:
#'   "pdf" (default, vector) or "png" (300 dpi). The combined PDF is always a PDF.
#' @param verbose logical; if TRUE prints the verdict, suggested settings and paths.
#'
#' @return (invisibly) a list with: \code{table} (the per-parameter data frame:
#'   mean, sd, median, hpd_low, hpd_high, ess, geweke_z, raftery_I, status),
#'   \code{converged} (overall logical), \code{suggested} (a list of run_gibbs
#'   gibbs_iter / gibbs_burn / gibbs_keep, or NULL), \code{heritability} (mean,
#'   median, 95\% HPD, MCSE, SD; or NULL), the raw \code{mcmc} object, the
#'   \code{effective_size} / \code{geweke} / \code{raftery} diagnostics, the
#'   \code{stats_file} / \code{plot_file} paths, and \code{plot_files} (the per-set
#'   plot paths when \code{separate_plots = TRUE}, otherwise empty).
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
                             separate_plots = FALSE,
                             plot_format = c("pdf", "png"),
                             verbose = TRUE) {

  if (!requireNamespace("coda", quietly = TRUE))
    stop("Package 'coda' is required for mcmc_diagnostics(). Install it with install.packages(\"coda\").")
  source      <- match.arg(source)
  plot_format <- match.arg(plot_format)
  if (is.null(output_dir)) output_dir <- input_files_dir
  ess_min <- 100
  
  ## ---- read samples --------------------------------------------------------
  if (source == "gibbs") {
    g <- .read_gibbs_matrix(file.path(input_files_dir, "gibbs_samples"))
    V <- g$values
    if (!is.null(random_effect)) {
      gi <- which(vapply(g$desc, function(d)
        length(d) >= 3 && !is.na(d[2]) && d[2] > 0 && d[2] == d[3] && d[2] == random_effect, logical(1)))
      if (length(gi)) colnames(V)[gi[1]] <- "genetic_var"
    }
  } else {
    g <- .read_postgibbs_matrix(file.path(input_files_dir, "postgibbs_samples"))
    V <- g$values
    if (genetic_col >= 1 && genetic_col <= ncol(V)) colnames(V)[genetic_col] <- "genetic_var"
  }
  rounds <- g$rounds
  if (nrow(V) < 10) stop("Too few samples (", nrow(V), ") for diagnostics.")
  
  ## ---- derived heritability + drop constant columns ------------------------
  is_var <- grepl("var", colnames(V)) & !grepl("^cov", colnames(V))
  if ("genetic_var" %in% colnames(V) && sum(is_var) >= 2)
    V <- cbind(V, heritability = V[, "genetic_var"] / rowSums(V[, is_var, drop = FALSE]))
  const <- apply(V, 2, function(z) !is.finite(stats::sd(z)) || stats::sd(z) == 0)
  if (any(const)) {
    if (verbose) base::message("mcmc_diagnostics: dropping constant parameter(s): ",
                               paste(colnames(V)[const], collapse = ", "))
    V <- V[, !const, drop = FALSE]
  }
  if (ncol(V) == 0) stop("No varying parameters left to diagnose.")
  
  ## ---- coda mcmc object ----------------------------------------------------
  thin  <- if (length(rounds) > 1) stats::median(diff(rounds)) else 1
  if (!is.finite(thin) || thin < 1) thin <- 1
  start <- if (length(rounds)) rounds[1] else 1
  S  <- coda::mcmc(V, start = start, thin = thin)
  nm <- colnames(V); np <- length(nm); n <- nrow(V)
  
  ## ---- diagnostics ---------------------------------------------------------
  smry <- summary(S)
  ess  <- coda::effectiveSize(S)
  gew  <- coda::geweke.diag(S)
  raf  <- try(coda::raftery.diag(S), silent = TRUE)
  st <- smry$statistics; if (is.null(dim(st))) { st <- t(as.matrix(st)); rownames(st) <- nm }
  qt <- smry$quantiles;  if (is.null(dim(qt))) { qt <- t(as.matrix(qt)); rownames(qt) <- nm }
  hpd <- coda::HPDinterval(S, prob = 0.95)
  gz  <- gew$z
  
  ## Raftery-Lewis burn-in (M), total (N), dependence factor (I)
  rI <- stats::setNames(rep(NA_real_, np), nm); rM <- NA_real_; rN <- NA_real_
  # resmatrix is a character message ("you need a sample size of at least ...") when the chain
  # is too short for Raftery-Lewis; require a numeric matrix so that case is skipped, not coerced.
  raf_ok <- !inherits(raf, "try-error") && !is.null(raf$resmatrix) && is.numeric(raf$resmatrix)
  if (!raf_ok)
    warning("Not enough samples for the Raftery-Lewis diagnostic; it was dropped ",
            "(raftery_I reported as NA). All other diagnostics still ran.", call. = FALSE)
  if (raf_ok) {
    rm <- raf$resmatrix; if (is.null(dim(rm))) rm <- matrix(rm, nrow = 1, dimnames = list(nm[1], names(rm)))
    rI[rownames(rm)] <- rm[, 4]; rM <- max(rm[, 1], na.rm = TRUE); rN <- max(rm[, 2], na.rm = TRUE)
  }
  
  ## ---- per-parameter verdict + table ---------------------------------------
  status <- vapply(seq_len(np), function(i) {
    z <- abs(gz[nm[i]]); I <- rI[nm[i]]; e <- ess[nm[i]]
    if ((is.finite(z) && z > 3) || (is.finite(I) && I > 5) || (is.finite(e) && e < ess_min / 2)) "failed"
    else if (is.finite(z) && z < 2 && (!is.finite(I) || I < 5) && is.finite(e) && e >= ess_min) "converged"
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
  if (raf_ok && is.finite(rN))
    suggested <- list(gibbs_iter = as.integer(ceiling(rN)),
                      gibbs_burn = as.integer(ceiling(rM)),
                      gibbs_keep = as.integer(max(ceiling(max(rI, na.rm = TRUE)), as.numeric(thin))))
  
  ## ---- heritability posterior with a 95% HPD interval ----------------------
  heritability <- NULL
  if ("heritability" %in% nm)
    heritability <- list(mean = st["heritability", "Mean"], median = qt["heritability", "50%"],
                         hpd_low = hpd["heritability", 1], hpd_high = hpd["heritability", 2],
                         mcse = st["heritability", "Time-series SE"], sd = st["heritability", "SD"])
  
  if (!converged)
    warning("MCMC not converged for: ",
            paste(tab$parameter[tab$status != "converged"], collapse = ", "),
            ". See $table and $suggested.", call. = FALSE)
  
  ## ---- text summary --------------------------------------------------------
  sp  <- file.path(output_dir, "mcmc_diagnostics.txt")
  con <- file(sp, "w")
  writeLines(c("BIGf90 MCMC diagnostics",
               sprintf("samples: %d | thin: %d | overall: %s",
                       n, thin, if (converged) "CONVERGED" else "NOT CONVERGED"),
               "", "== Per-parameter table =="), con)
  utils::capture.output(print(tab, row.names = FALSE), file = con)
  writeLines(c("",
               "== What to look for (per parameter) ==",
               "  ess        effective sample size  (want >= 100; higher is better)",
               "  geweke_z   start-vs-end z-score    (want |z| < 2; |z| > 3 = drift / not stationary)",
               "  raftery_I  dependence factor       (want < 5; large = slow mixing / a variance stuck near 0)",
               "  status     per-parameter verdict   (want 'converged'; 'borderline'/'failed' = run longer)",
               "  mean, median, 95% HPD [low, high]  = posterior estimate and credible interval"), con)
  if (!is.null(suggested))
    writeLines(c("", "== Suggested run_gibbs() settings ==",
                 sprintf("gibbs_iter = %d ; gibbs_burn = %d ; gibbs_keep = %d",
                         suggested$gibbs_iter, suggested$gibbs_burn, suggested$gibbs_keep)), con)
  if (!is.null(heritability))
    writeLines(c("", "== Heritability ==",
                 sprintf("mean = %.4f ; median = %.4f ; 95%% HPD = [%.4f, %.4f]",
                         heritability$mean, heritability$median,
                         heritability$hpd_low, heritability$hpd_high)), con)
  close(con)
  
  ## ---- diagnostic plots ----------------------------------------------------
  ## Each set is a self-contained helper, so it can go into the combined PDF and
  ## (when separate_plots = TRUE) into its own file for publication.
  ncp   <- ceiling(sqrt(np)); nrp <- ceiling(np / ncp)
  b_sig <- 2 / sqrt(n)                                   # +/- white-noise significance bound (2/sqrt(N))

  plt_trace <- function() graphics::plot(S)             # coda trace + posterior density

  plt_autocorr <- function() {
    op <- graphics::par(mfrow = c(nrp, ncp), mar = c(4, 4, 2, 1), oma = c(4, 0, 0, 0))
    on.exit(graphics::par(op))
    for (j in seq_len(np)) {
      coda::autocorr.plot(S[, j, drop = FALSE], auto.layout = FALSE)
      graphics::abline(h = c(-0.1, 0.1),     lty = 3, col = "blue")   # practical: |acf| < 0.1
      graphics::abline(h = c(-b_sig, b_sig), lty = 2, col = "red")    # significance: +/- 2/sqrt(N)
    }
    ## legend once, in the device's bottom-left, below all panels
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    graphics::legend("bottomleft", bty = "n", cex = 0.7, lty = c(3, 2), col = c("blue", "red"),
                     legend = c("|acf| < 0.1 (practical)",
                                sprintf("2/sqrt(N) = %.3f (significance)", b_sig)))
  }

  plt_geweke <- function() {
    op <- graphics::par(mfrow = c(nrp, ncp), mar = c(4, 4, 2, 1), oma = c(4, 0, 0, 0))
    on.exit(graphics::par(op))
    for (j in seq_len(np)) {
      try(coda::geweke.plot(S[, j, drop = FALSE], auto.layout = FALSE), silent = TRUE)
      graphics::abline(h = c(-2, 2), lty = 3, col = "blue")           # recommended: |z| < 2
    }
    ## legend once, in the device's bottom-left, below all panels
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    graphics::legend("bottomleft", bty = "n", cex = 0.7, lty = 3, col = "blue",
                     legend = "|z| < 2 (recommended)")
  }

  plt_qq <- function() {
    op <- graphics::par(mfrow = c(nrp, ncp), mar = c(4, 4, 2, 1)); on.exit(graphics::par(op))
    for (j in seq_len(np)) {
      stats::qqnorm(V[, j], main = nm[j], pch = 19, cex = 0.4)
      stats::qqline(V[, j], col = "red")
    }
  }

  sets <- list(trace = plt_trace, autocorr = plt_autocorr, geweke = plt_geweke, qq = plt_qq)

  ## combined multi-page PDF (all sets)
  pp <- file.path(output_dir, "mcmc_diagnostics.pdf")
  grDevices::pdf(pp, width = 8, height = 6)
  for (f in sets) f()
  grDevices::dev.off()

  ## optional: one file per diagnostic set (publication / other use)
  plot_files <- character(0)
  if (separate_plots) {
    for (s in names(sets)) {
      f <- file.path(output_dir, sprintf("mcmc_%s.%s", s, plot_format))
      if (plot_format == "png") grDevices::png(f, width = 8, height = 6, units = "in", res = 300)
      else                      grDevices::pdf(f, width = 8, height = 6)
      sets[[s]]()
      grDevices::dev.off()
      plot_files <- c(plot_files, f)
    }
  }

  ## QQ also rendered to the active device (RStudio Plots pane)
  op_screen <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op_screen), add = TRUE)
  plt_qq()
  
  if (verbose) {
    base::message(sprintf("mcmc_diagnostics: %s (%d/%d parameters). ESS(min) = %.0f.",
                          if (converged) "CONVERGED" else "NOT CONVERGED",
                          sum(status == "converged"), np, min(ess)))
    base::message("  targets: ess >= 100 (higher better) | |geweke_z| < 2 | raftery_I < 5")
    if (!is.null(suggested))
      base::message(sprintf("  suggested run_gibbs: gibbs_iter=%d, gibbs_burn=%d, gibbs_keep=%d",
                            suggested$gibbs_iter, suggested$gibbs_burn, suggested$gibbs_keep))
    base::message("  wrote ", sp, " and ", pp)
    if (length(plot_files))
      base::message("  wrote ", length(plot_files), " per-set plot file(s): ",
                    paste(basename(plot_files), collapse = ", "))
  }
  
  invisible(list(table = tab, converged = converged, suggested = suggested,
                 heritability = heritability, mcmc = S, effective_size = ess,
                 geweke = gew, raftery = raf, stats_file = sp, plot_file = pp,
                 plot_files = plot_files))
}
