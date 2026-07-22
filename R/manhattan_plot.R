#' Manhattan plot from a single-step GWAS
#'
#' Draws a Manhattan plot directly from the postGSf90 outputs produced by
#' \code{\link{run_postgs}} / \code{\link{run_gwas}}. Two y-axis statistics are
#' supported: the proportion of additive genetic variance explained (the ssGBLUP /
#' wssGBLUP standard, e.g. Vallejo et al. 2024), available from any genomic run, or
#' \eqn{-\log_{10}(p)} from the exact-GWAS p-values, which need \code{OPTION snp_p_value}
#' (set \code{snp_pvalue = TRUE} in \code{\link{run_gwas}}). Points are drawn per SNP
#' (\code{by = "snp"}, from \code{snp_sol}) or per non-overlapping window
#' (\code{by = "window"}, from \code{windows_variance}).
#'
#' Multiple-testing lines (Bonferroni, Benjamini-Hochberg FDR, or a Bonferroni on the
#' effective number of tests) are only meaningful for the p-value plot -- variance
#' explained is not a test statistic, so there use a numeric variance cutoff instead.
#'
#' @param x one of: a directory holding the outputs (default "."); a path to a
#'   \code{snp_sol} or \code{windows_variance} file; the list returned by
#'   \code{\link{run_gwas}}; or an already-read snp_sol data frame.
#' @param statistic y axis: "variance" (proportion of genetic variance explained,
#'   default) or "pvalue" (\eqn{-\log_{10}(p)}; needs OPTION snp_p_value). Falls back to
#'   variance when \code{by = "window"}.
#' @param by "snp" for per-SNP points (default) or "window" for per-window points.
#' @param chromosomes optional vector of chromosomes to keep (e.g. 1:28 to drop a
#'   rank-coded / unmapped chromosome). NULL keeps all.
#' @param snp_sol_file name of the per-SNP file when x is a directory. Defaults to "snp_sol".
#' @param windows_file name of the window file when x is a directory and by = "window".
#'   Defaults to "windows_variance".
#' @param threshold reference line(s). NULL for none. A number (or numeric vector) draws
#'   horizontal line(s) at those y values -- use this for the variance plot (e.g.
#'   \code{threshold = 0.02}). For the p-value plot a character vector picks computed
#'   lines: "bonferroni" (\eqn{-\log_{10}(\alpha/m)}), "fdr" (Benjamini-Hochberg at level
#'   \code{alpha}) and/or "meff" (Bonferroni on the effective number of tests, needs
#'   \code{meff}). TRUE is shorthand for "bonferroni".
#' @param alpha significance level for the "bonferroni" / "fdr" / "meff" lines. Defaults to 0.05.
#' @param meff effective number of independent tests for threshold = "meff" (e.g. from
#'   GEC, simpleM, or Li and Ji). NULL disables the meff line.
#' @param gc_correct logical (default FALSE); if TRUE apply genomic control to the
#'   p-value plot -- divide every chi-square by the genomic inflation factor lambda
#'   before recomputing p, recentring the null. The lambda shown is pre-correction.
#' @param highlight optional y value above which points are drawn in \code{highlight_col};
#'   NULL highlights nothing.
#' @param highlight_col colour for highlighted points. Defaults to "red".
#' @param colors two alternating colours for odd / even chromosomes.
#' @param point_cex point size. Defaults to 0.5.
#' @param main plot title. NULL builds a default.
#' @param save_to optional output file (.png or .pdf); NULL draws to the active device.
#' @param width,height size in inches when saving. Default 10 x 4.
#' @param res resolution in ppi for a .png. Defaults to 300.
#'
#' @return (invisibly) the plotted data frame with columns \code{chr}, \code{pos}, the
#'   plotted statistic \code{y} and the cumulative x coordinate \code{x_cum}.
#' @references Vallejo RL et al. (2024) Aquaculture 586:740819. Aguilar I et al. (2019)
#'   Front Genet 10:442. Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
#'   Devlin B, Roeder K (1999) Biometrics 55:997-1004.
#' @examples
#'
#' \donttest{
#'  # g <- run_gwas(execs, "gwas.par", "gwas_run", iterations = 2, windows_mbp = 1)
#'  # manhattan_plot(g, by = "window", threshold = 0.02)        # variance windows
#'  # manhattan_plot("gwas_run", statistic = "pvalue",          # needs snp_pvalue = TRUE
#'  #                threshold = c("bonferroni", "fdr"))
#' }
#'
#' @export
manhattan_plot <- function(x = ".",
                           statistic = c("variance", "pvalue"),
                           by = c("snp", "window"),
                           chromosomes = NULL,
                           snp_sol_file = "snp_sol",
                           windows_file = "windows_variance",
                           threshold = NULL,
                           alpha = 0.05,
                           meff = NULL,
                           gc_correct = FALSE,
                           highlight = NULL,
                           highlight_col = "red",
                           colors = c("grey35", "dodgerblue3"),
                           point_cex = 0.5,
                           main = NULL,
                           save_to = NULL,
                           width = 10, height = 4, res = 300) {

  statistic <- match.arg(statistic)
  by        <- match.arg(by)
  pv_mode   <- (by == "snp" && statistic == "pvalue")

  ## ---- resolve the data (chr, pos, and either variance or snp_sol cols) -----
  if(by == "window"){
    wv_path <- if(is.list(x) && !is.data.frame(x)){
      if(is.null(x$windows_variance))
        stop("This run_gwas() result has no windows_variance; re-run with windows_mbp or windows_snp.")
      x$windows_variance
    } else .resolve_path(x, windows_file)
    d <- .read_windows(wv_path)
  } else {
    d <- if(is.data.frame(x)) .name_snp_sol(x)
         else if(is.list(x) && !is.null(x$snp_sol)) x$snp_sol
         else .read_snp_sol(.resolve_path(x, snp_sol_file))
  }
  if(!is.null(chromosomes)) d <- d[d$chr %in% chromosomes, ]

  ## ---- y axis ---------------------------------------------------------------
  lambda <- NA_real_
  if(pv_mode){
    chisq  <- .snp_chisq(d)                        # errors clearly if no p-values
    lambda <- .lambda_gc(chisq)
    if(gc_correct) chisq <- chisq / lambda         # genomic control: recentre null median to 1
    d$p    <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
    d$y    <- -log10(d$p)
    ylab   <- expression(-log[10](italic(p)))
  } else {
    if(statistic == "pvalue")
      message("statistic = \"pvalue\" needs by = \"snp\"; plotting window variance instead.")
    d$y  <- d$variance_explained
    ylab <- "Proportion of genetic variance explained"
  }

  ## ---- order genomically and lay out a cumulative x axis --------------------
  d <- d[is.finite(d$chr) & is.finite(d$pos) & is.finite(d$y), ]
  d <- d[order(d$chr, d$pos), ]
  if(nrow(d) == 0) stop("Nothing to plot: no finite chr / pos / statistic values.")
  chrs   <- sort(unique(d$chr))
  gap    <- 0.02 * max(d$pos)
  offset <- 0
  centers <- numeric(length(chrs))
  d$x_cum <- NA_real_
  for(k in seq_along(chrs)){
    idx <- d$chr == chrs[k]
    cp  <- d$pos[idx]
    d$x_cum[idx] <- cp - min(cp) + offset
    centers[k]   <- offset + (max(cp) - min(cp)) / 2
    offset       <- offset + (max(cp) - min(cp)) + gap
  }
  pt_col <- ifelse(match(d$chr, chrs) %% 2 == 1, colors[1], colors[2])
  if(!is.null(highlight)) pt_col[d$y >= highlight] <- highlight_col

  ## ---- reference / significance lines --------------------------------------
  hlines <- .manhattan_lines(threshold, d, pv_mode, alpha, meff)

  ## ---- device + draw --------------------------------------------------------
  if(!is.null(save_to)){
    if(grepl("\\.pdf$", save_to, ignore.case = TRUE))
      grDevices::pdf(save_to, width = width, height = height)
    else
      grDevices::png(save_to, width = width, height = height, units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if(is.null(main))
    main <- if(pv_mode) sprintf("Manhattan plot (-log10 p%s)", if(gc_correct) ", GC-corrected" else "")
            else sprintf("Manhattan plot (variance explained, per %s)", by)

  graphics::plot(d$x_cum, d$y, col = pt_col, pch = 19, cex = point_cex,
                 xaxt = "n", xlab = "Chromosome", ylab = ylab, main = main, bty = "l")
  graphics::axis(1, at = centers, labels = chrs, las = 2, cex.axis = 0.7, tick = FALSE)
  if(length(hlines)){
    lc <- c("red", "blue", "darkgreen", "purple")[seq_along(hlines)]
    for(j in seq_along(hlines)) graphics::abline(h = hlines[[j]], col = lc[j], lty = 2)
    graphics::legend("topright", bty = "n", cex = 0.75, lty = 2, col = lc,
                     legend = sprintf("%s (%.2f)", names(hlines), unlist(hlines)))
  }
  if(pv_mode && is.finite(lambda))
    graphics::mtext(sprintf("lambda = %.3f%s", lambda, if(gc_correct) " (GC applied)" else ""),
                    side = 3, adj = 1, line = 0.1, cex = 0.7)

  invisible(d[, c("chr", "pos", "y", "x_cum")])
}
