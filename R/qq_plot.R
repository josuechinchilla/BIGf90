#' QQ plot of single-step GWAS p-values
#'
#' Draws a quantile-quantile plot of the exact-GWAS p-values from a postGSf90
#' \code{snp_sol} file -- overall, or one panel per chromosome -- and reports the
#' genomic inflation factor \code{lambda}. p-values come from the SNP solution and its
#' variance (snp_sol columns \code{snp_effect} and \code{var_a_hat}): the test
#' statistic is \eqn{a_i^2 / \mathrm{var}(a_i) \sim \chi^2_1} (Aguilar et al. 2019).
#' This needs \code{OPTION snp_p_value}; set \code{snp_pvalue = TRUE} in
#' \code{\link{run_gwas}}, otherwise \code{var_a_hat} is zero and no p-values exist.
#'
#' \code{lambda} is the median observed \eqn{\chi^2} over its null expectation
#' (\code{qchisq(0.5, 1) = 0.4549}); lambda near 1 is well calibrated, > 1 indicates
#' inflation (structure / relatedness / artifacts) and < 1 deflation (common in
#' single-step models, where the relationship matrix already absorbs structure). With
#' \code{gc_correct = TRUE} a single-parameter genomic-control correction divides every
#' statistic by the genome-wide lambda before recomputing p.
#'
#' @param x one of: a directory with \code{snp_sol} (default "."); a path to a snp_sol
#'   file; the list returned by \code{\link{run_gwas}}; or an already-read snp_sol
#'   data frame.
#' @param snp_sol_file name of the snp_sol file when x is a directory. Defaults to "snp_sol".
#' @param chromosomes optional vector of chromosomes to keep (e.g. 1:28 to drop a
#'   rank-coded / unmapped chromosome). NULL keeps all.
#' @param per_chromosome logical; if TRUE draw one QQ panel per chromosome. Defaults to FALSE.
#' @param gc_correct logical (default FALSE); if TRUE apply genomic control -- divide
#'   every chi-square by the genome-wide lambda before recomputing p.
#' @param point_cex point size. Defaults to 0.5.
#' @param col point colour. Defaults to "grey30".
#' @param main plot title (overall plot only). NULL builds a default.
#' @param save_to optional output file (.png or .pdf); NULL draws to the active device.
#' @param width,height size in inches when saving. Default 7 x 7.
#' @param res resolution in ppi for a .png. Defaults to 300.
#'
#' @return (invisibly) a data frame of the genomic inflation factor: one row
#'   (\code{chr = "all"}) for the overall plot, or one row per chromosome when
#'   \code{per_chromosome = TRUE}, with a \code{lambda} column.
#' @references Aguilar I et al. (2019) Front Genet 10:442. Devlin B, Roeder K (1999)
#'   Biometrics 55:997-1004.
#' @examples
#'
#' \donttest{
#'  # g <- run_gwas(execs, "gwas.par", "gwas_run", iterations = 1, snp_pvalue = TRUE)
#'  # qq_plot(g)                                 # overall
#'  # qq_plot("gwas_run", per_chromosome = TRUE) # one panel per chromosome
#' }
#'
#' @export
qq_plot <- function(x = ".",
                    snp_sol_file = "snp_sol",
                    chromosomes = NULL,
                    per_chromosome = FALSE,
                    gc_correct = FALSE,
                    point_cex = 0.5,
                    col = "grey30",
                    main = NULL,
                    save_to = NULL,
                    width = 7, height = 7, res = 300) {

  d <- if(is.data.frame(x)) .name_snp_sol(x)
       else if(is.list(x) && !is.null(x$snp_sol)) x$snp_sol
       else .read_snp_sol(.resolve_path(x, snp_sol_file))

  d$chisq <- .snp_chisq(d)                                  # errors clearly if no p-values
  if(!is.null(chromosomes)) d <- d[d$chr %in% chromosomes, ]
  d <- d[is.finite(d$chisq) & is.finite(d$chr), ]
  if(nrow(d) == 0) stop("No usable statistics to plot.")

  lambda_all <- .lambda_gc(d$chisq)                         # genome-wide inflation factor
  if(gc_correct) d$chisq <- d$chisq / lambda_all            # genomic control (single parameter)

  one_qq <- function(cs, ttl){
    lam <- .lambda_gc(cs)
    p   <- sort(stats::pchisq(cs, df = 1, lower.tail = FALSE))
    obs <- -log10(p)
    exp <- -log10(stats::ppoints(length(p)))
    graphics::plot(exp, obs, pch = 19, cex = point_cex, col = col,
                   xlab = expression(Expected ~ -log[10](italic(p))),
                   ylab = expression(Observed ~ -log[10](italic(p))),
                   main = ttl, bty = "l")
    graphics::abline(0, 1, col = "red")
    graphics::legend("topleft", bty = "n", cex = 0.85,
                     legend = sprintf("lambda = %.3f", lam))
    lam
  }

  ## ---- device --------------------------------------------------------------
  if(!is.null(save_to)){
    if(grepl("\\.pdf$", save_to, ignore.case = TRUE))
      grDevices::pdf(save_to, width = width, height = height)
    else
      grDevices::png(save_to, width = width, height = height, units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if(per_chromosome){
    chrs <- sort(unique(d$chr))
    nc <- ceiling(sqrt(length(chrs))); nr <- ceiling(length(chrs) / nc)
    op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(op), add = TRUE)
    lam <- vapply(chrs, function(k) one_qq(d$chisq[d$chr == k], paste("Chr", k)), numeric(1))
    out <- data.frame(chr = as.character(chrs), lambda = lam, stringsAsFactors = FALSE)
  } else {
    if(is.null(main)) main <- sprintf("QQ plot (single-step GWAS%s)", if(gc_correct) ", GC-corrected" else "")
    out <- data.frame(chr = "all", lambda = one_qq(d$chisq, main), stringsAsFactors = FALSE)
  }
  invisible(out)
}
