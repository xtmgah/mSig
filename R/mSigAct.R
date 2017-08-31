#
# mSigAct.R
#
# v 0.9
#
# 2017 08 11
#
# Copyright 2017 by Alvin Wei Tian Ng, Steven G. Rozen
#
# The code is released under GPL-3
# https://www.gnu.org/licenses/gpl-3.0.en.html
#
# This file contains functions for sparse maximum likelihood assignment of
# mutational signature activities and for testing whether a given signature is
# needed to explain a spectrum


#' Helper function, given signatures (sigs) and exposures (exp), return a
#' *proportional* reconstruction; in general, it is *not necessarily* scaled to
#' the actual spectrum counts.
#'
#' @param sigs signatures
#' @param exp  exposures to those signatures
#'
#' @return a matrix of proportions of reconstructed mutation counts
prop.reconstruct <- function(sigs, exp) {
  stopifnot(length(exp) == ncol(sigs))
  as.matrix(sigs) %*% exp
}


#' Define the objective function for use with nloptr (non-linear optimization)
#' Negative binomial maximum likelihood objective function
#' (nloptr minimizes the objective function.)
#' The lower the objective function, the better
#'
#' @param exp         matrix of exposures ("activities")
#' @param sigs        matrix of signatures
#' @param spectrum    spectrum to assess
#' @param nbinom.size dispersion parameter for the negative binomial
#'                    distribution: smaller is more dispersed
#'
#' @return -1 * the log(likelihood(spectrum | reconstruction))
#'
#' @importFrom stats dnbinom
obj.fun.nbinom.maxlh <- function(exp,
                                 spectrum,
                                 sigs,
                                 nbinom.size) {

  if (any(is.na(exp))) return(Inf)

  reconstruction <-  prop.reconstruct(sigs = sigs, exp = exp)

  ## catch errors with NA in the input or in reconstruction.
  if (any(is.na(reconstruction))) {
    save(reconstruction, spectrum, sigs, exp, nbinom.size,
         file='reconstruction.error.R')
  }
  stopifnot(!any(is.na(reconstruction)))

  loglh <- 0
  for (i in 1:length(spectrum)) { # Iterate over each channel in the
                                  # spectrum and sum the log
                                  # likelihoods.

    nbinom <- dnbinom(x=spectrum[i],
                      mu = reconstruction[i],
                      size=nbinom.size,
                      log=T)
    loglh <-loglh + nbinom
  }
  stopifnot(mode(loglh) == 'numeric' )
  -loglh
}


#' Use nloptr (numerical non-linear optimization) to find an assignmeent of
#' signature activites for one tumor. The nlpotr algorithm and the objective
#' function are arguments.
#'
#' @param spectrum    matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs        signatures used to reconstruct the tumor
#' @param algorithm   algorithm for nloptr, default is \code{NLOPT_LN_COBYLA}
#' @param maxeval     nloptr function stops when \code{maxeval} is reached
#' @param xtol_rel    stopping criterion for relative change reached (see ?nloptr)
#' @param obj.fun     objective function for evaluation
#' @param logger      logger object
#' @param ...         additional arguments passed to the function
#' @return nloptr object with the objective value, number of iterations and optimal solution
#'
#' @importFrom nloptr nloptr
nloptr.one.tumor <- function(spectrum,
                             sigs,
                             algorithm='NLOPT_LN_COBYLA', # ToDo: would it make sense to make this customisable?
                             maxeval=1000,
                             xtol_rel=0.001, # 0.0001,
                             obj.fun,
                             logger,
                             ... # additional arguments for obj.fun
                             ) {

  # In this case we got a numeric vector, not a matrix.
  # We assume the caller intended it as a single
  # column matrix.
  if (class(sigs) == 'numeric')
    sigs <- matrix(sigs, ncol = 1)

  # The if statement is an example of captuing state with
  # a call to save, which seems useful for debugging
  # in a call to mclapply (multi-core lapply).
  if (nrow(sigs)!= length(spectrum))
    save(spectrum, sigs, file='spectrum.sigs.debug.R')

  stopifnot(nrow(sigs) == length(spectrum))

  # x0 is uniform distribution of mutations across signatures
  # (Each signature gets the same number of mutations)
  x0 <- rep(sum(spectrum) / ncol(sigs), ncol(sigs))

  res <- nloptr( x0       = x0
               , eval_f   = obj.fun
               , lb       = rep(0, ncol(sigs))
               , opts     = list( algorithm = algorithm
                                , xtol_rel  = xtol_rel
                                , maxeval   = maxeval
                                )
               , spectrum = spectrum
               , sigs     = sigs
               , ...
               )
  names(res$solution) <- colnames(sigs)
  res
}


#' Calculates the log-likelihood and activities for 1 tumor
#'
#' @param spect       matrix of mutation counts, 96 rows and a column for each tumor
#' @param sigs        signatures used to reconstruct the tumor
#' @param algorithm   algorithm for nloptr, default is \code{NLOPT_LN_COBYLA}
#' @param obj.fun     objective function for evaluation
#' @param nbinom.size dispersion parameter for the negative binomial model
#' @param logger      logger object
#' @return Returns a list with elements:
#' loglh    - the log likelihood of the best solution (set of exposures) found
#' exposure - the vector of exposures that generate loglh, in this case
#' 'exposure' means the number of mutations ascribed to each signature
#'
#' @importFrom log4r debug warn
one.lh.and.exp <- function(spect,
                           sigs,
                           algorithm='NLOPT_LN_COBYLA',
                           obj.fun,
                           nbinom.size,
                           logger) {
  r <- nloptr.one.tumor( spect
                       , sigs
                       , maxeval     = 1e6
                       , xtol_rel    = 1e-7
                       , algorithm   = algorithm
                       , obj.fun     = obj.fun
                       , nbinom.size = nbinom.size
                       , logger      = logger
                       )
  debug(logger, paste(r$objective, r$iterations, sum(r$solution)))

  loglh <- r$objective

  if (loglh == Inf)
    warn(logger, "Got -Inf in one.lh.and.exp")

  # sum(recon) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spect) / sum(r$solution)) # sum(recon))

  list( loglh    = -loglh
      , exposure = exp
      )
}


# ToDo: what are the following to comments doing here?
## Assign activities and test whether signature activities can be removed

## Breadth-first first search down to fixed number of levels


#' Helper function: is the set 'probe' a superset of any set in 'background'?
#'
#' @param probe      set to probe
#' @param background background set of signatures
#'
#' @return TRUE if current set probe is a superset of any background set
#'
#' @importFrom sets set_is_proper_subset
is.superset.of.any <- function(probe, background) {
  for (b in background) {
    if (set_is_proper_subset(b, probe)) return(T)
  }
  return(F)
}


#' Assigns activities of mutational signatures on a per-tumor basis
#'
#' @param spect       matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs        signatures used to reconstruct the tumor
#' @param max.level   maximum number of active signatures
#' @param p.thresh    p-value threshold for significance
#' @param obj.fun     objective function for evaluation
#' @param nbinom.size dispersion parameter for the negative binomial model
#' @param logger      logger object
#'
#' @return matrix of assignments of activity to each mutational signature in \code{sigs}
#'
#' @importFrom log4r debug
#' @importFrom sets set set_combn set_union
#' @importFrom stats pchisq
sparse.assign.activity <- function(spect,
                                   sigs,
                                   max.level=5,
                                   p.thresh=0.05,
                                   obj.fun,
                                   nbinom.size,
                                   logger) {
  mode(spect) <- "numeric"
  start       <- one.lh.and.exp( spect
                               , sigs
                               , obj.fun     = obj.fun
                               , nbinom.size = nbinom.size
                               , logger      = logger
                               )
  lh.w.all <- start$loglh
  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp > 0.5)

  debug(logger, paste("Starting with", paste(names(start.exp)[non.0.exp.index], collapse = ",")))
  debug(logger, paste("max.level =", max.level))

  if (length(non.0.exp.index) < 2)
  {
    debug(logger, paste("returning, only", length(non.0.exp.index), "non-0 exposures"))
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)

  c.r.s <- set() # subsets of the signatures indices that cannot be removed -- ToDo: change the name to something meaningful!

  max.p <- numeric(max.level)
  best.subset <- non.0.exp.index
  best.try <- start
  best.try.is.start <- T

  for (df in 1:max.level)
  {
    debug(logger, paste("df =", df))
    subsets <- set_combn(non.0.exp.index, df)
    for (subset in subsets)
    {
      # subset is of class set (package sets)
      if (is.superset.of.any(subset, c.r.s)) next
      subset.to.remove.v <- as.numeric(subset)
      subset.name        <- paste(names(start.exp)[subset.to.remove.v], collapse = ",")
      tmp.set            <- setdiff(non.0.exp.index, subset.to.remove.v)
      try.sigs           <- sigs[ , tmp.set]

      if (length(tmp.set) == 1)
      {
        try.sigs <- as.matrix(try.sigs)
        rownames(try.sigs) <- rownames(sigs)
      }

      # Get the max lh for try.sig
      try <- one.lh.and.exp( spect
                           , try.sigs
                           , obj.fun     = obj.fun
                           , nbinom.size = nbinom.size
                           , logger      = logger
                           )
      # try contains maxlh and exposure
      statistic <- 2 * (lh.w.all - try$loglh)
      chisq.p   <- pchisq(statistic, df, lower.tail = F)
      debug(logger, paste("Trying", subset.name, "p =", chisq.p, ";"))
      info(logger, paste("loglh =", try$loglh, "; statistic  =", statistic, ";"))
      if (chisq.p > p.thresh)
      {
        # This an acceptable set of exposures
        debug(logger, "acceptable;")
        if (chisq.p > max.p[df])
        {
          # This is the best so far
          max.p[df]         <- chisq.p
          best.subset       <- tmp.set
          best.try          <- try
          best.try.is.start <- F
          debug(logger, "best")
        }
        else
          debug(logger, "not best")
      }
      else
      {
        c.r.s <- set_union(c.r.s, set(subset))
        debug(logger, "cannot remove")
      }
    } # end for (subset in subsets)
    if (max.p[df] == 0) break
  } # end for df in 1:max.level

  # Need to return the exposures in the context of the orginal signatures matrix
  out.exp <- numeric(ncol(sigs)) #all zeros
  names(out.exp) <- colnames(sigs)
  if (best.try.is.start)
    out.exp <- start.exp
  else
    out.exp[best.subset] <- best.try$exposure

  debug(logger, paste("max.p =", paste(max.p, collapse = ", ")))
  debug(logger, paste("Ending with", paste(names(start.exp)[best.subset], collapse = ",")))

  stopifnot(abs(sum(out.exp) - sum(spect)) < 2)
  out.exp
}


#' Likelihood ratio test for reconstruction of tumors with and without signature of interest
#'
#' @param spect       matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs        signatures used to reconstruct the tumor
#' @param sig.to.test maximum number of active signatures
#' @param algorithm   algorithm for nloptr, default is \code{NLOPT_LN_COBYLA}
#' @param obj.fun     objective function for evaluation
#' @param nbinom.size dispersion parameter for the negative binomial model
#' @param logger      logger object
#'
#' @return list object with log-likelihoods calculated with and without signature of interest, test statistic and chi-square p-value
#'
#' @importFrom stats pchisq
#' @importFrom log4r info
is.present.p.m.likelihood <- function(spect,
                                      sigs,
                                      sig.to.test,
                                      algorithm='NLOPT_LN_COBYLA',
                                      obj.fun,
                                      nbinom.size,
                                      logger) {

  loglh.with    <- one.lh.and.exp( spect
                                 , sigs
                                 , algorithm   = algorithm
                                 , obj.fun     = obj.fun
                                 , nbinom.size = nbinom.size
                                 , logger      = logger
                                 )$loglh
  loglh.without <- one.lh.and.exp( spect
                                 , sigs[, -sig.to.test]
                                 , algorithm   = algorithm
                                 , obj.fun     = obj.fun
                                 , nbinom.size = nbinom.size
                                 , logger      = logger
                                 )$loglh

  statistic <- 2 * (loglh.with - loglh.without)
  chisq.p   <- pchisq(statistic, 1, lower.tail=F)

  info(logger, paste0("statistic = ", statistic, ", chisq p = ", chisq.p))

  list( with      = loglh.with
      , without   = loglh.without
      , statistic = statistic
      , chisq.p   = chisq.p
      )
}


#' Calculates the chi-square p-pval of likelihoods of reconstructed mutations counts with and without the signature of interest
#'
#' @param spect            matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs             signatures used to reconstruct the tumor
#' @param target.sig.index index of target signature
#' @param obj.fun          objective function for evaluation
#' @param nbinom.size      dispersion parameter for the negative binomial model
#' @param logger           logger object
#'
#' @return a p-value for the likelihood of a target is required to reconstruct the mutation counts of a tumor
signature.presence.test <- function(spect,
                                    sigs,
                                    target.sig.index,
                                    obj.fun,
                                    nbinom.size,
                                    logger) {
  is.present.p.m.likelihood( spect
                           , sigs
                           , target.sig.index
                           , obj.fun     = obj.fun
                           , nbinom.size = nbinom.size
                           , logger      = logger
                           )$chisq.p
}


#' Calculates exposures for entire set of tumors.
#'
#' @param spect       matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs        signatures used to reconstruct the tumor
#' @param exp         exposure or activity used to reconstruct a tumor
#' @param obj.fun     objective function for evaluation
#' @param nbinom.size dispersion parameter for the negative binomial model
#'
#' @return a list of negative log-likelihoods of the reconstructed matrices
compute.all.neg.log.lh <- function(spect,
                                   sigs,
                                   exp,
                                   obj.fun,
                                   nbinom.size) {
  out <- numeric(ncol(spect))
  names(out) <- colnames(spect)
  for (i in 1:ncol(spect)) {
    recon <- prop.reconstruct(sigs, exp[ , i])
    out[i] <- obj.fun(exp[ ,1], spect[ ,i], sigs,
                      nbinom.size=nbinom.size)
  }
  out
}


#' Carries out sanity check on the signature and activity matrices from the mSigAct reconstructions
#'
#' @param spectrum matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs     signatures used to reconstruct the tumor
#' @param exposure exposure or activity used to reconstruct a tumor

sanity.check.ex <- function(spectrum, sigs, exposure) {
  ex.sums <- margin.table(exposure, 2)
  all.reconstruct <- as.matrix(sigs)  %*% exposure
  rec.sums <- margin.table(all.reconstruct, 2)
  stopifnot(abs(ex.sums - rec.sums) < 1)
  spect.sums <- margin.table(spectrum, 2)
  stopifnot(abs(spect.sums - ex.sums) < 1)
}


#' Plots reconstruction of log-likelihoods of reconstructed tumors
#'
#' @param spect       matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs        signatures used to reconstruct the tumor
#' @param ex          exposure or activity used to reconstruct a tumor
#' @param range       range of tumors to plot reconstruction
#' @param obj.fun     objective function for evaluation
#' @param nbinom.size dispersion parameter for the negative binomial model
#'
#' @return reconstruction plots with negative log-likelihoods
#'
#' @importFrom graphics abline axis plot
plot.recon.and.loglh <- function(spect,
                                 sigs,
                                 ex,
                                 range,
                                 obj.fun,
                                 nbinom.size) {
  plot.reconstruction(signatures      = sigs,
                      exposures.mat   = ex[ , range, drop=F],
                      input.genomes   = spect[ , range, drop=F],
                      normalize.recon = T)



  neg.ll <- compute.all.neg.log.lh(spect[ ,range,drop=F], sigs=sigs,
                                   exp = ex[ , range, drop=F],
                                   obj.fun=obj.fun,
                                   nbinom.size=nbinom.size)
  s.names <- names(neg.ll)
  l.range <- 1:length(s.names)
  plot(neg.ll, xaxt='n', ylab=('Neg ln likelihood'), xlab='', new=T)
  axis(side = 1, labels = s.names, at=l.range, las=2)
  abline(v=l.range, lty=3)
}


#' plot.recon.by.range calls the objective function as one of its analyses.
#' Plots reconstruction of log-likelihoods of for a range of reconstructed tumors
#'
#' @param path        path for pdfs to be saved
#' @param spect       matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs        signatures used to reconstruct the tumor
#' @param ex          exposure or activity used to reconstruct a tumor
#' @param range       range of tumors to plot reconstruction
#' @param obj.fun     objective function for evaluation
#' @param nbinom.size dispersion parameter for the negative binomial model
#'
#' @return T
#'
#' @importFrom grDevices cairo_pdf dev.off
plot.recon.by.range <- function(path,
                                spect,
                                sigs,
                                ex,
                                range,
                                obj.fun,
                                nbinom.size) {
  cairo_pdf(path,
            width=8.2677, height=11.6929, # for A4
            onefile=T)
  par(
    mfrow=c(3,1), # 3 graphs per page
    mar=c(1.5,1.1,4.6,1), # space between plot and figure, for axis labels etc
    oma=c(4,6,3,3) # outer margin between figure and device.
  )

  for (r in range)
    plot.recon.and.loglh(spect, sigs, ex, r,
                         obj.fun=obj.fun,
                         nbinom.size=nbinom.size)
  dev.off()
}


#' run.mSigAct
#'
#' Test and plot one group of spectra
#'
#' Generates 4 PDFs based on the input argument out.dir and out.prefix.
#' The generated files are:
#'  * <out.dir>/<out.prefix>.check.with.sig.pdf
#'  * <out.dir>/<out.prefix>.exposures.pdf
#'  * <out.dir>/<out.prefix>.pval.histogram.pdf
#'  * <out.dir>/<out.prefix>.reconstruction.err.pdf
#'
#' @param spectra         matrix of mutation counts, 96 rows and a column for a tumor
#' @param sigs            signatures used to reconstruct the tumor
#' @param target.sig.name name of signature to be tested
#' @param out.dir         Path to the output directory for PDFs. PDFs will be
#'                        generated if and only if an outdir is specified.
#' @param out.prefix      Prefix for the output files. Will be prepended to the
#'                        names of the generated PDFs if out.dir is specified.
#' @param nbinom.size     dispersion parameter for the negative binomial
#' @param obj.fun         objective function for evaluation
#' @param loglevel        loglevel to control how detailed logging should be done.
#'                        Valid choices: DEBUG, INFO, WARN, ERROR, default: WARN.
#' @param logfile         path to the file in which the logs should be saved
#' @param col             ## is this used?
#' @param mc.cores        Number of cores to use for computations.
#'
#' @return Output is an R a list with the elements: pval, exposure
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics hist
#' @importFrom log4r create.logger error warn
#' @importFrom parallel mclapply
#'
#' @export
run.mSigAct <- function( spectra
                       , sigs
                       , target.sig.name
                       , nbinom.size
                       , out.dir       = NULL
                       , out.prefix    = ""
                       , obj.fun       = obj.fun.nbinom.maxlh
                       , loglevel      = "WARN"
                       , logfile       = "./mSigAct.log"
                       , col           = NULL
                       , mc.cores      = 1
                       )
{
  # set up the logger.
  if (!(loglevel %in% c("DEBUG", "INFO", "WARN", "ERROR")))
  {
    logger <- create.logger(logfile = logfile, level = "WARN")
    warn(logger, paste("Invalid loglevel", loglevel, "using default loglevel WARN."))
  }
  else
    logger <- create.logger(logfile = logfile, level = loglevel)

  target.sig.index <- which(colnames(sigs) == target.sig.name)

  # If an output dir was given, it should exist
  if (!is.null(out.dir) && !dir.exists(out.dir))
  {
    error(logger, "The specified output directory does not exist!")
    stop("The specified output directory does not exist!", call. = FALSE)
  }

  # The target signature must be in the set of given signatures
  if (!(target.sig.name %in% colnames(sigs)))
  {
    error(logger, "The target signature is not in the set of given signatures!")
    stop("The target signature is not in the set of given signatures!", call. = FALSE)
  }

  # Need to match exactly one signature name
  if (length(target.sig.index) != 1)
  {
    error(logger, "The target signature must be exactly one signature!")
    stop("The target signature must be exactly one signature!", call. = FALSE)
  }

  s.spectra         <- spectra.columns.sort(spectra)
  s.spectra.to.list <- split(t(s.spectra), 1:ncol(s.spectra))

  out.pvals <- mclapply( X                = s.spectra.to.list
                       , FUN              = signature.presence.test
                       , sigs             = sigs
                       , target.sig.index = target.sig.index
                       , obj.fun          = obj.fun
                       , nbinom.size      = nbinom.size
                       , mc.cores         = mc.cores
                       , logger           = logger
                       )

  out.pvals        <- unlist(out.pvals)
  names(out.pvals) <- colnames(s.spectra)

  low.pval <- which(out.pvals < 0.05)
  if (!is.null(out.dir) && length(low.pval) > 0)
  {
    # Have to wrap column-wise index of s.spectra in as.matrix in case
    # length(low.pval) == 1, in which case indexing returns a vector

    check.w.sig <- s.spectra[, low.pval, drop=F]
    # The column names are lost if length(low.pval) == 1
    colnames(check.w.sig) = colnames(s.spectra)[low.pval]
    spec.path <- file.path(out.dir, paste(out.prefix, 'check.with.sig.pdf', sep = "."))
    pdf.mut.sig.profile(path = spec.path, check.w.sig)
  }

  out.exp <- mclapply( X           = s.spectra.to.list
                     , FUN         = sparse.assign.activity
                     , sigs        = sigs
                     , obj.fun     = obj.fun
                     , nbinom.size = nbinom.size
                     , mc.cores    = mc.cores
                     )

  out.exp  <-  do.call(cbind, out.exp)
  colnames(out.exp)  <-  colnames(s.spectra)
  sanity.check.ex(s.spectra, sigs, out.exp)

  # Generate PDFs if an outdir was given
  if (!is.null(out.dir))
  {
    # Histogram
    hist.path <- file.path(out.dir, paste(out.prefix, 'pval.histogram.pdf', sep = "."))
    pdf(hist.path, useDingbats = F)
    hist(out.pvals, breaks = seq(from = 0, to = 1, by = 0.01))
    dev.off()

    # Exposures
    approx.num.per.row <- 30
    starts <- seq(from = 1, to = ncol(s.spectra), by = approx.num.per.row)
    ranges <- lapply(starts, function(x) { x:(min(x + approx.num.per.row - 1, ncol(s.spectra))) } )
    exp.path <- file.path(out.dir, paste(out.prefix, 'exposures.pdf', sep="."))
    pdf.ex.by.range(exp.path, s.spectra, sigs, exp = out.exp, range = ranges, col = col)

    # Reconstruction
    recon.path <- file.path(out.dir, paste(out.prefix, 'reconstruction.err.pdf', sep = "."))
    plot.recon.by.range( recon.path
                       , s.spectra
                       , sigs
                       , out.exp
                       , range       = ranges
                       , obj.fun     = obj.fun
                       , nbinom.size = nbinom.size
                       )
  }

  list(pval = out.pvals, exposure = out.exp)
}
