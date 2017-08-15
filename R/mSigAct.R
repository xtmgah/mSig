#'
#' mSigAct.R
#'
#' v 0.9
#'
#' 2017 08 11
#'
#' Copyright 2017 by Alvin Wei Tian Ng, Steven G. Rozen
#'
#' The code is released under GPL-3
#' https://www.gnu.org/licenses/gpl-3.0.en.html
#'
#' This file contains functions for sparse maximum likelihood assignment of
#' mutational signature activities and for testing whether a given signature is
#' needed to explain a spectrum

# Check R version and dependencies
if (!(R.version$major >= "3")) stop("mSigAct only works with R 3.2 or newer")
if (!(R.version$minor >= "3.2")) stop("mSigAct only works with R 3.2 or newer")

# Dependencies -- ToDo: what's the proper way to define dependencies?
#library(nloptr)
#library(sets)
#library(parallel)


#' Helper function, given signatures (sigs) and exposures (exp), return a
#' *proportional* reconstruction; in general, it is *not necessarily* scaled to
#' the actual spectrum counts.
#'
#' @param sigs  signatures
#' @param exp   exposures to those signatures
#'
#' @return ToDo
prop.reconstruct <- function(sigs, exp) {
  stopifnot(length(exp) == ncol(sigs))
  as.matrix(sigs) %*% exp
}


#' Define the objective function for use with nloptr (non-linear optimization)
#' Negative binomial maximum likelihood objective function
#' (nloptr minimizes the objective function.)
#' The lower the objective function, the better
#'
#' @param exp         the matrix of exposures ("activities")
#' @param sigs        the matrix of signatures
#' @param spectrum    the spectrum to assess
#' @param nbinom.size the dispersion parameter for the negative binomial
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
#' @param spectrum    ToDo
#' @param sigs        ToDo
#' @param algorithm   ToDo
#' @param maxeval     ToDo
#' @param print_level ToDo
#' @param xtol_rel    ToDo
#' @param obj.fun     ToDo
#' @param ...         ToDo
#'
#' @return ToDo
#'
#' @importFrom nloptr nloptr
nloptr.one.tumor <- function(spectrum,
                             sigs,
                             algorithm='NLOPT_LN_COBYLA', # ToDo: would it make sense to make this customisable?
                             maxeval=1000,
                             print_level=0, # ToDo: is this the same as "trace" in other functions? why give it different names?
                             xtol_rel=0.001, # 0.0001,
                             obj.fun,
                             ... # additional arguments for obj.fun
                             ) {
  if (class(sigs) == 'numeric') {
    # In this case we got a numeric vector, not a matrix.
    # We assume the caller intended it as a single
    # column matrix.
    sigs <- matrix(sigs, ncol=1)
  }

  # The if statement is an example of captuing state with
  # a call to save, which seems useful for debugging
  # in a call to mclapply (multi-core lapply).
  if (nrow(sigs)!= length(spectrum)) {
    save(spectrum, sigs, file='spectrum.sigs.debug.R')
  }
  stopifnot(nrow(sigs) == length(spectrum))

  # x0 is uniform distribution of mutations across signatures
  # (Each signature gets the same number of mutations)
  x0 <- rep(sum(spectrum) / ncol(sigs), ncol(sigs))

  res <- nloptr(x0=x0,
                eval_f=obj.fun,
                lb=rep(0, ncol(sigs)),
                opts=list(algorithm=algorithm,
                          xtol_rel=xtol_rel,
                          print_level=print_level,
                          maxeval=maxeval),
                spectrum=spectrum,
                sigs=sigs,
                ...)
  names(res$solution) <- colnames(sigs)
  res
}


#' ToDo: document this function!
#'
#' @param spect       ToDo
#' @param sigs        ToDo
#' @param trace       ToDo
#' @param algorithm   ToDo
#' @param obj.fun     ToDo
#' @param nbinom.size ToDo
#'
#' @return Returns a list with elements:
#' loglh    - the log likelihood of the best solution (set of exposures) found
#' exposure - the vector of exposures that generate loglh, in this case
#' 'exposure' means the number of mutations ascribed to each signature
one.lh.and.exp <- function(spect,
                           sigs,
                           trace,
                           algorithm='NLOPT_LN_COBYLA',
                           obj.fun,
                           nbinom.size) {
  r <- nloptr.one.tumor(spect, sigs, maxeval=1e6,
                        xtol_rel = 1e-7,
                        algorithm=algorithm, obj.fun = obj.fun,
                        nbinom.size=nbinom.size)
  if (trace >  0) cat(r$objective, r$iterations, sum(r$solution), '\n')

  loglh <- r$objective

  if (loglh == Inf && trace > 0) cat("Got -Inf in one.lh.and.exp\n")

  # sum(recon) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spect) / sum(r$solution)) # sum(recon))

  list(loglh=-loglh, exposure=exp)
}

# ToDo: what are the following to comments doing here?
## Assign activities and test whether signature activities can be removed

## Breadth-first first search down to fixed number of levels

#' Helper function: is the set 'probe' a superset of any set in 'background'?
#'
#' @param probe      ToDo
#' @param background ToDo
#'
#' @return ToDo
#'
#' @importFrom sets set_is_proper_subset
is.superset.of.any <- function(probe, background) {
  for (b in background) {
    if (set_is_proper_subset(b, probe)) return(T)
  }
  return(F)
}


#' ToDo: document this function!
#'
#' @param spect       ToDo
#' @param sigs        ToDo
#' @param max.level   ToDo
#' @param p.thresh    ToDo
#' @param trace       ToDo
#' @param obj.fun     ToDo
#' @param nbinom.size ToDo
#'
#' @return ToDo
#'
#' @importFrom sets set set_combn set_union
#' @importFrom stats pchisq
sparse.assign.activity <- function(spect,
                                   sigs,
                                   max.level=5,
                                   p.thresh=0.05,
                                   trace=0,
                                   obj.fun,
                                   nbinom.size) {
  mode(spect) <-  'numeric'
  start <- one.lh.and.exp(spect, sigs, trace=0,
                          obj.fun=obj.fun,
                          nbinom.size=nbinom.size)
  lh.w.all <- start$loglh
  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp > 0.5)
  if (trace > 0) {
    cat('Starting with',
        paste(names(start.exp)[non.0.exp.index], collapse=','),
        '\n')
    cat('max.level =', max.level, '\n')
  }
  if (length(non.0.exp.index) < 2) {
    if (trace > 0) cat('returning, only', length(non.0.exp.index),
                       'non-0 exposures\n')
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)

  c.r.s <- set() # subsets of the signatures indices that cannot be removed -- ToDo: change the name to something meaningful!

  max.p <- numeric(max.level)
  best.subset <- non.0.exp.index
  best.try <- start
  best.try.is.start <- T

  for (df in 1:max.level) {
    if (trace > 0) cat('df =', df, '\n')
    subsets <- set_combn(non.0.exp.index, df)
    for (subset in subsets) {
      # subset is of class set (package sets)
      if (is.superset.of.any(subset, c.r.s)) next
      subset.to.remove.v <- as.numeric(subset)
      subset.name <- paste(names(start.exp)[subset.to.remove.v], collapse=',')
      tmp.set <- setdiff(non.0.exp.index, subset.to.remove.v)
      try.sigs <- sigs[ , tmp.set]

      if (length(tmp.set) == 1) {
        try.sigs <- as.matrix(try.sigs)
        rownames(try.sigs) <- rownames(sigs)
        if (trace > 0) cat("New code\n")
      }

      # Get the max lh for try.sig
      try <- one.lh.and.exp(spect, try.sigs, trace=0,
                            obj.fun=obj.fun,
                            nbinom.size=nbinom.size)
      # try contains maxlh and exposure
      statistic <- 2 * (lh.w.all - try$loglh)
      chisq.p <- pchisq(statistic, df, lower.tail=F)
      if (trace > 0) {
        cat('Trying', subset.name, 'p =', chisq.p, '; ')
        if (trace > 1) cat('loglh =', try$loglh, '; statistic  =', statistic, ';')
      }
      if (chisq.p > p.thresh) {
        # This an acceptable set of exposures
        if (trace > 0) cat (' acceptable;')
        if (chisq.p > max.p[df]) {
          # This is the best so far
          max.p[df] <- chisq.p
          best.subset <- tmp.set
          best.try <- try
          best.try.is.start <- F
          if (trace > 0) cat('best\n')
        } else {
          if (trace > 0) cat('not best\n')
        }
      } else {
        c.r.s <- set_union(c.r.s, set(subset))
        if (trace > 0)  {
          cat('cannot remove\n')
        }
      }
    } # end for (subset in subsets)
    if (max.p[df] == 0) break
  } # end for df in 1:max.level

  # Need to return the exposures in the context of the orginal signatures matrix
  out.exp <- numeric(ncol(sigs)) #all zeros
  names(out.exp) <- colnames(sigs)
  if (best.try.is.start) {
    out.exp <- start.exp
  } else {
    out.exp[best.subset] <- best.try$exposure
  }
  if (trace > 0) {
    cat('max.p =', paste(max.p, collapse = ', '), '\n')
    cat('Ending with',
        paste(names(start.exp)[best.subset], collapse=','),
        '\n')
  }
  stopifnot(abs(sum(out.exp) - sum(spect)) < 2)
  out.exp
}


#' Likelihood ratio test. ToDo: extend the documentation!
#'
#' @param spect       ToDo
#' @param sigs        ToDo
#' @param sig.to.test ToDo
#' @param trace       ToDo
#' @param algorithm   ToDo
#' @param obj.fun     ToDo
#' @param nbinom.size ToDo
#'
#' @return ToDo
#'
#' @importFrom stats pchisq
is.present.p.m.likelihood <- function(spect,
                                      sigs,
                                      sig.to.test,
                                      trace=0,
                                      algorithm='NLOPT_LN_COBYLA',
                                      obj.fun,
                                      nbinom.size) {

  loglh.with <- one.lh.and.exp(spect, sigs, trace=trace,
                               algorithm=algorithm, obj.fun=obj.fun,
                               nbinom.size=nbinom.size)$loglh
  loglh.without <- one.lh.and.exp(spect, sigs[ ,-sig.to.test], trace=trace,
                                  algorithm=algorithm,
                                  obj.fun=obj.fun,
                                  nbinom.size=nbinom.size)$loglh
  statistic <- 2 * (loglh.with - loglh.without)
  chisq.p <- pchisq(statistic, 1, lower.tail=F)
  if (trace > 0) cat('statistic  =', statistic, '\nchisq p =', chisq.p, '\n')

  list(with=loglh.with,
       without=loglh.without,
       statistic=statistic,
       chisq.p=chisq.p)
}


#' ToDo: document this function!
#'
#' @param spect            ToDo
#' @param sigs             ToDo
#' @param target.sig.index ToDo
#' @param trace            ToDo
#' @param obj.fun          ToDo
#' @param nbinom.size      ToDo
#'
#' @return ToDo
signature.presence.test <- function(spect,
                                    sigs,
                                    target.sig.index,
                                    trace=0,
                                    obj.fun,
                                    nbinom.size) {
  is.present.p.m.likelihood(spect,
                       sigs, target.sig.index,
                       trace=trace, obj.fun=obj.fun,
                       nbinom.size=nbinom.size)$chisq.p
}


#' Cacluate exposures for entire set of tumors.
#'
#' @param spect       ToDo
#' @param sigs        ToDo
#' @param exp         ToDo
#' @param obj.fun     ToDo
#' @param nbinom.size ToDo
#'
#' @return ToDo
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


#' ToDo: document this function! Is this some kind of invariant?
#'
#' @param spectrum ToDo
#' @param sigs     ToDo
#' @param exposure ToDo
#'
#' @return ToDo
sanity.check.ex <- function(spectrum, sigs, exposure) {
  ex.sums <- margin.table(exposure, 2)
  all.reconstruct <- as.matrix(sigs)  %*% exposure
  rec.sums <- margin.table(all.reconstruct, 2)
  stopifnot(abs(ex.sums - rec.sums) < 1)
  spect.sums <- margin.table(spectrum, 2)
  stopifnot(abs(spect.sums - ex.sums) < 1)
}


#' ToDo: document this function!
#'
#' @param spect       ToDo
#' @param sigs        ToDo
#' @param ex          ToDo
#' @param range       ToDo
#' @param obj.fun     ToDo
#' @param nbinom.size ToDo
#'
#' @return ToDo
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
#' ToDo: provide a useful documentation!
#'
#' @param path        ToDo
#' @param spect       ToDo
#' @param sigs        ToDo
#' @param ex          ToDo
#' @param range       ToDo
#' @param obj.fun     ToDo
#' @param nbinom.size ToDo
#'
#' @return ToDo
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

  for (r in range) {
    plot.recon.and.loglh(spect, sigs, ex, r,
                         obj.fun=obj.fun,
                         nbinom.size=nbinom.size)
  }
  dev.off()
}


#' run.mSigAct
#'
#' Test and plot one group of spectra
#'
#' ToDo: don't do this! let the user decide whether (s)he wants the side-effects
#'       or not. use a flag to turn this on/off!
#'
#' Side effects are to generate 4 pdfs based on the input argument
#' path.root. The names are:
#' <path.root>.check.with.sig.pdf
#' <path.root>.pos.with.sig.pdf ## FIX ME, probably removed
#' <path.root>.exposures.pdf
#' <path.root>.reconstruction.err.pdf
#'
#' @param spectra         ToDo
#' @param sigs            ToDo
#' @param target.sig.name ToDo
#' @param path.root       ToDo
#' @param obj.fun         ToDo
#' @param nbinom.size     ToDo
#' @param trace           ToDo
#' @param col             ToDo
#' @param mc.cores        ToDo -- what, 190?
#'
#' @return ToDo: Output is an R a list with the elements: pval, exposure
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics hist
#' @importFrom parallel mclapply
#'
#' @export
run.mSigAct <- function(spectra,
                        sigs,
                        target.sig.name,
                        path.root,
                        obj.fun,
                        nbinom.size,
                        trace=0,
                        col=NULL,
                        mc.cores=190) {

  target.sig.index <- which(colnames(sigs) == target.sig.name)

  # Need to match exactly one signature name
  stopifnot(length(target.sig.index) == 1)

  s.spectra <- sort.spectra.columns(spectra)
  s.spectra.to.list <- split(t(s.spectra), 1:ncol(s.spectra))

  out.pvals <-
      mclapply(
          X=s.spectra.to.list,
          FUN=signature.presence.test,
          sigs=sigs,
          target.sig.index=target.sig.index,
          trace=trace,
          obj.fun=obj.fun,
          nbinom.size=nbinom.size,
          mc.cores=mc.cores)

  out.pvals  <- unlist(out.pvals)
  names(out.pvals)  <- colnames(s.spectra)

  low.pval <- which(out.pvals < 0.05)
  if (length(low.pval) > 0) {
    # Have to wrap column-wise index of s.spectra in as.matrix in case
    # length(low.pval) == 1, in which case indexing returns a vector

    check.w.sig <- s.spectra[, low.pval, drop=F]
    # The column names are lost if length(low.pval) == 1
    colnames(check.w.sig) = colnames(s.spectra)[low.pval]
    spec.path <- paste(path.root, 'check.with.sig.pdf', sep='.')
    pdf.mut.sig.profile(path=spec.path, check.w.sig)
  }

  out.exp <-
      mclapply(
          X=s.spectra.to.list,
          FUN=sparse.assign.activity,
          sigs=sigs,
          obj.fun=obj.fun,
          nbinom.size=nbinom.size,
          mc.cores=mc.cores)

  out.exp  <-  do.call(cbind, out.exp)
  colnames(out.exp)  <-  colnames(s.spectra)
  sanity.check.ex(s.spectra, sigs, out.exp)

  # Plotting part
  hist.path <- paste(path.root, 'pval.histogram.pdf', sep='.')
  pdf(hist.path, useDingbats = F)
  hist(out.pvals, breaks=seq(from=0, to=1, by=0.01))
  dev.off()

  approx.num.per.row <- 30
  starts <- seq(from=1, to=ncol(s.spectra), by=approx.num.per.row)
  ranges <-
    lapply(starts,
           function(x) {
             x:(min(x+approx.num.per.row-1, ncol(s.spectra)))
             } )
  exp.path <- paste(path.root, 'exposures.pdf', sep='.')
  pdf.ex.by.range(exp.path, s.spectra, sigs, exp=out.exp,
                  range=ranges, col=col)
  recon.path <- paste(path.root, 'reconstruction.err.pdf', sep='.')

  plot.recon.by.range(recon.path, s.spectra,
                      sigs,
                      out.exp,
                      range = ranges,
                      obj.fun=obj.fun,
                      nbinom.size=nbinom.size)

  list(pval=out.pvals, exposure=out.exp)
}
