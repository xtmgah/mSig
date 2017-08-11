#
# mSigActTest
#
# 2017 05 30
#
# Tests for mSigAct

source('mSigTools.R')
source('mSigAct.v0.4.R')

# Get signatures

cosmic.sigs <- get.COSMIC.signatures()
cosmic.wes <- cosmic.sigs$exome
cosmic.wgs <- cosmic.sigs$genome

liver.sig.names <- paste('Signature', 
                         c(1, 4, 5, 6,
                           12, 16, 17, 'AA', 23, 
                           24),
                         sep = '.')

liver.wes.sigs <- cosmic.wes[, liver.sig.names]
liver.wgs.sigs <- cosmic.wgs[, liver.sig.names]

paper.liver.sig.names <- paste('Signature', c(1, 4, 5, 6,
                                        12, 16, 17, 23, 
                                        24, 'W6', 'AA'), sep = '.')

paper.liver.wes.sigs <- cosmic.wes[, paper.liver.sig.names]
paper.liver.wgs.sigs <- cosmic.wgs[, paper.liver.sig.names]

## Load a WES data set for testing

tmp.taiwan <- 
  read.and.prep.192.duke.nus.catalog('Poon_et_all_HCC_spectrum_counts.tsv')$channel.96
taiwan.hcc <- sort.spectra.columns(tmp.taiwan)
rm(tmp.taiwan)
short.taiwan.hcc <- taiwan.hcc[ , (ncol(taiwan.hcc)-10):ncol(taiwan.hcc)]
xtract.col.tw <- function(sample.name) {
  xtract.col(taiwan.hcc, sample.name)
}

# This is exactly the same test as in mSigAct.basic.test
# The only difference is that we read in the data from files.
short.taiwan.analysis <-
  process.one.group(spectra=short.taiwan.hcc,
                    sigs=liver.wes.sigs,
                    target.sig.name='Signature.AA',
                    path.root='mSigActTest-SHORT-WES',
                    obj.fun=obj.fun.nbinom.maxlh,
                    nbinom.size=5,
                    mc.cores=1)

expected.short.pval <-
  structure(c(0.036297099364187, 0.76185665842143, 0.000540732195433939, 
              0.999999809769497, 0.000672367936677276, 0.999997831038551, 0.120653696546247, 
              0.996613255999323, 0.801907129656769, 1.68741386815852e-05, 1.46775631766299e-05
  ), .Names = c("T68", "T41", "T74", "T16", "T46", "T82", "T50", 
                "T15", "T24", "T80", "T95"))

stopifnot(all.equal(short.taiwan.analysis$pval, expected.short.pval, tolerance = 0.005))

expected.exp <-
  structure(c(10.5883215430212, 0, 43.0167421150595, 0, 0, 0, 0, 
              6.39493634191927, 0, 0, 0, 0, 57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 37.5493947700166, 0, 11.6064867320358, 7.84411849794756, 
              0, 0, 0, 0, 0, 0, 37.8015117701765, 0, 0, 15.1984882298235, 0, 
              0, 0, 41.1408335720938, 0, 0, 0, 0, 10.8591664279062, 0, 0, 0, 
              0, 51, 0, 0, 0, 0, 0, 0, 0, 7.21731841509026, 14.5816564545762, 
              0, 0, 0, 25.2010251303335, 0, 0, 0, 0, 0, 0, 45, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 30.9066420035866, 13.0933579964134, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 29.2344487210762, 0, 14.7655512789238, 0, 
              0, 0, 0, 16.6275775292156, 0, 0, 0, 0, 12.3724224707844, 0, 0
  ), .Dim = 10:11,
  .Dimnames =
    list(c("Signature.1", "Signature.4", 
           "Signature.5", "Signature.6", "Signature.12", "Signature.16", 
           "Signature.17", "Signature.AA", "Signature.23", "Signature.24"
    ), c("T68", "T41", "T74", "T16", "T46", "T82", "T50", "T15", 
         "T24", "T80", "T95")))

stopifnot(all(all.equal(short.taiwan.analysis$exposure, expected.exp, tolerance = 0.005)))

## Load a WGS data set for testing

japan.wgs.hcc <-  
  read.96.duke.nus.format('LIRI-JP.WGS.96.mutation.counts.txt')

# Get a small subset of japan HCCs for testing
test.japan.wgs <- japan.wgs.hcc[ , c(
  # Arbitrary
  'RK001', 'RK002', 'RK003',
  # Probably not AA
  'RK119', 'RK030', 'RK227',
  # Probably AA
  'RK054', 'RK027'
)]

# STOP HERE

test.japan.analysis <- process.one.group(test.japan.wgs, liver.wgs.sigs,
                                         target.sig.name = 'Signature.AA',
                                         path.root = 'mSigActTest-WGS')
regress.japan.pval <-
  structure(c(2.9204079754646e-11, 4.88244331670293e-248, 4.77664635249271e-22, 
              1.76287439451299e-75, 0.0258984292206763, 7.29770301437411e-16, 
              0.0658417237333467, 2.92357348087274e-125),
            .Names = c("RK119", 
                       "RK027", "RK001", "RK030", "RK002", "RK227", "RK003", "RK054"
            ))
stopifnot(abs(test.japan.analysis$pval - regress.japan.pval) < 1e16)

regress.japan.exposure <-
  structure(c(4124.90032175754, 397.877885627809, 4078.65778077994, 
              13097.1580971145, 11656.3880302682, 3.86585833043579e-17, 1.47569710991701e-18, 
              300.278651619047, 285.739236077272, 1.74604540325637e-15, 146.606296037337, 
              3224.4622663913, 6601.99974757438, 3.26398491540463e-37, 1389.43018164366, 
              9572.1720523929, 5.8490609684051e-35, 2229.32945974629, 3.03970148826269e-20, 
              1.98853149735134e-20, 147.194154340277, 716.203955127181, 6486.0310025824, 
              5.40319383256908e-19, 5547.82711174046, 5702.02400782848, 1.56822068573014e-19, 
              501.670351230815, 264.04942057699, 1.55845876561982e-19, 0, 467.522368332369, 
              3674.4901195691, 0, 4922.39148270053, 3340.97878403383, 0, 838.195246231293, 
              125.422001158705, 0, 253.983431493705, 1238.96671223378, 2399.58109982024, 
              2.4124490409495e-35, 4617.14978671407, 4006.40185804102, 1.8602662604655e-33, 
              87.9171129821196, 2.56957039184572e-18, 4.50432544952547e-20, 
              392.557052573866, 285.875570093607, 1736.63583440389, 7.17006103065395e-21, 
              2893.41228139777, 2995.441652065, 6.00816650596711e-20, 288.208034811495, 
              147.225619185346, 197.643956533665, 227.228070221594, 612.050344936974, 
              2682.12811958984, 0, 408.34569825861, 2198.14413668746, 644.035709732559, 
              0, 62.0679215507593, 0, 169.152386251118, 462.863110668665, 1976.69876901154, 
              0, 442.445217451691, 2314.83702199208, 0, 845.317947014172, 91.6855488284509, 
              0), .Dim = c(10L, 8L), 
            .Dimnames = list(c("Signature.1", "Signature.4", 
                               "Signature.5", "Signature.6", "Signature.12", "Signature.16", 
                               "Signature.17", "Signature.AA", "Signature.23", "Signature.24"
            ), c("RK119", "RK027", "RK001", "RK030", "RK002", "RK227", "RK003", 
                 "RK054")))

stopifnot(abs(test.japan.analysis$exposure - regress.japan.exposure) < 1e-11)

### ==============================================================================
### FIX ME =======================================================================
### THIS IS FOR THE PAPER, REMOVE LATER
paper.japan.list.adj <- list()
paper.taiwan.list.adj <- list()
paper.japan.list.exp <- list()
paper.taiwan.list.exp <- list()

for (.nbinom.size in c(5,3)) {
  
  if (0) {
    
    file.out <- paste('new-mSigActtest-POIS-WGS', .nbinom.size, sep='-')
    
    
    p.test.japan.analysis <- process.one.group(test.japan.wgs, liver.wgs.sigs,
                                               target.sig.name = 'Signature.AA',
                                               path.root = file.out, 
                                               obj.fun=obj.fun.poisson.maxlh)
    
    data.frame(mult=test.japan.analysis$pval, pois=p.test.japan.analysis$pval)
  }
  
  p.t.file <- paste('new-taiwan-WES-POS', .nbinom.size, sep='-')
  
  p.taiwan.analysis <-  process.one.group(taiwan.hcc, paper.liver.wes.sigs,
                                          target.sig.name = 'Signature.AA',
                                          path.root = p.t.file,
                                          obj.fun=obj.fun.poisson.maxlh)
  
  cat('\n', .nbinom.size, '\n', p.taiwan.analysis$pval, '\n', 
      file='new-disp-params.txt', append=T)
  paper.taiwan.list.adj[[.nbinom.size]] <- p.taiwan.analysis$adj.p.val
  paper.taiwan.list.exp[[.nbinom.size]] <- p.taiwan.analysis$exposure
  
  p.j.file <- paste('new-japan-WGS-POS', .nbinom.size, sep='-')
  
  p.japan.analysis <- process.one.group(japan.wgs.hcc, paper.liver.wgs.sigs,
                                        target.sig.name = 'Signature.AA',
                                        path.root = p.j.file,
                                        obj.fun = obj.fun.poisson.maxlh)
  
  cat('\n', .nbinom.size, '\n', p.japan.analysis$pval, '\n', 
      file='new-disp-params.txt', append=T)
  paper.japan.list.adj[[.nbinom.size]] <- p.japan.analysis$adj.p.val
  paper.japan.list.exp[[.nbinom.size]] <- p.japan.analysis$exposure
}
### END FIX ME MATERIAL FOR THE PAPER ==================================
### ====================================================================

test.down.japan <- process.one.group(down.japan, liver.wgs.sigs,
                                         target.sig.name = 'Signature.AA',
                                         path.root = 'mSigActTest-down-WGS')



test.msi.analysis <- process.one.group(test.msi, liver.wgs.sigs,
                                         target.sig.name = 'Signature.AA',
                                         path.root = 'mSigActTest-MSI')

#### Long running test

taiwan.analysis <-  process.one.group(taiwan.hcc, liver.wes.sigs,
                                      target.sig.name = 'Signature.AA',
                                      path.root = 'mSigActTest-WES')

#### Various sanity checks and regression tests on the Taiwan analysis

taiwan.sanity.checks <- function() {
  stopifnot(rownames(cosmic.wes) == rownames((taiwan.hcc)))
  stopifnot(dim(taiwan.hcc) == c(96, 98))
  taiwan.exp.ml <- taiwan.analysis$exposure
  stopifnot(sum(taiwan.exp.ml  < 0.5) == 825) # The number of nearly 0 exposures
  taiwan.pvals <- taiwan.analysis$pval
  taiwan.fdrs <- taiwan.analysis$adj.p.val
  sums.aa <- taiwan.exp.ml['Signature.AA', ]
  high.fdrs <- taiwan.fdrs[taiwan.fdrs > 0.05]
  fdr.set <- names(high.fdrs)
  pvals.set <- names(taiwan.pvals[taiwan.pvals > 0.05])
  stopifnot(setdiff(fdr.set, pvals.set) == c("T57", "T47",  "T12" ))
  low.exp.set <- colnames(taiwan.exp.ml)[taiwan.exp.ml['Signature.AA', ] < 0.5]
  stopifnot(setdiff(low.exp.set, pvals.set) == c())
  # The next is a regression test
  stopifnot(round(sums.aa[setdiff(fdr.set, pvals.set)]) == c(8, 15, 13))
}
taiwan.sanity.checks()

# High AA
T43 <- xtract.col.tw('T43') # JH146
t43 <- T43

# Test compare.m.likelihood
tmp1 <- compare.m.likelihood(T43, liver.wes.sigs, aa.sig.index,
                             obj.fun=obj.fun.maxlh)
stopifnot(abs(tmp1$statitics - 3904.530383) < 0.0001)
stopifnot(tmp1$chisq.p < 1e-10)

tmp1 <- compare.m.likelihood(T43, liver.wes.sigs, 1,
                             obj.fun=obj.fun.maxlh)
stopifnot(tmp1$statistic < 1e-9)
stopifnot(1 - tmp1$chisq.p < 0.0001)

### test obj.fun.maxlh
# dmultinom returns -Inf if solution (exposures) contain 0
test.sigs <- cosmic.wes[, paste('Signature', c(1, 'AA'), sep = '.')]
t.out <- obj.fun.maxlh(x=c(0.00000000001, 1329.869),
                       spectrum=t43, sigs.as.matrix = test.sigs)
stopifnot(abs(t.out - 732.9628261) < 0.00001)

t.out <- obj.fun.poisson.maxlh(x=c(0.00000000001, 1329.869),
                       spectrum=t43, sigs.as.matrix = test.sigs)
stopifnot(abs(t.out - 737.425203) < 0.00001)


t.out <- obj.fun.maxlh(x=c(500, 1329.869-500),
                       spectrum=T43, sigs.as.matrix = as.matrix(test.sigs))
stopifnot(abs(t.out - 657.5251876) < 0.00001)

# Test nloptr.one.tumor
t43.out.ml <- nloptr.one.tumor(t43, liver.wes.sigs,
                               maxeval=1e6, xtol_rel=1e-5, obj.fun = obj.fun.maxlh)

t43.out.pois.ml <- nloptr.one.tumor(t43, liver.wes.sigs,
                               maxeval=1e6, xtol_rel=1e-5, obj.fun = obj.fun.poisson.maxlh)


# t43.out.ml
tmp.t43.sol <- T43.out.ml$solution
tmp.t43.sol <- tmp.t43.sol * (1 / sum(tmp.t43.sol))
stopifnot(abs(max(tmp.t43.sol) - 0.963542) < 0.0001)
stopifnot(names(which(tmp.t43.sol == max(tmp.t43.sol))) == 'Signature.AA')
stopifnot(abs(T43.out.ml$objective - 256.857753) < 0.0001)

# t43.out.ml
tmp.t43.pois.sol <- t43.out.pois.ml$solution
tmp.t43.pois.sol <- tmp.t43.pois.sol * (1 / sum(tmp.t43.pois.sol))
stopifnot(abs(max(tmp.t43.pois.sol) - 0.963545) < 0.0001)
stopifnot(names(which(tmp.t43.pois.sol == max(tmp.t43.pois.sol))) == 'Signature.AA')
stopifnot(abs(t43.out.pois.ml$objective - 261.32013) < 0.0001)


# Test reconstruct
T43.recon <- prop.reconstruct(liver.wes.sigs, T43.out.ml$solution)
scaled.recon <- scaled.reconstruct(T43, liver.wes.sigs, T43.out.ml$solution)
t43.abs.dif <- abs(T43 - scaled.recon)
t43.max.dif <- which(t43.abs.dif == max(t43.abs.dif))
stopifnot(T43[t43.max.dif, ] == 53)

test.olae <- one.lh.and.exp(spect=t43, 
                            sigs=liver.wes.sigs,
                            trace=0,
                            obj.fun=obj.fun.maxlh)
stopifnot(sum(test.olae$exposure) == sum(t43))
stopifnot(round(test.olae$exposure['Signature.AA'], digits=0) == 1152)

test.olae <- one.lh.and.exp(spect=t43, 
                            sigs=liver.wes.sigs,
                            trace=0,
                            obj.fun=obj.fun.poisson.maxlh)
stopifnot(abs(sum(test.olae$exposure) - sum(t43)) < 1e12)
stopifnot(round(test.olae$exposure['Signature.AA'], digits=0) == 1152)

aa.sig.index <- which(colnames(liver.wes.sigs) == 'Signature.AA')
sig.present.p(spect=t43, sigs=liver.wes.sigs, 
              target.sig.index = aa.sig.index, 
              trace=0, obj.fun = obj.fun.maxlh)

sig.present.p(spect=t43, sigs=liver.wes.sigs, 
              target.sig.index = aa.sig.index, 
              trace=0, obj.fun = obj.fun.poisson.maxlh)
rm(aa.sig.index)

####### OUT OF DATE BELOW HERE

# Low AA
T08 <- xtract.col.tw('T08') # CGMH553

# Almost no AA
T24 <- xtract.col.tw('T24') # JH104

# No AA by nsNMF
JH117 <- xtract.col.tw('JH117')

# Highest total number of mutatations
JH131 <- xtract.col.tw('JH131')

# ~90 mutations, mix of Sig 4 (lots) and Sig 24 (few)
JH170 <- xtract.col.tw('JH170')

# Led to error in sparse.assign.activity at one point
JH120 <- xtract.col.tw('JH120')


## Further testing of sparse.assign.activity (with tracing, need to visually inspect trace)

# Address differences between just removing AA from the best reconstruction and
# finding the reconstruction with the fewest signatures that can still
# reasonably explain the data.

# The affected tumors include JH33, JH170, and JH24


T68 <- xtract.col.tw('T68')
JH33.exp <- sparse.assign.activity(JH33, cosmic.liver, max.level = 5, trace=1)

JH170.exp <- sparse.assign.activity(JH170, cosmic.liver, max.level = 5, trace = 1)
pvals['JH170']

JH24 <- xtract.col.tw('JH24')
JH24.exp <- sparse.assign.activity(JH24, cosmic.liver, max.level= 5, trace = 1)
stopifnot(sum(JH24.exp) == 95)

JH131.out.ml <- nloptr.one.tumor(JH131, cosmic.liver,
                                 maxeval=1e6, obj.fun = obj.fun.maxlh)
JH131.out.ml$solution
stopifnot(abs(JH131.out.ml$objective - 265.8636) < 0.001)
stopifnot(abs(t.out - 657.5251875) < 0.0001)

## Test sparse.assign.activity

test.sparse.assign.activity <- function() {
  exp.out <- sparse.assign.activity(JH120, cosmic.liver)
  stopifnot(all(which(exp.out > 0.1) == c(3, 5, 7)))
  stopifnot(all(abs(exp.out[c(3, 5, 7)] - c(222.9107, 128.846, 441.243)) < 0.001))

  exp.out <- sparse.assign.activity(JH131, cosmic.liver)
  stopifnot(which(exp.out > 0.1) == c(3, 5, 7))
  stopifnot(all(abs(exp.out[c(3,5,7)] - c(14.2141, 85.1168, 1155.6690)) < 0.001))

  exp.out <- sparse.assign.activity(aa.hcc[, 'JH170'], cosmic.liver)
  stopifnot(which(exp.out > 0.1) == 3)
  stopifnot(exp.out[3] == 86)

  exp.out <- sparse.assign.activity(JH170, cosmic.liver, max.level=4)
  stopifnot(which(exp.out > 0.1) == c(3, 5))
  stopifnot(abs(exp.out[c(3,5)] - c(60.49395, 25.50605)) < 0.001)

  exp.out <- sparse.assign.activity(aa.hcc[, 'JH50'], cosmic.liver, max.level=2)
  stopifnot(which(exp.out > 0.1) == c(2, 3, 4, 7, 8))
  stopifnot(abs(exp.out - c(0, 26.29125, 50.00592, 11.09546, 0, 0, 63.55263, 20.05474)) < 0.001)

  exp.out <- sparse.assign.activity(aa.hcc[, 'JH50'], cosmic.liver)
  stopifnot(which(exp.out > 0.1) == c(2,3,7))
  stopifnot(abs(exp.out[c(2,3,7)] - c(45.84228, 64.42520, 60.73251)) < 0.0001)

  exp.out <- sparse.assign.activity(aa.hcc[, 'JH32'], cosmic.liver)
  stopifnot(which(exp.out > 0.1) == 3)
  stopifnot(exp.out[3] == 95)
}
test.sparse.assign.activity()

