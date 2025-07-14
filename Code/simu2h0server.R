# This code is for calculating type I error of simulation 2 (single path)
# The simulation is for Scenarios IV, V, and VI
# nohup R CMD BATCH --no-save --no-restore '--args scen="IV"' simu2h0remote.R logIV.out > lIV&

args = (commandArgs(T))
eval(parse(text = args[[1]]))
library(gtools) # used in nhppmle
library(dtt) # used in nhppmle
library(foreach)
library(doParallel)
library(doRNG)
cl = makeCluster(20)
registerDoParallel(cl)

workpath = "/home/xwwu/NHPPtest"
codepath = paste0(workpath, "/Code")
source(paste(codepath, "nhppfun.R", sep = "/"))
resultpath = paste0(workpath, "/Result")

interval = c(0, 1000)
nrep = 1000

lb = NA; ub = NA; a = NA
if (scen == "IV")
{
  set.seed(0)
  lambda1.fun = function(t) {(0.5 * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + 0.5)}
  dim11 = 64; dim12 = 2
  lambda2.fun = lambda1.fun
  dim21 = 64; dim22 = 2
}
if (scen == "V")
{
  set.seed(0)
  a = runif(4, -2, 2)
  lambda.fun0 = function(t)
  {
    a[1] + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))
  }
  lb = lambda.fun0(0)
  ub = optimize(lambda.fun0, interval, maximum = T)$objective
  lambda1.fun = function(t)
  {
    ((a[1] + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))) - lb) / (ub - lb)
  }
  dim11 = 32; dim12 = 4
  lambda2.fun = lambda1.fun
  dim21 = 32; dim22 = 4
}
if (scen == "VI")
{
  set.seed(0)
  lambda1.fun = function(t) {(0.5 * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + 0.5)}
  dim11 = 64; dim12 = 2
  a = runif(4, -2, 2)
  lambda.fun0 = function(t)
  {
    a[1] + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))
  }
  lb = lambda.fun0(0)
  ub = optimize(lambda.fun0, interval, maximum = T)$objective
  lambda2.fun = function(t)
  {
    ((a[1] + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))) - lb) / (ub - lb)
  }
  dim21 = 32; dim22 = 4
}

outfile = paste0(resultpath, "/simu2h0_scen", scen, ".RData")
res = foreach(i = 1 : nrep, .packages = c('gtools', 'dtt')) %dorng%
{
  temp1 = nhppSP(interval, lambda1.fun)
  arrtime1.vec = temp1$S.nhpp
  runt1 = system.time({lambda1.hat = nhppmle(interval, arrtime1.vec, dim11, dim12)})
  runt2 = system.time({pval1 = nhppGoF(interval, lambda1.hat, arrtime1.vec)})
  temp2 = nhppSP(interval, lambda2.fun)
  arrtime2.vec = temp2$S.nhpp
  runt3 = system.time({lambda2.hat = nhppmle(interval, arrtime2.vec, dim21, dim22)})
  runt4 = system.time({pval2 = nhppGoF(interval, lambda2.hat, arrtime2.vec)})
  arrtime.vec = sort(c(arrtime1.vec, arrtime2.vec))
  lambda.hat = function(t)
  {
    f1 = lambda1.hat
    f2 = lambda2.hat
    f1(t) + f2(t)
  }
  runt5 = system.time({pval = nhppGoF(interval, lambda.hat, arrtime.vec)})
  list(arrtime1.vec, arrtime2.vec, arrtime.vec, lambda1.hat, lambda2.hat, lambda.hat, runt1[3], runt2[3], runt3[3], runt4[3], runt5[3], pval1, pval2, pval)
}    
save(res, file = outfile)
