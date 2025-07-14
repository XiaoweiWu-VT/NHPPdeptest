# This code is for calculating power of simulation 2 (single path)
# The simulation is for Scenarios IV, V, and VI
# nohup R CMD BATCH --no-save --no-restore '--args scen="I" rho=0' simu2haserver.R logI0.out > lI0&

args = (commandArgs(T))
eval(parse(text = args[[1]]))
eval(parse(text = args[[2]]))
library(gtools) # used in nhppmle
library(dtt) # used in nhppmle
library(mvtnorm)
library(binaryLogic)
library(splines)
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
p = 0.5

lb = NA; ub = NA; a = NA
if (scen == "I")
{
  set.seed(0)
  lambda.fun = function(t) {(0.5 * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + 0.5)}
  dim1 = 64; dim2 = 2
}
if (scen == "II")
{
  set.seed(0)
  a = runif(4, -2, 2)
  lambda.fun0 = function(t)
  {
    a[1] + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))
  }
  lb = lambda.fun0(0)
  ub = optimize(lambda.fun0, interval, maximum = T)$objective
  lambda.fun = function(t)
  {
    ((a[1] + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))) - lb) / (ub - lb)
  }
  dim1 = 32; dim2 = 4
}
if (scen == "III")
{
  set.seed(2)
  a = c(runif(1, 0, 4), runif(1, -2, 0))
  lambda.fun0 = function(t) {a[1] * cos(2 * pi * 1.5 * t / (interval[2] - interval[1])) + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1]))}
  lb = lambda.fun0(1000)
  ub = optimize(lambda.fun0, interval, maximum = T)$objective
  lambda.fun1 = function(t)
  {
    ((a[1] * cos(2 * pi * 1.5 * t / (interval[2] - interval[1])) + a[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1]))) - lb) / (ub - lb)
  }
  degree = 6
  x.vec = seq(interval[1], interval[2], length.out = ngrid)
  basis.mat = bs(x.vec, degree = degree, intercept = T)
  y.vec = lambda.fun1(x.vec)
  bsm = lm(y.vec ~ 0 + basis.mat)
  c.vec = coef(bsm)
  lambda.fun = function(x) {vec2fun(x, interval, as.vector(basis.mat %*% c.vec))}
  dim1 = 32; dim2 = 3
}

outfile = paste0(resultpath, "/simu2ha_scen", scen, "_rho_", rho, ".RData")
res = foreach(i = 1 : nrep, .packages = c('gtools', 'dtt', 'mvtnorm', 'binaryLogic', 'splines')) %dorng%
{
  temp = nhpp2SPfilter(interval, lambda.fun, p, rho)
  arrtime1.vec = temp$S1
  runt1 = system.time({lambda1.hat = nhppmle(interval, arrtime1.vec, dim1, dim2)})
  runt2 = system.time({pval1 = nhppGoF(interval, lambda1.hat, arrtime1.vec)})
  arrtime2.vec = temp$S2
  runt3 = system.time({lambda2.hat = nhppmle(interval, arrtime2.vec, dim1, dim2)})
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
