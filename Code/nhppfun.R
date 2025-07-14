# This file contains all functions used for NHPP modeling

library(gtools)
library(splines)
library(dtt)
ngrid = 1000

hppSP = function(interval, lambda)
{
  # generate sample path for an HPP
  # Input:
  #         interval: a vector of length 2 specifying the starting and ending time of the HPP
  #         lambda: a positive scalar specifying the rate of the HPP
  # Output:
  #         a list of 2 components
  #           S: a vector of arrival times of the HPP, length npoints + 1
  #           X: a vector of interarrival times of the HPP, length npoints
  npoints = rpois(1, lambda * (interval[2] - interval[1]))
  S = sort(runif(npoints, interval[1], interval[2]))
  X = c(S[1], diff(S))
  return(list(S = S, X = X))
}

nhppSP = function(interval, lambda.fun)
{
  # generate sample path for an NHPP
  # Input:
  #         interval: a vector of length 2 specifying the starting and ending time of the NHPP
  #         lambda.fun: intensity function of the NHPP, which has been normalized into [0, 1]
  # Output:
  #         a list of 4 components
  #           X.hpp: a vector of interarrival times of the underlying HPP, length npoints
  #           S.hpp: a vector of arrival times of the underlying HPP, length npoints + 1
  #           X.nhpp: a vector of interarrival times of the NHPP, length npoints'
  #           S.nhpp: a vector of arrival times of the NHPP, length npoints' + 1
  delta.vec = seq(interval[1], interval[2], length.out = ngrid + 1)
  lambdafun.vec = lambda.fun(delta.vec)
  if (sum(lambdafun.vec < 0) > 0 | sum(lambdafun.vec > 1) > 0)
  {
    stop("Prescreening error: Intensity function has not been normalized into [0, 1]!")
  }
  lambda = 1
  hpp = hppSP(interval, lambda)
  S.hpp = hpp$S
  X.hpp = hpp$X
  prob.vec = lambda.fun(S.hpp) / lambda
  prob.vec[prob.vec < 0] = 0
  prob.vec[prob.vec > 1] = 1
  npoints = length(S.hpp)
  flag = (rbinom(npoints, 1, prob.vec) == 1)
  S.nhpp = S.hpp[flag]
  X.nhpp = c(S.nhpp[1], diff(S.nhpp))
  return(list(X.hpp = X.hpp, S.hpp = S.hpp, X.nhpp = X.nhpp, S.nhpp = S.nhpp))
}

demoscen = function(scenario, interval, ratio.vec)
{
  # demonstrate MLE/GoF based on single sample path, called by drawfig1&2.R
  # Input:
  #         scenario: an integer 1, 2, or 3
  #         interval: a vector of length 2 specifying the starting and ending locations
  #         ratio.vec: a vector of noise proportions in interarrival times
  # Output:
  #         a list of 5 components
  #           arrtime.vec: a vector of arrival times of the NHPP sample path
  #           lambda.fun: intensity function of the NHPP
  #           lambda.hat.fun: a function of the estimated intensity
  #           interarrtime.lst: a list of length nlength(ratio.vec) + 1
  #                             each contains a vector of interarrival times
  #           cdf.lst: a list of length nlength(ratio.vec) + 1
  #                    each contains a vector of cdf of interarrival times, the cdf is obtained by using estimated intensity
  n = length(ratio.vec)
  intarrtime.lst = vector("list", n + 1)
  cdf.lst = vector("list", n + 1)
  
  if (scenario == 1)
  {
    set.seed(0)
    lambda.fun = function(t) {(0.5 * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + 0.5)}
    dim1 = 64; dim2 = 2
  }
  if (scenario == 2)
  {
    set.seed(0)
    a.2 = runif(4, -2, 2)
    lambda.fun0 = function(t)
    {
      a.2[1] + a.2[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a.2[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a.2[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))
    }
    lb.2 = lambda.fun0(0)
    ub.2 = optimize(lambda.fun0, interval, maximum = T)$objective
    lambda.fun = function(t)
    {
      ((a.2[1] + a.2[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a.2[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a.2[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))) - lb.2) / (ub.2 - lb.2)
    }
    dim1 = 32; dim2 = 4
  }
  if (scenario == 3)
  {
    set.seed(2)
    a.3 = c(runif(1, 0, 4), runif(1, -2, 0))
    lambda.fun0 = function(t) {a.3[1] * cos(2 * pi * 1.5 * t / (interval[2] - interval[1])) + a.3[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1]))}
    lb.3 = lambda.fun0(1000)
    ub.3 = optimize(lambda.fun0, interval, maximum = T)$objective
    lambda.fun1 = function(t)
    {
      ((a.3[1] * cos(2 * pi * 1.5 * t / (interval[2] - interval[1])) + a.3[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1]))) - lb.3) / (ub.3 - lb.3)
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
  temp = nhppSP(interval, lambda.fun)
  arrtime.vec = temp$S.nhpp
  intarrtime.vec = sort(temp$X.nhpp)
  lambda.hat.fun = nhppmle(interval, arrtime.vec, dim1, dim2)
  intarrtime.lst[[1]] = intarrtime.vec
  cdf = Vectorize(function(x) {nhppcdf(x, interval[2], lambda.hat.fun)})
  cdf.vec = cdf(intarrtime.vec)
  cdf.lst[[1]] = cdf.vec
  nevents = length(intarrtime.vec)
  for (i in 1 : n)
  {
    ratio = ratio.vec[i]
    nsubstitute = round(nevents * ratio)
    intarrtime.part.vec = intarrtime.vec[-sample(nevents, nsubstitute)]
    intarrtime.mix.vec = sort(c(intarrtime.part.vec, runif(nsubstitute, min(intarrtime.part.vec), max(intarrtime.part.vec))))
    arrtime.mix.vec = cumsum(intarrtime.mix.vec)
    arrtime.mix.vec = arrtime.mix.vec[arrtime.mix.vec <= interval[2]]
    lambda.hat.mix.fun = nhppmle(interval, arrtime.mix.vec, dim1, dim2)
    intarrtime.lst[[i + 1]] = intarrtime.mix.vec
    cdf = Vectorize(function(x) {nhppcdf(x, interval[2], lambda.hat.mix.fun)})
    cdf.mix.vec = cdf(intarrtime.mix.vec)
    cdf.lst[[i + 1]] = cdf.mix.vec
  }

  return(list(arrtime.vec = arrtime.vec, lambda.fun = lambda.fun, lambda.hat.fun = lambda.hat.fun, intarrtime.lst = intarrtime.lst, cdf.lst = cdf.lst))
}

simu1scen = function(scenario, interval, nrep, npath)
{
  # simulation 1: calculate MISE and PNP
  # Input:
  #         scenario: an integer 1, 2, or 3
  #         interval: a vector of length 2 specifying the starting and ending locations
  #         nrep: a scalar, number of replicated simulations
  #         npath: a scalar, number of sample paths
  # Output:
  #         a list of 5 components, MISE and PNP
  if (scenario == 1)
  {
    set.seed(0)
    lambda.fun = function(t) {(0.5 * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + 0.5)}
    dim1 = 64; dim2 = 2
  }
  if (scenario == 2)
  {
    set.seed(0)
    a.2 = runif(4, -2, 2)
    lambda.fun0 = function(t)
    {
      a.2[1] + a.2[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a.2[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a.2[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))
    }
    lb.2 = lambda.fun0(0)
    ub.2 = optimize(lambda.fun0, interval, maximum = T)$objective
    lambda.fun = function(t)
    {
      ((a.2[1] + a.2[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1])) + a.2[3] * cos(2 * pi * 3 * t / (interval[2] - interval[1])) + a.2[4] * cos(2 * pi * 4 * t / (interval[2] - interval[1]))) - lb.2) / (ub.2 - lb.2)
    }
    dim1 = 32; dim2 = 4
  }
  if (scenario == 3)
  {
    set.seed(2)
    a.3 = c(runif(1, 0, 4), runif(1, -2, 0))
    lambda.fun0 = function(t) {a.3[1] * cos(2 * pi * 1.5 * t / (interval[2] - interval[1])) + a.3[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1]))}
    lb.3 = lambda.fun0(1000)
    ub.3 = optimize(lambda.fun0, interval, maximum = T)$objective
    lambda.fun1 = function(t)
    {
      ((a.3[1] * cos(2 * pi * 1.5 * t / (interval[2] - interval[1])) + a.3[2] * cos(2 * pi * 2 * t / (interval[2] - interval[1]))) - lb.3) / (ub.3 - lb.3)
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
  arrtime.lst = vector("list", nrep)
  lambda.hat.lst = vector("list", nrep)
  ISE.vec = rep(NA, nrep)
  pval.vec = rep(NA, nrep)
  pval.mat = matrix(NA, nrep, npath)
  set.seed(0)
  for (i in 1 : nrep)
  {
    if (i %% 10 == 0) print(i)
    if (npath == 1)
    {
      temp = nhppSP(interval, lambda.fun)
      arrtime.vec = temp$S.nhpp
      arrtime.lst[[i]] = arrtime.vec
      lambda.hat.fun = nhppmle(interval, arrtime.vec, dim1, dim2)
      diff.fun = function(x) {(lambda.hat.fun(x) - lambda.fun(x)) ^ 2}
      ISE.vec[i] = Lambda(lambda.fun = diff.fun, interval[1], interval[2], lambda.vec = NULL, usefun = T)
      pval.vec[i] = nhppGoF(interval, lambda.hat.fun, arrtime.vec)
    }
    if (npath > 1)
    {
      arrtimei.lst = vector("list", npath)
      for (j in 1 : npath)
      {
        temp = nhppSP(interval, lambda.fun)
        arrtime.vec = temp$S.nhpp
        arrtimei.lst[[j]] = arrtime.vec
      }
      arrtime.lst[[i]] = arrtimei.lst
      lambda.hat.fun = nhppmle(interval, arrtimei.lst, dim1, dim2)
      diff.fun = function(x) {(lambda.hat.fun(x) - lambda.fun(x)) ^ 2}
      ISE.vec[i] = Lambda(lambda.fun = diff.fun, interval[1], interval[2], lambda.vec = NULL, usefun = T)
      pval.mat[i, ] = nhppGoF(interval, lambda.hat.fun, arrtimei.lst)
    }
    lambda.hat.lst[[i]] = lambda.hat.fun
  }
  MISE = mean(ISE.vec)
  if (npath == 1)
  {
    PNP = mean(pval.vec > 0.05)
  }
  if (npath > 1)
  {
    PNP = mean(apply(matrix(pval.mat, ncol = 5), 1, median) > 0.05)
  }
  
  return(list(MISE, PNP, ISE.vec, pval.vec, pval.mat))
}

idct = function(t, dct.vec, interval)
{
  # perform inverse DCT for a vector in frequency domain
  # Input:
  #         t: a scalar, dummy variable
  #         dct.vec, a vector of DCT coefficients
  #         interval: a vector of length 2 specifying the starting and ending points
  # Output:
  #         a function of t which is the inverse DCT of dct.vec
  N = length(dct.vec)
  x = dct.vec[1] / 2
  ix.vec = setdiff(which(dct.vec != 0), 1)
  for (i in ix.vec)
  {
    x = x + dct.vec[i] * cos (pi * (i - 1) * t / (interval[2] - interval[1]))
  }
  return(x * 2 / N)
}

nhpplikeFD = function(interval, lambda.dct.vec, arrtime.vec)
{
  # calculate log-likelihood (in frequency domain) for a sample path of an NHPP
  # Input:
  #         interval: a vector of length 2 specifying the starting and ending time of the NHPP
  #         lambda.dct.vec: a vector of the intensity function in frequency domain
  #         arrtime.vec: a vector of the arrival times of the NHPP sample path
  # Output:
  #         a scalar of the log-likelihood
  lambda.vec = dct(lambda.dct.vec, inverted = T)
  if (sum(lambda.vec < 0) > 0)
  {
    ll = -1e6
  } else
  {
    dim1 = length(lambda.dct.vec)
    term1 = lambda.dct.vec[1] / dim1 * (interval[2] - interval[1])
    t.vec = seq(interval[1], interval[2], length.out = dim1)
    lambdaS.vec = approx(t.vec, lambda.vec, arrtime.vec)$y
    # lambdaS.vec[lambdaS.vec <= 0] = 1e-6
    term2 = sum(log(lambdaS.vec))
    # lambdaS.vec = sapply(arrtime.vec, idct, dct.vec = lambda.dct.vec, interval = interval)
    # lambda.fun = function(t) {idct(t, lambda.dct.vec, interval)}
    # lambdaS.vec = lambda.fun(arrtime.vec)
    # lambdaS.vec[lambdaS.vec <= 0] = 1e-6
    term2 = sum(log(lambdaS.vec))
    ll = term2 - term1
  }
  return(ll)
}

nhppmle = function(interval, arrtime, dim1, dim2)
{
  # calculate MLE for the intensity of NHPP based on the a single sample path
  # Input: 
  #         interval: a vector of length 2 specifying the starting and ending time of the NHPP
  #         arrtime: either a vector of the arrival times of the NHPP sample path, or a list of vectors
  #         dim1: a scalar, dimension of the original intensity function
  #         dim2: a scalar, reduced dimension of the intensity function
  # Output: 
  #         a function of the estimated intensity
  ix.mat = combinations(dim1 - 1, dim2 - 1, 2 : dim1)
  nmodel = nrow(ix.mat)
  nll.vec = rep(NA, nmodel)
  para.hat.mat = matrix(NA, nmodel, dim2)
  if (is.list(arrtime))
  {
    a0 = mean(sapply(arrtime, length)) / (interval[2] - interval[1]) * dim1
  } else
  {
    a0 = length(arrtime) / (interval[2] - interval[1]) * dim1
  }
  lambda.dct.hat.vec = c(a0, rep(0, dim1 - 1))#rep(0, dim1)
  for (i in 1 : nmodel)
  {
    pos.vec = c(1, ix.mat[i, ])
    lambda.dct.vec = c(a0, rep(0, dim1 - 1))#rep(0, dim1)
    obj = function(amp.vec, interval, lambda.dct.vec, pos.vec, arrtime)
    {
      lambda.dct.vec[pos.vec] = amp.vec
      # lambda.fun = function(t){idct(t, lambda.dct.vec, interval)}
      if (is.list(arrtime))
      {
        nsample = length(arrtime)
        ll.vec = rep(NA, nsample)
        for (j in 1 : nsample)
        {
          arrtime.vec = arrtime[[j]]
          # ll.vec[j] = nhpplike(interval, lambda.fun, arrtime.vec)
          ll.vec[j] = nhpplikeFD(interval, lambda.dct.vec, arrtime.vec)
        }
        nll = -sum(ll.vec)
      } else
      {
        # nll = -nhpplike(interval, lambda.fun, arrtime)
        nll = -nhpplikeFD(interval, lambda.dct.vec, arrtime)
      }
      return(nll)
    }
    res = optim(c(a0, rep(0, dim2 - 1)), obj, NULL, method = "L-BFGS-B", lower = rep(-100, dim2), upper = rep(100, dim2), interval = interval, lambda.dct.vec = lambda.dct.vec, pos.vec = pos.vec, arrtime = arrtime)
    para.hat.mat[i, ] = res$par
    nll.vec[i] = res$value
  }
  ix = which.min(nll.vec)
  pos.vec = c(1, ix.mat[ix, ])
  amp.vec = para.hat.mat[ix, ]
  lambda.dct.hat.vec[pos.vec] = amp.vec
  lambda.hat.fun = function(t){idct(t, lambda.dct.hat.vec, interval)}
  # return(list(lambda.hat.fun, lambda.dct.hat.vec, pos.vec, amp.vec))
  return(lambda.hat.fun)
}

calMISE = function(fun, funhat.lst, interval)
{
  # calculate MISE for estimating fun
  # Input: 
  #         fun: function to be estimated
  #         funhat.lst: a list of estimated functions in replicated simulations
  #         interval: a vector of length 2 specifying the starting and ending point of the function
  # Output: 
  #         a scalar of MISE
  nrep = length(funhat.lst)
  ISE.vec = rep(NA, nrep)
  for (i in 1 : nrep)
  {
    funhat = funhat.lst[[i]]
    diff.fun = function(x) {(funhat(x) - fun(x)) ^ 2}
    # integrate(diff.fun, interval[1], interval[2])
    integval = Lambda(lambda.fun = diff.fun, interval[1], interval[2], lambda.vec = NULL, usefun = T)
    ISE.vec[i] = integval
  }
  return(mean(ISE.vec))
}

vec2fun = function(x, interval, lambda.vec)
{
  # convert vector to function by interpolation
  # Input:
  #         x: a scalar/vector that is used to evaluate the function
  #         interval: a vector of length 2 specifying the support of the vector
  #         lambda.vec: a vector supported on interval
  # Output:
  #         a scalar/vector of the function value evaluated at x
  ndim = length(lambda.vec)
  t.vec = seq(interval[1], interval[2], length.out = ndim)
  fx = approx(t.vec, lambda.vec, x)$y
  return(fx)
}

Lambda = function(lambda.fun = NULL, t1, t2, lambda.vec = NULL, usefun = T)
{
  # calculate Lambda function (definite integral of lambda.fun)
  # Input:
  #         lambda.fun: intensity function of the NHPP
  #         t1: a scalar, starting point for numerical integration
  #         t2: a scalar, ending point for numerical integration
  #         lambda.vec: a vectorized intensity function, note, this argument is crucial for nhppmle() and nhpplike()
  #         usefun, a Boolean scalar, whether use lambda.fun or lambda.vec to do numerical integration
  # Output:
  #         a scalar of the integration
  if (usefun)
  {
    integ = integrate(Vectorize(lambda.fun), lower = t1, upper = t2, stop.on.error = F)
    if (integ$message == "maximum number of subdivisions reached") # when integrate() fails, do self calculation
    {
      delta.vec = seq(t1, t2, length.out = ngrid + 1)
      lambdafun.vec = lambda.fun(delta.vec)
      integval = (t2 - t1) / ngrid * (sum(lambdafun.vec[-c(1, ngrid + 1)]) + lambdafun.vec[1] / 2 + lambdafun.vec[ngrid + 1] / 2)
    } else
    {
      integval = integ$value
    }
  } else
  {
    t.vec = seq(t1, t2, length.out = length(lambda.vec))
    delta.vec = seq(t1, t2, length.out = ngrid + 1)
    lambdafun.vec = approx(t.vec, lambda.vec, delta.vec)$y
    integval = (t2 - t1) / ngrid * (sum(lambdafun.vec[-c(1, ngrid + 1)]) + lambdafun.vec[1] / 2 + lambdafun.vec[ngrid + 1] / 2)
  }
  return(integval)
}

nhppcdf = function(x, TT, lambda.fun)
{
  # evaluate interarrival time CDF on a given x for an NHPP
  # Input:
  #         x: a scalar within interval [0, TT]
  #         TT: a scalar specifying the ending time (of the NHPP)
  #         lambda.fun: intensity function of the NHPP
  # Output:
  #         a scalar of F(x)
  if (x == 0)
  {
    return(0)
  } else
  {
    denom = Lambda(lambda.fun = lambda.fun, 0, TT, lambda.vec = NULL, usefun = T)
    obj.fun = Vectorize(function(y)
    {
      expterm = Lambda(lambda.fun = lambda.fun, y, y + x, lambda.vec = NULL, usefun = T)
      # lambda.fun = function(x) {vec2fun(x, c(0, TT), lambda.vec)}
      return(lambda.fun(y + x) * exp(-expterm))
    })
    numer = Lambda(lambda.fun = obj.fun, 0, TT - x, lambda.vec = NULL, usefun = T)
    return(1 - numer / denom)
  }
}

nhppGoF = function(interval, lambda.fun, arrtime)
{
  # test GoF of an NHPP with vectorized intensity to sample path(s)
  # Input:
  #         interval: a vector of length 2 specifying the starting and ending time of the NHPP
  #         lambda.fun: intensity function of the NHPP
  #         arrtime: either a vector of the arrival times of the NHPP sample path, or a list of vectors
  # Output:
  #         a scalar/vector of p-value(s)
  TT = interval[2]
  cdf = Vectorize(function(x) {nhppcdf(x, TT, lambda.fun)})
  if (is.list(arrtime))
  {
    nsample = length(arrtime)
    pval.vec = rep(NA, nsample)
    for (j in 1 : nsample)
    {
      arrtime.vec = arrtime[[j]]
      intarrtime.vec = c(arrtime.vec[1], diff(arrtime.vec))
      intarrtime.vec = intarrtime.vec[!duplicated(intarrtime.vec)]
      test = ks.test(intarrtime.vec, cdf)
      pval.vec[j] = test$p.value
    }
    return(pval.vec)
  } else
  {
    intarrtime.vec = c(arrtime[1], diff(arrtime))
    intarrtime.vec = intarrtime.vec[!duplicated(intarrtime.vec)]
    test = ks.test(intarrtime.vec, cdf)
    pval = test$p.value
    return(pval)
  }
}

calPNP = function(arrtime.lst, lambda.hat.lst, interval)
{
  # calculate PNP for evaluating GoF
  # Input: 
  #         arrtime.lst: a list of the arrival times of the NHPP sample path in replicated simulations
  #         lambda.hat.lst: a list of estimated intensity functions in replicated simulations
  #         interval: a vector of length 2 specifying the starting and ending point of the function
  # Output: 
  #         a scalar of MISE
  nrep = length(arrtime.lst)
  if (length(lambda.hat.lst) != nrep)
  {
    stop("Input error!")
  }
  # pval.lst = vector("list", nrep)
  # runt.vec = rep(NA, nrep)
  pval.vec = rep(NA, nrep)
  for (i in 1 : nrep)
  {
    arrtime = arrtime.lst[[i]]
    lambda.hat = lambda.hat.lst[[i]]
    # runt = system.time({pval = nhppGoF(interval, lambda.hat, arrtime)})
    # pval = nhppGoF(interval, lambda.hat, arrtime)
    # pval.lst[[i]] = pval
    # runt.vec[i] = runt[3]
    pval.vec[i] = nhppGoF(interval, lambda.hat, arrtime)
  }
  # return(list(pval.lst, runt.vec))
  return(mean(pval.vec > 0.05))
}

calGoFmerge = function(arrtime1.lst, lambda1.hat.lst, arrtime2.lst, lambda2.hat.lst, arrtime.lst, lambda.hat.lst, interval)
{
  # calculate GoF p-values for N1(t), N2(t), and N(t)
  # Input: 
  #         arrtime1.lst: a list of the arrival times of the NHPP1 sample path
  #         lambda1.hat.lst: a list of estimated intensity functions for NHPP1
  #         arrtime2.lst: a list of the arrival times of the NHPP2 sample path
  #         lambda2.hat.lst: a list of estimated intensity functions for NHPP2
  #         arrtime.lst: a list of the arrival times of the merged NHPP sample path
  #         lambda.hat.lst: a list of estimated intensity functions for merged NHPP
  #         interval: a vector of length 2 specifying the starting and ending point of the function
  # Output: 
  #         a scalar of MISE
  nrep = length(arrtime1.lst)
  if (length(lambda1.hat.lst) != nrep || length(arrtime2.lst) != nrep || length(lambda2.hat.lst) != nrep || length(arrtime.lst) != nrep || length(lambda.hat.lst) != nrep)
  {
    stop("Input error!")
  }
  pval1.vec = rep(NA, nrep)
  pval2.vec = rep(NA, nrep)
  pval.vec = rep(NA, nrep)
  for (i in 1 : nrep)
  {
    arrtime1 = arrtime1.lst[[i]]
    lambda1.hat = lambda1.hat.lst[[i]]
    pval1.vec[i] = nhppGoF(interval, lambda1.hat, arrtime1)
    arrtime2 = arrtime2.lst[[i]]
    lambda2.hat = lambda2.hat.lst[[i]]
    pval2.vec[i] = nhppGoF(interval, lambda2.hat, arrtime2)
    arrtime = arrtime.lst[[i]]
    lambda.hat = lambda.hat.lst[[i]]
    pval.vec[i] = nhppGoF(interval, lambda.hat, arrtime)
    if (i %% 20 == 0) print(i)
  }
  return(list(pval1.vec = pval1.vec, pval2.vec = pval2.vec, pval.vec = pval.vec))
}

nhpp2SPfilter = function(interval, lambda.fun, p, rho)
{
  # generate one NHPP sample path, then another correlated NHPP sample path
  # Input:
  #         interval: a vector of length 2 specifying the starting and ending time of the NHPPs
  #         lambda.fun: the intensity function of base NHPP
  #         p: a scalar for CB samples
  #         rho: a scalar for CB samples
  # Output:
  #         a list of 4 components
  #           S1: a vector of arrival times of the NHPP1
  #           X1: a vector of interarrival times of the NHPP1
  #           S2: a vector of arrival times of the NHPP2
  #           X2: a vector of interarrival times of the NHPP2
  nhpp = nhppSP(interval, lambda.fun)
  S = nhpp$S.nhpp
  npoints = length(S)
  p.vec = rep(p, 2)
  r.mat = cbind(c(1, rho), c(rho, 1))
  para.lst = caliMVN(p.vec, r.mat)
  CB.mat = sampleCB1(npoints, p.vec, r.mat, para.lst)
  S1 = S[which(CB.mat[, 1] == 1)]
  X1 = c(S1[1], diff(S1))
  S2 = S[which(CB.mat[, 2] == 1)]
  X2 = c(S2[1], diff(S2))
  if ((sum(X1 < 0) > 0) || (sum(X2 < 0) > 0)) stop("attival times not sorted!")
  return(list(S1 = S1, X1 = X1, S2 = S2, X2 = X2))
}

caliMVN = function(p.vec, r.mat)
{
  # this function calculates the correlation matrix of the latent MVN and the
  # thresholding parameters
  # input: 
  #         p.vec: a length-m (m>=2) vector, the Bernoulli mean
  #         r.mat: an m bx m matrix, the correlation matrix
  # output:
  #         a list of 2 components: a matrix of size m by m, and a length-m
  #         vector
  m = length(p.vec)
  if (m < 2 | nrow(r.mat) != m | ncol(r.mat) != m)
  {
    stop("wrong dimension!")
  }else
  {
    d.vec = qnorm(p.vec, lower.tail = F)
    rho.mat = matrix(NA, m, m)
    diag(rho.mat) = 1
    for (i in 1 : (m - 1))
    {
      for (j in (i + 1) : m)
      {
        r = r.mat[i, j]
        rho = uniroot(objfun, r = r, p1 = p.vec[i], p2 = p.vec[j], lower = -1, upper = 1, tol = 1e-9)$root
        rho.mat[i, j] = rho
        rho.mat[j, i] = rho
      }
    }
    rho.mat = ngdef(rho.mat)
    
    return(list(rho.mat, d.vec))
  }
}

sampleCB1 = function(n, p.vec, r.mat, para.lst)
{
  # this function generates random samples for a multivariate cor-Bern RV 
  # (X1,...,Xm) where p(Xi=1)=pi and cor(Xi,Xj)=rij, by thresholding MVN
  # input: 
  #         n: a scalar, the sample size
  #         p.vec: a length-m (m>=2) vector, the Bernoulli mean
  #         r.mat: an m bx m matrix, the correlation matrix
  #         para.lst: a length-2 list, the output of caliMVN()
  # output:
  #         an n by m matrix of cor-Bern random samples
  m = length(p.vec)
  if (m < 2 | nrow(r.mat) != m | ncol(r.mat) != m)
  {
    stop("wrong dimension!")
  }else
  {
    rho.mat = para.lst[[1]]
    d.vec = para.lst[[2]]
    z.mat = rmvnorm(n, sigma = rho.mat)
    x.mat = matrix(as.integer(as.vector(z.mat) >= rep(d.vec, each = n)), ncol = m)
    
    return(x.mat)
  }
}
