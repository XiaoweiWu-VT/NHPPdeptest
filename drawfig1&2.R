# This code produces Fig 1 and Fig 2 in Section 3.1
# Figure 1, left panels show single sample path for Scenarios I, II, and III
# right panels show MLE results of the intensity functions
# Figure 2 shows CDF results of the interarrival times

library(latex2exp)
codepath = "C:/Users/admin-xw/Documents/CurrentProject/NHPPtest/Code"
source(paste(codepath, "nhppfun.R", sep = "/"))
resultpath = "C:/Users/admin-xw/Documents/CurrentProject/NHPPtest/Result"

nscen = 3
interval = c(0, 1000)
ratio.vec = c(0.05, 0.1)
arrtime.lst = vector("list", nscen)
lambda.lst = vector("list", nscen)
lambda.hat.lst = vector("list", nscen)
intarrtime.lst = vector("list", nscen * (length(ratio.vec) + 1))
cdf.lst = vector("list", nscen * (length(ratio.vec) + 1))
runt.vec = rep(NA, nscen)
for (scenario in 1 : nscen)
{
  runt = system.time({res = demoscen(scenario, interval, ratio.vec)})
  runt.vec[scenario] = runt[3]
  arrtime.lst[[scenario]] = res$arrtime.vec
  lambda.lst[[scenario]] = res$lambda.fun
  lambda.hat.lst[[scenario]] = res$lambda.hat.fun
  intarrtime.lst[[(scenario - 1) * (length(ratio.vec) + 1) + 1]] = res$intarrtime.lst[[1]]
  cdf.lst[[(scenario - 1) * (length(ratio.vec) + 1) + 1]] = res$cdf.lst[[1]]
  for (i in 1 : length(ratio.vec))
  {
    intarrtime.lst[[(scenario - 1) * (length(ratio.vec) + 1) + i + 1]] = res$intarrtime.lst[[i + 1]]
    cdf.lst[[(scenario - 1) * (length(ratio.vec) + 1) + i + 1]] = res$cdf.lst[[i + 1]]
  }
}
print(runt.vec) # 75.13 662.47  83.61

setEPS()
postscript(paste(resultpath, "/Fig1.eps", sep = ""))
# x11()
par(mfrow = c(3, 2), las = 1)
for (scenario in 1 : nscen)
{
  arrtime.vec = arrtime.lst[[scenario]]
  lambda.fun = lambda.lst[[scenario]]
  lambda.hat.fun = lambda.hat.lst[[scenario]]
  plot(arrtime.vec, rep(0.5, length(arrtime.vec)), type = "p", ylim = c(0, 1), xlab = "Arrival times in single sample path", ylab = "", yaxt = "n", main = chr(64 + 2 * scenario - 1))
  curve(lambda.fun, interval[1], interval[2], xlab = "True and estimated intensity functions", ylab = "", main = chr(64 + 2 * scenario))
  curve(lambda.hat.fun, interval[1], interval[2], col = "red", add = T)
  if (scenario == 1)
  {
    legend("bottom", TeX(c("$\\lambda(t)$", "$\\hat{\\lambda}(t)$")), lty = c(1, 1), col = c("black", "red"), cex = 0.9)
  }
}
dev.off()

setEPS()
postscript(paste(resultpath, "/Fig2.eps", sep = ""))
# x11()
par(mfrow = c(3, 3), las = 1)
for (scenario in 1 : nscen)
{
  plot.ecdf(intarrtime.lst[[(scenario - 1) * (length(ratio.vec) + 1) + 1]], xlab = TeX("$t$"), ylab = TeX("$F(t)$"), main = chr(64 + 3 * i - 2), cex.main = 1, cex.lab = 1)
  lines(intarrtime.lst[[(scenario - 1) * (length(ratio.vec) + 1) + 1]], cdf.lst[[(scenario - 1) * (length(ratio.vec) + 1) + 1]], type = "l", col ="red")
  if (scenario == 1)
  {
    legend("bottomright", legend = c("Empirical CDF", "CDF based on Eq. (4)"), col = c("black", "red"), lty = c(1, 1), pch = c(16, NA), cex = 0.9)
  }
  for (i in 1 : length(ratio.vec))
  {
    plot.ecdf(intarrtime.lst[[(scenario - 1) * (length(ratio.vec) + 1) + i + 1]], xlab = TeX("$t$"), ylab = TeX("$F(t)$"), main = chr(64 + 3 * i - 1), cex.main = 1, cex.lab = 1)
    lines(intarrtime.lst[[(scenario - 1) * (length(ratio.vec) + 1) + i + 1]], cdf.lst[[(scenario - 1) * (length(ratio.vec) + 1) + i + 1]], type = "l", col ="red")
  }
}
dev.off()
