# This code produces Tab 1 in Section 3.1 and Tab S1 in Supp

codepath = "C:/Users/admin-xw/Documents/CurrentProject/NHPPtest/Code"
source(paste(codepath, "nhppfun.R", sep = "/"))
resultpath = "C:/Users/admin-xw/Documents/CurrentProject/NHPPtest/Result"

nscen = 3
interval = c(0, 1000)
nrep = 1000
npath = 1 # 5
runt.vec = rep(NA, nscen)
# generate sample path
for (scenario in 2 : 2)#nscen)
{
  runt = system.time({res = simu1scen(scenario, interval, nrep, npath)})
  print(paste("MISE for Scenario ", scenario, ": ", res[[1]], sep = ""))
  print(paste("PNP for Scenario ", scenario, ": ", res[[2]], sep = ""))
  runt.vec[scenario] = runt[3]
}
print(paste("Running time: ", runt.vec, sep = ""))
