# This code produces Tab 1 and Tab 2 in Section 3.1

workpath = "C:/Users/admin-xw/Documents/CurrentProject/NHPPtest"
codepath = paste0(workpath, "/Code")
source(paste(codepath, "nhppfun.R", sep = "/"))
resultpath = paste0(workpath, "/Result")
datapath = paste0(workpath, "/Data")

nscen = 3
interval = c(0, 1000)
nrep = 1000
npath = 1 # 5
runt.vec = rep(NA, nscen)
for (scenario in 1 : nscen)
{
  runt = system.time({res = simu1scen(scenario, interval, nrep, npath)})
  # runt = system.time({arrtime.lst = simu1scen(scenario, interval, nrep, npath)})
  # save(arrtime.lst, file = paste0(datapath, "/simu_arrtime_scen", scenario, "_npath", npath, ".RData"))
  print(paste("MISE for Scenario ", scenario, ": ", res[[1]], sep = ""))
  print(paste("PNP for Scenario ", scenario, ": ", res[[2]], sep = ""))
  runt.vec[scenario] = runt[3]
}
print(paste("Running time: ", runt.vec, sep = ""))
