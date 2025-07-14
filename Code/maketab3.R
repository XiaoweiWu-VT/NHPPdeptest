# This code produces Tab 3 in Section 3.2

workpath = "C:/Users/admin-xw/Documents/CurrentProject/NHPPtest"
datapath = paste0(workpath, "/Data")

nrep = 1000
scen.vec = c("IV", "V", "VI")
FPR.vec = rep(NA, length(scen.vec))
for (j in 1 : length(scen.vec))
{
  scen = scen.vec[j]
  load(paste0(datapath, "/simu2h0_scen", scen, ".RData"))
  pval1.vec = rep(NA, nrep)
  pval2.vec = rep(NA, nrep)
  pval.vec = rep(NA, nrep)
  for (i in 1 : nrep)
  {
    pval1.vec[i] = unlist(res[[i]][12])
    pval2.vec[i] = unlist(res[[i]][13])
    pval.vec[i] = unlist(res[[i]][14])
  }
  FPR.vec[j] = sum((pval1.vec > 0.05) & (pval2.vec > 0.05) & (pval.vec <= 0.05)) / sum((pval1.vec > 0.05) & (pval2.vec > 0.05))
}

# IV
# 0.968/0.948/0.024/0.026/1965//1958//191943
# V
# 0.987/0.981/0.019/0.020/195907//195173//238599
# VI
# 0.944/0.981/0.019/0.020/1955/194505/57546/71271/188672

# arrtime1.lst = vector("list", nrep)
# arrtime2.lst = vector("list", nrep)
# arrtime.lst = vector("list", nrep)
# for (i in 1 : nrep)
# {
#   arrtime1.lst[[i]] = unlist(res[[i]][1])
#   arrtime2.lst[[i]] = unlist(res[[i]][2])
#   arrtime.lst[[i]] = unlist(res[[i]][3])
# }
# save(arrtime1.lst, arrtime2.lst, arrtime.lst, file = paste0(datapath, "/simu2h0_arrtimes_scen", scen, ".RData"))
