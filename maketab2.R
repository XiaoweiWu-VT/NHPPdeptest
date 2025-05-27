# This code produces Tab 2 in Section 3.2

scen = "VI"
load(paste("C:/Users/xwwu/Documents/NHPPtest/Result/simu2h0_scen_", scen, ".RData", sep = ""))
nrep = 1000
pval1.vec = rep(NA, nrep)
pval2.vec = rep(NA, nrep)
pval.vec = rep(NA, nrep)
runt1.vec = rep(NA, nrep)
runt2.vec = rep(NA, nrep)
runt3.vec = rep(NA, nrep)
runt4.vec = rep(NA, nrep)
runt5.vec = rep(NA, nrep)
for (i in 1 : nrep)
{
  pval1.vec[i] = unlist(res[[i]][12])
  pval2.vec[i] = unlist(res[[i]][13])
  pval.vec[i] = unlist(res[[i]][14])
  runt1.vec[i] = unlist(res[[i]][7])
  runt2.vec[i] = unlist(res[[i]][8])
  runt3.vec[i] = unlist(res[[i]][9])
  runt4.vec[i] = unlist(res[[i]][10])
  runt5.vec[i] = unlist(res[[i]][11])
}
print(paste("PNP for Scenario ", scen, " NHPP 1: ", mean(pval1.vec > 0.05), sep = ""))
print(paste("PNP for Scenario ", scen, " NHPP 2: ", mean(pval2.vec > 0.05), sep = ""))
print(paste("FPR for Scenario ", scen, ": ", mean((pval1.vec > 0.05) & (pval2.vec > 0.05) & (pval.vec <= 0.05)), sep = ""))
print(paste("Adjusted FPR for Scenario ", scen, ": ", sum((pval1.vec > 0.05) & (pval2.vec > 0.05) & (pval.vec <= 0.05)) / sum((pval1.vec > 0.05) & (pval2.vec > 0.05)), sep = ""))
print(paste("Running time for intensity estimation in Scenario ", scen, " NHPP 1: ", sum(runt1.vec), sep = ""))
print(paste("Running time for intensity estimation in Scenario ", scen, " NHPP 2: ", sum(runt3.vec), sep = ""))
print(paste("Running time for GoF test of NHPP 1 in Scenario ", scen, ": ", sum(runt2.vec), sep = ""))
print(paste("Running time for GoF test of NHPP 2 in Scenario ", scen, ": ", sum(runt4.vec), sep = ""))
print(paste("Running time for GoF test of pooled NHPP in Scenario ", scen, ": ", sum(runt5.vec), sep = ""))

# IV
# 0.968/0.948/0.024/0.026/1965//1958//191943
# V
# 0.987/0.981/0.019/0.020/195907//195173//238599
# VI
# 0.944/0.981/0.019/0.020/1955/194505/57546/71271/188672