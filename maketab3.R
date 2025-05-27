# This code produces Tab 3 in Section 3.2

scen = "III"; rho = 1
load(paste("C:/Users/xwwu/Documents/NHPPtest/Result/simu2ha_scen_", scen, "_rho_", rho, ".RData", sep = ""))
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
print(paste("TPR for Scenario ", scen, ": ", mean((pval1.vec > 0.05) & (pval2.vec > 0.05) & (pval.vec <= 0.05)), sep = ""))
print(paste("Adjusted TPR for Scenario ", scen, ": ", sum((pval1.vec > 0.05) & (pval2.vec > 0.05) & (pval.vec <= 0.05)) / sum((pval1.vec > 0.05) & (pval2.vec > 0.05)), sep = ""))
print(paste("Running time for intensity estimation in Scenario ", scen, " NHPP 1: ", sum(runt1.vec), sep = ""))
print(paste("Running time for intensity estimation in Scenario ", scen, " NHPP 2: ", sum(runt3.vec), sep = ""))
print(paste("Running time for GoF test of NHPP 1 in Scenario ", scen, ": ", sum(runt2.vec), sep = ""))
print(paste("Running time for GoF test of NHPP 2 in Scenario ", scen, ": ", sum(runt4.vec), sep = ""))
print(paste("Running time for GoF test of pooled NHPP in Scenario ", scen, ": ", sum(runt5.vec), sep = ""))

# I, 0
# 0.983/0.97/0.93/0.974/1920/1896/28184/28145/76586
# I, -1
# 0.979/0.979/0.019/0.020/1927/1914/28444/28439/94271
# I, 1
# 0.981/0.981/0.981/1/1924/1911/28622/28538/58233
# I, 0.1
# 0.982/0.975/0.952/0.993/1924/1900/28410/28431/75525
# I, -0.9
# 0.976/0.981/0.012/0.013/1919/1901/28316/28344/92297
# I, -0.8
# 0.973/0.984/0.026/0.027/1913/1901/28166/28203/90177
# I, -0.5
# 0.976/0.978/0.315/0.330/1918/1898/28274/28278/85727
# I, -0.1
# 0.981/0.979/0.902/0.939/1919/1900/28278/28306/78661
# I, -0.3
# 0.978/0.98/0.692/0.722/1914/1898/28242/28261/81984
# I, -0.4
# 0.977/0.975/0.508/0.534/1922/1907/28408/28378/84023
# II, -1
# 0.993/0.991/0.011/0.011/198728/198172/36879/36382/122678
# II, -0.5
# 0.99/0.995/0.249/0.253/199356/198686/37248/36528/111314
# II, -0.3
# 0.992/0.993/0.65/0.66/200889/199784/36965/36708/106067
# II, 0
# 0.987/0.986/0.952/0.977/200435/199113/36808/36697/98352
# II, 0.1
# 0.993/0.986/0.976/0.997/200582/199503/36824/36965/96196
# II, 1
# 0.984/0.984/0.984/1/200042/199373/36981/36846/71890
# III, -1
# 0.984/0.986/0.016/0.016/11291/11298/19316/19288/63909
# III, -0.5
# 0.988/0.988/0.365/0.374/11267/11251/19284/19423/58118
# III, -0.3
# 0.986/0.98/0.775/0.802/11285/11279/19431/19336/55553
# III, 0
# 0.987/0.99/0.974/0.997/11373/11361/19530/19436/52119
# III, 0.1
# 0.985/0.991/0.977/0.999/11339/11317/19392/19365/50461
# III, 1
# 0.991/0.991/0.991/1/11332/11301/19314/19263/39191
