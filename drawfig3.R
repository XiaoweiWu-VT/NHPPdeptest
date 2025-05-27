# This code produces Fig 3 in Section 4 (Real Data Application)

workpath = "C:/Users/admin-xw/Documents/CurrentProject/NHPPtest"
# read 13 ChIP-seq data (excluding CTCF, loc unordered)
data1.file = paste(workpath, "/RealData/mmc3.csv", sep = "")
data1.df = read.csv(data1.file)
nTF1 = ncol(data1.df)
bsmat1.lst = vector("list", nTF1)
names(bsmat1.lst) = unlist(strsplit(names(data1.df), split = ".bound.loci"))
for (i in 1 : nTF1)
{
  data.TF = data1.df[, i]
  data.TF.split = strsplit(data.TF, split = "[:-]")
  data.TF.mat = matrix(unlist(data.TF.split), nrow = 3)
  data.TF.posmat = matrix(as.numeric(data.TF.mat[2 : 3, ]), nrow = 2)
  data.TF.posvec = as.integer(round(colMeans(data.TF.posmat)))
  bsmat1.lst[[i]] = data.frame(data.TF.mat[1, ], data.TF.posvec)
  colnames(bsmat1.lst[[i]]) = c("Chromosome", "Coordinate")
}

# read Nr5a2 ChIP-seq data
data2.file = paste(workpath, "/RealData/mmc2.csv", sep = "")
data2.df = read.csv(data2.file, skip = 1)
data2.mat = data2.df[, 2 : 3]

# read Tcf3 ChIP-seq data
data3.file = paste(workpath, "/RealData/mES_OCT4_SOX2_NANOG_TCF3.WIG", sep = "")
data3.df = read.csv(data3.file, header = F, skip = 9381282)
data.TF.split = strsplit(data3.df[, 1], split = "\t")
data.TF.mat = matrix(unlist(data.TF.split), nrow = 3)
data.TF.posmat = matrix(as.numeric(data.TF.mat[2 : 3, ]), nrow = 2)
data.TF.posvec = as.integer(round(colMeans(data.TF.posmat)))
data3.mat = data.frame(data.TF.mat[1, ], data.TF.posvec)
colnames(data3.mat) = c("Chromosome", "Coordinate")

# combine ChIP-seq data
nTF = nTF1 + 1
bsmat.lst = vector("list", nTF)
TFname.vec = rep(NA, nTF)
TFname1.vec = names(bsmat1.lst)
j = 0
for (i in 1 : nTF1)
{
  if (TFname1.vec[i] != "CTCF")
  {
    j = j + 1
    bsmat.lst[[j]] = bsmat1.lst[[i]]
    TFname.vec[j] = TFname1.vec[i]
  }
}
bsmat.lst[[13]] = data2.mat
TFname.vec[13] = "Nr5a2"
bsmat.lst[[14]] = data3.mat
TFname.vec[14] = "Tcf3"
names(bsmat.lst) = TFname.vec

bsmat.Chr1.lst = vector("list", nTF)
for (i in 1 : nTF)
{
  bsmati = bsmat.lst[[i]]
  bsmat.Chr1.lst[[i]] = sort(bsmati[bsmati[, 1] == "chr1", 2])
}
l = min(sapply(bsmat.Chr1.lst, min))
u = max(sapply(bsmat.Chr1.lst, max))
nseg = 50
break.vec = seq(l, u, length.out = nseg + 1)
break.mat = rbind(break.vec[-(nseg + 1)], break.vec[-1])

source(paste(workpath, "/Code/nhppfun.R", sep = ""))
dim1 = 64; dim2 = 3
runt.vec = rep(NA, nseg)
npassTF.vec = rep(NA, nseg)
narr.mat = matrix(NA, nTF, nseg)
pvalgof.mat = matrix(NA, nTF, nseg)
dep.lst = vector("list", nseg)
interval = c(0, 1000)
for (k in 1 : nseg)
{
  runt = system.time({
    bslock.lst = vector("list", nTF)
    lambda.est.fun.lst = vector("list", nTF)
    for (i in 1 : nTF)
    {
      bsloci = bsmat.Chr1.lst[[i]]
      bslocik = bsloci[bsloci >= break.mat[1, k] & bsloci <= break.mat[2, k]]
      arrtime = interval[1] + (bslocik - break.mat[1, k]) / (break.mat[2, k] - break.mat[1, k]) * interval[2]
      bslock.lst[[i]] = arrtime
      narr.mat[i, k] = length(arrtime)
      runt1 = system.time({lambda.est.fun = nhppmle(interval, arrtime, dim1, dim2)})
      lambda.est.fun.lst[[i]] = lambda.est.fun
      runt2 = system.time({pval = nhppGoF(interval, lambda.est.fun, arrtime)})
      pvalgof.mat[i, k] = pval
    }
    
    ix.nhpp = which(pvalgof.mat[, k] >= 0.05)
    # print(paste("TFs that pass GoF test: ", TFname.vec[ix.nhpp], sep = ""))
    npassTF = length(ix.nhpp)
    npassTF.vec[k] = npassTF
    combmat = combinations(npassTF, 2, v = ix.nhpp)
    ncomb = nrow(combmat)
    pvaldep.mat = matrix(NA, ncomb, 3)
    for (j in 1 : ncomb)
    {
      ix1 = combmat[j, 1]
      ix2 = combmat[j, 2]
      pvaldep.mat[j, 1] = pvalgof.mat[ix1, k]
      pvaldep.mat[j, 2] = pvalgof.mat[ix2, k]
      arrtime1 = bslock.lst[[ix1]]
      arrtime2 = bslock.lst[[ix2]]
      arrtime = sort(c(arrtime1, arrtime2))
      lambda1.est.fun = lambda.est.fun.lst[[ix1]]
      lambda2.est.fun = lambda.est.fun.lst[[ix2]]
      lambdapool.est.fun = function(t)
      {
        f1 = lambda1.est.fun
        f2 = lambda2.est.fun
        f1(t) + f2(t)
      }
      runt3 = system.time({pvaldep.mat[j, 3] = nhppGoF(interval, lambdapool.est.fun, arrtime)})
    }
    ix.dep = which(pvaldep.mat[, 3] < 0.05) # which pairs are dependent?
    dep.mat = matrix(0, nTF, nTF) # diag(nTF)
    colnames(dep.mat) = TFname.vec
    rownames(dep.mat) = TFname.vec
    for (i in 1 : length(ix.dep))
    {
      ix = ix.dep[i]
      dep.mat[combmat[ix, 1], combmat[ix, 2]] = dep.mat[combmat[ix, 1], combmat[ix, 2]] + 1
      dep.mat[combmat[ix, 2], combmat[ix, 1]] = dep.mat[combmat[ix, 2], combmat[ix, 1]] + 1
    }
    # print(dep.mat)
  })
  print(paste("Segment ", k, " is done! Time = ", runt[3], sep = ""))
  dep.lst[[k]] = dep.mat
  runt.vec[k] = runt[3]
}

accdep.mat = matrix(0, nTF, nTF)
for (k in 1 : nseg)
{
  accdep.mat = accdep.mat + dep.lst[[k]]
}

save(dep.lst, narr.mat, pvalgof.mat, npassTF.vec, runt.vec, accdep.mat, file = paste(workpath, "/Result/dep.RData", sep = ""))
