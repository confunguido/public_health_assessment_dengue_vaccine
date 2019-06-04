
files = list.files()[
  which(sapply(strsplit(list.files(),'_'),function(ll)ll[1]=='brazil'&&ll[2]=='baseline'))]

infections = rep(0,10)
clinical = rep(0,10)
for(ff in files){
  x = read.csv(ff)
  infections = infections + table(c(x[,1],0:9)) - 1
  clinical = clinical + table(c(x[x[,3]>0,1],0:9)) - 1
}

infections.brazil.baseline = infections
clinical.brazil.baseline = clinical



files = list.files()[
  which(sapply(strsplit(list.files(),'_'),function(ll)ll[1]=='brazil'&&ll[2]=='nocatchup'))]

infections = rep(0,10)
clinical = rep(0,10)
for(ff in files){
  x = read.csv(ff)
  infections = infections + table(c(x[,1],0:9)) - 1
  clinical = clinical + table(c(x[x[,3]>0,1],0:9)) - 1
}

infections.brazil.nocatchup = infections
clinical.brazil.nocatchup = clinical



files = list.files()[
  which(sapply(strsplit(list.files(),'_'),function(ll)ll[1]=='brazil'&&ll[2]=='catchup'))]

infections = rep(0,10)
clinical = rep(0,10)
for(ff in files){
  x = read.csv(ff)
  infections = infections + table(c(x[,1],0:9)) - 1
  clinical = clinical + table(c(x[x[,3]>0,1],0:9)) - 1
}

infections.brazil.catchup = infections
clinical.brazil.catchup = clinical



files = list.files()[
  which(sapply(strsplit(list.files(),'_'),function(ll)ll[1]=='thailand'&&ll[2]=='baseline'))]

infections = rep(0,10)
clinical = rep(0,10)
for(ff in files){
  x = read.csv(ff)
  infections = infections + table(c(x[,1],0:9)) - 1
  clinical = clinical + table(c(x[x[,3]>0,1],0:9)) - 1
}

infections.thailand.baseline = infections
clinical.thailand.baseline = clinical



files = list.files()[
  which(sapply(strsplit(list.files(),'_'),function(ll)ll[1]=='thailand'&&ll[2]=='nocatchup'))]

infections = rep(0,10)
clinical = rep(0,10)
for(ff in files){
  x = read.csv(ff)
  infections = infections + table(c(x[,1],0:9)) - 1
  clinical = clinical + table(c(x[x[,3]>0,1],0:9)) - 1
}

infections.thailand.nocatchup = infections
clinical.thailand.nocatchup = clinical



files = list.files()[
  which(sapply(strsplit(list.files(),'_'),function(ll)ll[1]=='thailand'&&ll[2]=='catchup'))]

infections = rep(0,10)
clinical = rep(0,10)
for(ff in files){
  x = read.csv(ff)
  infections = infections + table(c(x[,1],0:9)) - 1
  clinical = clinical + table(c(x[x[,3]>0,1],0:9)) - 1
}

infections.thailand.catchup = infections
clinical.thailand.catchup = clinical



infections.brazil.baseline - infections.brazil.nocatchup
clinical.brazil.baseline - clinical.brazil.nocatchup

infections.thailand.baseline - infections.thailand.nocatchup
clinical.thailand.baseline - clinical.thailand.nocatchup


infections.brazil.baseline - infections.brazil.catchup
clinical.brazil.baseline - clinical.brazil.catchup

infections.thailand.baseline - infections.thailand.catchup
clinical.thailand.baseline - clinical.thailand.catchup



save(list=ls(),file='results.RData')
