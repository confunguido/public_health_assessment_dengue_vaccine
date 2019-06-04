# load packages
library(deSolve)
library(doParallel)
library(foreach)
library(pomp)


# define functions
# functions for tranforming parameters for optim
logit = function(x){log(x/(1-x))}
logit.inv = function(x){1 / (1 + exp(-x))}

# calculate the proportion of a certain age with a given number
# of previous infections assuming a fixed force of infection
calc_seropos = function(previnf,age,foi){
  choose(4,previnf) * (1 - dpois(0,foi * age)) ^ previnf * dpois(0, foi * age) ^ (4 - previnf)
}

# calculate the proportion of a given age group expected to fall into
# 1 of 9 categories once the vaccine has been administered
#   1 = protected, 0 previous infections
#   2 = unprotected, 0 previous infections
#   3 = protected, 1 previous infections
#   4 = unprotected, 1 previous infections
#   5 = protected, 2 previous infections
#   6 = unprotected, 2 previous infections
#   7 = protected, 3 previous infections
#   8 = unprotected, 3 previous infections
#   9 = protected, 4 previous infections
calc_ICs = function(age, foi, efficacy, progression){
  ICs = c(rep(c(efficacy, 1 - efficacy), 4) * calc_seropos(rep(0:3,each=2), age , foi),
    calc_seropos(4, age, foi))
  ICs.out = c(
    (1 - progression) * ICs[1:2],
    progression * ICs[1:2] + (1 - progression) * ICs[3:4],
    progression * ICs[3:4] + (1 - progression) * ICs[5:6],
    progression * ICs[5:6] + (1 - progression) * ICs[7:8],
    sum(progression * ICs[7:8]) + ICs[9])
  return(ICs.out)
}

# function to calculate derivatives for ODE solver
# transitions include waning of vaccine efficacy, infection,
# and conversion to disease of hospitalized status
func.transitions = function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    dP.0.p = -w * P.0.p
    dP.0.u = w * P.0.p - 4 * foi * P.0.u
    dP.1.p = -w * P.1.p
    dP.1.u = w * P.1.p - 3 * foi * P.1.u + 4 * foi * P.0.u
    dP.2.p = -w * P.2.p
    dP.2.u = w * P.2.p - 2 * foi * P.2.u + 3 * foi * P.1.u
    dP.3.p = -w * P.3.p
    dP.3.u = w * P.3.p - 1 * foi * P.3.u + 2 * foi * P.2.u
    dP.4 = foi * P.3.u
    dD = d0 * 4 * foi * P.0.u + d1 * 3 * foi * P.1.u + d2 * (2 * foi * P.2.u + foi * P.3.u)
    dH = h0 * d0 * 4 * foi * P.0.u + h1 * d1 * 3 * foi * P.1.u + h2 * d2 * (2 * foi * P.2.u + foi * P.3.u)
    list(c(dP.0.p,dP.0.u,dP.1.p,dP.1.u,dP.2.p,dP.2.u,dP.3.p,dP.3.u,dP.4,dD,dH))
  })
}

# calculate the proportion expected to have acquired disease or become hospitalized by
# the end of a trial of duration time.max
calc_disease_hosp_by_age = function(age,study.min,study.max,foi,efficacy,progression,w,d0,d1,d2,h0,h1,h2,time.inc=1/100){
  ICs = calc_ICs(age, foi, efficacy, progression)
  state = c(P.0.p = ICs[1],P.0.u = ICs[2],P.1.p = ICs[3],P.1.u = ICs[4],P.2.p = ICs[5],
            P.2.u = ICs[6],P.3.p = ICs[7],P.3.u = ICs[8],P.4 = ICs[9],D = 0,H = 0)
  parameters = c(w = w, foi = foi, d0 = d0, d1 = d1, d2 = d2, h0 = h0, h1 = h1, h2 = h2)
  times = seq(0,study.max,time.inc)
  out = ode(y = state, times = times, func = func.transitions, parms = parameters)
  return((out[which(times==study.max),] - out[which(times==study.min),])[c('D','H')])
}

# average the proportion expected to have acquired disease or become hospitalized
# by the end of a trial averaged over a specified age range between age.min and age.max
calc_disease_hosp_age_range = function(age.min,age.max,study.min,study.max,foi,efficacy,progression,w,d0,d1,d2,h0,h1,h2,age.inc=1/12){
  ages = seq(age.min,age.max+1,age.inc)
  results = matrix(0,length(ages),2)
  for(aa in 1:length(ages)){
    results[aa,] = calc_disease_hosp_by_age(ages[aa],study.min,study.max,foi,efficacy,progression,w,d0,d1,d2,h0,h1,h2)
  }
  return(colMeans(results))
}

# calculate the negative log likelihood of trial data
calc_NLL = function(foi.CYD14.tmp,foi.CYD15.tmp,foi.CYD57.tmp,efficacy.tmp,progression.tmp,waning.tmp){
  foi.tmp = c(foi.CYD14.tmp,foi.CYD15.tmp,foi.CYD57.tmp)[sapply(data$study,function(ss)which(unique(data$study)==ss))]
  efficacy.tmp = unlist(ifelse(data$group=='v',efficacy.tmp,0))
  progression.tmp = unlist(ifelse(data$group=='v',progression.tmp,0))
  attack.rate = rep(0,nrow(data))
  for(ii in 1:nrow(data)){
    attack.rate[ii] = calc_disease_hosp_age_range(
      age.min = data$age.min[ii],
      age.max = data$age.max[ii],
      study.min = data$study.min[ii],
      study.max = data$study.max[ii],
      foi = foi.tmp[ii],
      efficacy = efficacy.tmp[ii],
      progression = progression.tmp[ii],
      w = waning.tmp,
      d0,d1,d2,h0,h1,h2)[1]
  }
  return(-sum(dbinom(data$dengue,data$participants,attack.rate,log=T)))
}


# data from the CYD14 study
data = data.frame(
  study = c(
    'CYD14','CYD14','CYD14','CYD14','CYD14','CYD14',
    'CYD15','CYD15','CYD15','CYD15',
    'CYD57','CYD57','CYD57','CYD57','CYD57','CYD57','CYD57','CYD57'),
  age.min = c(
    2,6,12,2,6,12,
    9,12,9,12,
    4,4,6,6,4,4,6,6),
  age.max = c(
    5,11,14,5,11,14,
    11,16,11,16,
    5,5,11,11,5,5,11,11),
  study.min = c(
    0,0,0,0,0,0,
    0,0,0,0,
    2,2,3,3,2,2,3,3),
  study.max = c(
    1,1,1,1,1,1,
    1,1,1,1,
    3,3,4,4,3,3,4,4),
  dengue = c(
    15,10,2,1,8,4,
    10,6,9,6,
    5,5,17,11,1,3,10,14),
  participants = c(
    6778,1636,3598,813,1806,768,
    6029,7239,3005,3625,
    393,393,1738,1738,192,192,880,880),
  group = c(
    'v','v','v','c','c','c',
    'v','v','c','c',
    'v','v','v','v','c','c','c','c'))

# specify values for parameters that are fixed for model fitting
d0 = .4; d1 = .85; d2 = .1
h0 = .005 / d0; h1 = .1 / d1; h2 = .001 / d2


# set up details of parameter sweep, including number of different initial
# conditions to try and over what range of values for each parameter
num.ICs = 200
pars = sobol(vars = list(
  foi.CYD14.ic = c(-4,-1),
  foi.CYD15.ic = c(-4,-1),
  foi.CYD57.ic = c(-4,-1),
  efficacy.ic = c(0,1),
  progression.ic = c(0,1),
  waning.ic = c(3*7,3*360)),
  num.ICs)
pars$foi.CYD14.ic = 2 * 10 ^ pars$foi.CYD14
pars$foi.CYD15.ic = 2 * 10 ^ pars$foi.CYD14
pars$foi.CYD57.ic = 2 * 10 ^ pars$foi.CYD14
pars$waning.ic = 1 / pars$waning


# run the parameter sweep on multiple cores and save the output
cl = makeCluster(4)
registerDoParallel(cl)
sweep = foreach(ii = 1:num.ICs) %dopar% {
  library(deSolve)
  optim(
    par=c(
      log(pars$foi.CYD14.ic[ii]),
      log(pars$foi.CYD14.ic[ii]),
      log(pars$foi.CYD57.ic[ii]),
      logit(pars$efficacy.ic[ii]),
      logit(pars$progression.ic[ii]),
      log(pars$waning.ic[ii])),
    fn=function(par){calc_NLL(exp(par[1]),exp(par[2]),exp(par[3]),logit.inv(par[4]),logit.inv(par[5]),exp(par[6]))},control=list(maxit=50))
}
stopCluster(cl)
pars$foi.CYD14.mle = exp(sapply(sweep,function(ii)ii$par[1]))
pars$foi.CYD15.mle = exp(sapply(sweep,function(ii)ii$par[2]))
pars$foi.CYD57.mle = exp(sapply(sweep,function(ii)ii$par[3]))
pars$efficacy.mle = logit.inv(sapply(sweep,function(ii)ii$par[4]))
pars$progression.mle = logit.inv(sapply(sweep,function(ii)ii$par[5]))
pars$waning.mle = 1 / exp(sapply(sweep,function(ii)ii$par[6]))
pars$NLL = sapply(sweep,function(ii)ii$value)
save(list=ls(),file='trial_fit.RData')


# limit the parameters considered for further analyses to those
# that fit the data more or less equally well
pars = pars[which(pars$NLL<=min(pars$NLL)+5),]

