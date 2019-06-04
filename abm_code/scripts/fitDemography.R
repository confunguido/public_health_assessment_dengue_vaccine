# read in age distribution data and distribute over age in days
age = read.csv('../simulator/data/Input/census_peru.csv')
yrs = age$Total
age = age$Total.1 / sum(age$Total.1)
age.day = rep(age/365,each=365)
days = 1:length(age.day)

# function to calculate the goodness of fit of a given mortality model
mort.gof = function(mu){
  age.sim = 1
  age.mort = days * mu
  for(ii in 1:(length(days)-1)){
    age.sim = c(age.sim,age.sim[ii]*(1-age.mort[ii]))
  }
  age.sim = age.sim / sum(age.sim)
  return(sum((age.sim-age.day)^2))
}

# find the mortality model parameter that maximizes goodness of fit
mu.opt = optimize(mort.gof,interval=c(5e-9,1e-8))$minimum

# solve for age distribution to get necessary birth rate
age.sim = 1
age.mort = days * mu.opt
for(ii in 1:(length(days)-1)){
  age.sim = c(age.sim,age.sim[ii]*(1-age.mort[ii]))
}
age.sim = age.sim / sum(age.sim)
birth = age.sim[1]

# write birth and death rates to file
df = data.frame(birth = birth, death = mu.opt)
write.csv(df,'../simulator/data/Input/demo_rates.csv',row.names=FALSE)


# age.sim = 1
# age.mort = days * mu.opt
# for(ii in 1:(length(days)-1)){
#   age.sim = c(age.sim,age.sim[ii]*(1-age.mort[ii]))
# }
# age.sim = age.sim / sum(age.sim)
# 
# plot(age.day,type='l')
# lines(age.sim,col=2)
# 
# 
# age.model = age.day
# for(ii in 1:30000){
#   age.model[2:length(age.model)] = age.model[-length(age.model)] * (1 - age.mort[-1])
#   age.model[1] = 6.65032e-05
# }
