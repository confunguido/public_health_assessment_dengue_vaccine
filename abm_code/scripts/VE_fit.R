# 1 – AR_v / AR_p
# 1 – (Pr_dis|inf_v * Pr_inf|exp_v * Pr_exp_v) / (Pr_dis|inf_p * Pr_inf|exp_p * Pr_exp_p)
# 1 – (Pr_dis|inf_v * Pr_inf|exp_v) / (Pr_dis|inf_p * Pr_inf|exp_p)
# 1 – (Pr_dis|inf_v * Pr_inf|exp_v) / (Pr_dis|inf_p)
# 1 – (Pr_dis|inf_v / Pr_dis|inf_p) * Pr_inf|exp_v
# 1 - RR_dis|inf * RR_inf|exp

# vaccine efficacy at a given age
VE.age = function(a, b, c, age){
  1 - a / (1 + exp(b * (age - c)))
}

# calculate mean protection within a given age range
# assuming equal distribution within age range
VE = function(a, b, c, ages){
  mean(sapply(ages,function(aa) VE.age(a, b, c, aa)))
}

# specify age ranges from trial data
ages.young = seq(2, 9, .01)
ages.old = seq(9, 16, .01)

# VE estimates from trial data
VE.data = c(.701, .144, .819, .525)

# function to calculate goodness of fit between trial data and model
gof = function(a.pos, b.pos, c.pos, a.neg, b.neg, c.neg){
  (VE(a.pos,b.pos,c.pos,ages.young) - VE.data[1]) ^ 2 +
  (VE(a.neg,b.neg,c.neg,ages.young) - VE.data[2]) ^ 2 +
  (VE(a.pos,b.pos,c.pos,ages.old) - VE.data[3]) ^ 2 +
  (VE(a.neg,b.neg,c.neg,ages.old) - VE.data[4]) ^ 2
}

# find best-fit parameter values for age relationships
VE.fit = optim(par=c(1,0,9,1,0,9),fn=function(par){
  gof(par[1],par[2],par[3],par[4],par[5],par[6])})

# plot probability of protection against disease as a function of age
# for seronegative and seropositive individuals
ages.all = seq(0,40,.01)
plot(
  ages.all,
  pmax(0,VE.age(VE.fit$par[1], VE.fit$par[2], VE.fit$par[3], ages.all)),
  type='l', ylim=c(0,1), col = 2, lwd = 3, las=1,
  xlab = 'Age', ylab = 'Vaccine efficacy against virus-confirmed disease')
lines(
  ages.all,
  pmax(0,VE.age(VE.fit$par[4], VE.fit$par[5], VE.fit$par[6], ages.all)),
  col = 4, lwd = 3)
segments(c(2,2,9,9),VE.data,c(9,9,16,16),VE.data,col=c(2,4,2,4),lwd=2)
abline(v=c(2,9,16),lty=3)




# plot 3 different scenarios about the relative proportion of risk reduction
# attributable to reduction in risk or disease versus reduction in risk of infection
layout(matrix(1:3,1,3))

prop.inf = .1
RR_dis.pos = (1 - pmax(0,VE.age(VE.fit$par[1], VE.fit$par[2], VE.fit$par[3], ages.all))) ^ (1 - prop.inf)
RR_inf.pos = (1 - pmax(0,VE.age(VE.fit$par[1], VE.fit$par[2], VE.fit$par[3], ages.all))) ^ prop.inf
RR_dis.neg = (1 - pmax(0,VE.age(VE.fit$par[4], VE.fit$par[5], VE.fit$par[6], ages.all))) ^ (1 - prop.inf)
RR_inf.neg = (1 - pmax(0,VE.age(VE.fit$par[4], VE.fit$par[5], VE.fit$par[6], ages.all))) ^ prop.inf

RR_VCD_prot.fun(a, b, c, age)
ages.all = seq(0,40,.01)
plot(
  ages.all,
  RR_dis.pos,
  type='l', ylim=c(0,1), col = 2, lwd = 2, las=1,
  xlab = 'Age', ylab = 'Relative risk of infection or disease')
lines(
  ages.all,
  RR_dis.neg,
  col = 4, lwd = 2)
lines(
  ages.all,
  RR_inf.pos,
  col = 2, lwd = 2, lty = 3)
lines(
  ages.all,
  RR_inf.neg,
  col = 4, lwd = 2, lty = 3)


prop.inf = .49
RR_dis.pos = (1 - pmax(0,VE.age(VE.fit$par[1], VE.fit$par[2], VE.fit$par[3], ages.all))) ^ (1 - prop.inf)
RR_inf.pos = (1 - pmax(0,VE.age(VE.fit$par[1], VE.fit$par[2], VE.fit$par[3], ages.all))) ^ prop.inf
RR_dis.neg = (1 - pmax(0,VE.age(VE.fit$par[4], VE.fit$par[5], VE.fit$par[6], ages.all))) ^ (1 - prop.inf)
RR_inf.neg = (1 - pmax(0,VE.age(VE.fit$par[4], VE.fit$par[5], VE.fit$par[6], ages.all))) ^ prop.inf

RR_VCD_prot.fun(a, b, c, age)
ages.all = seq(0,40,.01)
plot(
  ages.all,
  RR_dis.pos,
  type='l', ylim=c(0,1), col = 2, lwd = 2, las=1,
  xlab = 'Age', ylab = 'Relative risk of infection or disease')
lines(
  ages.all,
  RR_dis.neg,
  col = 4, lwd = 2)
lines(
  ages.all,
  RR_inf.pos,
  col = 2, lwd = 2, lty = 3)
lines(
  ages.all,
  RR_inf.neg,
  col = 4, lwd = 2, lty = 3)


prop.inf = .9
RR_dis.pos = (1 - pmax(0,VE.age(VE.fit$par[1], VE.fit$par[2], VE.fit$par[3], ages.all))) ^ (1 - prop.inf)
RR_inf.pos = (1 - pmax(0,VE.age(VE.fit$par[1], VE.fit$par[2], VE.fit$par[3], ages.all))) ^ prop.inf
RR_dis.neg = (1 - pmax(0,VE.age(VE.fit$par[4], VE.fit$par[5], VE.fit$par[6], ages.all))) ^ (1 - prop.inf)
RR_inf.neg = (1 - pmax(0,VE.age(VE.fit$par[4], VE.fit$par[5], VE.fit$par[6], ages.all))) ^ prop.inf

RR_VCD_prot.fun(a, b, c, age)
ages.all = seq(0,40,.01)
plot(
  ages.all,
  RR_dis.pos,
  type='l', ylim=c(0,1), col = 2, lwd = 2, las=1,
  xlab = 'Age', ylab = 'Relative risk of infection or disease')
lines(
  ages.all,
  RR_dis.neg,
  col = 4, lwd = 2)
lines(
  ages.all,
  RR_inf.pos,
  col = 2, lwd = 2, lty = 3)
lines(
  ages.all,
  RR_inf.neg,
  col = 4, lwd = 2, lty = 3)

