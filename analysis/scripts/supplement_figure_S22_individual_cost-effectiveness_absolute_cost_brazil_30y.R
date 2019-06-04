#=============================================================================#
# Author: Guido Espana & Yutong Yao & Alex Perkins
#=============================================================================#
# user input ---------------
#=============================================================================#
rm(list = ls())
library(RColorBrewer)
library(tidyverse)
library(randomForest)
library(grDevices)
library(here)
library(mgcv)
setwd(here())
set.seed(123)
#=============================================================================#
# Economic parameters ---------------
#=============================================================================#
country = "Brazil"
econ = list()
econ$c.vax = 23
if(country == "Brazil"){
  econ$c.amb = 200 - 60
  econ$c.hosp = 500 - 200
  econ$c.death = 11000
  econ$gdp = 8649
}else if (country == "Philippines"){
  econ$c.amb = 40 - 20
  econ$c.hosp = 500 - 400
  econ$c.death = 3000
  econ$gdp = 2951
}else{
  print("Country not found")
}  

econ$c.test = 10
econ$d.amb = 0.545 * 4/365
econ$d.hosp = 0.545 * 14/365
econ$d.death = 1
econ$d.vax = 0

p.death = 0.005
discounting.rate = 0.03
discounting.daly.rate = 0.03
daly.threshold = 3 * econ$gdp

save.fig = T
#=============================================================================#
# Functions ---------------
#=============================================================================#
plot_color_bar_horizontal = function(myColStrip_in, myColBreaks_in,myColBar_in, myTicks_in, 
                                     myBarLims_in,
                                     label_in = "",...){
  image(myColStrip_in,
        col = myColBar_in, breaks = myColBreaks_in,
        xlim = c(0,1.0), ylim = c(0,1.0),axes = F, add = F)
  xlims = par("usr")[c(1,2)]
  slope_x = diff(xlims) / diff(myBarLims_in) 
  b_x = xlims[2] - slope_x * myBarLims_in[2]
  myBarTicks = myTicks_in * slope_x + b_x
  myBarTickLabels = myTicks_in
  axis(1,labels = rep("",length(myBarTicks)), at = myBarTicks,cex.axis = 0.8,tck = -0.4)
  axis(1, at =  myBarTicks, labels = myBarTickLabels, cex.axis = 0.8,line = -1.0,lwd = 0)
  mtext(label_in, side = 1, outer = F, lwd = 2,...)
}

calculate_individual_costs = function(dis_table.in, econ.in, p.death.in, 
                                         disc.rate.in, disc.daly.in, daly.threshold.in)
{
  n.years = nrow(dis_table.in)
  disc.rate = 1/((1+disc.rate.in)^(0:(n.years - 1)))
  disc.daly = 1/((1+disc.daly.in)^(0:(n.years - 1)))
  
  # Disease cases in the dataframe dis_table.in do not include hospitalizations...
  # no need to substract
  cost.treatment = sum(
    ( (dis_table.in$DisAverted) * econ.in$c.amb + 
        (dis_table.in$HospAverted) * econ.in$c.hosp + 
        (dis_table.in$HospAverted) * p.death.in * econ.in$c.death
    ) * disc.rate, na.rm = T
  )
  
  daly.treat = sum(
    ((dis_table.in$DisAverted ) * econ$d.amb + 
       (dis_table.in$HospAverted) * econ$d.hosp 
    ) * disc.daly, na.rm = T
  )
  
  daly.death = sum(
    ((dis_table.in$DaysLostDeathAverted) / 365) * econ.in$d.death * disc.daly, na.rm = T
  )
  
  daly.vax = sum(dis_table.in$Vaccinated * econ.in$d.vax * disc.daly, na.rm = T)
  
  net.daly = daly.treat + daly.death - daly.vax
  
  threshold.cost = (cost.treatment + net.daly * daly.threshold.in ) 
  return(list(threshold.cost = threshold.cost, net.treat.cost = cost.treatment , net.daly.factor = net.daly))
}

flip.matrix = function(x.in){
  return(t(x.in[nrow(x.in):1,]))
}
image.fxn = function(x.in, ...){
  # image(t(x.in[nrow(x.in):1,]), xaxt = 'n', yaxt = 'n', ...)
  image(t(x.in), xaxt = 'n', yaxt = 'n', ...)
}

#=============================================================================#
# Read data and set economic variables ---------------
#=============================================================================#
summary_list = readRDS('../data/output/20190218_output/summary_cohort_economics_test.RDS')

cc = 1

cohort_dataframe_test = data.frame(
  Specificity = rep(0,length(summary_list)), Sensitivity = 0, SP9 = 0, 
  PropVax = 0,PropPos = 0,PropNeg = 0,
  Dis = 0, DisVax = 0, Hosp = 0, HospVax = 0,
  net.daly.factor = 0, net.treat.cost = 0, threshold.cost = 0
)

for(ff in 1:length(summary_list)){
  seroneg_table = summary_list[[ff]][[cc]]$Seronegative %>% filter(Year <= 100)
  seropos_table = summary_list[[ff]][[cc]]$Seropositive %>% filter(Year <= 100)
  cohort_table = select(seroneg_table, c('Year'))
  cohort_dataframe_test[ff,c('Sensitivity', 'Specificity', 'SP9')] = 
    seroneg_table[1,c('Sensitivity', 'Specificity', 'SP9')]
  
  cohort_table$DisAverted = 
    (seroneg_table$Dis + seropos_table$Dis) / (seroneg_table$PopSize + seropos_table$PopSize) - 
    (seroneg_table$DisVax + seropos_table$DisVax) / (seroneg_table$PopSizeVax + seropos_table$PopSizeVax)
  cohort_table$Dis = 
    (seroneg_table$Dis + seropos_table$Dis) / (seroneg_table$PopSize + seropos_table$PopSize) 
  cohort_table$DisVax = 
    (seroneg_table$DisVax + seropos_table$DisVax) / (seroneg_table$PopSizeVax + seropos_table$PopSizeVax) 
  cohort_table$Hosp = 
    (seroneg_table$Hosp + seropos_table$Hosp) / (seroneg_table$PopSize + seropos_table$PopSize) 
  cohort_table$HospVax = 
    (seroneg_table$HospVax + seropos_table$HospVax) / (seroneg_table$PopSizeVax + seropos_table$PopSizeVax)
  
  cohort_table$HospAverted = 
    (seroneg_table$Hosp + seropos_table$Hosp) / (seroneg_table$PopSize + seropos_table$PopSize) - 
    (seroneg_table$HospVax + seropos_table$HospVax) / (seroneg_table$PopSizeVax + seropos_table$PopSizeVax)
  
  cohort_table$Vaccinated = (seroneg_table$PopSizeVax + seropos_table$PopSizeVax) / 
    (seroneg_table$PopSizeVax + seropos_table$PopSizeVax + seroneg_table$PopSize + seropos_table$PopSize)
  
  cohort_table$daysLostAverted = 
    (seroneg_table$discountedLifeExpectancy  +  seropos_table$discountedLifeExpectancy) /  (seroneg_table$PopSize + seropos_table$PopSize) - 
    (seroneg_table$discountedLifeExpectancyVax + seropos_table$discountedLifeExpectancyVax) / (seroneg_table$PopSizeVax + seropos_table$PopSizeVax)

    cohort_table$Tested = seroneg_table$Tested
  
  threshold_cost_cohort_test = calculate_individual_costs(
    cohort_table, econ.in = econ, p.death.in = p.death, 
    disc.rate.in = discounting.rate, disc.daly.in = discounting.daly.rate, 
    daly.threshold.in = daly.threshold)
  
  cohort_dataframe_test$net.treat.cost[ff] =  threshold_cost_cohort_test$net.treat.cost
  cohort_dataframe_test$net.daly.factor[ff] =  threshold_cost_cohort_test$net.daly.factor
  cohort_dataframe_test$threshold.cost[ff] = threshold_cost_cohort_test$threshold.cost
  
  cohort_dataframe_test$PropNeg[ff] = 
    (seroneg_table$PopSize + seroneg_table$PopSizeVax)[cc] /
    (seroneg_table$PopSizeVax + seropos_table$PopSizeVax + seroneg_table$PopSize + seropos_table$PopSize)[cc]
  cohort_dataframe_test$PropPos[ff] = 
    (seropos_table$PopSize + seropos_table$PopSizeVax)[cc] /
    (seroneg_table$PopSizeVax + seropos_table$PopSizeVax + seroneg_table$PopSize + seropos_table$PopSize)[cc]
  cohort_dataframe_test$PropVax[ff] = cohort_table$Vaccinated[cc] / cohort_table$Tested[cc]
  cohort_dataframe_test$Dis[ff] = sum(cohort_table$Dis, na.rm = T)
  cohort_dataframe_test$DisVax[ff] = sum(cohort_table$DisVax, na.rm = T)
  cohort_dataframe_test$Hosp[ff] = sum(cohort_table$Hosp, na.rm = T)
  cohort_dataframe_test$HospVax[ff] = sum(cohort_table$HospVax, na.rm = T)

  if(ff %% 100 == 0){cat("\r",ff)}
}

cohort_dataframe_test = cohort_dataframe_test %>% drop_na()
cohort_dataframe_test = cohort_dataframe_test %>% drop_na()
cohort_dataframe_test$PropVax[cohort_dataframe_test$PropVax > 1] = 1
cohort_dataframe_test$PropVax[cohort_dataframe_test$PropVax > 1] = 1

#=============================================================================#
# Fit Random Forest ---------------
#=============================================================================#
m = "gam"
if(m == "randomforest"){
  model_cost_treat_test = randomForest(net.treat.cost ~ (
      Specificity + Sensitivity + SP9 )^3,
      data = cohort_dataframe_test,ntree = 1000, importance = T)
  
  model_daly_test = randomForest(net.daly.factor ~ (
    Specificity + Sensitivity + SP9 )^3,
    data = cohort_dataframe_test,ntree = 1000, importance = T)
  
  model_threshold_test = randomForest(threshold.cost ~ (
    Specificity + Sensitivity + SP9 )^3,
    data = cohort_dataframe_test,ntree = 1000, importance = T)
  
  model_propvax_test = randomForest(PropVax ~ s(
    Specificity + Sensitivity + SP9 )^3,
    data = cohort_dataframe_test,ntree = 1000, importance = T)
}else if (m == "gam"){
  model_cost_treat_test = gam(net.treat.cost ~ s(
    Specificity,Sensitivity,SP9, bs = "gp"),
    data = cohort_dataframe_test,
    family = "gaussian")
  
  model_daly_test = gam(net.daly.factor ~ s(
    Specificity,Sensitivity,SP9, bs = "gp"),
    data = cohort_dataframe_test,
    family = "gaussian")
  
  model_threshold_test = gam(threshold.cost ~ s(
    Specificity,Sensitivity,SP9, bs = "gp"),
    data = cohort_dataframe_test,
    family = "gaussian")
  
  model_propvax_test = gam(PropVax ~s (
    Specificity,Sensitivity,SP9, bs = "gp"),
    data = cohort_dataframe_test,
    family = "gaussian")
}

#=============================================================================#
# Plot heatmaps of Threshold cost of the test when vaccine is $69 v2.0---------
#=============================================================================#
sensitivity_array = seq(from = 0, by = 0.02, to = 1.0)
specificity_array = seq(from = 0, by = 0.02, to = 1.0)
cost_effective_matrix_test_screening = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
n.breaks = 2
gdp.vec = c(1,3)
daly.thresholds = econ$gdp*gdp.vec

lb = -1;ub = 1

myColBar = c("#E0E0E0FF", "#00FF00FF")
myColBreaks = seq(lb,ub,by = (ub - lb)/n.breaks)
myColBreaksTicks = c(-1,0,1)
myColStrip = as.matrix((myColBreaks + diff(myColBreaks)[1] / 2)[-length(myColBreaks)])


if(save.fig == T){
  jpeg(sprintf('../figures/supplement_figure_S22_Individual_cost_test_%d_vax_%d_%s_30y.jpeg',econ$c.test,econ$c.vax,country),
       width=6,height=2.5,units='in',res=400)
}

layout(
  matrix(1:10,2,5,byrow = T),
  widths = rep(2,5),heights = c(2,2)
)
par(mar = c(0.5,0.2,0.1,0.5), oma = c(2.5,2.5,1,1))

sensitivity.specificity.grid = expand.grid(Specificity = specificity_array,Sensitivity = sensitivity_array)

for(dd in 1:length(daly.thresholds)){
  daly.threshold.tmp = daly.thresholds[dd]
  print(daly.thresholds[dd])
  for(SP9_tmp in c(0.1,0.3,0.5,0.7,0.9)){
    print(SP9_tmp)
    test_data = sensitivity.specificity.grid
    test_data$SP9 = SP9_tmp
    prop.vax = SP9_tmp * sensitivity.specificity.grid$Sensitivity + (1- sensitivity.specificity.grid$Specificity)*(1-SP9_tmp)
    prediction.grid =  predict(model_cost_treat_test,test_data,predict.all = F) + 
      predict(model_daly_test, test_data,predict.all = F) * daly.threshold.tmp - econ$c.vax * prop.vax - econ$c.test
    k = 1
    for(sen in 1:length(sensitivity_array)){
      for(sp in 1:length(specificity_array)){
        cost_effective_matrix_test_screening[sen,sp] =  prediction.grid[k]
        k = k + 1
      }
    }
    print(max(cost_effective_matrix_test_screening))
    #correct bounds
    cost_matrix = cost_effective_matrix_test_screening
    cost_matrix[cost_matrix < lb] = lb
    cost_matrix[cost_matrix > ub] = ub
    image(
      specificity_array,sensitivity_array, t(cost_matrix), 
      col = myColBar,
      breaks = seq(lb,ub,by = (ub - lb)/n.breaks),axes = F
    )
    breaks.array = c(-1,0,1)
    contour(
      x=specificity_array,
      y=sensitivity_array,
      z=t(cost_matrix),
      level=breaks.array,
      lwd= c(rep(0.0,length(which(breaks.array < 0))), 1, rep(0.0,length(which(breaks.array > 0)))),
      c(rep(0.1,n.breaks/2),2,rep(0.1,n.breaks/2)),
      add=T,drawlabels=F
    )
    if(SP9_tmp == 0.9){
      mtext(text = sprintf("%.0fx GDP",gdp.vec[dd]), side = 4, line = 0, cex = 0.6)    
    }
    if(dd == 1){
      sp9str = sprintf("%.1f",SP9_tmp)
      mtext(text = bquote(PE[9] ~ " = " ~ .(sp9str)), side = 3, line = 0,cex = 0.7)
    }
    if(SP9_tmp == 0.1){
      axis(2,labels = rep("",6), at = seq(from=0,to=1.0,by=0.2),cex.axis = 0.8,tck = -0.03)
      axis(2,at = seq(from=0,to=0.8,by=0.2),cex.axis = 0.7,line = -0.7,lwd = 0)
    }else{
      axis(2,labels = F,tick = F)
    }
    box()
    if(dd == 2){
      axis(1,labels = rep("",6), at = seq(from=0,to=1.0,by=0.2),cex.axis = 0.8,tck = -0.03)
      axis(1,at = seq(from=0,to=0.8,by=0.2),cex.axis = 0.7,line = -0.7,lwd = 0)
    }else{
      axis(1,labels = F,tick = F)
    }
  }
}
mtext(text = "Specificity", side = 1, line = 1, cex = 0.8,outer = T)  
mtext(text = "Sensitivity", side = 2, line = 1, cex = 0.8,outer = T)        

if(save.fig == T){dev.off()}
