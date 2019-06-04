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
save.fig = T
#=============================================================================#
# Functions---------------
#=============================================================================#
plot_color_bar = function(myColStrip_in, myColBreaks_in,myColBar_in, myTicks_in, plot_line = TRUE, label_in = "",...){
  image(t(myColStrip_in),
        col = myColBar_in, breaks = myColBreaks_in,
        xlim = c(0.0,1.0), axes = F)
  ylims = par("usr")[c(3,4)]
  myBarTicks = seq(from = ylims[1], to = ylims[2], by = diff(ylims)/(length(myTicks_in) - 1))
  myBarTickLabels = myTicks_in
  axis(4,labels = rep("",length(myBarTicks)), at = myBarTicks,cex.axis = 0.8,tck = -0.1)
  axis(4, at =  myBarTicks, labels = myBarTickLabels, cex.axis = 0.8,line = -1.0,lwd = 0)
  if(plot_line){
    abline(h = myBarTicks[2], col = "black", lwd = 2) 
  }
  mtext(label_in, side = 4, outer = F, lwd = 2,...)
}

#=============================================================================#
# Load and pre-process data from Test results---------------
#=============================================================================#
summary_list = readRDS('../data/output/20190218_output/summary_cohort_economics_test.RDS')

cc = 1

cohort_dataframe_seroneg_test = data.frame(
  Specificity = rep(0,length(summary_list)), Sensitivity = 0, SP9 = 0,
  DisAverted = 0, HospAverted = 0
)

cohort_dataframe_seropos_test = data.frame(
  Specificity = rep(0,length(summary_list)), Sensitivity = 0, SP9 = 0, 
  DisAverted = 0, HospAverted = 0
)

for(ff in 1:length(summary_list)){
  seroneg_table = summary_list[[ff]][[cc]]$Seronegative %>% filter(Year <= 50)
  seropos_table = summary_list[[ff]][[cc]]$Seropositive %>% filter(Year <= 50)
  
  cohort_dataframe_seropos_test[ff,c('Sensitivity', 'Specificity', 'SP9')] = 
    seroneg_table[1,c('Sensitivity', 'Specificity', 'SP9')]

  cohort_dataframe_seroneg_test[ff,c('Sensitivity', 'Specificity', 'SP9')] = 
    seroneg_table[1,c('Sensitivity', 'Specificity', 'SP9')]
     
  cohort_dataframe_seroneg_test$DisAverted[ff] = sum(seroneg_table$DisVax / seroneg_table$PopSizeVax, na.rm = T) / 
    sum(seroneg_table$Dis / seroneg_table$PopSize, na.rm = T)
  cohort_dataframe_seropos_test$DisAverted[ff] = sum(seropos_table$DisVax / seropos_table$PopSizeVax, na.rm = T) / 
    sum(seropos_table$Dis / seropos_table$PopSize, na.rm = T)

  cohort_dataframe_seroneg_test$HospAverted[ff] = sum(seroneg_table$HospVax / seroneg_table$PopSizeVax, na.rm = T) / 
    sum(seroneg_table$Hosp / seroneg_table$PopSize, na.rm = T)
  cohort_dataframe_seropos_test$HospAverted[ff] = sum(seropos_table$HospVax / seropos_table$PopSizeVax, na.rm = T) / 
    sum(seropos_table$Hosp / seropos_table$PopSize, na.rm = T)
  
  if(ff %% 100 == 0){cat("\r",ff)}
}

cohort_dataframe_seroneg_test = cohort_dataframe_seroneg_test[!is.infinite(cohort_dataframe_seroneg_test$HospAverted),]
cohort_dataframe_seropos_test = cohort_dataframe_seropos_test[!is.infinite(cohort_dataframe_seropos_test$HospAverted),]

cohort_dataframe_seroneg_test = cohort_dataframe_seroneg_test %>% drop_na()
cohort_dataframe_seropos_test = cohort_dataframe_seropos_test %>% drop_na()

#=============================================================================#
# Fit random forest to Average Cases by Cohort ---------------
#=============================================================================#
m = "randomforest"
if(m == "randomforest"){
  n.tree = 1000    
  model_hosp_seroneg_test = randomForest(
    HospAverted ~
      (Sensitivity + Specificity +  SP9 )^3,
    data = as.data.frame(cohort_dataframe_seroneg_test),
    importance = T,
    ntree = n.tree)
  
  model_hosp_seropos_test = randomForest(
    HospAverted ~
      (Sensitivity + Specificity +  SP9 )^3,
    data = as.data.frame(cohort_dataframe_seropos_test),
    importance = T,
    ntree = n.tree)
}else if (m == "gam"){
  cohort_dataframe_seroneg_test$HospAverted = cohort_dataframe_seroneg_test$HospAverted + 1e-5
  cohort_dataframe_seropos_test$HospAverted = cohort_dataframe_seropos_test$HospAverted + 1e-5  
  model_hosp_seroneg_test = gam(
    log(HospAverted) ~
      s(Specificity,Sensitivity,SP9, bs = "gp"),
    data = as.data.frame(cohort_dataframe_seroneg_test),
    family = "gaussian")
  
  model_hosp_seropos_test = gam(
    log(HospAverted) ~
      s(Specificity,Sensitivity,SP9, bs = "gp"),
    data = as.data.frame(cohort_dataframe_seropos_test),
    family = "gaussian")
}

#=============================================================================#
# Plot heatmaps of cases averted Test vs Only Vax-----------
#=============================================================================#
sp9_array = c(0.1,0.3,0.5,0.7,0.9)
sensitivity_array = seq(from = 0, by = 0.02, to  = 1.0)
specificity_array = seq(from = 0, by = 0.02, to  = 1.0)
lb = 0;ub = 2;n.breaks = 50
# myColBar = rainbow(n.breaks + 10)[-((n.breaks+1):(n.breaks+10))]
myColBar = c(rainbow(n.breaks/2,start=0,end=2/6),rainbow(n.breaks/2,start=3/6,end=5/6))
myColBreaks = seq(lb,ub,by = (ub - lb) / n.breaks)
myColBreaksTicks = seq(lb,ub,by = (ub - lb)/2)
myColStrip = as.matrix((myColBreaks + diff(myColBreaks)[1] / 2)[-length(myColBreaks)])

if(save.fig == T){
  jpeg(sprintf('../figures/supplement_figure_S3_per_capita_risk_serostatus_10y.jpeg'),
       width=6.35,height=2.5,units='in',res=400)
}
layout(
  cbind(matrix(1:10,2,5),c(11,11)), widths = c(rep(2,5),1), heights = rep(1,2)
)

par(mar = c(0.5,0.2,0.1,0.5), oma = c(2.5,2.5,1,2.5))
sensitivity.specificity.grid = expand.grid(Specificity = specificity_array,Sensitivity = sensitivity_array)
for(sp9_tmp in sp9_array){
  print(sp9_tmp)
  cases_averted_prediction_test = data.frame(
    'Specificity' = sensitivity.specificity.grid$Specificity, 'Sensitivity' = sensitivity.specificity.grid$Sensitivity,
    'SP9' = sp9_tmp)
  if(m == "randomforest"){  
    hosp.prediction.seroneg.grid =  predict(
      model_hosp_seroneg_test, cases_averted_prediction_test, predict.all = F)
    hosp.prediction.seropos.grid =  predict(
      model_hosp_seropos_test, cases_averted_prediction_test, predict.all = F)
  }else if (m == "gam"){
    hosp.prediction.seroneg.grid =  exp(predict(
      model_hosp_seroneg_test, cases_averted_prediction_test))
    hosp.prediction.seropos.grid =  exp(predict(
      model_hosp_seropos_test, cases_averted_prediction_test))
  }    
  hosp_averted_matrix_seroneg_test_novax = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
  hosp_averted_matrix_seropos_test_novax = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
  k = 1
  for(sen in 1:length(sensitivity_array)){
      for(spec in 1:length(specificity_array)){
          hosp_averted_matrix_seroneg_test_novax[sen,spec] = hosp.prediction.seroneg.grid[k]
          hosp_averted_matrix_seropos_test_novax[sen,spec] = hosp.prediction.seropos.grid[k]
          k = k + 1
      }
  }
  # correct bounds - seroneg
  threshold_matrix = hosp_averted_matrix_seroneg_test_novax
  threshold_matrix[threshold_matrix < lb] = lb
  threshold_matrix[threshold_matrix > ub] = ub

  image(
      specificity_array,sensitivity_array, t(threshold_matrix), 
      col = myColBar,
      breaks = seq(lb,ub,by = (ub - lb)/n.breaks),axes = F
  )
  contour(
      x=specificity_array,
      y=sensitivity_array,
      z=t(threshold_matrix),
      level=seq(lb,ub,by = (ub - lb)/n.breaks),
      lwd=c(rep(0.1,n.breaks/2),2,rep(0.1,n.breaks/2)),
      add=T,drawlabels=T
  )
  
  sp9str = sprintf("%.1f",sp9_tmp)
  mtext(text = bquote(PE[9] ~ " = " ~ .(sp9str)), side = 3, line = 0,cex = 0.7)
  if(sp9_tmp == 0.1){
    axis(2,labels = rep("",6), at = seq(from=0,to=1.0,by=0.2),cex.axis = 0.8,tck = -0.03)
    axis(2,at = seq(from=0,to=0.8,by=0.2),cex.axis = 0.7,line = -0.7,lwd = 0)
  }else{
    axis(2,labels = F,tick  = F)
  }
  box()
  if(sp9_tmp == 0.9){
    mtext(text = "Hospitalized\nseronegative", side = 4, line = 1, cex = 0.6)    
  }
  # correct bounds - seropos
  threshold_matrix = hosp_averted_matrix_seropos_test_novax
  threshold_matrix[threshold_matrix < lb] = lb
  threshold_matrix[threshold_matrix > ub] = ub
  image(
    specificity_array,sensitivity_array, t(threshold_matrix), 
    col = myColBar,
    breaks = seq(lb,ub,by = (ub - lb)/n.breaks),axes = F
  )
  contour(
    x=specificity_array,
    y=sensitivity_array,
    z=t(threshold_matrix),
    level=seq(lb,ub,by = (ub - lb)/n.breaks),
    lwd=c(rep(0.1,n.breaks/2),2,rep(0.1,n.breaks/2)),
    add=T,drawlabels=T
  )
  if(sp9_tmp == 0.1){
    axis(2,labels = rep("",6), at = seq(from=0,to=1.0,by=0.2),cex.axis = 0.8,tck = -0.03)
    axis(2,at = seq(from=0,to=0.8,by=0.2),cex.axis = 0.7,line = -0.7,lwd = 0)
  }else{
    axis(2,labels = F,tick  = F)
  }
  box()
  axis(1,labels = rep("",6), at = seq(from=0,to=1.0,by=0.2),cex.axis = 0.8,tck = -0.03)
  axis(1,at = seq(from=0,to=0.8,by=0.2),cex.axis = 0.7,line = -0.7,lwd = 0)
  if(sp9_tmp == 0.9){
    mtext(text = "Hospitalized\nseropositive", side = 4, line = 1, cex = 0.6)    
  }
}
mtext(text = "Sensitivity", side = 2, line = 1,cex = 0.8, outer = T)    
mtext(text = "Specificity", side = 1, line = 1, cex = 0.8, outer =T)    

# Plot color bar
par(mar = c(0.5,2,0.1,0.5))
plot_color_bar(myColStrip_in = myColStrip, 
               myColBreaks_in = myColBreaks, myColBar_in = myColBar, 
               myTicks_in = myColBreaksTicks, label_in = "Relative risk", 
               line = 1.5,cex = 0.8)
if(save.fig == T){dev.off()}

