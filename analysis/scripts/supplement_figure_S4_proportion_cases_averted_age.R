##=============================================================================#
## Author: Guido Espana & Yutong Yao & Alex Perkins
##=============================================================================#
## user input ---------------
##=============================================================================#
rm(list = ls())
library(RColorBrewer)
library(tidyverse)
library(randomForest)
library(grDevices)
library(here)
library(RColorBrewer)
library(mgcv)
setwd(here())

set.seed(123)
save.fig = T

##=============================================================================#
## Functions---------------
##=============================================================================#
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

calculate_ci_randomForest = function(model.in, data.in){
  pred_tmp = predict(model.in, data.in, predict.all = T)
  pred_ci = apply(pred_tmp$individual,1,sd,na.rm = T) * 1.96
  return(list(aggregate = pred_tmp$aggregate, CI = pred_ci))
}
##=============================================================================#
## Load and pre-process data ---------------
##=============================================================================#
summary_list_test = readRDS('../data/output/20190218_output/summary_routine_vaccination_life_exp_age_sweep.RDS')

test_novax_hosp_per = data.frame(
  Sensitivity = rep(0,length(summary_list_test)),
  Specificity = 0, SP9 = 0, Age = 0, Cases_averted = 0)
test_novax_symp_per = test_novax_hosp_per

for(ff in 1:length(summary_list_test)){
  tmp_summary_test = filter(summary_list_test[[ff]],Year <= 50)
  test_novax_symp_per$Specificity[ff] =tmp_summary_test$Specificity[1]
  test_novax_symp_per$Sensitivity[ff] = tmp_summary_test$Sensitivity[1]
  test_novax_symp_per$Age[ff] = tmp_summary_test$AgeVax[1]
  test_novax_symp_per$SP9[ff] =tmp_summary_test$SP9Prevax[1]
  test_novax_symp_per$Cases_averted[ff] = 
    sum(tmp_summary_test$DisAverted) / sum(tmp_summary_test$Dis)
  test_novax_hosp_per$Cases_averted[ff] = 
    sum(tmp_summary_test$HospAverted) / sum(tmp_summary_test$Hosp)
}
test_novax_hosp_per[,1:4] = test_novax_symp_per[,1:4]


##=============================================================================#
## Fit random forest to cases averted ---------------
##=============================================================================#
m = "gam"
if(m == "randomforest"){
    n.trees = 1000
    model_cases_novax = randomForest(
        Cases_averted ~
          (Sensitivity + Specificity + SP9 + Age) ^ 4,
      data = test_novax_symp_per,
      importance = T,
      ntree = n.trees)
  
  model_hosp_novax = randomForest(
      Cases_averted ~
          (Sensitivity + Specificity + SP9 + Age) ^ 4,
      data = test_novax_hosp_per,
      importance = T,
  ntree = n.trees)
}else if (m == "gam"){
    model_cases_novax = gam(
        Cases_averted ~
            s(Sensitivity,Specificity,SP9,Age),
        data = test_novax_symp_per,
        family = "gaussian")
  
  
    model_hosp_novax = gam(
        Cases_averted ~
            s(Sensitivity,Specificity,SP9,Age),
        data = test_novax_hosp_per,
        family = "gaussian")
}


##=============================================================================#
## Plot heatmaps of hospitalization cases averted Test vs No Vax-----------
##=============================================================================#
sp9_array = c(0.1,0.3,0.5,0.7,0.9)
sensitivity_array = seq(from = 0, by = 0.01, to  = 1.0)
specificity_array = seq(from = 0, by = 0.01, to  = 1.0)
lb = -0.15;ub = 0.15;n.breaks = 50

myColBar = c(rainbow(n.breaks/2,start=0,end=2/6),rainbow(n.breaks/2,start=3/6,end=5/6))
myColBreaks = seq(lb,ub,by = (ub - lb)/n.breaks)
myColBreaksTicks = seq(lb,ub,by = (ub - lb)/2)
myColStrip = as.matrix((myColBreaks + diff(myColBreaks)[1] / 2)[-length(myColBreaks)])

if(save.fig == T){
  jpeg(sprintf('../figures/supplement_figure_S4_cases_averted_heatmap_age.jpeg'),
       width=6.35,height=10,units='in',res=400)
}
layout(
  cbind(matrix(1:40,8,5,byrow = T),rep(41,8)), widths = c(rep(2,5),1), heights = rep(1,8)
)

par(mar = c(0.5,0.2,0.1,0.5), oma = c(2.5,2.5,1,2.5))
for(age in 9:16){
  sensitivity.specificity.grid = expand.grid(Specificity = specificity_array,Sensitivity = sensitivity_array)
  print(age)
  for(sp9_tmp in sp9_array){
    cases_averted_matrix_test_novax = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
    hosp_averted_matrix_test_novax = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
    cases_averted_prediction_test_novax = data.frame(
      'Specificity' = sensitivity.specificity.grid$Specificity, 'Sensitivity' = sensitivity.specificity.grid$Sensitivity,
      'SP9' = sp9_tmp, 'Age' = age)
    cases.prediction.grid =  predict(
      model_cases_novax, cases_averted_prediction_test_novax, predict.all = F)
    hosp.prediction.grid =  predict(
      model_hosp_novax, cases_averted_prediction_test_novax, predict.all = F)
    k = 1
    for(sen in 1:length(sensitivity_array)){
      for(spec in 1:length(specificity_array)){
        cases_averted_matrix_test_novax[sen,spec] = cases.prediction.grid[k]
        hosp_averted_matrix_test_novax[sen,spec] = hosp.prediction.grid[k]
        k = k + 1
      }
    }
    # correct bounds
    threshold_matrix = hosp_averted_matrix_test_novax
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
      add=T,drawlabels=F
    )
    
    if(sp9_tmp == 0.1){
      axis(2,labels = rep("",6), at = seq(from=0,to=1.0,by=0.2),cex.axis = 0.8,tck = -0.03)
      axis(2,at = seq(from=0,to=0.8,by=0.2),cex.axis = 0.7,line = -0.7,lwd = 0)
    }else{
      axis(2,labels = F,tick = F)
    }
    axis(1,labels = F,tick = F)
    box()
    if(sp9_tmp == 0.9){
      mtext(text = sprintf("Age:%.0f",age), side = 4, line = 0, cex = 0.6)    
    }
    sp9str = sprintf("%.1f",sp9_tmp)
    if(age == 9){
      mtext(text = bquote(PE[9] ~ " = " ~ .(sp9str)), side = 3, line = 0,cex = 0.7)      
    }
  }
}

mtext(text = "Sensitivity", side = 2, line = 1.2,cex = 0.8, outer = T)    
mtext(text = "Specificity", side = 1, line = 1.2, cex = 0.8, outer = T)    

# Plot color bar
par(mar = c(0.5,2,0.1,0.5))
plot_color_bar(myColStrip_in = myColStrip, 
               myColBreaks_in = myColBreaks, myColBar_in = myColBar, 
               myTicks_in = myColBreaksTicks, label_in = "Proportion of cases averted", 
               cex = 0.6, line = 1.5)

if(save.fig == T){dev.off()}

# ##=============================================================================#
# ## Plot optimal age of vaccination -----------
# ##=============================================================================#
# sensitivity_array = c(0.9,0.9,0.9,0.7)
# specificity_array = c(0.9,0.95,0.7,0.9)
# sp9_array = c(0.1,0.3,0.5,0.7,0.9)
# age_array = 9:16
# optimal_age_matrix = matrix(0,nrow = length(sp9_array), ncol = length(sensitivity_array))
# 
# for(ss in 1:length(sensitivity_array)){
#   sp9.age.grid = expand.grid(SP9 = sp9_array,Age = age_array)
#   hosp_averted_table = data.frame(
#     'Specificity' = specificity_array[ss], 'Sensitivity' = sensitivity_array[ss],
#     'SP9' = sp9.age.grid$SP9, 'Age' = sp9.age.grid$Age)
#   
#   hosp_averted_table$PredAverted =  predict(
#     model_hosp_novax, hosp_averted_table, predict.all = F)
#   
#   for(sp in 1:length(sp9_array)){
#     sp9_tmp = sp9_array[sp]
#     tmp_table = filter(hosp_averted_table, SP9 == sp9_tmp)
#     optimal_age_matrix[sp,ss] = tmp_table[which.max(tmp_table$PredAverted),'Age']
#   }
# }
# 
# col_bar = RColorBrewer::brewer.pal(4,"Dark2")
# col_bar = rep("black",4)
# 
# point_types = c(8,9,10,11)
# point_types = rep(8,4)
# 
# if(save.fig == T){
#   jpeg(sprintf('../figures/experimental_test_1b_optimal_age.jpeg'),
#        width=6,height=5.5,units='in',res=400)
# }
# layout(
#   matrix(1:4,2,2,byrow = T)
# )
# 
# par(mar = c(1,1,1,1), oma = c(3,3,3,3))
# plot(sp9_array,optimal_age_matrix[,1], type = "l", col = col_bar[1], lwd = 2, axes = F,
#      ylim = c(9,16), ylab = "Optimal Age", xlab = "PE9")
# grid(0,NULL)
# points(sp9_array,optimal_age_matrix[,1], pch = point_types[1], col = col_bar[1])
# box();axis(2)
# 
# mtext(side = 3, line = 0, sprintf("sens:%.2f,spec:%.2f",sensitivity_array[1], specificity_array[1]))
# 
# plot(sp9_array,optimal_age_matrix[,2], type = "l",  col = col_bar[2], lwd = 2, axes = F,
#      ylim = c(9,16))
# points(sp9_array,optimal_age_matrix[,2], pch = point_types[2], col = col_bar[2])
# box()
# grid(0,NULL)
# mtext(side = 3, line = 0, sprintf("sens:%.2f,spec:%.2f",sensitivity_array[2], specificity_array[2]))
# 
# plot(sp9_array,optimal_age_matrix[,3], type = "l", col = col_bar[3], lwd = 2, ylim = c(9,16))
# points(sp9_array,optimal_age_matrix[,3], pch = point_types[3],  col = col_bar[3])
# grid(0,NULL)
# mtext(side = 3, line = 0, sprintf("sens:%.2f,spec:%.2f",sensitivity_array[3], specificity_array[3]))
# 
# plot(sp9_array,optimal_age_matrix[,4], type = "l", col = col_bar[4], lwd = 2, axes = F, 
#      ylim = c(9,16))
# points(sp9_array,optimal_age_matrix[,4], pch = point_types[4], col = col_bar[4])
# box();axis(1)
# grid(0,NULL)
# mtext(side = 3, line = 0, sprintf("sens:%.2f,spec:%.2f",sensitivity_array[4], specificity_array[4]))
# 
# mtext(side = 1, line = 2, text = "PE9", outer = T)
# mtext(side = 2, line = 2, text = "Age", outer = T)
# if(save.fig == T){dev.off()}


