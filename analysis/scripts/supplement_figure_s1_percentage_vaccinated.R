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
  axis(4, at =  myBarTicks, labels = myBarTickLabels, cex.axis = 0.7,line = -1.0,lwd = 0)
  if(plot_line){
    abline(h = myBarTicks[2], col = "black", lwd = 2) 
  }
  mtext(label_in, side = 4, outer = F, lwd = 2,...)
}

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
  axis(1,labels = rep("",length(myBarTicks)), at = myBarTicks,cex.axis = 0.8,tck = -0.1)
  axis(1, at =  myBarTicks, labels = myBarTickLabels, cex.axis = 0.8,line = -1.0,lwd = 0)
  # axis(1, at = myBarTicks, labels = myBarTickLabels,las = 1)
  mtext(label_in, side = 1,  outer = F, lwd = 2,...)
}
#=============================================================================#
# Plot heatmaps of cases averted Test vs No Vax-----------
#=============================================================================#
sp9_array = c(0.1,0.3,0.5,0.7,0.9)
sensitivity_array = seq(from = 0, by = 0.05, to  = 1.0)
specificity_array = seq(from = 0, by = 0.05, to  = 1.0)
lb = 0.0;ub = 1.0;n.breaks = 20

myColBar = rainbow(n.breaks + 10)[-((n.breaks+1):(n.breaks+10))]
myColBreaks = seq(lb,ub,by = (ub - lb)/n.breaks)
myColBreaksTicks = seq(lb,ub,by = (ub - lb)/2)
myColStrip = as.matrix((myColBreaks + diff(myColBreaks)[1] / 2)[-length(myColBreaks)])
fticks = seq(from=0,to=1.0,by=0.2)

if(save.fig == T){
  jpeg(sprintf('../figures/supplement_results_exp_1_fig_s1_percentage_vaccinated.jpeg'),
       width=6.2,height=1.5,units='in',res=400)
}
layout(
  matrix(1:6,1,6),widths = c(rep(4,5),1),heights = c(2,2,2,2,2,3)
)

par(mar = c(0.5,0.2,0.1,0.5), oma = c(2.5,2.5,1,2.5))
sensitivity.specificity.grid = expand.grid(Specificity = specificity_array,Sensitivity = sensitivity_array)
for(sp9_tmp in sp9_array){
  print(sp9_tmp)
  proportion_vaccinated_matrix = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
  proportion_vaccinated_pos_matrix = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
  proportion_vaccinated_neg_matrix = matrix(0,nrow = length(sensitivity_array), ncol = length(specificity_array))
  k = 1
  for(sen in 1:length(sensitivity_array)){
    for(spec in 1:length(specificity_array)){
      proportion_vaccinated_matrix[sen,spec] = sp9_tmp * sensitivity_array[sen] + (1-specificity_array[spec])*(1-sp9_tmp)
      k = k + 1
    }
  }
  # correct bounds
  threshold_matrix = proportion_vaccinated_matrix
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
    lwd=c(rep(0.1,n.breaks/2),0.1,rep(0.1,n.breaks/2)),
    add=T,drawlabels=T,labcex = 0.6
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
  sp9str = sprintf("%.1f",sp9_tmp)
  mtext(text = bquote(PE[9] ~ " = " ~ .(sp9str)), side = 3, line = 0,cex = 0.7)
}
mtext(text = "Sensitivity", side = 2, line = 1,cex = 0.8, outer = T)    
mtext(text = "Specificity", side = 1, line = 1, cex = 0.8, outer = T)    

# Plot color bar
plot_color_bar(myColStrip_in = myColStrip, 
               myColBreaks_in = myColBreaks, myColBar_in = myColBar, plot_line = FALSE,
               myTicks_in = myColBreaksTicks, label_in = "Proportion vaccinated", 
               cex = 0.6,line = 0.8)
if(save.fig == T){dev.off()}

