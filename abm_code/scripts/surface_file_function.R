loc = read.csv('../simulator/data/Input/locations_20140807.csv')

if(!('mosquitoes' %in% names(loc))){
  if(!require(mgcv)){install.packages('mgcv');library(mgcv)}
  load('Surface')
  
  s = exp(predict(Surface,data.frame(
    Date=as.Date('01-01-2004','%m-%d-%Y'),Xcor=loc$xcor,Ycor=loc$ycor)))
  
  s = s / mean(s)
  
  write.table(
    cbind(loc,mosquitoes=s),
    file='../simulator/data/Input/locations_20150801.csv',
    quote=FALSE,sep=',',row.names=FALSE)
}
