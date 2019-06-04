if(!require(mgcv)){install.packages('mgcv');library(mgcv)}
load('Surface')

# Plot time and surface if I had used GAMs instaed of what I use (gives relatively fine output)

par(mfrow=c(1,2))
plot(Surface)

# Predict what the mean number of mosquitoes are at the point (690000,9582000) on day 10619 (which corresponds to 1999-01-28 (the first observation) (as.numeric(as.Date("1999-01-28"))))

s = exp(predict(Surface,data.frame(Date=10619,Xcor=loc[,1],Ycor=loc[,2])))

loc = read.csv('../simulator/data/Input/locations_20140807.csv')
loc = cbind(loc$xcor,loc$ycor)
