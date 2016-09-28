######################################################################
### Quantile Regeression Forest and Regression Krigeage comparison###
#####################################################################

### Author : K. Vaysse (UMR LISAH - INRA Montpellier - SIG L-R)
### last update : 27/09/2016


##### library loading #####
library(raster)
library(randomForest)
library(gstat)
library(rgdal)
library(maptools)
library(plyr)
library(ggplot2)
library(grid)
library(shapefiles)
library(maps)
library(maptools)
library(spatstat)
library(sp)
library(lattice)
library(latticeExtra)
library(colorspace)
library(fields)
library(Hmisc)
library(GISTools)
library(quantregForest)
library(plyr)

###########################

#### Soil data and model loading ###
load()#add the location of the .RDATA files with the soil data for pH, clay and organic carbon and associated model for Quantile Regression Forest and QRF.

####################################

#### Cross-validation for models ###

##case of regression krigeage##
#pH#
data <- nadftph2_cal#choice of the table use to define the data using in the model
data<-data[!is.na(data$MNT),]#clean some NA values
data<-data[!is.na(data$EMBERGER),]#clean some NA values

transformation="no"
case="ph"
fold = 100 #Number of time where the data will be randomly separate
size_validation = 0.75*nrow(data)# the pourcentage, the ratio of calibration profil will be use for validation
remove(rf_result)
remove(model_result)
remove(result)
remove(r2)

#Creating a progress bar to know the status of CV
progress.bar <- create_progress_bar("text")
progress.bar$init(fold)

for(iteration in 1:fold){
 
  train <- sample(nrow(data),size_validation) ## Selection of the training set
  ## Model calibration
  mymodel <- randomForest(ph2~MNT+SLOPE+PLANCURV+PROCURV+TPI+TWI+MRVBF25+MRVBF50+MRVBF100+MRRTF25+MRRTF50+MRRTF100+MINERALOGY+HARDNESS+TEXTURE+EMBERGER+MARTONNE+TMAX+TMIN+PRECIP+CLASS,data=data[train,],mtry=7,nodesize=10, ntree=1000)
model_prediction<-predict(mymodel,data[train,])#run validation of the model
model_result <- as.data.frame(cbind(prediction, data[train, c(11,85:86)]))
model_result$residu<-model_result$prediction-model_result$ph2
coordinates(model_result)=~x_rgf93.1+y_rgf93.1
rk.ev = variogram(residu~1, model_result, cutoff=60000, width=5000)
plot(rk.ev, xlab="distance", ylab="semivariance",main="experimental variogram of residual pH [5-15 cm]")
rk.vgm = fit.variogram(rk.ev, vgm(nugget = 0.5,psill=0.2, "Sph",range = 40000, cutoff=50000, width=5000))
plot(rk.ev, rk.vgm, xlab="distance (m)", ylab="semivariance",main="residual variogram of pH [5-15 cm]",pch=19,col="black")

## model validation
  
  prediction <-predict(mymodel,data[-train,])#run validation of the model
  ##Preparation of the randomForest Result for kriging
rf_result <- as.data.frame(cbind(prediction, data[-train, c(11,85:86)]))
names(rf_result)[2]<-"Actual"
rf_result$residu<-rf_result$prediction-rf_result$Actual
##Some kriging ....
coordinates(rf_result)=~x_rgf93.1+y_rgf93.1
rk.cv <- krige.cv(residu~1, rf_result, rk.vgm,verbose=F)
quantile.cv <- as.data.frame(setNames(replicate(22,numeric(nrow(rk.cv)), simplify = F),c("q005","q025","q05","q10","q15","q20","q25","q30","q35","q40","q45","q55","q60","q65","q70","q75","q80","q85","q90","q95","q975","q995") ))
rk.cv<-cbind(rk.cv,quantile.cv)
pi_coef<-as.vector(c(-2.58,-1.96,-1.645,-1.282,-1.04,-0.83,-0.675,-0.51,-0.38,-0.29,-0.12,0.12,0.29,0.38,0.51,0.675,0.83,1.04,1.282,1.645,1.96,2.58))# coefficient a appliquer sur le standard deviation de krigeage pour obtenir les intervalles de prédictions respectivement à 10,20,30,40,50 , 60 , 70 , 80 ,90, 95 ,99 %
for(iterationbis in 1 : length(pi_coef)){
  rk.cv[,iterationbis+8]=rk.cv$var1.pred+sqrt(rk.cv$var1.var)*pi_coef[iterationbis]
}
rk_result <- as.data.frame(setNames(replicate(23,numeric(nrow(rk.cv)), simplify = F),c("Prediction","q005","q025","q05","q10","q15","q20","q25","q30","q35","q40","q45","q55","q60","q65","q70","q75","q80","q85","q90","q95","q975","q995") ))
rf_result@data<-cbind(rf_result@data,rk_result)
rf_result@data$Prediction<-(rf_result@data$prediction-rk.cv$var1.pred)
for(iterationbis in 1:length(pi_coef)){
  rf_result@data[,iterationbis+4]=rf_result@data$prediction+rk.cv[,(8+iterationbis)]
}
rf_result@data$residub<-(rf_result@data$Prediction-rf_result@data$Actual)
if(iteration==1){
    result <- as.data.frame(rf_result@data)
  }else{
    result<-rbind(result,rf_result@data)
  }
  progress.bar$step()#add one step in the progress bar
}
r2<-1-(var(result$residub)/var(result$Actual))
me<-mean(result$residub)
mse<-mean(result$residub^2)
rmse<-sqrt(mean(result$residub^2))
r2
me
mse
rmse

##Accuracy plot##
####Accuracy plot####
result$cond99<-NA
result$cond95<-NA
result$cond90<-NA
result$cond80<-NA
result$cond70<-NA
result$cond60<-NA
result$cond50<-NA
result$cond40<-NA
result$cond30<-NA
result$cond20<-NA
result$cond10<-NA
  ##Calculation of the indicators
 if(transformation=="yes"&case=="oc"){
    result$Actual<-(log(result$Actual+1))
 }
 if(transformation=="yes"&case=="clay"){
    result$Actual<-(sqrt(result$Actual))
 }
total<-length(result$Actual)
for(i in 1:total){
    if(result$Actual[i]<=result$q995[i]&result$Actual[i]>=result$q005[i])
    {result$cond99[i]<-1}else{result$cond99[i]<-0}
    if(result$Actual[i]<=result$q975[i]&result$Actual[i]>=result$q025[i])
    {result$cond95[i]<-1}else{result$cond95[i]<-0}
    if(result$Actual[i]<=result$q95[i]&result$Actual[i]>=result$q05[i])
    {result$cond90[i]<-1}else{result$cond90[i]<-0}
    if(result$Actual[i]<=result$q90[i]&result$Actual[i]>=result$q10[i])
    {result$cond80[i]<-1}else{result$cond80[i]<-0}
    if(result$Actual[i]<=result$q85[i]&result$Actual[i]>=result$q15[i])
    {result$cond70[i]<-1}else{result$cond70[i]<-0}
    if(result$Actual[i]<=result$q80[i]&result$Actual[i]>=result$q20[i])
    {result$cond60[i]<-1}else{result$cond60[i]<-0}
    if(result$Actual[i]<=result$q75[i]&result$Actual[i]>=result$q25[i])
    {result$cond50[i]<-1}else{result$cond50[i]<-0}
  if(result$Actual[i]<=result$q70[i]&result$Actual[i]>=result$q30[i])
  {result$cond40[i]<-1}else{result$cond40[i]<-0}
  if(result$Actual[i]<=result$q65[i]&result$Actual[i]>=result$q35[i])
  {result$cond30[i]<-1}else{result$cond30[i]<-0}
  if(result$Actual[i]<=result$q60[i]&result$Actual[i]>=result$q40[i])
  {result$cond20[i]<-1}else{result$cond20[i]<-0}
    if(result$Actual[i]<=result$q55[i]&result$Actual[i]>=result$q45[i])
    {result$cond10[i]<-1}else{result$cond10[i]<-0}
}
  number99<-sum(result$cond99)
  number95<-sum(result$cond95)
  number90<-sum(result$cond90)
  number80<-sum(result$cond80)
  number70<-sum(result$cond70)
  number60<-sum(result$cond60)
  number50<-sum(result$cond50)
  number40<-sum(result$cond40)
  number30<-sum(result$cond30)
  number20<-sum(result$cond20)
  number10<-sum(result$cond10)
  pourcent99<-number99/total*100
  pourcent95<-number95/total*100
  pourcent90<-number90/total*100
  pourcent80<-number80/total*100
  pourcent70<-number70/total*100
  pourcent60<-number60/total*100
  pourcent50<-number50/total*100
  pourcent40<-number40/total*100
  pourcent30<-number30/total*100
  pourcent20<-number20/total*100
  pourcent10<-number10/total*100
  y<-as.data.frame(cbind(c(99,95,90,80,70,60,50,40,30,20,10),c(pourcent99,pourcent95,pourcent90,pourcent80,pourcent70,pourcent60,pourcent50,pourcent40,pourcent30,pourcent20,pourcent10)))
  names(y)<-c("PI","Pourcentage")
  y<-y/100
  ap<-qplot(PI,Pourcentage,data=y)
  ap + geom_abline(intercept = 0,colour="red",size=2)+geom_point(colour = "black", size = 3.5) +labs(title="a")+theme_bw()+scale_y_continuous(limits = c(0, 0.99))+scale_x_continuous(limits = c(0, 0.99))+ylab(label = "Proportion within interval")

write.table(y, "D://CD/Projet_DSM/Quantile_Regression_Forest/results/results_dataframe/rkpi_result_ph.csv", row.names=FALSE)
write.table(result, "D://CD/Projet_DSM/Quantile_Regression_Forest/results/results_dataframe/rk_result_ph.csv", row.names=FALSE)