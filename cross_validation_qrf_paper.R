######################################################################
### Quantile Regeression Forest and Regression Krigeage comparison###
#####################################################################

### Author : K. Vaysse (UMR LISAH - INRA Montpellier - SIG L-R)
### last update : 14/09/2016


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

##case of quantile regression forest##
#pH#
#data <- nadftph2_cal#choice of the table use to define the data using in the model
data<-nadftarg2_cal
data<-data[!is.na(data$MNT),]#clean some NA values
data<-data[!is.na(data$EMBERGER),]#clean some NA values

# in this cross validation example, the data represent the all the soil profile available to calibrate the
# model.  
transformation="yes"
case="clay"
fold = 100 #Number of time where the data will be randomly separate
size_validation = 0.75*nrow(data)# the pourcentage, the ratio of calibration profil will be use for validation
remove(result)

#Creating a progress bar to know the status of CV
progress.bar <- create_progress_bar("text")
progress.bar$init(fold)

for(iteration in 1:fold){
 
  train <- sample(nrow(data),size_validation) ## Selection of the training set
  mymodel <- quantregForest(x=data[train,c(58:59,65:83)],y=data[train,23],mtry=7,nodesize=10,ntree=1000)#run calibration of the model
  prediction <-predict(mymodel,data[-train,c(58:59,65:83)],quantile=c(0.005,0.025,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.975,0.995))#run validation of the model
  ##Preparation of the result table to calculate indicators
  if(iteration==1){
    result <- as.data.frame(cbind(prediction, data[-train, 23]))
  }else{
    new_result <- as.data.frame(cbind(prediction, data[-train, 23]))
    result<-rbind(result,new_result)
  }
  progress.bar$step()#add one step in the progress bar
}
names(result)[1:24]<-c("q005","q025","q05","q10","q15","q20","q25","q30","q35","q40","q45","Predicted","q55","q60","q65","q70","q75","q80","q85","q90","q95","q975","q995","Actual")

  ##Calculation of the indicators
 if(transformation=="yes"&case=="oc"){
   result$Predicted<-(exp(result$Predicted)-1)
    result$Actual<-(exp(result$Actual)-1)
 }
 if(transformation=="yes"&case=="clay"){
   result$Predicted<-(result$Predicted)*(result$Predicted)
    result$Actual<-(result$Actual)*(result$Actual)
    }

result$residu <-(result$Predicted) - (result$Actual)
r2<-1-(var(result$residu)/var(result$Actual))
me<-mean(result$residu)
mse<-mean(result$residu^2)
rmse<-sqrt(mean(result$residu^2))
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

write.table(y, "D://CD/Projet_DSM/Quantile_Regression_Forest/results/results_dataframe/qrfpi_result_clay.csv", row.names=FALSE)
write.table(result, "D://CD/Projet_DSM/Quantile_Regression_Forest/results/results_dataframe/qrf_result_clay.csv", row.names=FALSE)
  
  
