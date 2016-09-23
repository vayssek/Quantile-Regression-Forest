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
data <- nadftph2_cal#choice of the table use to define the data using in the model
data<-data[!is.na(data$MNT),]#clean some NA values
data<-data[!is.na(data$EMBERGER),]#clean some NA values

# in this cross validation example, the data represent the all the soil profile available to calibrate the
# model.  

fold = 100 #Number of time where the data will be randomly separate
size_validation = 0.75*nrow(data)# the pourcentage, the ratio of calibration profil will be use for validation
vector_r2 <- data.frame()
vector_me <- data.frame()
vector_mse <- data.frame()
vector_rmse <- data.frame()


#Creating a progress bar to know the status of CV
progress.bar <- create_progress_bar("text")
progress.bar$init(fold)

for(iteration in 1:fold){
 
  train <- sample(nrow(data),size_validation) ## SELECTION DE 25% POUR LE JEU DE VALIDATION EXTERNE
  mymodel <- quantregForest(x=data[train,c(58:59,65:83)],y=data[train,11],mtry=7,nodesize=10,ntree=1000)
  prediction <-predict(mymodel,data[-train,c(58:59,65:83)],quantile=c(0.50))
  
  result <- as.data.frame(cbind(prediction, data[-train, 11]))
  names(result) <- c("Predicted", "Actual")
  result$residu <-(result$Predicted) - (result$Actual)
  r2<-1-(var(result$residu)/var(result$Actual))
  me<-mean(result$residu)
  mse<-mean(result$residu^2)
  rmse<-sqrt(mean(result$residu^2))
  vector_r2<-rbind(vector_r2,r2)
  vector_me<-rbind(vector_me,me)
  vector_mse<-rbind(vector_mse,mse)
  vector_rmse<-rbind(vector_rmse,rmse)
  # append this iteration's test set to the test set copy data frame
  # keep only the Sepal Length Column
  
  progress.bar$step()
}

mean(vector_r2[,1])
mean(vector_me[,1])
mean(vector_mse[,1])
mean(vector_rmse[,1])

  
