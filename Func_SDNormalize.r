#############################################################################
#File Name: Func_SDNormalize.r
#Author: Yanbo Xia
#Start Date: 2015.09.03
#Start Notes: 
#Description: Given a data set (already mean zero), normalize (by dividing 
#   by sample standard deviation).
#
#############################Edit Notes######################################
#Date: 
#Type: 
#Description: 

#############################################################################
SDNormalize<-function(dataset) {
  RecordLength<-ncol(dataset)
  Dim<-nrow(dataset)
  NYears<-RecordLength%/%T
  DaysRes<-RecordLength%%T
  TimeAxis<-1:T
  
  #Compute sample sd for each season (day of year)
  SDAll<-c()
  for (i in 1:Dim) {
    row<-c()
    Samples<-c()
    for (k in 1:NYears) {
      Samples<-cbind(Samples,dataset[i,((k-1)*T+1):(k*T)])
    }
    for (nu in 1:T) {
      if (nu<=DaysRes) {
        CrntSample<-c(Samples[nu,],dataset[i,NYears*T+nu])
      } else {
        CrntSample<-Samples[nu,]
      }
      row<-c(row,sd(CrntSample))
    }
    SDAll<-rbind(SDAll,row)
  }
  row.names(SDAll)<-city.names
  #SDAll<<-SDAll
  SDEst<-PlotNFit(SDAll,"RPlot_FitSD.pdf","Standard Deviation")
  
  #All seasons devide by estimated standard deviation
  dataNorm<-dataset
  for (t in 1:RecordLength) {
    nu<-(t-1)%%T+1
    for (i in 1:Dim) {
      dataNorm[i,t]<-dataNorm[i,t]/SDEst[i,nu]
    }
  }
  
  return(dataNorm)
}