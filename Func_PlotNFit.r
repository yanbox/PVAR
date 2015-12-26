#############################################################################
#File Name: Func_PlotNFit.r
#Author: Yanbo Xia
#Start Date: 2015.09.03
#Description: Defines the PlotNFit function.
#
#############################Edit Notes######################################
#Date: 
#Type: 
#Description: 

#############################################################################
#Spline degree of freedom
#SplineDF<-13
#For a given matrix d*T, denoting daily measurements for every station,
#Plot graphs, fit and plot, return fitted set.
PlotNFit<-function(dataset,file.name,y.var) {
  Dim<-nrow(dataset)
  color.set<-c("Black","Red","Blue","Green")
  TimeAxis<-1:T
  
  pdf(file=file.name,width=16,height=10)
  #Seperate plots for each city
  for (i in 1:Dim) {
    plot(TimeAxis,dataset[i,],type="l",
         xlab="Day of year", ylab=paste(y.var,"(*0.1 Degrees Centigrade)"),
         main=paste("Daily ",y.var,": ",city.names[i]))
  }
  #Combined plot
  plot(TimeAxis,dataset[1,],type="l",xlim=c(0,T+1),ylim=c(min(dataset)-3,max(dataset)+3),
       xlab="Day of year", ylab=paste(y.var,"(*0.1 Degrees Centigrade)"),
       main=paste("Daily ",y.var,": all cities"))
  for (i in 2:Dim) {
    lines(TimeAxis,dataset[i,],col=color.set[i])
  }
  legend(x=350,legend=city.names,text.col=color.set)
  
  #Shift sd and fit with spline
  dataShift<-cbind(dataset[,282:365],dataset[,1:281])
  dataEstShift<-c()
  for (i in 1:Dim) {
    SS.fit<-smooth.spline(TimeAxis,dataShift[i,],df=SplineDF,all.knots=TRUE)
    row<-SS.fit$y
    dataEstShift<-rbind(dataEstShift,row)
  }
  row.names(dataEstShift)<-city.names
  #Shift back
  dataEst<-cbind(dataEstShift[,85:365],dataEstShift[,1:84])
  
  #Plot fitted data
  for (i in 1:Dim) {
    plot(TimeAxis,dataset[i,],type="l",col="Blue",xlim=c(0,T+1),ylim=c(min(dataset)-3,max(dataset)+3),
         xlab="Day of year", ylab=paste(y.var,"(*0.1 Degrees Centigrade)"),
         main=paste("Fitted Daily ",y.var,": ",city.names[i]))
    lines(TimeAxis,dataEst[i,],col="Red")
  }
  
  dev.off()
  return(dataEst)
}
