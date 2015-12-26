#############################################################################
#File Name: DataProcessing.r
#Author: Yanbo Xia
#Start Date: 2015.09.02
#Description: Prepare and fit real temperature data.
#
#############################Edit Notes######################################
#Date: 2015.09.27
#Type: Correction
#Description: In function MeanZero, "#Substract fitted mean from data",
#   I substracted mu instead of muEst. Corrected this.

#Date: 2015.12.26
#Type: Migration
#Description: Migrate to Palmetto. Start using GitHub?

#############################################################################
#Current time
Time.Start<-Sys.time()
#Spline degree of freedom
SplineDF<-21

# WorkSPC<-"/home/yanbox/MVPARMA/WorkSpace"
# RepoDir<-"/home/yanbox/MVPARMA/Repo"
DBDrive<-"C"
WorkSPC<-"C:\\WorkSpace\\MVPARMA_WorkSpace"
RepoDir<-paste(DBDrive,":\\Dropbox\\[Graduate Courses]\\[Research]1403 PARMA\\Programs\\Repo")
setwd(WorkSPC)
city.names<-c("Atlanta","Charlotte","Columbia","Raleigh")
StationData<-list()
T<-365    #Period=365days(1year)
#Location of Dropbox
DBDrive<-"D"
#Define the PlotNFit function:
source(paste(RepoDir,"\\Func_PlotNFit.r",sep=''))
#Standardizing function: SDNormalize
source(paste(RepoDir,"\\Func_SDNormalize.r",sep=''))
#VAR Likelihood function:
source(paste(RepoDir,"\\Func_VAR_Likelihood.r",sep=''))
#Sample AutoCovariance function:
source(paste(RepoDir,"\\Repo\\Func_SampleCov.r",sep=''))

Tmax_All<-c()
longitude<-c()
latitude<-c()
for (city in city.names) {
  file.name<-paste("Station_",city,".csv",sep='')
  StationData[[city]]<-read.csv(file.name,header=TRUE)
  tempSet<-StationData[[city]]$TMAX/10
  Tmax_All<-rbind(Tmax_All,tempSet)
  longitude<-c(longitude,mean(StationData[[city]]$LONGITUDE))
  latitude<-c(latitude,mean(StationData[[city]]$LATITUDE))
  
  #Plot raw Tmax series
  TimeAxis<-(1:length(tempSet))/T
  file.name<-paste("RawTmax_Station_",city,".png",sep='')
  png(filename=file.name,width=800,height=400)
  plot(TimeAxis,tempSet,type="l",col="Black",
       xlab="Year",ylab=expression(paste("Tmax(",degree,"C)",sep='')),
       #main=paste("Daily Maximum Temperature: ",city,sep='')
       )
  dev.off()
}

MeanZero<-function(dataset) {
  RecordLength<-ncol(dataset)
  Dim<-nrow(dataset)
  NYears<-RecordLength%/%T
  DaysRes<-RecordLength%%T
  TimeAxis<-1:T
  mu<-matrix(0,nrow=Dim,ncol=T)
  for (t in 1:RecordLength) {
    mu[,(t-1)%%T+1]<-mu[,(t-1)%%T+1]+dataset[,t]
  }
  if (DaysRes>0) {
    for (nu in 1:DaysRes) {
      mu[,nu]<-mu[,nu]/(NYears+1)
    }
  }
  for (nu in (DaysRes+1):T) {
    mu[,nu]<-mu[,nu]/NYears
  }
  muShift<-cbind(mu[,282:365],mu[,1:281])
  muEstShift<-c()
  for (i in 1:Dim) {
    SS.fit<-smooth.spline(TimeAxis,muShift[i,],df=SplineDF,all.knots=TRUE)
    row<-SS.fit$y
    muEstShift<-rbind(muEstShift,row)
  }
  row.names(muEstShift)<-city.names
  #Shift back
  muEst<-cbind(muEstShift[,85:365],muEstShift[,1:84])
  
  #Plots
  res<-mu-muEst
  for (i in 1:Dim) {
    #Plot fitted daily mean
    file.name<-paste("RPlot_FittedDailyMean_",city.names[i],".png",sep='')
    png(filename=file.name,width=800,height=400)
    plot(TimeAxis,mu[i,],type="l",col="Blue",
         xlab="Day of year",ylab=expression(paste("Daily Average( ",degree,"C)",sep='')),
         #main=paste("Fitted daily mean: ",city.names[i],sep='')
         )
    lines(TimeAxis,muEst[i,],col="Red")
    dev.off()
    #Plot residuals of fitted daily mean
    file.name<-paste("RPlot_DailyRes_",city.names[i],".png",sep='')
    png(filename=file.name,width=800,height=400)
    plot(TimeAxis,res[i,],type="l",
          xlab="Day of year",ylab=expression(paste("Residual( ",degree,"C)",sep='')),
         #main=paste("Residuals of fitted daily mean: ",city.names[i],sep='')
         )
    dev.off()
  }
  mu<<-mu
  muEst<<-muEst
    
  #Substract fitted mean from data
  NewData<-c()
  for (i in 1:Dim) {
    if (DaysRes==0) {
      row<-dataset[i,]-muEst[i,]
    } else {
      row<-c(rep(muEst[i,],NYears),muEst[1:DaysRes,])
      row<-dataset[i,]-row
    }
    NewData<-rbind(NewData,row)
  }
  row.names(NewData)<-city.names
  return(NewData)
}

Tmax_All_0<-MeanZero(Tmax_All)

Tmax_All_1<-SDNormalize(Tmax_All_0)

#Plot normalized data
RecordLength<-ncol(Tmax_All_1)
Dim<-nrow(Tmax_All_1)
pdf(file="Normalized_data.pdf",width=160,height=10)
for (i in 1:Dim) {
  plot((1:RecordLength)/365,Tmax_All_1[i,],type="l",
       xlab="Year",ylab="Tmax(*0.1 Degrees Centigrade)",
       main=paste("Normalized Tmax: ",city.names[i],sep=''))
}
dev.off()
#Plot only 3 years
RecordLength<-365*3
pdf(file="Normalized_data_3yrs.pdf",width=16,height=10)
for (i in 1:Dim) {
  plot((1:RecordLength)/365,Tmax_All_1[i,1:RecordLength],type="l",
       xlab="Year",ylab="Tmax(*0.1 Degrees Centigrade)",
       main=paste("Normalized Tmax: ",city.names[i],sep=''))
}
dev.off()

#Fit MVPAR model
#Vector dimension: number of stations
d<-nrow(Tmax_All_1)
#Auto-Regressive degree
p<-8
#Total length of time in practice data.
N<-ncol(Tmax_All_1)
#Coordinates of stations
x<-longitude; y<-latitude
Iter<-0
par<-c(0,0,1,1)
fit<-nlminb(start=par,objective=nLikelihood,X=Tmax_All_1)
nlik<-nLikelihood(fit$par,Tmax_All_1)

#Current time
Time.End<-Sys.time()



#Plot residuals of residuals
# file.name<-"TMP_RPlot_Residuals_Atlanta_long.png"
# png(filename=file.name,width=16000,height=1000)
# plot((1:21900)/365,Tmax_All_0[1,],type="l",
#      xlab="Year",ylab="Residual(*0.1 Degrees Centigrade)",
#      main="Residuals of all 60 years: Atlanta")
# dev.off()
# file.name<-"TMP_RPlot_Residuals_Atlanta_3yrs.png"
# png(filename=file.name,width=1600,height=1000)
# plot((1:(3*365))/365,Tmax_All_0[1,1:(3*365)],type='l',
#      xlab="Year",ylab="Residual(*0.1 Degrees Centigrade)",
#      main="Residuals of first 3 years: Atlanta")
# dev.off()










# TimeAxis<-1:T
# plot(TimeAxis,mu[1,],type="l",xlim=c(0,366),ylim=c(85,345),
#      xlab="Day of year",ylab="Tmax")
# lines(TimeAxis,mu[2,],col="Red")
# lines(TimeAxis,mu[3,],col="Blue")
# lines(TimeAxis,mu[4,],col="Green")
# legend(x=350,legend=c("Atlanta","Charlotte","Columbia","Raleigh"),
#         text.col=c("Black","Red","Blue","Green"))
# abline(v=281,col="Purple")
# 
# mu1<-cbind(mu[,282:365],mu[,1:281])
# plot(TimeAxis,mu1[1,],type="l",xlim=c(0,366),ylim=c(85,345),
#      xlab="Day of year",ylab="Tmax")
# lines(TimeAxis,mu1[2,],col="Red")
# lines(TimeAxis,mu1[3,],col="Blue")
# lines(TimeAxis,mu1[4,],col="Green")
# legend(x=350,legend=c("Atlanta","Charlotte","Columbia","Raleigh"),
#        text.col=c("Black","Red","Blue","Green"))
# #abline(v=281,col="Purple")
# 
# x<-mu1[1,]
# SS.fit1<-smooth.spline(TimeAxis,x,df=9,all.knots=TRUE)
# plot(TimeAxis,x,type="l",main="Daily average w/ fitted smoothing spline",ylab="Tmax (degrees Centigrade)",xlab="nu",col=3)
# lines(SS.fit1,col=2)
# #Plot residuals
# plot(SS.fit1$y-mu1.ATL,type="l",col=4)
