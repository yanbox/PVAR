MeanZero_tmp<-function(dataset) {
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
  
  
  mu<<-mu
  muEst<<-muEst
  
  #Substract fitted mean from data
  NewData<-c()
  for (i in 1:Dim) {
    if (DaysRes==0) {
      row<-dataset[i,]-mu[i,]
    } else {
      row<-c(rep(mu[i,],NYears),mu[1:DaysRes,])
      row<-dataset[i,]-row
    }
    NewData<-rbind(NewData,row)
  }
  row.names(NewData)<-city.names
  return(NewData)
}

tempdata<-MeanZero(Tmax_All)

##########################################################
##########################################################
##########################################################


atl<-tempdata[1,]
clt<-tempdata[2,]

h<-3
N<-length(atl)



AutoCov.Atl<-matrix(0,nrow=h+1,ncol=365)
for (i in 0:h) {
  for (t in 1:(N-i)) {
    nu<-(t-1)%%365+1
    AutoCov.Atl[i+1,nu]<-AutoCov.Atl[i+1,nu]+atl[t]*atl[t+i]
  }
}
AutoCov.Atl<-AutoCov.Atl/60

TimeAxis<-1:365
file.name<-"RPlot_DailySampleVar_ATL.png"
png(filename=file.name,width=800,height=400)
plot(TimeAxis,AutoCov.Atl[1,],type="l",
     xlab="Day of year",ylab="Daily Sample Variance",
     cex.lab=1.3
     #main=paste("Fitted daily mean: ",city.names[i],sep='')
)
dev.off()

file.name<-"RPlot_DailyAutoCov3_ATL.png"
png(filename=file.name,width=800,height=400)
plot(TimeAxis,AutoCov.Atl[4,],type="l",
     xlab="Day of year",ylab="Daily Sample AutoCov (lag 3)",
     cex.lab=1.3
     #main=paste("Fitted daily mean: ",city.names[i],sep='')
)
dev.off()

#AutoCorrelation
AutoCor.Atl<-AutoCov.Atl
for (i in 0:h) {
  for (t in 1:365) {
    t1<-(t+i-1)%%365+1
    AutoCor.Atl[i+1,t]<-AutoCor.Atl[i+1,t]/sqrt(AutoCov.Atl[1,t]*AutoCov.Atl[1,t1])
  }
}
#Plot AutoCorrelation
file.name<-"RPlot_DailyAutoCor3_ATL.png"
png(filename=file.name,width=800,height=400)
plot(TimeAxis,AutoCor.Atl[4,],type="l",
     xlab="Day of year",ylab="Daily Sample AutoCor (lag 3)",
     cex.lab=1.3
     #main=paste("Fitted daily mean: ",city.names[i],sep='')
)
dev.off()

#Cross autocovariance ATL, CLT
Cov.Atl<-AutoCov.Atl[1,]
Cov.Clt<-rep(0,365)
CrsCov<-matrix(0,nrow=7,ncol=365)
CrsCor<-CrsCov
for (t in 1:N) {
  nu<-(t-1)%%365+1
  Cov.Clt[nu]<-Cov.Clt[nu]+clt[t]^2
}
Cov.Clt<-Cov.Clt/60
for (h in -3:3) {
  if (h<0) {
    timespand<-(1-h):N
  } else {
    timespand<-1:(N-h)
  }
  for (t in timespand) {
    nu<-(t-1)%%365+1
    CrsCov[4+h,nu]<-CrsCov[h+4,nu]+atl[t+h]*clt[t]
  }
}
CrsCov<-CrsCov/60
for (h in -3:3) {
  for (nu in 1:365) {
    nu1<-(nu+h-1)%%365+1
    CrsCor[4+h,nu]<-CrsCov[4+h,nu]/sqrt(Cov.Atl[nu1]*Cov.Clt[nu])
  }
}
#Plots
TimeAxis<-1:365
for (h in -3:3) {
  file.name<-paste("RPlot_CrossCov",as.character(h),".png",sep='')
  png(filename=file.name,width=800,height=400)
  plot(TimeAxis,CrsCov[4+h,],type="l",
       xlab="Day of year",ylab=paste("Daily Sample Cross Cov (lag ",as.character(h),")"),
       cex.lab=1.3
       #main=paste("Fitted daily mean: ",city.names[i],sep='')
  )
  dev.off()
  file.name<-paste("RPlot_CrossCor",as.character(h),".png",sep='')
  png(filename=file.name,width=800,height=400)
  plot(TimeAxis,CrsCor[4+h,],type="l",
       xlab="Day of year",ylab=paste("Daily Sample Cross Cor (lag ",as.character(h),")"),
       cex.lab=1.3
       #main=paste("Fitted daily mean: ",city.names[i],sep='')
  )
  dev.off()
}

