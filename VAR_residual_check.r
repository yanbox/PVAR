#############################################################################
#File Name: VAR_residual_Check.r
#Author: Yanbo Xia
#Start Date: 2015.10.14
#Start Notes: 
#Description: Test residuals from the VAR model.
#
#############################Edit Notes######################################
#Date: 
#Type: 
#Description: 

#############################################################################

RecordLength<-ncol(Tmax_All_1)
Dim<-length(city.names)

#Estimated parameter values
par.est<-fit$par
wdir<-par.est[1:3]
rho<-par.est[4]

#Construct Gamma's (for solving all Phi's)
#Matrix G: "the big matrix"
#(for the equation GX=b)
G<-c()
b<-c()
for (j in 1:p) {
  col<-c()
  for (i in 1:p) {
    col<-rbind(col,gmm(j-i,rho,wdir))
  }
  G<-cbind(G,col)
  b<-cbind(b,gmm(j,rho,wdir))   #Correction(2015.08.05)
}
#print(G)
#print(b)

#Now solve for P: the "stacked" matrix of all Phi's.
P<-b%*%solve(G)
#print(P)
#Phi(1), ..., Phi(p) as a list
Phi<-list()
for (i in 1:p) {
  next_Phi<-P[,((i-1)*d+1):(i*d)]
  Phi[[i]]<-next_Phi
}
#   print(Phi)

#Solve for Sigma
Sigma<-gmm(0,rho,wdir)
for (i in 1:p) {
  for (j in 1:p) {
    Sigma<-Sigma-Phi[[i]]%*%gmm(j-i,rho,wdir)%*%t(Phi[[j]])
  }
}
invSigma<-solve(Sigma)

#Compute residuals
VAR.residuals<-matrix(0,nrow=Dim,ncol=RecordLength-p)
for (i in 1:p) {
  VAR.residuals<-VAR.residuals+Phi[[i]]%*%Tmax_All_1[,(p+1-i):(RecordLength-i)]
}


VAR.var<-var(t(VAR.residuals))
VAR.var
#Compared to estimated Sigma of the parametric model
Sigma
#Residuals' lag h auto-covariance
h<-1
cov(t(VAR.residuals[,1:(RecordLength-p-h)]),t(VAR.residuals[,(1+h):(RecordLength-p)]))
#Q-Q plot of single entries
i<-1
qqnorm(VAR.residuals[i,])
#Plot residual against fitted values
i<-1

png(filename="1111_ATL_Res_vs_Fitted.png",width=800,height=400)
plot(Tmax_All_1[i,(p+1):RecordLength],VAR.residuals[i,])
dev.off()

#Fit a scalar AR(8) model
ARMA.fit<-arima(Tmax_All_1[1,],order=c(8,0,1))
ARMA.fitted<-Tmax_All_1[1,]-ARMA.fit$residuals
png(filename="1111_ATL_Res_vs_Fitted_ARMA.png",width=800,height=400)
plot(ARMA.fitted,ARMA.fit$residuals)
dev.off()

#Plot residuals against nu
season<-rep(1:365,60)[1:21892]
png(filename="1119_ATL_Res_vs_season.png",width=800,height=500)
plot(season, VAR.residuals[1,])
dev.off()

#20151118
#Compare MLE of auto-correlation with sample auto-correlation
#MLE of auto-correlations
MaxLag<-50
AutoCorr.VAR<-c()
for (h in 0:MaxLag) {
  AutoCorr.VAR<-c(AutoCorr.VAR,gmm(h,rho,wdir)[1,1])
}
AutoCorr.VAR
#Sample auto-correlations
AutoCorr.Smpl<-c()
for (h in 0:MaxLag) {
  ACh<-cov(t(Tmax_All_1[,(h+1):RecordLength]),t(Tmax_All_1[,1:(RecordLength-h)]))
  AutoCorr.Smpl<-cbind(AutoCorr.Smpl,diag(ACh))
}
pdf(file="1119_Compare_AutoCorr.pdf",width=11,height=8.5)
for (city in city.names) {
  plot(0:MaxLag,AutoCorr.VAR,type="l",col="Red"
       ,main=paste("Compare Autocorrelations: ",city,sep=''))
  lines(0:MaxLag,AutoCorr.Smpl[city,])
}
dev.off()
#Compare cross-auto-correlations
#MLE of cross-auto-correlations
CrossAutoCorr.VAR<-c()
for (h in 0:MaxLag) {
  CrossAutoCorr.VAR<-cbind(CrossAutoCorr.VAR,gmm(h,rho,wdir)[1,2:4])
}
row.names(CrossAutoCorr.VAR)<-city.names[2:4]
#Sample cross-auto-correlations
CrossAutoCorr.Smpl<-c()
for (h in 0:MaxLag) {
  ACh<-cov(t(Tmax_All_1[,(h+1):RecordLength]),t(Tmax_All_1[,1:(RecordLength-h)]))
  CrossAutoCorr.Smpl<-cbind(CrossAutoCorr.Smpl,ACh[1,2:4])
}
pdf(file="1119_Compare_CrossAutoCorr.pdf",width=11,height=8.5)
for (city in city.names[2:4]) {
  plot(0:MaxLag,CrossAutoCorr.VAR[city,],type="l",col="Red"
       ,main=paste("Compare Cross-Autocorrelations: Atlanta with ",city,sep=''))
  lines(0:MaxLag,CrossAutoCorr.Smpl[city,])
}
dev.off()


#Multivariate Q-test
library('MTS')
mq
