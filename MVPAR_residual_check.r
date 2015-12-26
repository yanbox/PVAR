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
par.est<-MVPAR.fit$par
wdir<-matrix(par.est[1:9],nrow=3,ncol=3)
g_rho<-par.est[10]

##############################################
#Solve for all Phi's and all Sigma's
##############################################
Phi<-list()
Sigma<-list()
invSigma<-list()
#Phi_nu(h) for nu goes from p+1 to T, and then from 1 to p,
#so that the last matrix G we construct is what we'll use later in the likelihood.
for (next_nu in p:(p+T-1)) {
  nu<-next_nu%%T+1
  G<-c()
  b<-c()
  for (j in 1:p) {
    col<-c()
    for (i in 1:p) {
      col<-rbind(col,gmm(j-i,nu-j,g_rho,wdir))
    }
    G<-cbind(G,col)
    b<-cbind(b,gmm(j,nu-j,g_rho,wdir))   #Correction(2015.08.05)
  }
  #Now solve for P: the "stacked" matrix of all Phi_{p+1}'s.
  invG<-solve(G)
  P<-b%*%invG
  #print(P)
  #Phi(1), ..., Phi(p) as a list
  Phi_nu<-list()
  for (i in 1:p) {
    Phi_nu[[i]]<-P[,((i-1)*d+1):(i*d)]
  }
  Phi[[nu]]<-Phi_nu
  
  #Now solve for Sigma_nu
  Sigma_nu<-gmm(0,nu,g_rho,wdir)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma_nu<-Sigma_nu-Phi[[nu]][[i]]%*%gmm(i-j,nu-j,g_rho,wdir)%*%t(Phi[[nu]][[j]])
    }
  }
  Sigma[[nu]]<-Sigma_nu
  invSigma[[nu]]<-solve(Sigma_nu)
}
detSigma<-c()
for (nu in 1:T) {
  detSigma<-c(detSigma,det(Sigma[[nu]]))
}

#Compute residuals
MVPAR.residuals<-c()
for (t in (p+1):RecordLength) {
  nu<-(t-1)%%T+1
  res.t<-Tmax_All_1[,t]
  for (i in 1:p) {
    res.t<-res.t-Phi[[nu]][[i]]%*%Tmax_All_1[,t-i]
  }
  MVPAR.residuals<-cbind(MVPAR.residuals,res.t)
}





MVPAR.var<-var(t(MVPAR.residuals))
MVPAR.var
#Compared to estimated Sigma of the parametric model
Sigma
#Residuals' lag h auto-covariance
h<-1
cov(t(VAR.residuals[,1:(RecordLength-p-h)]),t(VAR.residuals[,(1+h):(RecordLength-p)]))
#Q-Q plot of single entries
i<-1
qqnorm(VAR.residuals[i,])
#Plot residual against fitted values
i<-4
plot(Tmax_All_1[i,(p+1):RecordLength],MVPAR.residuals[i,])
Box.test(MVPAR.residuals[i,],lag=10)

png(filename="1112_ATL_Res_vs_Fitted_VPAR.png",width=800,height=400)
plot(Tmax_All_1[i,(p+1):RecordLength],MVPAR.residuals[i,])
dev.off()

xx<-Tmax_All_1[1,]
acf(xx,lag.max=1000)
pacf(xx,lag.max=1000)
arima(xx,order=c(8,0,1))
