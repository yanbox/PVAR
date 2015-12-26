#############################################################################
#File Name: Func_VAR_Likelihood_wsigma2.r
#Author: Yanbo Xia
#Start Date: 2015.09.03
#Start Notes: Copied from ./Repo/VAR_Try_more_values.r on 2015.09.03
#Description: Defines the likelihood function of the VAR model.
#
#Note: Before each run, change output file "fileName" 
#    to save to a different file.
#############################Edit Notes######################################
#Date: 2015.09.03
#Type: Modification
#Description: Added on variable to model: sigma^2
#   Changed "g_rho" into "rho"
#Note: This is a subversion of the original file "Func_VAR_Likelihood.r".
#Note: Currently we'll still use the model without sigma^2, since the data
#   has already been normalized.

#############################################################################
library(MASS)
#library("MASS", lib.loc="C:/Program Files/R/library")
#library("MASS", lib.loc="C:/Program Files/[Schoolwork]/R/library")
#Installing "stringr" 
# install.packages("stringr")
#Loading "stringr" 
library(stringr)

#Write results to this file.
#Output file name is based on current date.
today<-gsub("-","",Sys.Date())
fileName<-str_c("E:\\R_WorkSpace\\Output_",today,".txt")
# fileName<-"E:\\R_WorkSpace\\Output_150805.txt"
logFile<-"E:\\R_WorkSpace\\Log.txt"
cat("", file=fileName,sep="",append=FALSE)
cat("", file=logFile,sep="",append=FALSE)

#Matrices
Mtx0<-matrix(0,nrow=d,ncol=d)
Mtx1<-diag(1,d)

#Modified time-distance at lag h=t_i-t_j, given "wind direction" wdir
dst<-function(i,j,h,wdir) {
  dx<-x[i]-x[j]+wdir[1]*h
  dy<-y[i]-y[j]+wdir[2]*h
  dt<-wdir[3]*h
  return(sqrt(dx^2+dy^2+dt^2))
  #   return(dx^2+dy^2+dt^2)
}

#Construct Gamma(h)
gmm<-function(h,sig2,rho,wdir) {
  gmm<-matrix(,nrow=d,ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      gmm[i,j]<-sig2*exp(-dst(i,j,h,wdir)/rho)
    }
  }
  return(gmm)
}

#The likelihood:
nLikelihood<-function(par,X) {
  #   print(par,digits=12)
  #   ParHist<<-cbind(ParHist,as.vector(par))
  wdir<-par[1:3]
  sig2<-par[4]
  rho<-par[5]
  N<-ncol(X)
  
  #Construct Gamma's (for solving all Phi's)
  #Matrix G: "the big matrix"
  #(for the equation GX=b)
  G<-c()
  b<-c()
  for (j in 1:p) {
    col<-c()
    for (i in 1:p) {
      col<-rbind(col,gmm(j-i,sig2,rho,wdir))
    }
    G<-cbind(G,col)
    b<-cbind(b,gmm(j,sig2,rho,wdir))   #Correction(2015.08.05)
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
  Sigma<-gmm(0,sig2,rho,wdir)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma<-Sigma-Phi[[i]]%*%gmm(j-i,sig2,rho,wdir)%*%t(Phi[[j]])
    }
  }
  invSigma<-solve(Sigma)
  
  #Series {W_t}
  W<-c()
  for (t in (p+1):N) {
    nextW<-X[,t]
    for (k in 1:p) {
      nextW<-nextW-Phi[[k]]%*%X[,t-k]
    }
    W<-cbind(W,nextW)
  }
  #print(W)
  
  #The covariance matrix G_p:
  #Note(20150805): G_p is essentially the same as G.
  Gp<-c()
  for (j in 1:p) {
    col<-c()
    for (i in 1:p) {
      col<-rbind(col,gmm(j-i,sig2,rho,wdir))
    }
    Gp<-cbind(Gp,col)
  }
  #print(Gp)
  invGp<-solve(Gp)
  #print(invGp)
  
  #The vector of the first p observations: XXp
  XXp<-c()
  for (t in 1:p) {
    XXp<-rbind(as.matrix(X[,t]),XXp)
  }
  #print(XXp)
  
  
  #Likelihood
  #"Le" denotes the exponential part of the likelihood
  Le<-t(XXp)%*%invGp%*%XXp
  for (t in (p+1):N) {
    Wt<-as.matrix(W[,t-p])
    Le<-Le+t(Wt)%*%invSigma%*%Wt
  }
  nlogL<-Le/2+log(det(Gp))/2+log(det(Sigma))*(N-p)/2
  #   print("negative log likelihood:")
  #   print(c(nlogL,Le/2,log(det(Gp))/2,log(det(Sigma))*(N-p)/2),digit=15)
  
  Iter<<-Iter+1
  cat(Iter, file=logFile,sep=" ",append=TRUE)
  cat(" :: ", file=logFile,sep=" ",append=TRUE)
  cat(sprintf("%.12f",par), file=logFile,sep=", ",append=TRUE)
  cat(" :: ", file=logFile,sep=" ",append=TRUE)
  cat(sprintf("%.20f",nlogL), file=logFile,sep="\n",append=TRUE)
  return(nlogL)
}
