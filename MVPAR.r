#############################################################################
#File Name: MVPAR.r
#Author: Yanbo Xia
#Start Date: 2015.08.08
#Description: Given different settings of parameters (listed in par_list),
#   generate random samples based on corresponding MVPAR(p) model, and 
#   compute MLE.
#
#Note: Before each run, change output file "fileName" 
#    to save to a different file.
#############################Edit Notes######################################
#Date: 2015.08.08
#Type: 
#Description: 

#############################################################################

library(MASS)
library(stringr)

#Vector dimension: number of stations
d<-4
#Auto-Regressive degree
p<-3
#Time period
T<-365
#Total length of time in practice data.
#N<-5000
#Coordinates of stations
#x<-c(1,0,3,0,1); y<-c(0,1,2,0,2)
x<-longitude; y<-latitude
#Write results to this file.
#Output file name is based on current date.
today<-gsub("-","",Sys.Date())
fileName<-str_c("E:\\R_WorkSpace\\Output_",today,".txt")
# fileName<-"E:\\R_WorkSpace\\Output_150805.txt"
logFile<-"E:\\R_WorkSpace\\Log.txt"
cat("", file=fileName,sep="",append=FALSE)
cat("", file=logFile,sep="",append=FALSE)

par_list<-list(
  c(0.3,-0.1,0.15,-0.2,0.05,0,0.2,0,0.2,1)
  )

#The modified distance
dst<-function(i,j,h,nu,wdir) {
  dx<-(x[i]-x[j])+wdir[1,1]*h+
    wdir[2,1]*(cos(2*pi*(nu+h)/T)-cos(2*pi*nu/T))+
    wdir[3,1]*(sin(2*pi*(nu+h)/T)-sin(2*pi*nu/T))
  dy<-(y[i]-y[j])+wdir[1,2]*h+
    wdir[2,2]*(cos(2*pi*(nu+h)/T)-cos(2*pi*nu/T))+
    wdir[3,2]*(sin(2*pi*(nu+h)/T)-sin(2*pi*nu/T))
  dt<-wdir[1,3]*h+
    wdir[2,3]*(cos(2*pi*(nu+h)/T)-cos(2*pi*nu/T))+
    wdir[3,3]*(sin(2*pi*(nu+h)/T)-sin(2*pi*nu/T))
  return(sqrt(dx^2+dy^2+dt^2))
}

#Gamma_nu(h)
gmm<-function(h,nu,g_rho,wdir) {
  gmm<-matrix(,nrow=d,ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      gmm[i,j]<-exp(-dst(i,j,h,nu,wdir)/g_rho)
    }
  }
  return(gmm)
}

nL.MVPAR<-function(par,X) {
  wdir<-matrix(par[1:9],nrow=3,ncol=3)
  g_rho<-par[10]
  N<-ncol(X)
  
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
  
  ##############################################
  #Prepare the series
  ##############################################
  #Series W_t
  W<-c()
  for (t in (p+1):N) {
    nextW<-X[,t]
    nu<-(t-1)%%T+1
    for (k in 1:p) {
      nextW<-nextW-Phi[[nu]][[k]]%*%X[,t-k]
    }
    W<-cbind(W,nextW)
  }
  #The vector of the first p observations: XXp
  XXp<-c()
  for (t in 1:p) {
    XXp<-rbind(as.matrix(X[,t]),XXp)
  }
  
  ##############################################
  #Likelihood
  ##############################################
  #"Le" denotes the exponential part of the likelihood
  Le<-t(XXp)%*%invG%*%XXp
  for (t in (p+1):N) {
    nu<-(t-1)%%T+1
    Wt<-as.matrix(W[,t-p])
    Le<-Le+t(Wt)%*%invSigma[[nu]]%*%Wt
  }
  nlogL<-Le/2+log(det(G))/2
  for (nu in 1:T) {
    if ((nu-p-1)%%T<=(N-p-1)%%T) {
      nlogL<-nlogL+log(detSigma[nu])*((N-p)%/%T+1)
    } else {
      nlogL<-nlogL+log(detSigma[nu])*((N-p)%/%T)
    }
  }
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
# 
# est_list<-list()
# 
# #####################################################
# #Simulate with each of the above parameter sets.
# #####################################################
# for (trial in 1:length(par_list)) {
#   #####################################################
#   #Generate data using the given parameters.
#   #####################################################
#   wdir<-matrix(par_list[[trial]][1:9],nrow=3,ncol=3)
#   g_rho<-par_list[[trial]][10]
#   crntEst<-c()
#   
#   ##############################################
#   #Solve for all Phi's and all Sigma's
#   ##############################################
#   Phi<-list()
#   Sigma<-list()
#   invSigma<-list()
#   #Phi_nu(h) for nu goes from p+1 to T, and then from 1 to p,
#   #so that the last matrix G we construct is what we'll use later in the likelihood.
#   for (next_nu in p:(p+T-1)) {
#     nu<-next_nu%%T+1
#     G<-c()
#     b<-c()
#     for (j in 1:p) {
#       col<-c()
#       for (i in 1:p) {
#         col<-rbind(col,gmm(j-i,nu-j,g_rho,wdir))
#       }
#       G<-cbind(G,col)
#       b<-cbind(b,gmm(j,nu-j,g_rho,wdir))   #Correction(2015.08.05)
#     }
#     #Now solve for P: the "stacked" matrix of all Phi_{p+1}'s.
#     invG<-solve(G)
#     P<-b%*%invG
#     #print(P)
#     #Phi(1), ..., Phi(p) as a list
#     Phi_nu<-list()
#     for (i in 1:p) {
#       Phi_nu[[i]]<-P[,((i-1)*d+1):(i*d)]
#     }
#     Phi[[nu]]<-Phi_nu
#     
#     #Now solve for Sigma_nu
#     Sigma_nu<-gmm(0,nu,g_rho,wdir)
#     for (i in 1:p) {
#       for (j in 1:p) {
#         Sigma_nu<-Sigma_nu-Phi[[nu]][[i]]%*%gmm(i-j,nu-j,g_rho,wdir)%*%t(Phi[[nu]][[j]])
#       }
#     }
#     Sigma[[nu]]<-Sigma_nu
#     invSigma[[nu]]<-solve(Sigma_nu)
#   }
#   detSigma<-c()
#   for (nu in 1:T) {
#     detSigma<-c(detSigma,det(Sigma[[nu]]))
#   }
#   
#   #Generate data
#   All_X<-c()
#   mup<-matrix(0,nrow<-d*p,ncol<-1)
#   Xp<-as.matrix(mvrnorm(n=1,mup,G))
#   for (t in 1:p) {
#     All_X<-cbind(All_X,Xp[((p-t)*d+1):((p-t+1)*d),])
#     #   All_X<-cbind(Xp[((t-1)*p+1):(t*p),],All_X)
#   }
#   mu<-matrix(0,nrow<-d,ncol<-1)
#   for (t in (p+1):N) {
#     nu<-(t-1)%%T+1
#     nextX<-mvrnorm(n=1,mu,Sigma[[nu]])
#     for (i in 1:p) {
#       nextX<-nextX+Phi[[nu]][[i]]%*%All_X[,t-i]
#     }
#     All_X<-cbind(All_X,nextX)
#   }
#   
#   #####################################################
#   #Starting from different values, compute MLE
#   #####################################################
#   nChangeStrt<-3
#   true_par<-par_list[[trial]]
#   stepSize<-sum(true_par)/10
#   for (estTrial in 1:20) {
#     par<-true_par
#     changePos<-sample(1:9,size=nChangeStrt,replace=TRUE)
#     for (iChangeStrt in 1:nChangeStrt) {
#       if (rbinom(1,2,0.5)) {
#         par[changePos[iChangeStrt]]<-par[changePos[iChangeStrt]]+stepSize
#       } else {
#         par[changePos[iChangeStrt]]<-par[changePos[iChangeStrt]]-stepSize
#       }
#     }
#     Iter<-0
#     fit<-nlminb(start=par,objective=nLikelihood,X=All_X)
#     nlik<-nLikelihood(fit$par,All_X)
#     
#     #Write to file:
#     cat("Starting value: ", file=fileName,sep="\t",append=TRUE)
#     cat(par, file=fileName,sep=", ",append=TRUE)
#     cat("", file=fileName,sep="\n",append=TRUE)
#     cat("Estimated value: ", file=fileName,sep="\t",append=TRUE)
#     cat(fit$par, file=fileName,sep=", ",append=TRUE)
#     cat("", file=fileName,sep="\n",append=TRUE)
#     cat("Likelihood: ", file=fileName,sep="\t",append=TRUE)
#     cat(sprintf("%.20f",nlik), file=fileName,sep=", ",append=TRUE)
#     cat("", file=fileName,sep="\n",append=TRUE)
#     cat("# of iterations: ", file=fileName,sep="\t",append=TRUE)
#     cat(Iter, file=fileName,sep=", ",append=TRUE)
#     cat("\n", file=fileName,sep="\n",append=TRUE)
#     
#     #Save estimated parameters to a list
#     crntEst<-cbind(crntEst,fit$par)
#   }
# 
# }

Time.Strt<-Sys.time()
MVPAR.fit<-nlminb(start=c(0.3,-0.1,0.15,-0.2,0.05,0,0.2,0,0.2,1),objective=nL.MVPAR,X=Tmax_All_1)
Time.End<-Sys.time()













