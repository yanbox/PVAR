#############################################################################
#File Name: VAR_Try_more_values.r
#Author: Yanbo Xia
#Start Date: 2015.08.05
#Start Notes: Copied from ./Archive/150722_Try_more_values.r
#Description: Given different settings of parameters (listed in par_list),
#   generate random samples based on corresponding VAR(p) model.
#   For each sample, start from different initial estimations,
#   and find MLE using nlminb (requires package "MASS"?).
#
#Note: Before each run, change output file "fileName" 
#    to save to a different file.
#############################Edit Notes######################################
#Date: 2015.08.03
#Type: Bug detected
#Description: Tried a model without wdir[1] and wdir[2]
#   (simulation and MLE both using model with only (wdir[3],rho)),
#   estimated correctly.
#   In contrast, the full parameterized model yields different MLE 
#   using different initial estimations.

#Date: 2015.08.05
#Type: Bug correction
#Description: In defining nLikelihood(), when defining G and b,
#   changed b<-cbind(b,gmm(i,g_rho,wdir))
#   to      b<-cbind(b,gmm(j,g_rho,wdir))
#   ("i" changed to "j")

#Date: 2015.08.05
#Type: Automation improvement
#Description: Automatically set output file name using string manipulation
#   functions. The function str_c requires loading package "stringr".

#Date: 2015.08.07
#Type: Comment
#Description: (Final commit of the VAR model?) The model turns out to be 
#   identifiable as long as the position of the cities span the 2-D space.
#   The former problem occurs when d<3 (or when cities are situated on the
#   same straight line), however, I noticed that the values of
#   \hat{wdir[1]}^2+\hat{wdir[2]}^2+\hat{wdir[3]}^2
#   are consistent, and close to the true value.

#Date: 2015.08.09
#Type: Bug
#Description: invSigma are misused as InvSigma in some places.

#############################################################################
library(MASS)
#library("MASS", lib.loc="C:/Program Files/R/library")
#library("MASS", lib.loc="C:/Program Files/[Schoolwork]/R/library")
#Installing "stringr" 
# install.packages("stringr")
#Loading "stringr" 
library(stringr)

#Vector dimension: number of stations
d<-4
#Auto-Regressive degree
p<-1
#Total length of time in practice data.
N<-5000
#Coordinates of stations
x<-c(1,0,3,0,1); y<-c(0,1,2,0,2)
#Write results to this file.
#Output file name is based on current date.
today<-gsub("-","",Sys.Date())
fileName<-str_c("E:\\R_WorkSpace\\Output_",today,".txt")
# fileName<-"E:\\R_WorkSpace\\Output_150805.txt"
logFile<-"E:\\R_WorkSpace\\Log.txt"
cat("", file=fileName,sep="",append=FALSE)
cat("", file=logFile,sep="",append=FALSE)

#List of parameter sets to be tested.
par_list<-list(
#   c(7,3,5,1)
  c(0,0,0.5,1)
  ,c(-0.05,0.03,0.2,0.7)
  ,c(1,1.5,4,1)
  ,c(4,5,3,1)
  ,c(-0.15,0.2,0.5,1.2)
  )
# par_list<-list(c(-0.05,0.03,0.2,0.7))

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
gmm<-function(h,g_rho,wdir) {
  gmm<-matrix(,nrow=d,ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      gmm[i,j]<-exp(-dst(i,j,h,wdir)/g_rho)
    }
  }
  return(gmm)
}

#The likelihood:
nLikelihood<-function(par,X) {
#   print(par,digits=12)
  #   ParHist<<-cbind(ParHist,as.vector(par))
  wdir<-par[1:3]
  g_rho<-par[4]
  N<-ncol(X)
  
  #Construct Gamma's (for solving all Phi's)
  #Matrix G: "the big matrix"
  #(for the equation GX=b)
  G<-c()
  b<-c()
  for (j in 1:p) {
    col<-c()
    for (i in 1:p) {
      col<-rbind(col,gmm(j-i,g_rho,wdir))
    }
    G<-cbind(G,col)
    b<-cbind(b,gmm(j,g_rho,wdir))   #Correction(2015.08.05)
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
  Sigma<-gmm(0,g_rho,wdir)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma<-Sigma-Phi[[i]]%*%gmm(j-i,g_rho,wdir)%*%t(Phi[[j]])
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
      col<-rbind(col,gmm(j-i,g_rho,wdir))
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

est_list<-list()

#####################################################
#Simulate with each of the above parameter sets.
#####################################################
for (trial in 1:length(par_list)) {
  #####################################################
  #Generate data using the given parameters.
  #####################################################
  wdir<-par_list[[trial]][1:3]
  g_rho<-par_list[[trial]][4]
  crntEst<-c()
  
  #Matrix G: "the big matrix"
  #(for the equation GX=b)
  G<-c()
  b<-c()
  for (j in 1:p) {
    col<-c()
    for (i in 1:p) {
      col<-rbind(col,gmm(j-i,g_rho,wdir))
    }
    G<-cbind(G,col)
    b<-cbind(b,gmm(j,g_rho,wdir))
  }
#   print(G)
#   print(b)
  
  #Now solve for P: the "stacked" matrix of all Phi's.
  P<-b%*%solve(G)
  print(P)
  Phi<-list()
  for (i in 1:p) {
    next_Phi<-P[,((i-1)*d+1):(i*d)]
    Phi[[i]]<-next_Phi
  }
#   print(Phi)
  
  #Solve for Sigma
  Sigma<-gmm(0,g_rho,wdir)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma<-Sigma-Phi[[i]]%*%gmm(j-i,g_rho,wdir)%*%t(Phi[[j]])
    }
  }
  #Generate data
  All_X<-c()
  mup<-matrix(0,nrow<-d*p,ncol<-1)
  Xp<-as.matrix(mvrnorm(n=1,mup,G))
  for (t in 1:p) {
    All_X<-cbind(All_X,Xp[((p-t)*d+1):((p-t+1)*d),])
    #   All_X<-cbind(Xp[((t-1)*p+1):(t*p),],All_X)
  }
  mu<-matrix(0,nrow<-d,ncol<-1)
  Eps<-t(mvrnorm(n=N,mu,Sigma))
  for (t in (p+1):N) {
    #   Eps<-mvrnorm(n=1,mu,Sigma)
    #   nextX<-Eps
    nextX<-Eps[,t]
    for (i in 1:p) {
      nextX<-nextX+Phi[[i]]%*%All_X[,t-i]
    }
    All_X<-cbind(All_X,nextX)  
  }
#   print(All_X)
#   print("Sigma:")
#   print(Sigma)
#   print("The gammas are:")
#   for (i in 0:p) {
#     print(i)
#     print(gmm(i,g_rho,wdir))
#   }
#   print(Phi)

  #####################################################
  #Write to file: Summary of simulated data.
  #####################################################
  Y<-All_X
  mu1<-mean(Y[1,])
  mu2<-mean(Y[2,])
  muY<-c(mu1,mu2)
  Y<-Y-muY
  
  S<-list(Mtx0,Mtx0,Mtx0)
  for (h in 0:2) {
    for (i in 1:(N-h)) {
      y1<-as.matrix(Y[,i+h])
      y2<-as.matrix(Y[,i])
      S[[h+1]]<-S[[h+1]]+y1%*%t(y2)
    }
    S[[h+1]]<-S[[h+1]]/N
  }

  cat("#############################################################################", file=fileName,sep="\n",append=TRUE)
  cat("\tData generated using parameter set:", file=fileName,sep="\t",append=TRUE)
  cat(par_list[[trial]], file=fileName,sep=", ",append=TRUE)
  cat("", file=fileName,sep="\n",append=TRUE)
  cat("\tSeries length:", file=fileName,sep="\t",append=TRUE)
  cat(N, file=fileName,sep="\n",append=TRUE)
  cat("\tTheoretical Gamma(0):", file=fileName,sep="\t",append=TRUE)
  cat(gmm(0,g_rho,wdir), file=fileName,sep=", ",append=TRUE)
  cat("", file=fileName,sep="\n",append=TRUE)
  cat("\tSample Gamma(0):", file=fileName,sep="\t",append=TRUE)
  cat(S[[1]], file=fileName,sep=", ",append=TRUE)
  cat("", file=fileName,sep="\n",append=TRUE)
  cat("\tTheoretical Gamma(1):", file=fileName,sep="\t",append=TRUE)
  cat(gmm(1,g_rho,wdir), file=fileName,sep=", ",append=TRUE)
  cat("", file=fileName,sep="\n",append=TRUE)
  cat("\tSample Gamma(1):", file=fileName,sep="\t",append=TRUE)
  cat(S[[2]], file=fileName,sep=", ",append=TRUE)
  cat("", file=fileName,sep="\n",append=TRUE)
  cat("\tTheoretical Gamma(2):", file=fileName,sep="\t",append=TRUE)
  cat(gmm(2,g_rho,wdir), file=fileName,sep=", ",append=TRUE)
  cat("", file=fileName,sep="\n",append=TRUE)
  cat("\tSample Gamma(2):", file=fileName,sep="\t",append=TRUE)
  cat(S[[3]], file=fileName,sep=", ",append=TRUE)
  cat("", file=fileName,sep="\n",append=TRUE)
  cat("#############################################################################", file=fileName,sep="\n\n",append=TRUE)

  #####################################################
  #Starting from different values, compute MLE
  #####################################################
  stepSize<-sum(abs(par_list[[trial]][1:3])/3)
  for (st1 in (-2:2)*stepSize) {
    for (st2 in (-2:2)*stepSize) {
      for (st3 in (-2:2)*stepSize) {
        par<-par_list[[trial]]+c(st1,st2,st3,0)
        Iter<-0
        fit<-nlminb(start=par,objective=nLikelihood,X=All_X)
        nlik<-nLikelihood(fit$par,All_X)
        
        #Write to file:
        cat("Starting value: ", file=fileName,sep="\t",append=TRUE)
        cat(par, file=fileName,sep=", ",append=TRUE)
        cat("", file=fileName,sep="\n",append=TRUE)
        cat("Estimated value: ", file=fileName,sep="\t",append=TRUE)
        cat(fit$par, file=fileName,sep=", ",append=TRUE)
        cat("", file=fileName,sep="\n",append=TRUE)
        cat("Likelihood: ", file=fileName,sep="\t",append=TRUE)
        cat(sprintf("%.20f",nlik), file=fileName,sep=", ",append=TRUE)
        cat("", file=fileName,sep="\n",append=TRUE)
        cat("# of iterations: ", file=fileName,sep="\t",append=TRUE)
        cat(Iter, file=fileName,sep=", ",append=TRUE)
        cat("\n", file=fileName,sep="\n",append=TRUE)
        
        #Save estimated parameters to a list
        crntEst<-cbind(crntEst,fit$par)
      }
    }
  }
#   stepSize<-sum(abs(par_list[[trial]][1:3])/20)
#   for (st1 in (-10:10)*stepSize) {
#     par<-par_list[[trial]]+c(st1,st2,st3,0)
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
#     cat(nlik, file=fileName,sep=", ",append=TRUE)
#     cat("", file=fileName,sep="\n",append=TRUE)
#     cat("# of iterations: ", file=fileName,sep="\t",append=TRUE)
#     cat(Iter, file=fileName,sep=", ",append=TRUE)
#     cat("\n", file=fileName,sep="\n",append=TRUE)
#     
#     #Save estimated parameters to a list
#     est_list[[trial]]<-cbind(est_list[[trial]],fit$par)
  }

  est_list[[trial]]<-crntEst

}






