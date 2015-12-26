#############################################################################
#File Name: Func_SampleCov.r
#Author: Yanbo Xia
#Start Date: 2015.09.03
#Description: Defines the SampleCov function.
#
#############################Edit Notes######################################
#Date: 
#Type: 
#Description: 

#############################################################################
#Given a data set containing (normalized) temperature data of all cities,
#the SampleCov function computes \Gamma(h) for h=0,1,...,Maxh.
#The results are returned in a list, where the [[k]]th element is \Gamma(k-1)
SampleCov<-function(dataset,Maxh) {
  RecordLength<-ncol(dataset)
  Dim<-nrow(dataset)
  SampleCov<-matrix(,nrow=Dim,ncol=Dim)
  
  #Substract mean
  Y<-c()
  mu<-c()
  for (i in 1:Dim) {
    Y<-rbind(Y,dataset[i,]-mean(dataset[i,]))
    mu<-c(mu,mean(dataset[i,]))
  }
  
  Gmm<-list()
  for (h in 0:Maxh) {
    S<-matrix(0,nrow=Dim,ncol=Dim)
    for (t in 1:(RecordLength-h)) {
      y1<-as.matrix(Y[,t+h])
      y2<-as.matrix(Y[,t])
      S<-S+y1%*%t(y2)
    }
    Gmm[[h+1]]<-S/RecordLength
  }
  
  #Record the mean vector.
  Gmm[[Maxh+2]]<-mu
  
  return(Gmm)
}




