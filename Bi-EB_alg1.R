#####################
###Algorhitm 2 Bi-EB function
#####################

#################


alg2<-function(data,mu1,mu2,var1,var2,p1){
  mysum <- function(x) {
    sum(x[is.finite(x)])
  }
  
  beta1<-beta2<-apply(data,1,sd,na.rm=TRUE)
  alpha1<-alpha2<-apply(data,2,sd,na.rm=TRUE)
  
  p2<-1-p1
  loglik<- rep(NA, 200)
  loglik[1]<-0
  loglik[2]<-mysum((p1*p2*(((data-mu1-beta1-alpha1)^2)/(-2*var1))-log(sqrt(2*pi*var1)))+
                     ((1-p1*p2)*(((data-mu2-beta2-alpha2)^2)/(-2*var2))-log(sqrt(2*pi*var2)))+
                     p1*log(p1)+((1-p1)*log(1-p1))+
                     p2*log(p2)+((1-p2)*log(p2)))
  
  
  k=2
  
  while(abs(loglik[k]-loglik[k-1])>=0.00001){
  
    #E step
    
    phat1<-matrix(ncol=ncol(data),nrow=nrow(data))
    
    for(i in 1:nrow(data)){
      for(j in 1:ncol(data)){
        phat1[i,j]<-(p1 *p2*(1/sqrt(2*pi*var1))*exp(((data[i,j]-mu1-beta1[i]-alpha1[j])^2)/(-2*var1)))/
          ((p1*p2*(1/(sqrt(2*pi*var1))*exp(((data[i,j]-mu1-beta1[i]-alpha1[j])^2)/(-2*var1))))
           +(((1-p2)*(1-p1))*(1/(sqrt(2*pi*var2))*exp(((data[i,j]-mu2-beta2[i]-alpha2[j])^2)/(-2*var2))))
           +(((1-p2)*(p1))*(1/(sqrt(2*pi*var2))*exp(((data[i,j]-mu2-beta2[i]-alpha2[j])^2)/(-2*var2))))
           +(((p2)*(1-p1))*(1/(sqrt(2*pi*var2))*exp(((data[i,j]-mu2-beta2[i]-alpha2[j])^2)/(-2*var2)))))
        
      }
    }
    phat1[is.na(phat1)] <- 0.5
    
    
    phat2<-1-phat1
    
    #M step
    
    #beta1<-mysum(phat1*(apply(data,1,function(x) x-mu1))^2)/mysum(phat1)   
    
    mu1<-mysum(as.matrix(phat1*(data-beta1-alpha1)))/mysum(phat1)
    mu2<-mysum(as.matrix(phat2*(data-beta2-alpha2)))/mysum(phat2)
    
    var1<-mysum(as.matrix(phat1*(data-mu1-beta1-alpha1)^2))/mysum(phat1)  
    var2<-mysum(as.matrix(phat2*(data-mu2-beta2-alpha2)^2))/mysum(phat2) 
    
    beta1k<-apply(phat1*(data-mu1-alpha1),1,sum)/mysum(phat1) 
    
    beta2k<-apply(phat2*(data-mu2-alpha2),1,sum)/mysum(phat2) 
    
    alpha1k<-apply(phat1*(data-mu1-beta1),2,sum)/mysum(phat1) 
    
    alpha2k<-apply(phat2*(data-mu2-beta2),2,sum)/mysum(phat2) 
    beta1<-beta1k
    beta2<-beta2k
    alpha1<-alpha1k
    alpha2<-alpha2k
    
    p1<-mean(phat1)
    
    p2<-1-p1
    
    loglik[k+1]<- mysum((p1*p2*(((data-mu1-beta1-alpha1)^2)/(-2*var1))-log(sqrt(2*pi*var1)))+
                          ((1-p1*p2)*(((data-mu2-beta2-alpha2)^2)/(-2*var2))-log(sqrt(2*pi*var2)))+
                          p1*log(p1)+((1-p1)*log(1-p1))+
                          p2*log(p2)+((1-p2)*log(p2)))
    k<-k+1
    
  }
  list("mu" = c(mu1, mu2),
       "var" = c(var1, var2),
       "p" = c(p1, p2),
       "phat1"=phat1,
       "phat2"=phat2,
       "beta1mean"=mean(beta1),
       "beta1var"=var(beta1),
       "beta2mean"=mean(beta2),
       "beta2var"=var(beta2),
       "alpha1mean"=mean(alpha1),
       "alpha1var"=var(alpha1),
       "alpha2mean"=mean(alpha2),
       "alpha2var"=var(alpha2),
       "k"=k,
       "loglik"=loglik,
       "alpha1"=c(alpha1),
       "alpha2"=c(alpha2),
       "beta1"=c(beta1),
       "beta2"=c(beta2))
}
