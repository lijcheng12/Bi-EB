############################################################
###    Step 3. Performance evaluation of Bi-clustering   ###
###    Recovery and Relevance Function                   ###
############################################################

##########
#bicres1>> TRUE set  ; bicres2>>RESULT set

#T : bicluster elements in the true set
#R : bicluster elements in the true set

#function parameters: biclusters dimensions, bicluster output of algorithm
#################

Recovery<-function(row.bic11,row.bic21,col.bic11,col.bic21,
                   row.bic12,row.bic22,col.bic12,col.bic22,
                   row.bic13,row.bic23,col.bic13,col.bic23,
                   bicres2,n){
  le2 <- bicres2@Number
  jacvec <- c()
  number<-c()
  Jmax<-c()

  bic.dim<-c(row.bic11,row.bic21,col.bic11,col.bic21,
             row.bic12,row.bic22,col.bic12,col.bic22,
             row.bic13,row.bic23,col.bic13,col.bic23)
  for (i in seq(from=1,to=9,by=4)) {
    jacvec2 <- c()
    alle1<-matrix(0,nrow = 200,ncol=300)
    
    alle1[c(bic.dim[i]:bic.dim[i+1]),c(bic.dim[i+2]:bic.dim[i+3])]<-1
    
    for (j in 1:le2) {
      
      # Taking into account the case when NUMBER=1 and the output is in WRONG FORMAT
      #	Note on 'Wrong Format': This occurs only when the expected BC for the algorithm is >1, but only 1 BC is discovered. (NOT when you are asking an algorithm for 1 BC)
      #		Why? as.matrix() of a vector results into a matrix with 1 column always
      
      if(le2==1 & (dim(bicres2@NumberxCol)[2]==1)){
        res2.NumberxCol <- bicres2@NumberxCol[,1]
        
      }
      else{
        res2.NumberxCol <- bicres2@NumberxCol[j,]
      }
     
      alle2 <- bicres2@RowxNumber[, j] %*% t(res2.NumberxCol)
      alle <- alle1 + alle2
      loalle <- alle > 0
      loalle1 <- alle1 > 0
      loalle2 <- alle2 > 0
      jacvec2 <- c(jacvec2,((sum(loalle1) + sum(loalle2) -sum(loalle))/sum(loalle))) 
      ##up to here this will be comparision of each biclusters in result with one bicluster (i) in the true set.
      ##So jacvec2 is the vector of J index in recovery for i=1
      #then finding max of jacvec2 for each i
      #next is to find vectors of J for all i in TRUE set 
      #and finding max of them
      #then sumation of max over i
      #then devide by number of elements in True set
      #sum(loalle2) is the number of elements in bicluster j Result set
      #sum(loalle1) is the number of elements in bicluster i True set
    }
    number<-c(number,sum(loalle1))
    Jmax<-c(Jmax,max(jacvec2))
  }
  #T.elements<-sum(number)
  jacvec<-sum(Jmax)
  rec <- jacvec/n
  return(rec)
}



##Relevance function

Relevance<-function(row.bic11,row.bic21,col.bic11,col.bic21,
                    row.bic12,row.bic22,col.bic12,col.bic22,
                    row.bic13,row.bic23,col.bic13,col.bic23,
                    bicres2){
  le2 <- bicres2@Number
  jacvec <- c()
  number<-c()
  Jmax<-c()
  #  bic.dim<-c(1,25,1,25,150,200,250,300,0,0,0,0)
  
  bic.dim<-c(row.bic11,row.bic21,col.bic11,col.bic21,
             row.bic12,row.bic22,col.bic12,col.bic22,
             row.bic13,row.bic23,col.bic13,col.bic23)
  for (j in 1:le2)
  {
    jacvec2 <- c()
    
    if(le2==1 & (dim(bicres2@NumberxCol)[2]==1)){
      res2.NumberxCol <- bicres2@NumberxCol[,1]
      
    }
    else{
      res2.NumberxCol <- bicres2@NumberxCol[j,]
    }
    
    for (i in seq(from=1,to=9,by=4)) {
      
      
      alle1<-matrix(0,nrow = 200,ncol=300)
      
      alle1[c(bic.dim[i]:bic.dim[i+1]),c(bic.dim[i+2]:bic.dim[i+3])]<-1
      
      
      alle2 <- bicres2@RowxNumber[, j] %*% t(res2.NumberxCol)
      alle <- alle1 + alle2
      loalle <- alle > 0
      loalle1 <- alle1 > 0
      loalle2 <- alle2 > 0
      jacvec2 <- c(jacvec2,((sum(loalle1) + sum(loalle2) -sum(loalle))/sum(loalle))) 

    }
    number<-c(number,sum(loalle2))
    Jmax<-c(Jmax,max(jacvec2))
  }
  #R.elements<-sum(number)
  jacvec<-sum(Jmax)
  rel <- jacvec/le2
  return(rel)
}

