#############################################################################
###################### Step 3. Searching bicluster in Bi-EB    ##############
#############################################################################

#bic final result

#alg1.phat1<-data.frame(resultalg1$phat1)

##putting rows together as one box

#sort(apply(resultalg1$phat1,1,function(x) sum(x>=.8)))

list1<-apply(resultalg1$phat1,1,function(x) which(x>=.8,arr.ind = TRUE))


#taking genes with those samples

library (plyr)

df <- ldply (list1, rbind)
for(i in 1:nrow(df)){
  row.names(df)[i]<-paste("gene", i, sep = "")   
}


#apply(df,1,function(x) alg1.phat1[i,names(alg1.phat1) %in% paste("X",x,sep="")])

list2<-list()

for(i in 1: nrow(alg1.phat1)){
  list2[[i]]<-alg1.phat1[i,names(alg1.phat1) %in% paste("X",df[i,],sep="")]
}


#delete rows with all NA

list2.full<-list2[(which(sapply(list2, function(x) length(is.na(x) %in% "FALSE")>0)))]

phat1.df <- ldply(sapply(list2,data.frame), rbind)
colnames(phat1.df)
row.names(phat1.df)<-paste(rep("gene",NROW(phat1.df)),which(sapply(list2, function(x) length(is.na(x) %in% "FALSE")>0)))


detach(package:plyr)
#######sorting manually o get bicluster
#sort by genes with large numbr of samples (rows)

s.s.phat1.df<-rbind(apply(phat1.df,2,function(x) length(x[!is.na(x)])),phat1.df)
s.s.phat1.df<-s.s.phat1.df[,order(s.s.phat1.df[1,],decreasing = T)]

#sort by sample based on large number of genes (column)
s.g.phat1.df<-cbind(apply(s.s.phat1.df,1,function(x) length(x[!is.na(x)])),s.s.phat1.df)

s.g.phat1.df<-s.g.phat1.df[order(s.g.phat1.df[,1],decreasing = T),]

heatmap.2(as.matrix(s.g.phat1.df[-1,-1]),Colv = F,Rowv = F, distfun=dist,keysize=1,dendrogram="none",col =
            colorRampPalette(c("green","red"))(100),offsetCol=1,3,scale="none",trace="none",key.ylab="")


##new searching method

bicscore<-function(data)
{
  score<-mean(data,na.rm=T)
  score
}

mat<-data.matrix(s.g.phat1.df[-1,-1])
mat[is.na(mat)]<-0


   box<-function(data,i,m,alpha){
     bic<-list()
     score<-rep(NA,nrow(data))
     logr<-i-1
    logc<-i+1
    score[1]<-bicscore(data[1,c(logr:logc)])
    j=1
    while(score[j]>alpha){
      bic[[j]]<-data[c(1:j),c(c(logr-1):c(logc+m))]
      score[j+1]<-bicscore(bic[[j]])
      j<-j+1
      ifelse(logr>2, logr<-logr-1,logr<-2)
      ifelse(logc<ncol(data)-m, logc<-logc+m,logc<-ncol(data)-m)
     }
    return(bic)
   }
   
   #for luminal
   bicluster<-box(mat,11,20,.9)
   
   #for basal
   bicluster<-box(mat,16,6,.9)
  
##getting bicluster names
   bic.final<-bicluster[[10]]
   bic.final<-cbind(rownames(bic.final),bic.final)
   colnames(bic.final)[1]<-"genesID"
   
   rownames(data3)<-mydata.ba[-1,1]
   genes<-data.frame(paste("gene",c(1:44),sep=" "),rownames(data3))
   names(genes)<-c("genesID","genenames")
   bic.final.genes<-merge(genes,bic.final,by.x="genesID",by.y="genesID")
   
   bic.genes<-bic.final.genes$genenames
   samples<-t(data.frame(colnames(data3)))
   colnames(samples)<-paste("X",c(1:ncol(data3)),sep="")
   bic.samples<-samples[colnames(samples) %in% names(bic.final.genes)]
   
   write.xlsx(bic.genes,"H:\\research\\Bicluster\\algorithm development\\results\\real data\\TCGA&cell\\alg1\\basal\\new search algorithm\\bic final genes.xlsx",append = T,sheetName ="basal.genes" )
   write.xlsx(bic.samples,"H:\\research\\Bicluster\\algorithm development\\results\\real data\\TCGA&cell\\alg1\\basal\\new search algorithm\\bic final genes.xlsx",append = T,sheetName ="basal.samples" )
   

  
   
