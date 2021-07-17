
###################################
###algorithm 1, function
 #get source codes and truth tables (it is on git repo for V2_studies 'https://git.illumina.com/ClinicalGenomics/V2_studies.git')
    gitpath = "...\\Bi-EB\\"
ssource(file.path(gitpath,"Bi-EB\\src\\Bi_EB_alg1.R"))

###Define inputpath to data directory
inputpath<- ".../data"
output <- ".../result"
########################
##Load  TCGA data

##reading log(protein/gene) of TCGA

data<-read.csv(file.path(inputpath,"data\\tissue.log.ratio.overlap.all-no_NA.csv",stringsAsFactors=F))

##whole data
data2<-t(data)

data3<-data2[-c(1,2),-c(1,2)]
data3<-as.matrix(apply(data3,2,as.numeric))


##normalizing
datanorm<-apply(data3,1,function(x) (exp(x)-mean(exp(x)))/sd(exp(x)))

datanorm<-t(datanorm)

##by  subtype...add subtype to data
class<-read.csv(file.path(inputpath,"\\For_Correlation_TCGA_gene.protein&PT.Matched.csv"))[,c(2,6)]

data.wclass<-merge(class,data,by.x="Complete.TCGA.ID",by.y="pro..1..1.",all.y=T)


#TNBC
data.b<-data.wclass[data.wclass$Our.Classification %in% "Basal-like" ,]
data.b<-rbind(data.wclass[c(358,359),],data.b)

data.la<-data.wclass[data.wclass$Our.Classification %in% "Luminal A" ,]
data.la<-rbind(data.wclass[c(358,359),],data.la)

data2<-t(data.la)
data3<-data2[-c(1,2,3),-c(1,2)]
data3<-as.matrix(apply(data3,2,as.numeric))

#plot data
library(gplots)

heatmap.2(datanorm, distfun=dist,keysize=1,Rowv = FALSE,dendrogram="none",col =
            colorRampPalette(c("green","white","red"))(100),offsetCol=1,3,scale="none",trace="none",key.ylab="")


dev.off()
heatmap.2(data3, distfun=dist,keysize=1,dendrogram="none",col =
            colorRampPalette(c("green","white","red"))(100),offsetCol=1,3,scale="none",trace="none",key.ylab="")

##################
##Algorithm 1 implementation

#for normalized data
mu1<-12
mu2<-mean(datanorm,na.rm=TRUE)

#for regular log ratio data
 mu1<-max(data3,na.rm=TRUE)
 mu2<-min(data3,na.rm=TRUE)
####

var1<-var2<-sd(datanorm,na.rm=TRUE)^2

var1<-var2<-sd(data3,na.rm=TRUE)^2

##apply Bi_EB function
resultalg1<-alg1(data3,mu1,mu2,var1,var2,.5)

resultalg1<-alg1(datanorm,mu1,mu2,var1,var2,.5)


resultalg1$mu
resultalg1$var
resultalg1$k
resultalg1$p
resultalg1$beta1mean
resultalg1$beta1var
resultalg1$beta2mean
resultalg1$beta2var
plot(resultalg1$loglik)

alg1.phat1<-data.frame(resultalg1$phat1)

# following code limits the lowest and highest color to 5%, and 95% of your range, respectively
quantile.range <- quantile(alg1.phat1, probs = seq(0, 1, 0.01),na.rm=T)
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)


library(gplots)
heatmap.2(as.matrix(resultalg1$phat1), distfun=dist,keysize=1,Rowv = FALSE,dendrogram="none",breaks = palette.breaks,
          col =colorRampPalette(c("green","white","red"))(6-1),offsetCol=1,3,scale="none",trace="none",key.ylab="")

