#install.packages('clusterGeneration')
library(clusterGeneration)
library(devtools) # Tools to Make Developing R Packages Easier
#library(annotater) # Annotate Package Load Calls
#library(RcppArmadillo) # not installed on this machine
library(data.table) # Extension of `data.frame` # Extension of `data.frame`
library(Matrix) # Sparse and Dense Matrix Classes and Methods
library(tidyverse) # Easily Install and Load the 'Tidyverse'
# library(lsgl) # not installed on this machine
# library(lmridge) # not installed on this machine
library(data.table) # Extension of `data.frame` # Extension of `data.frame`
library(glmnet) # Lasso and Elastic-Net Regularized Generalized Linear Models
#library(KRIS) # not installed on this machine
library(monomvn) # Estimation for MVN and Student-t Data with Monotone Missingness
#library(DMwR) # Functions and data for "Data Mining with R"
library(MASS)
library(dplyr)
#library(ggpmisc)
library(corrplot)

#sample size
set.seed(1)
n=10000
#rep

#generate LD matrix----
#two blocks
d=c(50,50)
LD=list()
for(i in 1:length(d)){
  LD[[i]]=rcorrmatrix(d[i], alphad = 0.1)
  #LD[[i]]=diag(1,d[i], d[1])
}


#whole LD matrix
LDall=LD[[1]]
for(i in 2:length(d)){
  M1=LDall;M2=LD[[i]]
  LDall=rbind(cbind(M1,matrix(0,nrow=nrow(M1),ncol=ncol(M2))),cbind(matrix(0,nrow=nrow(M2),ncol=ncol(M1)),M2))
}

#LD score matrix
LDscore=t(LDall)%*%LDall
ldscore=diag(LDscore)
ldscore1=ldscore[1:50]
ldscore2=ldscore[51:100]
pos1=which(ldscore1==max(ldscore1))[1]
pos2=which(ldscore2==max(ldscore2))[1]
ldscore[pos1]
ldscore[pos2]

#effect size generation----
set.seed(1)
rep=200
h2=0.8 #herit
p=sum(d)
mu=rep(0,p)
h2m=h2/p
d1=d[1]
d2=d[2]
Rsb=diag(1,p,p)
rhos=-0.9
Rsb[pos1,pos2]=rhos
Rsb[pos2,pos1]=rhos
Rsa=h2m*Rsb
b=mvrnorm(n = rep, mu,Rsa, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
cor(b[,pos1],b[,pos2])

#zscore generation----
Zscore=list()
for(i in 1:nrow(b)){
mean=sqrt(n)*LDall%*%b[i,]
z=mvrnorm(n = 1, mean,LDall, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
Zscore[[i]]=z
}
tmp=t(as.data.frame(Zscore))
rownames(tmp)<-paste0('GWAS',1:rep)
colnames(tmp)<-paste0('SNP',1:p)
Zscore=tmp

cor(b[,pos1],b[,pos2])
cor(Zscore[,pos1],Zscore[,pos2])
#cor(Zscore[,1],Zscore[,2])

#regression----

data=list()
# ii=1
# for(i in 1:rep){
# for(j in 1:50){
#   for(k in 51:100){
# y=Zscore[i,j]*Zscore[i,k]
# x1=n*h2/p*LDscore[j,k]
# x2=LDall[j,k]
# g=2*n*h2/p*LDall[j,pos1]*LDall[k,pos2]
# data[[ii]]=c(j,k,y,g)
# ii=ii+1
#   }}}

data=matrix(0,nrow=rep*50*50,ncol=4)
ii=1
for(i in 1:rep){
  for(j in 1:50){
    for(k in 51:100){
      y=Zscore[i,j]*Zscore[i,k]
      x1=n*h2/p*LDscore[j,k]
      x2=LDall[j,k]
      g=2*n*h2/p*LDall[j,pos1]*LDall[k,pos2]
      data[ii,]=c(j,k,y,g)
      ii=ii+1
    }}}

tmp=as.data.frame(data)
#rownames(tmp)<-paste0('comb',1:nrow(tmp))
colnames(tmp)<-c('j','k','y','g')
data2=tmp
#data2=as.data.frame(data2)
out=lm(y~g+1,data2)
tmp<-summary(out)
result_out=c(cor(b[,pos1],b[,pos2]),cor(Zscore[,pos1],Zscore[,pos2]),tmp$coefficients[2,])

#regression 2

head(data2)
data2$jk=factor(data2$j):factor(data2$k)
data3=aggregate(data2[,3:4],by=list(data2$jk),mean)
out=lm(y~g+1,data3)
tmp<-summary(out)
result_out=c(cor(b[,pos1],b[,pos2]),cor(Zscore[,pos1],Zscore[,pos2]),tmp$coefficients[2,])
result_out
tmp
plot(data3$g,data3$y)

#plot----

LDscore[pos2,pos1]

data3=data2[1:2500,]
corrplot(LDall)
corrplot(LDscore/max(ldscore))
heatmap(LDscore,Rowv=NULL,Colv=NULL)
plot(data=data2,y,g)
plot(data2$g,data2$y)
plot(data3$g,data3$y)




