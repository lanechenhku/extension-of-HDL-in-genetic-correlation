#new correlation method
library(CorBin)
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

#X <- cBern(50, rep(0.5,10), 0.5, type="exchange")
#cor(X)
set.seed(1)
down=0.5
up1=0.95
up2=0.85
ldnext1=c(seq(down,up1,length.out = d[1]/2),seq(up2,down,length.out = d[1]/2-1))
X <- cBern(200, rep(0.5,50), ldnext1, type="DCP")
a=cor(X)
set.seed(2)
ldnext2=c(seq(down,up1,length.out = d[2]/2),seq(up2,down,length.out = d[2]/2-1))
Y <- cBern(200, rep(0.5,50), ldnext2, type="DCP")
b=cor(Y)
#X <- cBern(5, c(0.4,0.5,0.6), c(0.2,0.3), type="1-dependent")
LD[[1]]=a
LD[[2]]=b

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
pos2=which(ldscore2==max(ldscore2))[1]+50
pos1
pos2
ldscore[pos1]
ldscore[pos2]

#plot
corrplot(LDall)
corrplot(LDscore/max(ldscore))

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

rhos
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

data=matrix(0,nrow=rep*50*50,ncol=5)
ii=1
for(i in 1:rep){
  for(j in 1:50){
    for(k in 51:100){
      y=Zscore[i,j]*Zscore[i,k]
      x1=n*h2/p*LDscore[j,k]
      x2=LDall[j,k]
      g=n*h2/p*LDall[j,pos1]*LDall[k,pos2]
      data[ii,]=c(j,k,y,g,i)
      ii=ii+1
    }}}

tmp=as.data.frame(data)
#rownames(tmp)<-paste0('comb',1:nrow(tmp))
colnames(tmp)<-c('j','k','y','g','i')
data2=tmp
#data2=as.data.frame(data2)
head(data2)
data2$jk=factor(data2$j):factor(data2$k)

#OLS 1
out=lm(y~g+0,data2)
tmp<-summary(out)
tmp
#regression 2


data3=aggregate(data2[,3:4],by=list(data2$jk),mean)
out=lm(y~g+0,data3)
tmp<-summary(out)
tmp
out2=lm(y~g+1,data3)
tmp.2<-summary(out2)
tmp.2
result_out=c(rhos,cor(b[,pos1],b[,pos2]),cor(Zscore[,pos1],Zscore[,pos2]),tmp$coefficients[1,1:2],tmp.2$coefficients[2,1:2])
tmp2=data.frame(result_out) %>% t() %>% as.data.frame()
colnames(tmp2)<-c('rhos','true beta emp.corr','naive est.','reg(+0).slope Est.','Std. Error','reg(+1).slope Est.','Std. Error')
tmp2
plot(data3$g,data3$y,xlab = "regressor_jk", ylab = "EZjZk",
pch = 15, frame = FALSE)
abline(lm(y ~ g+1, data = data3), col = "blue")
abline(lm(y ~ g+0, data = data3), col = "red")

#plot----

LDscore[pos2,pos1]

corrplot(LDall)
corrplot(LDscore/max(ldscore))
# heatmap(LDscore,Rowv=NULL,Colv=NULL)
# plot(data=data2,y,g)
# plot(data2$g,data2$y)
# plot(data3$g,data3$y)




