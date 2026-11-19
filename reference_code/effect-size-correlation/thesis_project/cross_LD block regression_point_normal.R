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
library(DMwR) # Functions and data for "Data Mining with R"
library(MASS)

#sample size
set.seed(1)
n=10000
#rep
rep=50
#generate LD matrix
#two blocks
#d1=50;d2=50;d=900
d=c(4,5,50,50,50,50)
LD=list()
for(i in 1:length(d)){
  LD[[i]]=rcorrmatrix(d[i], alphad = 5)
  #LD[[i]]=diag(1,d[i], d[1])
}
#generate genotype
p=sum(d)
mu=rep(0,p)
#whole LD matrix
LDall=LD[[1]]
for(i in 2:length(d)){
  M1=LDall;M2=LD[[i]]
  LDall=rbind(cbind(M1,matrix(0,nrow=nrow(M1),ncol=ncol(M2))),cbind(matrix(0,nrow=nrow(M2),ncol=ncol(M1)),M2))
}

#genotype
X=mvrnorm(n , mu, LDall, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
apply(X,2,mean)
X=scale(X)
#average differece Er-r
#sum(cor(X)[1:d1,1:d1]-LD[[1]])/(d1*d1)


#effect size correlation
h2=0.8 #herit
h2m=h2/p
d1=d[1]
d2=d[2]
R0=rcorrmatrix(d1+d2, alphad = 5) #block one and block 2
#R0[1:d1,1:d1]=diag(1,nrow=d1,ncol=d1)
#R0[(d1+1):(d1+d2),(d1+1):(d1+d2)]=diag(1,nrow=d2,ncol=d2)
#Rsa=R0 #effect size correlation: full connection --not PD
Rsb=diag(1,p,p)
l=c(1:2)
Rs=R0[l,d1+l]
Rsb[l,d1+l]=Rs
Rsb[d1+l,l]=t(Rs) #a few connection---3*3=9 pairs

#new way
# Rsb=rcorrmatrix(p, alphad = 3)
Rsa=h2m*Rsb
#b=mvrnorm(n = 1, mu,Rsb, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

b=mvrnorm(n = rep, mu,Rsa, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
b[,c(3,4,7,8,9)]=0#point normal
apply(b^2,1,sum)
#bb=cor(b)
#mean(Rs)
#mean(bb[l,d1+l]-(Rs))
#phenotype


#crossprod<-function(k){
z12rep<-list()
for(k in 1:nrow(b)){
  beta=b[k,]
  Xb=X%*%beta
  sum(beta^2)
  var(Xb)
  var(X[,1]) #chec variance=1
  error <- rnorm(n, sd=sqrt(var(Xb)/h2*(1-h2)))
  Y=Xb+error #phenotype
  mean(Y)
  var(Y)
  var(Xb)/var(Y) #check h2
  #GWAS
  z<-c()
  for(i in 1:ncol(X)){
    a=summary(lm(Y~X[,i]))$coefficients
    z[i]=a[2,3]}
  chi2=z^2
  var(chi2)
  mean(chi2)
  #LD score regression----
  # lds=apply(LDall^2,2,sum)
  # mean(chi2)-mean(lds)*n*h2m+1
  # a<-lm(chi2-1~lds+0) %>% summary() %>% coefficients()
  # h2hat=a[1,1]/n*p
  # h2hat
  #cross-block LDSC----
  z1=z[1:d1]
  z2=z[(d1+1):(d1+d2)]
  
  z12<-c()
  ii=0
  for(i in 1:d1){
    for(j in 1:d2){
      ii=ii+1
      z12[ii]=z1[i]*z2[j]
    }
  }
  z12rep[[k]]<-z12
}

z12all=z12rep %>% data.frame() %>% t()
#z score for 100 reps
colnames(z12all)=paste0('col',1:ncol(z12all))
rownames(z12all)=paste0('row',1:nrow(z12all))


Rs12=Rsb[1:d1,(d1+1):(d1+d2)]
rhos12<-c()
reg<-list()
ii=0
#cross R matrix
for(i in 1:d1){
  for(j in 1:d2){
    ii=ii+1
    # z12[ii]=z1[i]*z2[j]
    rhos12[ii]=Rs12[i,j]
    crossr=c()
    for(k in 1:d1){
      for(l in 1:d2){
        crossr=c(crossr,LD[[1]][i,k]*LD[[2]][j,l])
      }}
    reg[[ii]]=crossr
  }
}
crossLD=as.data.frame(reg) 
colnames(crossLD)=paste0('col',1:ncol(crossLD))
rownames(crossLD)=paste0('row',1:nrow(crossLD))

#evaluate and estimate

Lassosol<-function(crossLD,z12){
  xx=as.matrix(crossLD) 
  y=z12 
  fit<-cv.glmnet(x=xx,y)
  plot(fit)
  summary(fit)
  fit$lambda.1se
  fit$glmnet.fit
  est=coef(fit,s='lambda.min') %>% as.matrix()
  #est2<-glmnet(x=xx,y=z12,lambda=40*fit$lambda.1se) %>% coef() %>% as.matrix()
  est=est[-1,1]
  #est2=est2[-1,1]
  co=n*h2m
  rhos12hat2=est/co
  return(rhos12hat2)}

co=n*h2m
fit.rs<-function(crossLD,z12){
  # direct solve
  dir.est=rep(NA,length(rhos12))
  if(nrow(crossLD)==ncol(crossLD)){
    rhoshat=solve(crossLD,z12)
    dir.est=rhoshat/co}
  
  #lasso solve
  lasso.est=Lassosol(crossLD,z12)
  
  #lm
  yy=scale(z12);xx=scale(crossLD)
  ols.est=rep(NA,length(rhos12))
  ols.p=rep(NA,length(rhos12))
  
  if(nrow(crossLD)>ncol(crossLD)){
    ols.fit<-lm(yy~xx+0) %>% summary()
    bbb=ols.fit$coefficients[,1]
    bbc=ols.fit$coefficients[,4]
    ols.est[1:length(bbb)]=bbb
    ols.p[1:length(bbb)]=bbc}
  #marginal lm
  mar.fit<-function(crossLD,z12){
    z<-c()
    p<-c()
    y=scale(z12);x=scale(crossLD)
    for(i in 1:ncol(crossLD)){
      a=summary(lm(y~x[,i]))$coefficients
      z[i]=a[2,1]
      p[i]=a[2,4]}
    return(list(z,p))}
  
  margin.est<-mar.fit(crossLD,z12)[[1]]
  margin.p<-mar.fit(crossLD,z12)[[2]]
  
  #true
  true.value=rhos12
  
  compare=data.frame(dir.est,lasso.est,ols.est,ols.p,margin.est,margin.p,true.value)
  return(compare)}

#in one simulation
results.1sim=fit.rs(crossLD,z12all[1,])
#in k simulation
k=50
results.mean=fit.rs(crossLD,apply(z12all[1:k,],2,mean)) #use mean of all simu


z12long=flatten(list(z12all[1:k,])) %>% unlist()
crossLDlong=crossLD
for(i in 2:nrow(z12all[1:k,])){crossLDlong=rbind(crossLDlong,crossLD)}

results.all=fit.rs(crossLDlong,z12long) #use mean of all simu

results.1sim
results.mean
results.all