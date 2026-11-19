#WLS/GLS/MLE data
library(mvtnorm)
library(stats4)
#data2$jk=factor(data2$j):factor(data2$k)
data3=aggregate(data2[,3:4],by=list(data2$jk),mean)
data4=aggregate(data2[,3],by=list(data2$jk),var)
colnames(data4)[2]<-'var_y'
data4$g=data3$g

#diagnostic outlier
data2.outlier=data2[which(data2$y<(-2000)),]
table(data2.outlier$i) 
table(data2.outlier$j) 
table(data2.outlier$k) 
hist(Zscore[26,],breaks = 30)

# by GWAS
slope.cal<-function(dt){
  out<-lm(y~g+1,data=dt)
  slope=coefficients(out)[2]
  return(slope)
}
d2.tmp=split( data2 , f = data2$i )
slope.by.gwas=lapply(d2.tmp,slope.cal)
slope.tmp=unlist(slope.by.gwas)
mean(slope.tmp)
hist(slope.tmp,breaks = 30,xlab='per-GWAS slope estimate',main=paste0("Hist: mean=",round(mean(slope.tmp),2)))
abline(v=mean(slope.tmp),col='red')
#three plot comparison
par(mfrow=c(2,2))
plot(data2$g,data2$y,xlab='predictor X_kl',ylab="raw response Z_kZ_l")
abline(h=0,col='red')
plot(data3$g,data3$y,xlab='predictor X_kl',ylab="GWAS-wide average(Z_kZ_l)")
abline(lm(y ~ g+1, data = data3), col = "blue")
abline(lm(y ~ g+0, data = data3), col = "red")
abline(a=0,b=rhos,col='green')
plot(data4$g,data4$var_y,xlab='predictor X_kl',ylab="GWAS-wide variance(Z_kZ_l)")
plot(data3$y,data4$var_y,xlab='GWAS-wide average(Z_kZ_k)',ylab="GWAS-wide variance(Z_kZ_k)")

par(mfrow=c(1,1))
#plot(data2$g,data2$y,xlab='predictor X_kl',ylab="raw response Z_kZ_l")

#mean undo the weight
plot(data3$g,data3$y/(n*h2m),main='Average Y after undoing the Nh2/m weight')
abline(v=0,h=-0.9,col='red')
#variance calculation----

#data3: y & x
#data4: vary & x
tmp.vec=c()
#data3$Group.1
k=0
for(i in 1:length(ldscore1)){
  for(j in 1:length(ldscore2)){
    k=k+1
    tmp.vec[k]=(n*h2m)^2*ldscore1[i]*ldscore2[j]+n*h2m*ldscore1[i]+n*h2m*ldscore2[j]+1
  }}
data4$True_var_y=tmp.vec+rhos*data3$g

lm.fit<-lm(var_y~True_var_y,data4)
summary(lm.fit)

p <- ggplot(data4, aes(x=True_var_y, y=var_y)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab('Theorectical variance(Z_kZ_k)')+ylab("GWAS-wide variance(Z_kZ_k)")
#回归线
#添加回归曲线
p1<-p + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)

p1+geom_abline(slope=1,intercept = 0,col='red')

plot(data4$True_var_y,data4$var_y,xlab='Theorectical variance(Z_kZ_k)',ylab="GWAS-wide variance(Z_kZ_k)")
abline(lm.fit,col='red')
abline(a=0,b=1,col='blue')

#wls 1----
#define weights to use
model=lm(y~g+1,data3)
wt <- 1 / (lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2)
#perform weighted least squares regression
wls_model1 <- lm(y ~ g+1, data = data3, weights=wt)
coeff=coefficients(wls_model1)
# equation of the line : 
eq = paste0("y = ", round(coeff[2],2), "*x ", round(coeff[1],2))

par(mfrow=c(2,2))
plot(data3$g,data3$y,xlab='X: predictor (weighted) Nh2/m*r_ki*r_lj',ylab="Y: GWAS-wide average(Z_kZ_l)",main=paste0('red: ',eq))
abline(lm(y ~ g+1, data = data3), col = "blue")
abline(wls_model1, col = "red")
abline(a=0,b=rhos,col='green')
plot(data3$g,wt,xlab='X: predictor (weighted) Nh2/m*r_ki*r_lj',ylab='weight')
plot(data4$var_y,1/wt,xlab='GWAS-wide empirical variance',ylab='1/weight')
plot(data4$True_var_y,1/wt,xlab='theoretical variance',ylab='1/weight')

#view summary of model
summary(wls_model1)
summary(model)

#WLS 2: theoretical var

#define weights to use
model=lm(y~g+1,data3)
wt <- 1 / data4$True_var_y
#perform weighted least squares regression
wls_model2 <- lm(y ~ g+1, data = data3, weights=wt)
coeff=coefficients(wls_model2)
# equation of the line : 
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))

par(mfrow=c(2,2))
plot(data3$g,data3$y,xlab='X: predictor (weighted) Nh2/m*r_ki*r_lj',ylab="Y: GWAS-wide average(Z_kZ_l)",main=paste0('red: ',eq))
abline(lm(y ~ g+1, data = data3), col = "blue")
abline(wls_model2, col = "red")
abline(a=0,b=rhos,col='green')
plot(data3$g,wt,xlab='X: predictor (weighted) Nh2/m*r_ki*r_lj',ylab='weight')
plot(data4$var_y,1/wt,xlab='GWAS-wide empirical variance',ylab='1/weight')
plot(data4$True_var_y,1/wt,xlab='theoretical variance',ylab='1/weight')

#view summary of model
summary(wls_model2)

##WLS 3&4

#define weights to use
model=lm(y~g+1,data3)
wt <- 1 / data4$var_y
wt1=wt
#perform weighted least squares regression
wls_model3 <- lm(y ~ g+1, data = data3, weights=wt)
coeff=coefficients(wls_model3)
# equation of the line : 
eq1 = paste0("y = ", round(coeff[2],2), "*x ", round(coeff[1],2))
summary(wls_model3)

#IRLS
rhos.tmp=coefficients(model)[2]
wt <- 1 / (tmp.vec+rhos.tmp*data3$g)
wt2=wt
#perform weighted least squares regression
wls_model4 <- lm(y ~ g+1, data = data3, weights=wt)
coeff=coefficients(wls_model4)
# equation of the line : 
eq2 = paste0("y = ", round(coeff[2],2), "*x ", round(coeff[1],2))
summary(wls_model4)

par(mfrow=c(1,1))
nf <- layout( matrix(c(1,1,2,3), nrow=2, byrow=TRUE) )
plot(data3$g,data3$y,xlab='X: predictor Nh2/m*r_ki*r_lj',ylab="Y: GWAS-wide average(Z_kZ_l)",main=paste0('red: WLS1 ',eq1,'\n','~ WLS2 ',eq2))
abline(wls_model3, col = "blue")
abline(wls_model4, col = "red")
abline(a=0,b=rhos,col='green')
#plot(data3$g,wt,xlab='X: predictor (weighted) Nh2/m*r_ki*r_lj',ylab='weight')
plot(data4$True_var_y,1/wt1,xlab='theoretical variance',ylab='1/weight WLS1')
plot(data4$True_var_y,1/wt2,xlab='theoretical variance',ylab='1/weight WLS2')

#view summary of model
summary(wls_model2)

#MLE likelihood

#negative log likelihood

P1=n*h2m*LDscore+LDall
pos1
pos2
P2.tmp=matrix(LDall[,pos1])%*%t(matrix(LDall[,pos2]))+matrix(LDall[,pos2])%*%t(matrix(LDall[,pos1]))
P2.tmp[pos1,pos2]
P2=rhos*n*h2m*P2.tmp
#sigma covariance matrix
C1=P1+P2
#trial
x <- rmvnorm(1, sigma = C1)
loglik <- dmvnorm(x, sigma = C1, log = T) 

# manual cal
#y <- t(x)
# deter2pi <- function(mat){
#   determinant(2*pi*mat,log=FALSE)[[1]][1]
# }
# ll <- (-0.5*t(y)%*%solve(C1,y)-0.5*log(deter2pi(C1)))[1]
# loglik/ll

nll <- function(rhos.est){
  P2=rhos.est*n*h2m*P2.tmp
  C1=P1+P2
  loglik <-function(x){dmvnorm(x, sigma = C1, log = T)}
  ll=apply(Zscore,1,loglik)
  return(-sum(ll))
}
#plot
range=seq(-0.95,-0.7,0.005)
ll.tmp=c()
for(i in 1:length(range)){ll.tmp[i]=nll(range[i])}
plot(range,-ll.tmp,xlab='rhos_ij est.',ylab='LogLik')
abline(v=rhos,col='red')
abline(v=range[which(ll.tmp==min(ll.tmp))],col='blue')

est <- stats4::mle(minuslog=nll, start=list(rhos.est=0))
fit=est
#summary(est)
summary(fit)
logLik(fit)
vcov(fit)

plot(profile(fit), absVal = FALSE)
confint(fit)


