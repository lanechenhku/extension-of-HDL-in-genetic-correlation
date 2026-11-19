big.rep=50
rhos=0.9
date='20220912'

kk=1
set.seed(kk)
n=10000#sample size
#rep

#generate LD matrix----
#two blocks
len.block=10
d=c(len.block,len.block)
p=sum(d) #no. of SNP
LD=list()

#X <- cBern(50, rep(0.5,10), 0.5, type="exchange")
#cor(X)
set.seed(1)
down=0.5
up1=0.95
up2=0.85
ldnext1=c(seq(down,up1,length.out = d[1]/2),seq(up2,down,length.out = d[1]/2-1))
X <- cBern(200, rep(0.5,d[1]), ldnext1, type="DCP")
a=cor(X)
set.seed(2)
ldnext2=c(seq(down,up1,length.out = d[2]/2),seq(up2,down,length.out = d[2]/2-1))
Y <- cBern(200, rep(0.5,d[2]), ldnext2, type="DCP")
b=cor(Y)
#X <- cBern(5, c(0.4,0.5,0.6), c(0.2,0.3), type="1-dependent")
LD[[1]]=a
LD[[2]]=b

length(LD[[1]])
length(LD[[2]])

#whole LD matrix
LDall=LD[[1]]
for(i in 2:length(d)){
  M1=LDall;M2=LD[[i]]
  LDall=rbind(cbind(M1,matrix(0,nrow=nrow(M1),ncol=ncol(M2))),cbind(matrix(0,nrow=nrow(M2),ncol=ncol(M1)),M2))
}

LDscore=t(LDall)%*%LDall
ldscore=diag(LDscore)
ldscore1=ldscore[1:d[1]]
ldscore2=ldscore[(d[1]+1):p]
pos1=which(ldscore1==max(ldscore1))[1]
pos2=which(ldscore2==max(ldscore2))[1]+d[1]
pos1
pos2
ldscore[pos1]
ldscore[pos2]

#GWAS
set.seed(kk)
rep=200 #200 GWAS
#repeat effect size generation
h2=0.8 #mean herit

mu=rep(0,p)
d1=d[1]
d2=d[2]
Rsb=diag(1,p,p)
#add rhos
Rsb[pos1,pos2]=rhos
Rsb[pos2,pos1]=rhos


#f distribution
v1=100
v2=100
#p=100
set.seed(kk+1)
w=rf(p,v1,v2)
sc=100/(100-2)
w=w/sc
mean(w)
var(w)
hist(w,breaks=20)

h2all=rnorm(rep,h2,0.05)
hist(h2all,breaks=15)
max(h2all)<1

h2mall=h2all/p

b=list()
for(i in 1:rep){
h2p=rep(h2mall[i],p)
h2p=w*h2p

# h2p*Rsb
# sqrt(h2p)%*%Rsb%*%sqrt(h2p)
Rsa=h2p*Rsb
# head(diag(Rsa))
# head(h2p)
b.tmp=mvrnorm(n = 1, mu,Rsa, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
b[[i]]=b.tmp
}
b=as.data.frame(b)
b=t(b)
cor(b[,pos1],b[,pos2])
w[pos1]
w[pos2]

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

#generate data
Z2=Zscore^2
wLD2=n*h2mall[1]*LDall^2

wLD2=data.frame(wLD2)
colnames(wLD2)=paste0('wLD',1:p)

data3=wLD2
data3$z2=Z2[1,]

for(i in 2:nrow(Z2)){
  
  dt.tmp=data.frame(n*h2mall[i]*LDall^2)
  colnames(dt.tmp)=paste0('wLD',1:p)
  dt.tmp$z2=Z2[i,]
  data3=rbind(data3,dt.tmp)
  
}

data3$z2=data3$z2-1
#ols----

model=lm(z2~.+0,data3)
fit.ols=summary(model)
coeff.ols=coefficients(model)
slope=coeff.ols
slope
dt=data.frame(w.true=w,w.mle=slope)
p0 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "OLS estmates of w")+ ylab( "true w parameters")
#回归线
#添加回归曲线
p0 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)
mean(slope,na.rm=T)

tmp=data.frame(y=data3$z2,y_hat=model$fitted.values)
# plot(fit.ols$residuals)
# plot(model)
# plot(data3$z2)
# plot(tmp$y_hat)
# 

#wls----

#define weights to use
slope
e.var=2*(model$fitted.values+1)^2
#wLDall=data3[,-ncol(data3)]
#e.var0=2*(apply(wLDall,1,sum))^2
#e.var=e.var0
wt<-1/e.var


#wt <- 1 / (lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2)
#wt<-c(rep(10,10),rep(0.1,3980),rep(10,10)) ----check


#perform weighted least squares regression
hist(wt)

wls_model1 <- lm(z2~.+0, data = data3, weights=wt)
fit.wls=summary(wls_model1)
coeff.wls=coefficients(wls_model1)
coeff.wls
slope=coeff.wls
dt=data.frame(w.true=w,w.mle=slope)
p0 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "WLS estmates of w")+ ylab( "true w parameters")
#回归线
#添加回归曲线
p0 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)
mean(slope,na.rm=T)
#comp
fit.ols$coefficients
fit.wls$coefficients

summary(fit.ols$coefficients[,2])
summary(fit.wls$coefficients[,2])
summary(model)
summary(wls_model1)

#wls new----

#define weights to use
slope
#e.var=2*(model$fitted.values+1)^2
wLDall=data3[,-ncol(data3)]
e.var0=2*(apply(wLDall,1,sum))^2
#e.var=e.var0
wt<-1/e.var


#wt <- 1 / (lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2)
#wt<-c(rep(10,10),rep(0.1,3980),rep(10,10)) ----check


#perform weighted least squares regression
hist(wt)

wls_model1 <- lm(z2~.+0, data = data3, weights=wt)
fit.wls=summary(wls_model1)
coeff.wls=coefficients(wls_model1)
coeff.wls
slope=coeff.wls
dt=data.frame(w.true=w,w.mle=slope)
p0 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "WLS estmates of w")+ ylab( "true w parameters")
#回归线
#添加回归曲线
p0 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)
mean(slope,na.rm=T)
#comp
fit.ols$coefficients
fit.wls$coefficients

summary(fit.ols$coefficients[,2])
summary(fit.wls$coefficients[,2])
summary(model)
summary(wls_model1)

#IWLS----
test <- vector('list', 50)
test[[1]] <- list(coef0 = coef(model), res0 = model$residuals)

wls_func_lm <- function(obj){
  dt=data3
  weight <- obj[['res0']]^2
  lm_reg <- lm(z2~.+0, data = dt, weights=1/weight)
  dt$z2=lm_reg$residuals^2
  lm_reg_res <- lm(z2 ~ .+0, data = dt)
  return(list(coef0 = coef(lm_reg), res0 = lm_reg_res$fitted.values)) 
}

for (i in 1:10){
  test[[i+1]] <- wls_func_lm(test[[i]])
  dif <- sum(test[[i+1]][['coef0']] - test[[i]][['coef0']])^2
  cat(i, dif , '\n')
}

coeff.wls10=test[[10]][['coef0']]
slope=coeff.wls10
dt=data.frame(w.true=w,w.mle=slope)
p0 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "IRLS estmates of w")+ ylab( "true w parameters")
#回归线
#添加回归曲线
p0 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)
mean(slope,na.rm=T)

#GLS----
LDscore=t(LDall)%*%LDall
V0=matrix(0,nrow=length(Z2),ncol=length(Z2))
ngwas=rep
for(i in 1:ngwas){
  wLDsc2=n*h2mall[i]*LDscore
  a=1+nrow(wLDsc2)*(i-1)
  b=nrow(wLDsc2)*i
  V0[a:b,a:b]=wLDsc2}
V=2*V0^2
diag(V) %>% head()

#Z2=Zscore^2
wLDall=data3[,-ncol(data3)]
e.var0=2*(apply(wLDall,1,sum))^2
head(e.var0)

V.diag=diag(diag(V))
X<-as.matrix(wLDall)
V.inv<-solve(V)
V[1:p,1:p]
V.inv[1:p,1:p]

#gls
tmp=solve(t(X)%*%V.inv%*%X)
tmp2=tmp%*%t(X)%*%V.inv
beta.gls<-tmp2%*%data3$z2
coeff.gls=beta.gls
#se
res<-matrix(Z2)-X%*%beta

sig<-sum(res^2)/model$df
se.gls=sqrt(diag(solve(t(X)%*%V.inv%*%X))*sig)

# #wls
# 
# V.diag.inv=solve(V.diag)
# V.diag.inv[1:20,1:20]
# tmp=solve(t(X)%*%V.diag.inv%*%X)
# tmp2=tmp%*%t(X)%*%V.diag.inv
# beta.wls<-tmp2%*%data3$z2
# beta.wls
# 
# #ols
# tmp=solve(t(X)%*%X)
# tmp2=tmp%*%t(X)
# beta.ols<-tmp2%*%data3$z2
# beta.ols
# 
# beta.gls-beta.wls


slope=coeff.gls
dt=data.frame(w.true=w,w.mle=slope)
p0 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "GLS estmates of w")+ ylab( "true w parameters")
#回归线
#添加回归曲线
p0 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)

#check
mean(slope,na.rm=T)

# #gls another way
# library(nlme)
# tmp=diag(V)^(-(1/2))
# corV=diag(tmp)%*%V%*%t(diag(tmp))
# corV[1:3,1:3]
# gls.model<-gls(z2~.,correlation=corV,data=data3)
# summary(gls.model)
# intervals(gls.model)

#gls new----
LDscore=t(LDall)%*%LDall
V0=matrix(0,nrow=length(Z2),ncol=length(Z2))
ngwas=rep
for(i in 1:ngwas){
  wLDsc2=n*h2mall[i]*diag(w)%*%LDscore
  a=1+nrow(wLDsc2)*(i-1)
  b=nrow(wLDsc2)*i
  V0[a:b,a:b]=wLDsc2}
V=2*V0^2
#check V property

diag(V) %>% head()
hist(diag(V),breaks=50)
diagV.scaled=diag(V)/mean(diag(V))
hist(diagV.scaled,breaks=50)
log10diagV=log(diag(V),10)
hist(log10diagV,breaks=50)
summary(diag(V))
summary(V[row(V)!=col(V)])
sum(V[row(V)!=col(V)]==0)
sum(V[row(V)!=col(V)]!=0)

#Z2=Zscore^2
wLDall=data3[,-ncol(data3)]

V.diag=diag(diag(V))
X<-as.matrix(wLDall)
V.inv<-solve(V)
V[1:p,1:p]
V.inv[1:p,1:p]

#gls
tmp=solve(t(X)%*%V.inv%*%X)
tmp2=tmp%*%t(X)%*%V.inv
beta.gls<-tmp2%*%data3$z2
coeff.gls=beta.gls
#se
res<-matrix(Z2)-X%*%beta

sig<-sum(res^2)/model$df
se.gls=sqrt(diag(solve(t(X)%*%V.inv%*%X))*sig)

slope=coeff.gls
dt=data.frame(w.true=w,w.mle=slope)
p0 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "GLS estmates of w")+ ylab( "true w parameters")
#回归线
#添加回归曲线
p0 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)

#MLE for w----
dim(Zscore)
Zscore.tmp=Zscore %>% data.frame()
Zscore.tmp$h2mall=h2mall
nll <- function(wst1,wst2,wst3,wst4,wst5,wst6,wst7,wst8,wst9,wst10,wst11,wst12,wst13,wst14,wst15,wst16,wst17,wst18,wst19,wst20){
  wst=c(wst1,wst2,wst3,wst4,wst5,wst6,wst7,wst8,wst9,wst10,wst11,wst12,wst13,wst14,wst15,wst16,wst17,wst18,wst19,wst20)
  wLDscore=t(LDall)%*%diag(wst)%*%LDall
 
  loglik <-function(x){
    P1=n*x[p+1]*wLDscore+LDall
    C1=P1
    y=x[-(p+1)]
    f=dmvnorm(y, sigma = C1, log = T)
    return(f)}
  ll=apply(Zscore.tmp,1,loglik)
  return(-sum(ll))
}

#MLE evaluation

w.est <- stats4::mle(minuslog=nll, start=list(wst1=1,wst2=1,wst3=1,wst4=1,wst5=1,wst6=1,wst7=1,wst8=1,wst9=1,wst10=1,wst11=1,wst12=1,wst13=1,wst14=1,wst15=1,wst16=1,wst17=1,wst18=1,wst19=1,wst20=1))
fit=w.est
#summary(est)
w.mle=summary(fit)

wst=w.mle@coef

#plot w
coeff.mle=wst[,1]
slope=coeff.mle
dt=data.frame(w.true=w,w.mle=coeff.mle)
p0 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "MLE estimates of w")+ ylab( "true w parameters")
#回归线
#添加回归曲线
p0 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)

#GLS----

plot(coeff.ols,coeff.wls)
abline(a=0,b=1)

