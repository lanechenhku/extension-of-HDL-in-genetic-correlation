big.rep=50
rhos=0.9
date='20220819'

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
h2=0.8 #herit

mu=rep(0,p)
h2m=h2/p
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
h2p=rep(h2m,p)
h2p=w*h2p

# h2p*Rsb
# sqrt(h2p)%*%Rsb%*%sqrt(h2p)
Rsa=h2p*Rsb
# head(diag(Rsa))
# head(h2p)
b=mvrnorm(n = rep, mu,Rsa, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
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

#OLS
Z2=Zscore^2
wLD2=n*h2m*LDall^2
wLD2=data.frame(wLD2)
colnames(wLD2)=paste0('wLD',1:p)

data3=wLD2
data3$z2=Z2[1,]

for(i in 2:nrow(Z2)){
  
  dt.tmp=wLD2
  dt.tmp$z2=Z2[i,]
  data3=rbind(data3,dt.tmp)

}

data3$z2=data3$z2-1
#ols

model=lm(z2~.+0,data3)
fit.ols=summary(model)
coeff=coefficients(model)
slope=coeff
slope
dt=data.frame(w.true=w,w.mle=slope)
p <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "OLS estmates of weights")+ ylab( "true weights")
#回归线
#添加回归曲线
p + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)
mean(slope,na.rm=T)

plot(fit.ols$residuals)
plot(model)
plot(data3$z2)
plot(tmp$y_hat)

tmp=data.frame(y=data3$z2,y_hat=model$fitted.values)
#wls----

#define weights to use
slope
e.var=model$fitted.values+1
wt<-1/(2*e.var^2)

#perform weighted least squares regression
hist(wt)
#wt<-c(rep(10,10),rep(0.1,3980),rep(10,10))
wls_model1 <- lm(z2~.+0, data = data3, weights=wt)
fit.wls=summary(wls_model1)
coeff=coefficients(wls_model1)
coeff
slope=coeff
dt=data.frame(w.true=w,w.mle=slope)
p <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "WLS1 estmates of weights")+ ylab( "true weights")
#回归线
#添加回归曲线
p + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
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
