big.rep=10
n=10000#sample size
date='20221012'
rep=50 #No GWAS
perf.rep.ols=list()
perf.rep.wls=list()
perf.rep.gls=list()
h2=0.4 #mean herit
#generate LD matrix----
#two blocks
len.block=10
d=c(len.block,len.block)
p=sum(d) #no. of SNP
nsnp=p
LD=list()


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


for(kk in 1:big.rep){
#GWAS
set.seed(kk)
#repeat effect size generation


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
set.seed(2)
w=rf(p,v1,v2)
sc=100/(100-2)
w=w/sc
mean(w)
var(w)
hist(w,breaks=20)

set.seed(2)
h2all=rnorm(rep,h2,0.05)
hist(h2all,breaks=15)
max(h2all)<1

h2mall=h2all/p

b=list()
set.seed(kk)
for(i in 1:rep){
  h2p=rep(h2mall[i],p)
  h2p=w*h2p
  Rsa=h2p*Rsb
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
SSD=sum((slope-w)^2)
perf.ols=c(R2=cor(w,slope)^2,SSD=SSD,ASSD=SSD/length(w))

#wls----

#define weights to use
slope
e.var=2*(model$fitted.values+1)^2
wt<-1/e.var

wls_model1 <- lm(z2~.+0, data = data3, weights=wt)
fit.wls=summary(wls_model1)
coeff.wls=coefficients(wls_model1)
coeff.wls
slope=coeff.wls
SSD=sum((slope-w)^2)
perf.wls=c(R2=cor(w,slope)^2,SSD=SSD,ASSD=SSD/length(w))

# #GLS----
# wLDall=data3[,-ncol(data3)]
# X<-as.matrix(wLDall)
# V0=matrix(0,nrow=length(Z2),ncol=length(Z2))
# ngwas=rep
# for(i in 1:ngwas){
#   wLDsc2=n*h2mall[i]*LDscore
#   a=1+nrow(wLDsc2)*(i-1)
#   b=nrow(wLDsc2)*i
#   V0[a:b,a:b]=wLDsc2}
# V=2*V0^2
# V.inv<-solve(V)
# tmp=solve(t(X)%*%V.inv%*%X)
# tmp2=tmp%*%t(X)%*%V.inv
# beta.gls<-tmp2%*%data3$z2
# coeff.gls=beta.gls
# slope=coeff.gls
# SSD=sum((slope-w)^2)
# perf.gls=c(R2=cor(w,slope)^2,SSD=SSD,ASSD=SSD/length(w))

perf.rep.ols[[kk]]=perf.ols
# perf.rep.gls[[kk]]=perf.gls
perf.rep.wls[[kk]]=perf.wls
print(kk)
}

tmp=perf.rep.ols
perf.rep.ols1=as.data.frame(t(as.data.frame(tmp)))
apply(perf.rep.ols1,2,mean)
tmp=perf.rep.wls
perf.rep.wls1=as.data.frame(t(as.data.frame(tmp)))
apply(perf.rep.wls1,2,mean)
# tmp=perf.rep.gls
# perf.rep.gls1=as.data.frame(t(as.data.frame(tmp)))
# apply(perf.rep.gls1,2,mean)
