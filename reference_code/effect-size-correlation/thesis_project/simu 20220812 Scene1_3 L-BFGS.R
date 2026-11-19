big.rep=50
rhos=0.9
date='20220812'
kk=200
# result_rep=list()
# for(kk in 1:big.rep){

set.seed(kk)
n=10000#sample size
#rep

#generate LD matrix----
#two blocks
len.block=20
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
v1=20
v2=20
#p=100
set.seed(kk+2)
w=rf(p,v1,v2)
sc=v2/(v2-2)
w=w/sc
mean(w)
var(w)
hist(w,breaks=20)

# heritability matrix 
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
#cor(Zscore[,1],Zscore[,2])

# P1=n*w*h2m*LDscore+LDall
# P2.tmp=matrix(LDall[,pos1])%*%t(matrix(LDall[,pos2]))+matrix(LDall[,pos2])%*%t(matrix(LDall[,pos1]))
# P2=rhos*n*w*h2m*P2.tmp

# C1=P1+P2

#MLE for w
nll <- function(wst1,wst2,wst3,wst4,wst5,wst6,wst7,wst8,wst9,wst10,wst11,wst12,wst13,wst14,wst15,wst16,wst17,wst18,wst19,wst20,wst21,wst22,wst23,wst24,wst25,wst26,wst27,wst28,wst29,wst30,wst31,wst32,wst33,wst34,wst35,wst36,wst37,wst38,wst39,wst40){
  wst=c(wst1,wst2,wst3,wst4,wst5,wst6,wst7,wst8,wst9,wst10,wst11,wst12,wst13,wst14,wst15,wst16,wst17,wst18,wst19,wst20,wst21,wst22,wst23,wst24,wst25,wst26,wst27,wst28,wst29,wst30,wst31,wst32,wst33,wst34,wst35,wst36,wst37,wst38,wst39,wst40)
  wLDscore=t(LDall)%*%diag(wst)%*%LDall
  P1=n*h2m*wLDscore+LDall
  # P2.tmp=matrix(LDall[,pos1])%*%t(matrix(LDall[,pos2]))+matrix(LDall[,pos2])%*%t(matrix(LDall[,pos1]))
  # P2=rhos.est*n*h2m*P2.tmp
  # #sigma covariance matrix
  C1=P1
  loglik <-function(x){dmvnorm(x, sigma = C1, log = T)}
  ll=apply(Zscore,1,loglik)
  return(-sum(ll))
}

start_time <- Sys.time()
#MLE
w.est <- stats4::mle(minuslog=nll, method='L-BFGS-B',start=list(wst1=1,wst2=1,wst3=1,wst4=1,wst5=1,wst6=1,wst7=1,wst8=1,wst9=1,wst10=1,wst11=1,wst12=1,wst13=1,wst14=1,wst15=1,wst16=1,wst17=1,wst18=1,wst19=1,wst20=1,wst21=1,wst22=1,wst23=1,wst24=1,wst25=1,wst26=1,wst27=1,wst28=1,wst29=1,wst30=1,wst31=1,wst32=1,wst33=1,wst34=1,wst35=1,wst36=1,wst37=1,wst38=1,wst39=1,wst40=1),lower=rep(0,p),upper=rep(p,p))
fit=w.est
#summary(est)
w.mle=summary(fit)

wst=w.mle@coef
#time gap
end_time <- Sys.time()
end_time - start_time

#plot w
dt=data.frame(w.true=w,w.mle=wst[,1])
p1 <- ggplot(dt, aes(w.mle, w.true)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "maximum likelihood estmates of weights")+ ylab( "true weights")
#回归线
#添加回归曲线
p1 + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)



