big.rep=50
rhos=0.9
date='20220616'
result_rep=list()
for(kk in 1:big.rep){
  set.seed(kk)
  rep=200 #200 GWAS
  #repeat effect size generation
  h2=0.8 #herit
  p=sum(d) #no. of SNP
  mu=rep(0,p)
  h2m=h2/p
  d1=d[1]
  d2=d[2]
  Rsb=diag(1,p,p)
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
  
  
  data=matrix(0,nrow=rep*50*50,ncol=4)
  ii=1
  for(i in 1:rep){
    for(j in 1:50){
      for(k in 51:100){
        y=Zscore[i,j]*Zscore[i,k]
        x1=n*h2/p*LDscore[j,k]
        x2=LDall[j,k]
        g=n*h2/p*LDall[j,pos1]*LDall[k,pos2]
        data[ii,]=c(j,k,y,g)
        ii=ii+1
      }}}
  
  tmp=as.data.frame(data)
  #rownames(tmp)<-paste0('comb',1:nrow(tmp))
  colnames(tmp)<-c('j','k','y','g')
  data2=tmp
  #data2=as.data.frame(data2)
  
  data2$jk=factor(data2$j):factor(data2$k)
  data3=aggregate(data2[,3:4],by=list(data2$jk),mean)
  
  P1=n*h2m*LDscore+LDall
  P2.tmp=matrix(LDall[,pos1])%*%t(matrix(LDall[,pos2]))+matrix(LDall[,pos2])%*%t(matrix(LDall[,pos1]))
  P2=rhos*n*h2m*P2.tmp
  #sigma covariance matrix
  C1=P1+P2
  nll <- function(rhos.est){
    P2=rhos.est*n*h2m*P2.tmp
    C1=P1+P2
    loglik <-function(x){dmvnorm(x, sigma = C1, log = T)}
    ll=apply(Zscore,1,loglik)
    return(-sum(ll))
  }
  
  est <- stats4::mle(minuslog=nll, start=list(rhos.est=0))
  fit=est
  #summary(est)
  out=summary(fit)
  tmp1=out
  
  model=lm(y~g+1,data3)
  tmp2<-summary(model)
  
  #IRLS
  rhos.tmp=coefficients(model)[2]
  wt <- 1 / (tmp.vec+rhos.tmp*data3$g)
  wt2=wt
  #perform weighted least squares regression
  wls_model4 <- lm(y ~ g+1, data = data3, weights=wt)
  coeff=coefficients(wls_model4)
  # equation of the line : 
  eq2 = paste0("y = ", round(coeff[2],2), "*x ", round(coeff[1],2))
  tmp3=summary(wls_model4)
  
  result_out=c(cor(b[,pos1],b[,pos2]),cor(Zscore[,pos1],Zscore[,pos2]),tmp1@coef,tmp2$coefficients[2,1:2],tmp3$coefficients[2,1:2])
  result_rep[[kk]]=result_out
  print(kk)
}


tmp=as.data.frame(t(as.data.frame(result_rep)))
colnames(tmp)<-c('true beta emp.corr','naive est.','MLE.slope Est.','Std. Error','OLS(+1).slope Est.','Std. Error','WLS4(+1).slope Est.','Std. Error')
# tmp$`WLS4(+1).slope Est.`=2*tmp$`WLS4(+1).slope Est.`
# tmp$`OLS(+1).slope Est.`=2*tmp$`OLS(+1).slope Est.`
out_table=round(tmp,3)
out_table$theo.rhos=rhos
rownames(out_table)<-paste0('rep',1:nrow(out_table))
mean=apply(out_table,2,FUN=mean)
out_table=rbind(out_table,mean)
rownames(out_table)[nrow(out_table)]='Mean'

write.csv(out_table,paste0('out_table_compare_',date,rhos,'.csv'))
#MSE

tmp.out=out_table[-nrow(out_table),]
colnames(tmp.out)
tmp.out=(tmp.out[,c(1,2,3,5,7)]-tmp.out$theo.rhos)^2
mean=apply(tmp.out,2,FUN=mean)
tmp.out=rbind(tmp.out,mean,sqrt(mean),sqrt(mean)/abs(rhos)*100)
rownames(tmp.out)[(nrow(tmp.out)-2):nrow(tmp.out)]=c('MSE','RMSE','RMSE/|rhos|(%)')
write.csv(tmp.out,paste0('out_table_MSE_',date,rhos,'.csv'))


