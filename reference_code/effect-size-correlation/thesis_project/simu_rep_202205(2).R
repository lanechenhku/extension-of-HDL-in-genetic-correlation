big.rep=50
rhos=-0.9
date='20220516'
result_rep=list()
for(kk in 1:big.rep){
  set.seed(kk)
  rep=200
  #repeat effect size generation
  h2=0.8 #herit
  p=sum(d)
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
        g=2*n*h2/p*LDall[j,pos1]*LDall[k,pos2]
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
  out=lm(y~g+0,data3)
  tmp<-summary(out)
  out2=lm(y~g+1,data3)
  tmp2<-summary(out2)
  result_out=c(cor(b[,pos1],b[,pos2]),cor(Zscore[,pos1],Zscore[,pos2]),tmp$coefficients[1,],tmp2$coefficients[2,])
  result_rep[[kk]]=result_out
}

tmp=as.data.frame(t(as.data.frame(result_rep)))
colnames(tmp)<-c('true beta emp.corr','naive est.','reg(+0).slope Est.','Std. Error','t value','Pr(>|t|)','reg(+1).slope Est.','Std. Error','t value','Pr(>|t|)')
tmp$`2*reg(+0).Est`=2*tmp$`reg(+0).slope Est.`
tmp$`2*reg(+1).Est`=2*tmp$`reg(+1).slope Est.`
#out_table=cbind(round(tmp[,1:4],3),tmp[,5:6])
out_table=cbind(round(tmp[,1:4],3),tmp[,5:6],round(tmp[,7:8],3),tmp[,9:10],round(tmp[,11:12],3))
out_table$theo.rhos=rhos
rownames(out_table)<-paste0('rep',1:nrow(out_table))
mean=apply(out_table,2,FUN=mean)
out_table=rbind(out_table,mean)
rownames(out_table)[nrow(out_table)]='Mean'
write.csv(out_table,paste0('out_table_compare_',date,rhos,'.csv'))
#MSE
tmp.out=out_table[-nrow(out_table),]
tmp.out=(tmp.out[,c(1,2,3,7,11,12)]-tmp.out$theo.rhos)^2
mean=apply(tmp.out,2,FUN=mean)
tmp.out=rbind(tmp.out,mean,sqrt(mean),sqrt(mean)/abs(rhos)*100)
rownames(tmp.out)[(nrow(tmp.out)-2):nrow(tmp.out)]=c('MSE','RMSE','RMSE/|rhos|(%)')
write.csv(tmp.out,paste0('out_table_MSE_',date,rhos,'.csv'))


