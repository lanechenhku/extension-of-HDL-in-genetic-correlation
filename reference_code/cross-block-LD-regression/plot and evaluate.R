corrplot::corrplot(LDall[1:9,1:9])
corrplot::corrplot(Rsb[1:9,1:9])
result.1sim <- read.csv("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/result.1sim.csv")
result.100mean <- read.csv("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/result.100sim.csv")
q=c(1,seq(10,200,10))
result=list()
for(i in 1:length(q)){
  k=q[i]
result[[i]] <- read.csv(paste0("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/result.",k,"sim.csv"))}

plot.perform<-function(result.1sim){
y=result.1sim
a<-abs(y$dir.est)>mean(abs(y$dir.est))
b1<-y$lasso.est!=0
b2<-abs(y$lasso.est)>mean(abs(y$lasso.est))
c1<-y$margin.p<0.05/20
c2<-abs(y$margin.est)>mean(abs(y$margin.est))
d<-y$true.value!=0
z=data.frame(a,b1,b2,c1,c2,d) 
zz=apply(z,2,as.numeric)
corrplot::corrplot(t(zz),method='square',col = 'black')
zz1=data.frame(y$dir.est,y$lasso.est,y$margin.est,y$true.value)
print(cor(zz1))
corrplot::corrplot(t(zz1),method='circle')
}
plot.perform(result.100mean)
plot.perform(result.1sim)

perform<-function(result.1sim){
  y=result.1sim
  a<-abs(y$dir.est)>mean(abs(y$dir.est))
  b1<-y$lasso.est!=0
  b2<-abs(y$lasso.est)>mean(abs(y$lasso.est))
  c1<-y$margin.p<0.05/20
  c2<-abs(y$margin.est)>mean(abs(y$margin.est))
  d<-y$true.value!=0
  cal<-function(a){
    TP=sum((a==1)+(d==1)>1)
    TN=sum((a==0)+(d==0)>1)
   return(c(TP,TN)) 
  }
  z=data.frame(a,b1,b2,c1,c2) 
  zz=apply(z,2,as.numeric)
  out=apply(zz,2,cal) %>% data.frame()
  rownames(out)=c('TP','TN')
  #corrplot::corrplot(t(zz),method='square',col = 'black')
  zz1=data.frame(y$dir.est,y$lasso.est,y$margin.est,y$true.value)
  out2=cor(zz1)[4,1:3]
  #corrplot::corrplot(t(zz1),method='circle')
  return(list(out,out2))
}
cor_result1<-list()
test_result1<-list()

for(i in 1:length(q)){
  cor_result1[[i]]<-perform(result[[i]])[2]
  test_result1[[i]]<-perform(result[[i]])[1]
}

cor_result<-cor_result1 %>% data.frame() %>% t() %>% data.frame()
colnames(cor_result)<-c('Direct solution','Lasso','Marginal LR')
rownames(cor_result)<-1:length(q)
cor_result$K<-q

test_result<-test_result1 %>% data.frame() %>% t() %>% data.frame()
test_result$K=rep(q,each=5)
test_result$method=rep(c('Direct Sol.','Lasso','Lasso+T','LR+pval','LR+T'),21)

#tidy(test_result)
data=gather(cor_result,key ='method', value= 'correlation',-K)
data1=gather(test_result,key ='type', value= 'count',-K,-method)
ggplot(data=data)+geom_point(aes(x=K,y=correlation,color=method))+theme_bw()+
  ylab('correlation with true value')+
  geom_line(aes(x=K,y=correlation,color=method))

ggplot(data=data1,aes(x=K,y=count,fill=type))+geom_bar(position = 'stack',stat='identity',alpha=0.8)+theme_bw()+
  ylab('TP+TN')+facet_grid(.~method)+abline(h=4,col='black')


