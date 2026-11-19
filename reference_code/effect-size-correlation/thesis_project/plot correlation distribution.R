a=seq(-0.5,0.9,0.1)
b=c(1,2,3,4,5,8,10,9,8,7,5,3,3,2,1)
d=c()
for(i in 1:length(a)){
  c=rep(a[i],b[i])
  d=c(d,c)
}
library(ggplot2)
dataset=data.frame(X=d)
ggplot(dataset, aes(x = X)) + xlim(c(-1,1))+geom_vline(xintercept=mean(d),col='blue')+
  geom_density()+geom_vline(xintercept=0,col='red')+xlab('GWAS effect size estimate correlation')
lines(density(d), col="blue", lwd=2) # add a density estimate with defaults
lines(density(d, adjust=2), lty="dotted", col="darkgreen", lwd=2) 


d1=rnorm(1000,0,0.5)
d2=rnorm(1000,0,0.2)
d3=rnorm(1000,0,0.05)
hist(d1)
dataset=data.frame(d1,d2,d3)
ggplot(dataset) + xlim(c(-1,1))+
  geom_density(aes(x = d1),col='red')+ geom_density(aes(x = d2),col='green')+ geom_density(aes(x = d3),col='blue')+geom_vline(xintercept=0,col='black')+xlab('GWAS effect size estimate correlation')
