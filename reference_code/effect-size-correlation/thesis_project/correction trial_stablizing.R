#install.packages(file.choose(),repos=NULL)
#install.packages("C:/Users/22367/AppData/Local/Temp/Rtmpotonfs/file360878d6388f/RandomFields_3.3.14.tar.gz",repos=NULL)
library("geoR")
library('tidyverse')
library('devtools')
library(MASS)
library(car)
# install_version("geoR", version = "1.7-4.1", repos = "http://cran.r-project.org")
# options(gsubfn.engine = "R")
# install_github('cran/RandomFields')


data3.tmp=data3

#transform----
hist(log(data3.tmp$g))
hist(scale(data3.tmp$g))
hist(log(data3.tmp$y))

#box cox
#BC transform
model=boxcox(data3.tmp$y~1)
data3.tmp$yy=-data3.tmp$y
min(data3.tmp$yy)

tmp=with(data3.tmp,boxCox(yy~g,data=data3.tmp, param=c("lambda"),
       lambda = seq(-3, 3, 1/10), gamma=seq(0.5, 100, 1/2), family="bcnPower"))
lambda=tmp$x[which(tmp$y==max(tmp$y))]

tmp2=with(data3.tmp,boxCox(yy~g,data=data3.tmp, param=c("gamma"),
                      lambda = lambda, gamma=seq(0.5, 50,1/2), family="bcnPower"))

gamma=tmp2$x[which(tmp2$y==max(tmp2$y))]
data3.tmp$yyy=bcnPower(data3.tmp$yy, lambda, jacobian.adjusted = FALSE, gamma=-min(data3.tmp$yy))
# lambda.tmp=1/3
# gamma.tmp=-min(data3.tmp$yy)

ggplot(data3.tmp, aes(yyy)) +            # ggplot2 histogram with default bins
  geom_histogram()+xlab('BC-transformed response')

out2.tmp=lm(yyy~g+1,data3.tmp)
summary(out2.tmp)
plot(out2.tmp)

res.tmp = resid(out2.tmp)
#We now plot the residual against the observed values of the variable waiting.
res.tmp=data.frame(res.tmp)
ggplot(res.tmp, aes(res.tmp)) +            # ggplot2 histogram with default bins
  geom_histogram()+xlab('residual')

#simple cutoff----

bbb=abs(data3$g)
#removing small
ddd=c(seq(0,0.99,by=0.01))
thres=c()
slope=c()
for(i in 1:length(ddd)){
  thres[i]=ddd[i]
  aa=ddd[i]
  out2.tmp=lm(y~g+1,data3.tmp[which(bbb>=quantile(bbb,aa)),])
  tmp=summary(out2.tmp)
  slope[i]=tmp$coefficients[2,1]}
cutoff.result=round(t(data.frame(thres,slope)),2) %>% as.data.frame()

plot(100*thres,slope,ylab='slope estimate',xlab='percentile threshold (% of small data removed)',ylim=c(min(min(slope),rhos)-0.1,max(max(slope),rhos)+0.1))
abline(h=-0.9,col='red')


#removing big

ddd=c(seq(0.01,0.99,by=0.01))
thres=c()
slope=c()
for(i in 1:length(ddd)){
  thres[i]=ddd[i]
  aa=ddd[i]
  out2.tmp=lm(y~g+1,data3.tmp[which(bbb<=quantile(bbb,aa)),])
  tmp=summary(out2.tmp)
  slope[i]=tmp$coefficients[2,1]}
cutoff.result=round(t(data.frame(thres,slope)),2) %>% as.data.frame()

#data.frame(thres,slope)
plot(100*thres,slope,ylab='slope estimate',xlab='percentile threshold (% of small data kept)',ylim=c(min(min(slope),rhos)-0.1,max(max(slope),rhos)+0.1))
abline(h=-0.9,col='red')

