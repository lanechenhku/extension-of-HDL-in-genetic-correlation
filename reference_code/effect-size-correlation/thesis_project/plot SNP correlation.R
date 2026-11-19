#install.packages('tidyverse')
library(tidyverse)
SNP17 <- read.csv("C:/Users/admin/Desktop/effect size correlation/PheWAS_associations SNP17.csv", header=FALSE)
SNP12<- read.csv("C:/Users/admin/Desktop/effect size correlation/PheWAS_associations SNP12.csv", header=FALSE)
SNP17<-SNP17[-1,]
SNP12<-SNP12[-1,]
d1=SNP17[,c(6,8:13)]
name=c('group','p','beta','sd','n','dataset','pheno')
colnames(d1)=name
d2=SNP12[,c(6,8:13)]
colnames(d2)=name
d=left_join(d1,d2,by='pheno')
dd=na.omit(d)
dd$z.x=as.numeric(dd$beta.x)/as.numeric(dd$sd.x)
dd$z.y=as.numeric(dd$beta.y)/as.numeric(dd$sd.y)
par(mfrow=c(1,3))
a=cor(as.numeric(dd$beta.x),as.numeric(dd$beta.y))
plot(as.numeric(dd$beta.x),as.numeric(dd$beta.y),xlab='rs78378222',ylab='rs3730556',main=paste0('beta: ','r=',round(a,3)))
a=cor(as.numeric(dd$beta.x)/as.numeric(dd$sd.x),as.numeric(dd$beta.y)/as.numeric(dd$sd.y))
plot(as.numeric(dd$beta.x)/as.numeric(dd$sd.x),as.numeric(dd$beta.y)/as.numeric(dd$sd.y),xlab='rs78378222',ylab='rs3730556',main=paste0('Z score: ','r=',round(a,3)))
a=cor(-log(as.numeric(dd$p.x)),-log(as.numeric(dd$p.y)))
plot(-log(as.numeric(dd$p.x),10),-log(as.numeric(dd$p.y),10),xlab='rs78378222',ylab='rs3730556',main=paste0('-log(p-val): ','r=',round(a,3)))

dd$dataset=dd$dataset.x
dd$group=dd$group.x
dd$dataset[grep('GWAS_UKB',dd$dataset.x)]='GWAS_UKB'
dd$dataset[grep('GWAS_CKDGenConsortium',dd$dataset.x)]='GWAS_CKDGenConsortium'
dd$dataset[grep('GWAS_GEFOS',dd$dataset.x)]='GWAS_GEFOS'
dd$dataset[grep('GWAS_GODARTS',dd$dataset.x)]='GWAS_GODARTS'
dd$dataset[grep('GWAS_HOORN',dd$dataset.x)]='GWAS_HOORN'
dd$dataset[grep('TOPMED',dd$dataset.x)]='TOPMED'

library(ggpubr)
library(ggpmisc)
#install.packages('ggpmisc')
ggplot(data=dd,aes(x=z.x,y=z.y,col=dataset))+geom_point()+xlab('rs78378222')+ylab('rs3730556')+title('Z score')+
  geom_smooth(method = "lm",
              se = F,
              color = "red")
  # stat_cor(data = dd,
  #          method = "spearman")

ggplot(data=dd,aes(x=z.x,y=z.y,col=group))+geom_point()+xlab('rs78378222')+ylab('rs3730556')+title('Z score')+
geom_line(stat='smooth',method = 'lm', formula = y ~ x,se=F,alpha=0.5,size=0.5) + #绘制回归直线
  stat_cor(data = dd, method = "pearson")
           #stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + #添加回归方程和调整R方
  # stat_fit_tb(tb.type = 'fit.anova') + #添加方差分析表
  #theme_classic() #设置主题为classic

aaa=function(x){cor(x$z.x,x$z.y)}
aaa(dd)
aggregate(dd[,c('z.x','z.y')],list(factor(dd$dataset)),aaa)
summary(factor(dd$group))
