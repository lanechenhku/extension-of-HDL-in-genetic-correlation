library(ggplot2) #加载ggplot2包
library(dplyr) #加载dplyr包
library(ggpmisc) #加载ggpmisc包
library(car)
#hist
ggplot(data3, aes(g)) +            # ggplot2 histogram with default bins
  geom_histogram()+xlab('regressor')
ggplot(data3, aes(y)) +            # ggplot2 histogram with default bins
  geom_histogram()+xlab('response')

#plot by gwas
d.tmp=data2[which(data2$i<=10),]
d.tmp$`GWAS_No.`=as.factor(d.tmp$i)
ggplot(d.tmp, aes(g, y,color=GWAS_No.)) +
geom_point(size = 3, alpha = 0.6)+xlab( "regressor_jk")+theme_bw()+ ylab( "(Z_jZ_k)_i")

#residual
#out2=lm(y~g+1,data3)
# tmp.2<-summary(out2)
# tmp.2
out2
res = resid(out2)
#We now plot the residual against the observed values of the variable waiting.

res.tmp=data.frame(res)
ggplot(res.tmp, aes(res)) +            # ggplot2 histogram with default bins
  geom_histogram()+xlab('residual')

plot(out2)

plot(data3$g, res, ylab="Residuals", xlab="Regressor",  main="Residual Plot") 
abline(0, 0)  

plot(out2,1)
mean(res)
plot(out2, 3)

plot(1:length(res),res)
#scatter plot
plot(data3$g,data3$y,xlab = "regressor_jk", ylab = "EZjZk",
     pch = 15, frame = FALSE)
abline(lm(y ~ g+1, data = data3), col = "blue")
abline(lm(y ~ g+0, data = data3), col = "red")
abline(a=0,b=rhos,col='green')

#regression line

p <- ggplot(data3, aes(g, y)) +
  geom_point(color = "grey50",size = 3, alpha = 0.6)+xlab( "regressor_jk")+ ylab( "E(Z_jZ_k)")
#回归线
#添加回归曲线
p + stat_smooth(color = "skyblue", fill = "skyblue", method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 5, #公式字体大小
    label.x = 0.1,  #位置 ，0-1之间的比例
    label.y = 0.95)

#tests

durbinWatsonTest(out2)
ncvTest(out2)

#collinearity
head(data3)
