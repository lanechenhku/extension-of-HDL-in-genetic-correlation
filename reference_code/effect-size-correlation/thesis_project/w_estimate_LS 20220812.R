P=n*h2m*LDall^2
P.mean=apply(P,2,mean)
z2.mean=apply(Zscore^2,2,mean)
Z2=Zscore^2
dt.tmp=data.frame(P.mean,z2.mean,w)
dt.tmp=data.frame(P1=P[,1],Z2.1=Z2[,1],w)
fit.tmp=lm(Z2.1~P1+1,data=dt.tmp)
res_sq=(fit.tmp$residuals)^2
mean(res_sq)

dt.tmp=data.frame(P1=P[,2],Z2.1=Z2[,2],w)
fit.tmp=lm(Z2.1~P1+1,data=dt.tmp)
res_sq=(fit.tmp$residuals)^2
mean(res_sq)

dt.tmp=data.frame(P1=P[,3],Z2.1=Z2[,3],w)
fit.tmp=lm(Z2.1~P1+1,data=dt.tmp)
res_sq=(fit.tmp$residuals)^2
mean(res_sq)

dt.tmp=data.frame(P1=P[,4],Z2.1=Z2[,4],w)
fit.tmp=lm(Z2.1~P1+1,data=dt.tmp)
res_sq=(fit.tmp$residuals)^2
mean(res_sq)

plot(z2.mean,w)
plot(z2.mean,P.mean)

z2.mean.sd=scale(z2.mean)+1
#MLE
start_time <- Sys.time()
w.est <- stats4::mle(minuslog=nll, method='L-BFGS-B',start=list(wst1=1,wst2=1,wst3=1,wst4=1,wst5=1,wst6=1,wst7=1,wst8=1,wst9=1,wst10=1,wst11=1,wst12=1,wst13=1,wst14=1,wst15=1,wst16=1,wst17=1,wst18=1,wst19=1,wst20=1,wst21=1,wst22=1,wst23=1,wst24=1,wst25=1,wst26=1,wst27=1,wst28=1,wst29=1,wst30=1,wst31=1,wst32=1,wst33=1,wst34=1,wst35=1,wst36=1,wst37=1,wst38=1,wst39=1,wst40=1),lower=rep(0,p),upper=rep(p,p))
fit=w.est
#summary(est)
w.mle=summary(fit)

wst=w.mle@coef
#time gap
end_time <- Sys.time()
end_time - start_time
