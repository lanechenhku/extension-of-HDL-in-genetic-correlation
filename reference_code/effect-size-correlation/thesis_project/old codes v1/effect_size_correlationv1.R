#20201009#
library(tidyverse)
library(data.table)
GWAS.altas.genes <- fread("/home/unfated/realdata/202010gene/magma.P.r3.txt",stringsAsFactors = F)
#GWAS.altas.genes <- fread("C:/Users/22367/Desktop/fine mapping real data/magma.P.r3.txt",stringsAsFactors = F)
GWAS.altas.genes<-GWAS.altas.genes %>% as.data.frame()
a=ncol(GWAS.altas.genes)
colnames(GWAS.altas.genes)[2:a]<-paste0('GWAS',1:(a-1))
colnames(GWAS.altas.genes)[1:5]
colnames(GWAS.altas.genes)[1]<-'gene'
GWAS.altas.genes.all=GWAS.altas.genes
GWAS.altas.genes=GWAS.altas.genes[-1,-1]
nrow(GWAS.altas.genes)
GWAS.altas.genes.all[,1][1:100]
#GWAS.altas.genes=GWAS.altas.genes %>% t()
#col.names(GWAS.altas.genes)=GWAS.altas.genes.all$gene[]
# missing.gene=apply(GWAS.altas.genes,1,function(x)sum(is.na(x)))
# missing.gwas=apply(GWAS.altas.genes,2,function(x)sum(is.na(x)))
# head(missing.gene)

#----
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

## load the library 
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
## list the contents that are loaded into memory 
ls('package:TxDb.Hsapiens.UCSC.hg19.knownGene') 
## show the db object that is loaded by calling it's name 
TxDb.Hsapiens.UCSC.hg19.knownGene
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

#write.table(GWAS.altas.genes.all$gene[-1],'genename.txt',row.names=F,col.names=F,quote = F)

#map gene name----
library(org.Hs.eg.db)
ls('package:org.Hs.eg.db') 
org.Hs.egENSEMBL
## Bimap interface:
genename=GWAS.altas.genes.all$gene[-1]
x <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID 
mapped_genes <- mappedkeys(x) 
# Convert to a list 
ENS.name<- as.list(x[mapped_genes]) 
#if(length(xx) > 0) 
#{ # Get the Ensembl gene IDs for the first five genes 
 # xx[1:5] 
  # Get the first one 
  #xx[[1]] } 
ENS.name[1:5]

EG.name <- as.list(org.Hs.egENSEMBL2EG) 

EG.name[1:5]

EG.name['ENSG00000175899']
genename[1:5]
EG.name[genename[1:5]]
genename.EG<-EG.name[genename] 
head(genename.EG)
length(genename.EG)

genename.EG[[3]]

listtovec <- function(genename.EG) {
  genename.EG.ID<-NA
  for(i in 1:length(genename.EG)){
    a=genename.EG[[i]]
    if(length(a)==0){a=NA}
    genename.EG.ID[i]=a
  }
  return(genename.EG.ID)
}
genename.EG.ID=listtovec(genename.EG)
genename.EG.ID[1:100]

## Bimap interface: chromosome
x <- org.Hs.egCHR # Get the entrez gene identifiers that are mapped to a chromosome 
mapped_genes <- mappedkeys(x) 
# Convert to a list 
xx <- as.list(x[mapped_genes]) 
xx[1:5]
gene.chr=xx[genename.EG.ID]
gene.chr.ID=listtovec(gene.chr)
gene.chr.ID[20000:20100]
genename[1:5]
genename.EG[1:5]
genename.EG.ID[1:5]
gene.chr.ID[1:5]
gene.chr[1:5]
GWAS.altas.genes$chr=gene.chr.ID
rownames(GWAS.altas.genes)<-genename

#effect size correlation dist----
genes.chr1<-GWAS.altas.genes %>% subset(.,chr==1)
genes.chr1<-genes.chr1 %>% as.matrix() %>% t()
genes.chr1.1<-genes.chr1 

genes.chr1= data.frame(genes.chr1.1)
genes.chr2<-GWAS.altas.genes %>% subset(.,chr==2)
genes.chr2<-genes.chr2 %>% as.matrix() %>% t() %>%data.frame()

genes.chr1[1:3,1:3]
genes.chr2[1:3,1:3]
#genes.chr1 %>% as.data.frame.table() %>% 
#fwrite(.,'/home/unfated/realdata/202010gene/genes.chr1.txt',row.names=T,col.names=T)

write.table(genes.chr1,'/home/unfated/realdata/202010gene/genes.chr1.txt',row.names=T,col.names=T,quote=F)
#genes.chr2 %>% as.data.frame.table() %>% 
#fwrite(.,'/home/unfated/realdata/202010gene/genes.chr2.txt',row.names=T,col.names=T)
write.table(genes.chr2,'/home/unfated/realdata/202010gene/genes.chr2.txt',row.names=T,col.names=T,quote=F)


genes.chr1<-apply(genes.chr1,2,function(x) as.numeric(x))
genes.chr2<-apply(genes.chr2,2,function(x) as.numeric(x))

# new 20201112
library(tidyverse)

genes.chr1=read.table('/home/unfated/realdata/202010gene/genes.chr1.txt')
head(genes.chr1)
genes.chr2=read.table('/home/unfated/realdata/202010gene/genes.chr2.txt')
ncol(genes.chr1)
ncol(genes.chr2)
#k=10
library(combinat, help, pos = 2, lib.loc = NULL)

count<-function(k) {
 # k=10
s1=sample(1:ncol(genes.chr1),k)
s1
s2=sample(1:ncol(genes.chr2),k)
s2
#-log((genes.chr1[,2]),10)%>%head()
#head(genes.chr1[,2])
#log(10,10)
genes.chr1.1<-apply(genes.chr1[,s1],2,function(x) -log(as.numeric(x),10))
genes.chr2.2<-apply(genes.chr2[,s2],2,function(x) -log(as.numeric(x),10))
genes.chr1.1[1:3,1:3]


corr.chr1.chr2=cor(genes.chr1.1,genes.chr2.2,use='pairwise.complete')
corr.chr1.chr2%>%nrow()
corr.chr1.chr2%>%ncol()
corr12=unlist(corr.chr1.chr2)
length(corr12)
head(corr12)
corr12[100]

#sum(corr12>0.01)
breakpoint=c(-1,-0.5,-0.25,-0.1,-0.05,-0.025,-0.01,-0.005,-0.001,0,0.001,0.005,0.01,0.025,0.05,0.1,0.25,0.5,1)
freq=cut(corr12,breaks=breakpoint)%>%table(.)
print('freq')
print(freq)
proportion=freq/length(corr12)
print('proportion')
print(proportion)
b=list(freq,proportion)
return(b)}

kk=c(10,20,50,100,200,300,400,500)

for(k in 1:10){
freq.n=list()
prop.n=list()
for(i in 1:length(kk)){
freq.n[[i]]=count(kk[i])[[1]]
prop.n[[i]]=count(kk[i])[[2]]
}


freq=freq.n%>%data.frame()
prop=prop.n%>%data.frame()
prop
freq.final=freq[,c(1,2*c(1:length(kk)))]
colnames(freq.final)=c('bin',kk)
prop.final=prop[,c(1,2*c(1:length(kk)))]
colnames(prop.final)=c('bin',kk)
write.csv(freq.final,paste0('/home/unfated/realdata/202010gene/freq.logcor20201113.real',k,'.csv'),row.names=F)
write.csv(prop.final,paste0('/home/unfated/realdata/202010gene/prop.logcor20201113.real',k,'.csv'),row.names=F)
print(k)
}



#sim----
n.size=ceiling(sum(!is.na(genes.chr1))/ncol(genes.chr1))

count2<-function(genes.chr1.1,genes.chr2.2){
  corr.chr1.chr2=cor(genes.chr1.1,genes.chr2.2,use='pairwise.complete')
corr12=unlist(corr.chr1.chr2)
 breakpoint=c(-1,-0.5,-0.25,-0.1,-0.05,-0.025,-0.01,-0.005,-0.001,0,0.001,0.005,0.01,0.025,0.05,0.1,0.25,0.5,1)
freq=cut(corr12,breaks=breakpoint)%>%table(.)
print('freq')
print(freq)
proportion=freq/length(corr12)
print('proportion')
print(proportion)
b=list(freq,proportion)
return(b)}

#simulation

 kk=c(10,20,50,100,200,300,400,500)

for(k in 1:10){
freq.sim=list()
prop.sim=list()
for(j in 1:length(kk)){
g1=matrix(0,nrow=n.size,ncol=kk[j])
g2=matrix(0,nrow=n.size,ncol=kk[j])
for(i in 1:kk[j]){
g1[,i]=-log(runif(n.size,-0,1),10)
g2[,i]=-log(runif(n.size,-0,1),10)}

freq.sim[[j]]=count2(g1,g2)[[1]]
prop.sim[[j]]=count2(g1,g2)[[2]]}

freq=freq.sim%>%data.frame()
prop=prop.sim%>%data.frame()
freq.final=freq[,c(1,2*c(1:length(kk)))]
colnames(freq.final)=c('bin',kk)
prop.final=prop[,c(1,2*c(1:length(kk)))]
colnames(prop.final)=c('bin',kk)
freq.final
prop.final
write.csv(freq.final,paste0('/home/unfated/realdata/202010gene/freq.logcor20201113.sim',k,'.csv'),row.names=F)
write.csv(prop.final,paste0('/home/unfated/realdata/202010gene/prop.logcor20201113.sim',k,'.csv'),row.names=F)
print(k)
}

#perm

count3<-function(k) {

 #k=20
s1=sample(1:ncol(genes.chr1),k)
s1
s2=sample(1:ncol(genes.chr2),k)
s2
#-log((genes.chr1[,2]),10)%>%head()
#head(genes.chr1[,2])
#log(10,10)
#k=10
a=apply(genes.chr1[,s1],2,function(x) -log(as.numeric(x),10))
b=apply(genes.chr2[,s2],2,function(x) -log(as.numeric(x),10))
resampling<-function(x) {
  aa=x
  aa[!is.na(x)]<-sample(na.omit(x),length(na.omit(x)),replace=F)
  return(aa)}
genes.chr1.1<-t(apply(a,1,resampling))
genes.chr2.2<-t(apply(b,1,resampling))



corr.chr1.chr2=cor(genes.chr1.1,genes.chr2.2,use='pairwise.complete')
corr.chr1.chr2%>%nrow()
corr.chr1.chr2%>%ncol()
corr12=unlist(corr.chr1.chr2)
length(corr12)

#sum(corr12>0.01)
breakpoint=c(-1,-0.5,-0.25,-0.1,-0.05,-0.025,-0.01,-0.005,-0.001,0,0.001,0.005,0.01,0.025,0.05,0.1,0.25,0.5,1)
freq=cut(corr12,breaks=breakpoint)%>%table(.)
print('freq')
print(freq)
proportion=freq/length(corr12)
print('proportion')
print(proportion)
b=list(freq,proportion)
return(b)}

#permu
kk=c(10,20,50,100,200,300,400,500)

for(k in 1:10){
freq.n=list()
prop.n=list()
for(i in 1:length(kk)){
freq.n[[i]]=count3(kk[i])[[1]]
prop.n[[i]]=count3(kk[i])[[2]]
}


freq=freq.n%>%data.frame()
prop=prop.n%>%data.frame()
prop
freq.final=freq[,c(1,2*c(1:length(kk)))]
colnames(freq.final)=c('bin',kk)
prop.final=prop[,c(1,2*c(1:length(kk)))]
colnames(prop.final)=c('bin',kk)
write.csv(freq.final,paste0('/home/unfated/realdata/202010gene/freq.logcor20201113.perm',k,'.csv'),row.names=F)
write.csv(prop.final,paste0('/home/unfated/realdata/202010gene/prop.logcor20201113.perm',k,'.csv'),row.names=F)
print(k)
}

# 20210205----
library(tidyverse)

genes.chr1=read.table('/home/unfated/realdata/202010gene/genes.chr1.txt')
head(genes.chr1)
genes.chr2=read.table('/home/unfated/realdata/202010gene/genes.chr2.txt')

genes.chr1=genes.chr1[-nrow(genes.chr1),]
genes.chr2=genes.chr2[-nrow(genes.chr2),]

str_gene<-function(s1,s2,threshold=0.5){
genes.chr1.1<-apply(genes.chr1[,s1],2,function(x) -log(as.numeric(x),10))
#2377
aa1=apply(genes.chr1.1,2,function(x) sum(is.na(x)))
genes.chr1.1=genes.chr1.1[,aa1<2377]
genes.chr2.2<-apply(genes.chr2[,s2],2,function(x) -log(as.numeric(x),10))
aa2=apply(genes.chr2.2,2,function(x) sum(is.na(x)))
genes.chr2.2=genes.chr2.2[,aa2<2377]
#genes.chr1.1[1:3,1:3]
corr.chr1.chr2=cor(genes.chr1.1,genes.chr2.2,use='pairwise.complete')
#corr.chr1.chr2%>%nrow()
#corr.chr1.chr2%>%ncol()
#corr.chr1.chr2[which(corr.chr1.chr2>0.1)]
corr12=unlist(corr.chr1.chr2)
#which(corr12>0.1)
#
x1=c()
x2=c()
x3=c()
a=rownames(corr12)
b=colnames(corr12)
for(i in 1:nrow(corr12)){
  #row chr1
  for(j in 1:ncol(corr12)){
    #col chr2
    if(corr12[i,j]>threshold){
      x1=c(x1,corr12[i,j])
      x2=c(x2,a[i])
      x3=c(x3,b[j])
    }
    
  }
}
s.gene = data.frame(size=x1,chr1=x2,chr2=x3)
return(s.gene)}

ss1=list(1:500,501:1000,1001:1500,1501:ncol(genes.chr1))
ss2=list(1:400,401:800,801:ncol(genes.chr2))
ncol(genes.chr2)
length(ss1)

gene.out1=list()
for(i in 1:4){
gene.out1[[i]]=str_gene(ss1[[i]],ss2[[1]])
print(i)
}

gene.out2=list()
for(i in 1:4){
gene.out2[[i]]=str_gene(ss1[[i]],ss2[[2]])
print(i)
}

gene.out3=list()
for(i in 1:4){
gene.out3[[i]]=str_gene(ss1[[i]],ss2[[3]])
print(i)
}

head(gene.out1[[1]])

xx=list()
for(i in 1:4){
xx[[i]]=rbind(gene.out1[[i]],gene.out2[[i]])%>% rbind(.,gene.out3[[i]])
write.csv(gene.out1[[i]],paste0('/home/unfated/realdata/202010gene/v2.genes.larger0.5.chr1.',i,'.chr2.1.csv'),row.names=F,quote=F)
write.csv(gene.out2[[i]],paste0('/home/unfated/realdata/202010gene/v2.genes.larger0.5.chr1.',i,'.chr2.2.csv'),row.names=F,quote=F)
write.csv(gene.out3[[i]],paste0('/home/unfated/realdata/202010gene/v2.genes.larger0.5.chr1.',i,'.chr2.3.csv'),row.names=F,quote=F)
}
gene.out.final=xx[[1]]
for(i in 2:4){
gene.out.final=rbind(gene.out.final,xx[[i]])}

head(gene.out.final)
nrow(gene.out.final)
write.csv(gene.out.final,'/home/unfated/realdata/202010gene/v2.genes.larger0.5.chr1.chr2.all.csv',row.names=F,quote=F)

gene.out.final[which(gene.out.final$size==max(gene.out.final$size)),]

out1=genes.chr1[,which(colnames(genes.chr1)=='ENSG00000120341')]
sum(na.omit(out1<10^(-6)))
write.table(out1,'out.pval.txt')
out2=genes.chr2[,which(colnames(genes.chr2)=='ENSG00000151353')]
sum(na.omit(out2<10^(-6)))
a=out1<10^(-6)
b=out2<10^(-6)
out=a+b
sum(na.omit(out)>1)
write.table(out1,'out.pval.txt')

top=read.csv('/home/unfated/realdata/202010gene/out2.csv')
top

t=0.05/20188
x1=c()
x2=c()
x3=c()
x4=list()
x5=list()
for(i in 1:nrow(top)){
out1=genes.chr1[,which(colnames(genes.chr1)==top$chr1[i])]
out2=genes.chr2[,which(colnames(genes.chr2)==top$chr2[i])]
a=out1<t
b=out2<t
out=a+b
x1[i]=sum(na.omit(out1<t))
x2[i]=sum(na.omit(out2<t))
x3[i]=sum(na.omit(out)>1)
x4[[i]]=apply(gwas[which(a==T),],2,function(x)length(unique(x)))
x5[[i]]=apply(gwas[which(b==T),],2,function(x)length(unique(x)))
}

y=data.frame(x4)
nrow(y)
ncol(y)
y=t(y)
colnames(y)=colnames(gwas)
rownames(y)=1:nrow(y)
y1=y
y=data.frame(x5)
nrow(y)
ncol(y)
y=t(y)
colnames(y)=colnames(gwas)
rownames(y)=1:nrow(y)
y2=y

top$gwas1=x1
top$gwas2=x2
top$gwas.overlap=x3
top$percent=x3/min(x1,x2)
top.final=cbind(top[,c(2,5,7)],y1)%>%cbind(.,top[,c(6,8)])%>% cbind(.,y2)%>%cbind(.,top[,c(9)])
write.csv(top.final,'top.gene.final.csv',row.names=F,quote=F)
out2=genes.chr2[,which(colnames(genes.chr2)=='ENSG00000151353')]
sum(na.omit(out2<10^(-6)))

gwas.info=read.delim('/home/unfated/realdata/202010gene/gwasATLAS_v20191115.txt',header=T)
gwas=gwas.info[,7:10]
colnames(gwas.info)[7:10]

#20210208----


genes.chr1=read.table('/home/unfated/realdata/202010gene/genes.chr1.txt')
head(genes.chr1)
genes.chr2=read.table('/home/unfated/realdata/202010gene/genes.chr2.txt')

genes.chr1=genes.chr1[-nrow(genes.chr1),]
genes.chr2=genes.chr2[-nrow(genes.chr2),]

colnames(genes.chr1)[1]
x1=apply(genes.chr1,2,function(x)sum(is.na(x)))
genes.chr11=genes.chr1[,which(x1/nrow(genes.chr1)<=0.5)]
ncol(genes.chr11)
x2=apply(genes.chr2,2,function(x)sum(is.na(x)))
sum(x2/nrow(genes.chr2)>0.5)
genes.chr22=genes.chr2[,which(x2/nrow(genes.chr2)<=0.5)]

y1=apply(genes.chr11,2,function(x)mean(x,na.rm=T))
z1=apply(genes.chr11,1,function(x)mean(x,na.rm=T))
y11=apply(genes.chr11,2,function(x)median(x,na.rm=T))
z11=apply(genes.chr11,1,function(x)median(x,na.rm=T))
y2=apply(genes.chr22,2,function(x)mean(x,na.rm=T))
z2=apply(genes.chr22,1,function(x)mean(x,na.rm=T))
y22=apply(genes.chr22,2,function(x)median(x,na.rm=T))
z22=apply(genes.chr22,1,function(x)median(x,na.rm=T))
med1=median(unlist(genes.chr11),na.rm=T)
mean1=mean(unlist(genes.chr11),na.rm=T)
med2=median(unlist(genes.chr22),na.rm=T)
mean2=mean(unlist(genes.chr22),na.rm=T)


df1=data.frame(y1,y11)
df2=data.frame(y2,y22)
df11=data.frame(z1,z11)
df22=data.frame(z2,z22)


cols=c('mean'='red','median'='blue')

p1<-ggplot(df1) + 
geom_histogram(aes(x=y1,color='mean'),fill="red",alpha=0.5)+ 
geom_histogram(aes(x=y11,color="median"), fill="blue",alpha=0.5)+ 
geom_vline(aes(xintercept=mean1),
              color="blue", linetype="dashed", size=1)+
geom_vline(aes(xintercept=med1),
             color="red", linetype="dotted", size=1)+
             xlab('mean/median p value')+
                   scale_colour_manual(name="",values=cols)+
ggtitle('a) Distribution of across-GWASs mean/median p value  of genes from Chromosome 1')+
  theme(plot.title = element_text(size=5),legend.position=NULL)+
  theme_bw()


p2<-ggplot(df2) + 
geom_histogram(aes(x=y2,color='mean'),fill="red",alpha=0.5)+ 
geom_histogram(aes(x=y22,color="median"), fill="blue",alpha=0.5)+ 
geom_vline(aes(xintercept=mean2),
              color="blue", linetype="dashed", size=1)+
geom_vline(aes(xintercept=med2),
             color="red", linetype="dotted", size=1)+
              scale_colour_manual(name="",values=cols)+
             xlab('mean/median p value')+ 
ggtitle('b) That of genes from Chr 2')+
  theme(plot.title = element_text(size=5))+
    theme_bw()


p11<-ggplot(df11) + 
geom_histogram(aes(x=z1,color='mean'),fill="red",alpha=0.5)+ 
geom_histogram(aes(x=z11,color="median"), fill="blue",alpha=0.5)+ 
geom_vline(aes(xintercept=mean1),
              color="blue", linetype="dashed", size=1)+
geom_vline(aes(xintercept=med1),
             color="red", linetype="dotted", size=1)+
             xlab('mean/median p value')+
                   scale_colour_manual(name="",values=cols)+
ggtitle('c) Distribution of across-chr1-gene mean/median p value of GWASs')+
  theme(plot.title = element_text(size=5),legend.position=NULL)+
    theme_bw()


p22<-ggplot(df22) + 
geom_histogram(aes(x=z2,color='mean'),fill="red",alpha=0.5)+ 
geom_histogram(aes(x=z22,color="median"), fill="blue",alpha=0.5)+ 
geom_vline(aes(xintercept=mean2),
              color="blue", linetype="dashed", size=1)+
geom_vline(aes(xintercept=med2),
             color="red", linetype="dotted", size=1)+
             xlab('mean/median p value')+
              scale_colour_manual(name="",values=cols)+
ggtitle('d) That of across-chr2-gene mean/median p value of GWAS')+
  theme(plot.title = element_text(size=5))+
    theme_bw()


library(Rmisc)
p<-multiplot(p1,p11,p2,p22,cols=2)


png('/home/unfated/realdata/202010gene/mean1.png',width=1200,height=1100)
print(multiplot(p1,p11,p2,p22,cols=2))
dev.off()

gene.chr.all=cbind(genes.chr11,genes.chr22)
ncol(genes.chr11)
gene.chr.all=apply(gene.chr.all,2,function(x)-log(as.numeric(x),10))
ncol(gene.chr.all)
t=200
trial=gene.chr.all[,c((ncol(genes.chr11)-t+1):(ncol(genes.chr11)+t))]
head(trial)
ncol(trial)

library(RColorBrewer)

png('/home/unfated/realdata/202010gene/heatmap2.png',width=1000,height=1000)
print(heatmap(as.matrix(trial[1:400,]),Rowv=NA,Colv=NA,labRow=NA,labCol=NA,xlab='last 200 Chr1 genes + first 200 Chr genes',ylab='The first 400 GWASs',
col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
)
dev.off()



png('/home/unfated/realdata/202010gene/heatmap3.png',width=1000,height=1000)
print(heatmap(as.matrix(trial[1:400,]),Rowv=NA,Colv=NA,labRow=NA,labCol=NA,xlab='last 200 Chr1 genes + first 200 Chr genes',ylab='The first 400 GWASs'

))
dev.off()


cor400=cor(trial,use='pairwise.complete')
ncol(trial)
head(cor400)
ncol(cor400)


library(corrplot)
mycols <- colorRampPalette(colors = c("blue",'white',"red"))(100)
png('/home/unfated/realdata/202010gene/corrmap3.png',width=1000,height=1000)
print(heatmap(as.matrix(cor400),Rowv=NA,Colv=NA,labRow=NA,labCol=NA,xlab='last 200 Chr1 genes + first 200 Chr genes',ylab='last 200 Chr1 genes + first 200 Chr genes'
,main='correlation matrix'
,scale="none"
,col=mycols

))
dev.off()

png('/home/unfated/realdata/202010gene/corrmap2.png',width=1000,height=1000)
print(
  corrplot(cor400,method='color',type='upper',
  ,tl.pos='n',
,main='correlation matrix of last 200 genes at CHR 1 + first 200 genes on CHR 2')
)
dev.off()














