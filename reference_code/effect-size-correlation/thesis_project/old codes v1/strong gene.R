library(tidyverse)
data=read.csv(file.choose())
df=data
ggplot(df, aes(x=size)) + geom_histogram()
# Change the width of bins
# ggplot(df, aes(x=size)) + 
#   geom_histogram(binwidth=1)
# Change colors
p<-ggplot(df, aes(x=size)) + 
  geom_histogram(color="black", fill="grey")
p

# Add mean line
p+ geom_vline(aes(xintercept=mean(size)),
              color="blue", linetype="dashed", size=1)+
geom_vline(aes(xintercept=median(size)),
             color="red", linetype="dotted", size=1)+xlab('effect size correlation (r)')+
  ggtitle('9859 cross-chr1-chr2 gene pairs with r>0.5 ')+theme_bw()


x1=data.frame(df[,2])
write.table(x1,'chr1.str.gene.txt',row.names = F,col.names = F,quote = F)
x2=data.frame(df[,3])
write.table(x2,'chr2.str.gene.txt',row.names = F,col.names = F,quote = F)
unique(df$chr1) %>% length()
unique(df$chr2) %>% length()
#ID
library(openxlsx)
#library(tidyverse)
library(readxl)
chr1 <- read_excel("chr1_ID/idmap.xlsx")
chr2 <- read_excel("chr2_ID/idmap.xlsx")

chr11=read.table('chr1_ID.txt')
chr22=read.table('chr2_ID.txt')
df2=read.csv(file.choose())
aa=which(df2$size%in%df$size)
df$chr1.sym=chr11$V1[aa]
df$chr2.sym=chr22$V1[aa]

out1=subset(df,size>0.8)
a=out1$chr1.sym %>% as.character() %>% unique
b=out1$chr2.sym %>% as.character() %>% unique
c=c(a,b)
#write.table(data.frame(c),'topgene0.85.txt',row.names = F,col.names = F,quote = F)
colnames(chr1)[4]='chr1.sym'
colnames(chr2)[4]='chr2.sym'
out1$chr1.sym=as.character(out1$chr1.sym)
out1$chr2.sym=as.character(out1$chr2.sym)
chr1$chr1.sym=as.character(chr1$chr1.sym)
chr2$chr2.sym=as.character(chr2$chr2.sym)
out1=left_join(out1,chr1,by='chr1.sym')
out1=left_join(out1,chr2,by='chr2.sym')
colnames(out1)
out2=out1[,c(1,4,8,9,5,14,15)]
#write.csv(out2,'out.csv')
write.csv(out1[,1:5],'out2.csv')

#gwas info----

