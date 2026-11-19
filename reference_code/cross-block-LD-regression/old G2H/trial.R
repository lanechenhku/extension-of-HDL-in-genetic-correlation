# cd C:\Users\22367\Downloads\selenium-java-4.0.0
#java -jar selenium-server-standalone-3.11.0.jar

library(reticulate)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
SNPid1=c('rs78378222')
SNPid2=c('rs3730556')
SNPpair=c(SNPid1,SNPid2)

effect_size_plot<-function(source='GWAS ATLAS',
                           SNPid1=NA, 
                           SNPid2=NA,
                           SNPid3=NA,
                           SNPid4=NA,
                           SNPid5=NA,
                           path = "C:/Users/22367/Downloads",
                           SNPpair=SNPpair,
                           chosen='-log Pval',
                           logbase=10){
  #py_module_available("pandas")#检查“pandas”是否安装
  #path <- "C:/Users/22367/Downloads"
  if(source=='GWAS ATLAS'){
    file_default<-'Genome wide association study ATLAS'}
  
  fileNames <- dir(path) 
  fileN=fileNames[grep(file_default,fileNames)]
  file.remove(paste0(path,'/',fileN))
  
  #input----
  # SNPid1=NA
  # SNPid2=NA
  # SNPid3=NA
  # SNPid4=NA
  # SNPid5=NA
  
  
  #logbase=10
  SNPlist=na.omit(c(SNPid1,SNPid2,SNPid3,SNPid4,SNPid5))
  SNPno=length(SNPlist)
  #xpath='/html/body/div[1]/div[2]/div/div[4]/div[2]/div/div/div[1]/a[2]'
  if(source=='GWAS ATLAS'){
    URL=c("https://atlas.ctglab.nl/PheWAS")
  }
  write.table(c(SNPpair,URL),"C:/Users/22367/Desktop/thesis/cross-block-LD-regression/old G2H/r.txt",row.names = F,col.names = F,quote = F)
  
  #get data----
  
  source_python("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/getdata.py")
 # py_run_file("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/old G2H/getdata.py")
  
  #read file----
  path <- "C:/Users/22367/Downloads"
  #aaa=paste0('Genome wide association study ATLAS (',c(1:(SNPno-1)),').csv')
  fileNames <- c('Genome wide association study ATLAS.csv','Genome wide association study ATLAS (1).csv')
  filePath <- sapply(fileNames, function(x){ 
    paste(path,x,sep='/')})   
  data <- lapply(filePath, function(x){
    read.csv(x)}) 
  
  #----
  
  # id1=which(SNPlist==SNPpair[1])
  # id2=which(SNPlist==SNPpair[2])
  # df1=data[[id1]]
  # df2=data[[id2]]
  df1=data[[1]]
  df2=data[[2]]
  df=left_join(df1,df2,by='atlas.ID')
  df=na.omit(df)
  variable=c("atlas.ID",  "PMID.x" ,"Domain.x" , "Trait.x", "P.value.x","P.value.y")
  dd=df[,variable]
  # dd$z.x=as.numeric(dd$beta.x)/as.numeric(dd$sd.x)
  # dd$z.y=as.numeric(dd$beta.y)/as.numeric(dd$sd.y)
  if(chosen=='-log Pval'){
    dd$z.x=-log(as.numeric(dd$P.value.x),logbase)
    dd$z.y=-log(as.numeric(dd$P.value.y),logbase)
  }
  dd$dataset=as.character(dd$PMID.x)
  dd$group=dd$Domain.x
  dd$dataset[grep('BioRxiv',dd$PMID.x)]='BioRxiv'
  # dd$dataset[grep('GWAS_CKDGenConsortium',dd$dataset.x)]='GWAS_CKDGenConsortium'
  # dd$dataset[grep('GWAS_GEFOS',dd$dataset.x)]='GWAS_GEFOS'
  # dd$dataset[grep('GWAS_GODARTS',dd$dataset.x)]='GWAS_GODARTS'
  # dd$dataset[grep('GWAS_HOORN',dd$dataset.x)]='GWAS_HOORN'
  # dd$dataset[grep('TOPMED',dd$dataset.x)]='TOPMED'
  p1<-ggplot(data=dd,aes(x=z.x,y=z.y,col=dataset))+geom_point()+xlab(SNPpair[1])+ylab(SNPpair[2])+ggtitle(paste0(chosen,' scatterplot'))+
    geom_line(stat='smooth',method = 'lm', formula = y ~ x+1,se=F,alpha=0.5,size=0.5) + #绘制回归直线
    stat_cor(data = dd, method = "pearson")
  
  p2<-ggplot(data=dd,aes(x=z.x,y=z.y,col=group))+geom_point()+xlab(SNPpair[1])+ylab(SNPpair[2])+ggtitle(paste0(chosen,' scatterplot'))+
    geom_line(stat='smooth',method = 'lm', formula = y ~ x+1,se=F,alpha=0.5,size=0.5) + #绘制回归直线
    stat_cor(data = dd, method = "pearson")
  p<-list(p1,p2)
  return(p)
}

p<-effect_size_plot(source='GWAS ATLAS',
                    SNPid1=SNPid1,
                    SNPid2=SNPid2,
                           SNPid3=SNPid1,
                           SNPid4=SNPid2,
                           SNPid5=NA,
                           path = "C:/Users/22367/Downloads",
                           SNPpair=c('rs78378222','rs3730556'),
                           chosen='-log Pval',
                           logbase=10)
