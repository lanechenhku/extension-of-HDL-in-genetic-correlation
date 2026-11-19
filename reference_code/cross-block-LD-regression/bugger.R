install.packages('RCurl')
#利用RCurl包批量下载（抓取）文件

library('RCurl')

html=getURL("http://rfunction.com/code/1202/")

#查看网页源码，之后确定抓取信息的“代码字串”特征

temp=strsplit(html,"<li><a href=\"")[[]]   #列表

files=strsplit(temp,"\"")   #转义字符“\”

files=lapply(files,function(x){x[]}) 

#可通过命令“?lapply”查看该函数具体用法

files

class(files)

files=unlist(files)    #二维转一维

files=files[-(:)]     #剔除特定元素

#----------------

#设置下载目录，循环遍历

setwd('G:\\R_Project\\RCurl抓取的文件')     #注意‘\\’转义

dir()

i=
  
  base="http://rfunction.com/code/1202/"

for(i in :length(files)){
  
  url=paste(base,files[i],sep='')    #拼接url
  
  temp=getBinaryURL(url)    #获取网页内容
  
  note=file(paste("",files[i],sep='.'),open="wb")   #文件属性
  
  writeBin(temp,note)     #文件写入内容
  
  close(note)         #关闭文件
  
}