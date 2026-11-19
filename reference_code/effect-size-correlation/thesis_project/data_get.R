install.packages('RCurl')
install.packages('devtools')
install.packages('downloader')
install.packages('rvest')
install.packages('RSelenium')
#利用RCurl包批量下载（抓取）文件
library(devtools)
library(tidyverse)

install_github(repo= "crubba/Rwebdriver", username = "crubba")
library('RCurl')

library('Rwebdriver')
library(downloader)
library(rvest)
library(RSelenium)
#爬取网页图片的连接----
url = 'https://list.tmall.com/search_product.htm?spm=a220m.1000858.1000723.1.1c5a7c5awATvmK&&from=rs_1_key-top-s&q=%BB%AA%CE%AA%CA%D6%BB%FA'
remDr <-remoteDriver(  
  
  browserName = "chrome",  
  
  #remoteServerAddr = "localhost",  
  
  port = 4444L  
  
) 

remDr$open() 
remDr$navigate(url) 
picture<- remDr$getPageSource()[[1]]%>%read_html()%>% html_nodes("a.productImg img")%>%html_attr("src")
#网页爬取时没有刷新，所以只爬取到前5条连接，所以先以此为例
picture = picture[1:5]
#构造正确的连接
picture_1 = paste("http:",picture,collapse = "")
#下载图片
for(i in 1:length(picture)) 
{
  download(picture[i],paste("C:\\Users\\Administrator\\Desktop\\picture",i,".jpg",sep = ""), mode = "wb")
}

#批量下载----
url = 'https://elections.maryland.gov/elections/2012/election_data/index.html'
url="https://msk.hugeamp.org/variant.html?variant=rs78378222"
remDr$open() 
remDr$navigate(url) 
#爬取网页的CSV文件连接,选择前5条连接作为例子
links = remDr$getPageSource()[[1]]%>%read_html()%>% html_nodes("tr td a")%>%html_attr("href")
links =  links[1:5]
#构造正确的连接
links_1 = as.list(links)
#构造下载CSV函数
downCSV = function(filename,url,folder){
  dir.create(folder,showWarnings = FALSE)
  Furl = str_c(url,filename)
  if(!file.exists(str_c(folder,"/",filename))){
    download.file(Furl,destfile = str_c(folder,"/",filename))
    Sys.sleep(1)}
}
#批量下载
library(plyr)
l_ply(links_1,downCSV,url="https://www.elections.maryland.gov/elections/2012/election_data/",folder = 'TEST')
