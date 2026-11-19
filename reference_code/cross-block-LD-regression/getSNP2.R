install.packages('RCurl')
install.packages('devtools')
install.packages('downloader')
install.packages('rvest')
install.packages('RSelenium')
install.packages('bitops')
#利用RCurl包批量下载（抓取）文件
library(devtools)
library(tidyverse)
library(bitops)

install_github(repo= "crubba/Rwebdriver", username = "crubba")
library(RCurl)

library(Rwebdriver)
library(downloader)
library(rvest)
library(RSelenium)

remDr <-remoteDriver(  
  
  browserName = "chrome",  
  
  #remoteServerAddr = "localhost",  
  
  port = 4444L  
  
) 

url="https://atlas.ctglab.nl/PheWAS"
remDr$open() 
remDr$maxWindowSize()
remDr$navigate(url) 
remDr$getTitle()
#remDr$screenshot()
remDr$screenshot(file = 'test.png')

# png(filename = 'fig.png')
# print(fig)
# dev.off()
page = remDr$getPageSource()
#webEle = remDr$findElement(using='partial link text',value="Explore Region")
remDr$mouseMoveToLocation(x=10,y=100)
xpath='/html/body/div[1]/div[2]/div[2]/div/div[4]/button'
xpath2='/html/body/div[1]/div[2]/div[5]/div/div/div[2]/div/div[1]/div[2]/button'
xpath3='/html/body/div[1]/div[2]/div[5]/div/div/div[2]/div/div[1]/div[2]/button'
xpath4='/html/body/div[1]/div[2]/div[6]/div/span/div[4]/div[2]/div[2]/div[6]/button'
button= remDr$findElement(using='xpath',value=xpath4)
button$getElementText()
button$highlightElement()
remDr$mouseMoveToLocation(webElement = button)
remDr$screenshot(file = 'test2.png')
remDr$click(buttontId = 0)
button$clickElement()
button$isElementSelected()
button$isElementDisplayed()
button$describeElement()
Sys.sleep(4)
xls_link <- button$getElementAttribute("href")
curl_download(url = xls_link, destfile = 'test.xls')  

