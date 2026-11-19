#install.packages('reticulate')
library(reticulate)
# py_available()#[1] FALSE  #检查您的系统是否安装过Python
# use_condaenv("C:/Users/22367/AppData/Roaming/Microsoft/Windows/Start Menu/Programs/Anaconda3 (64-bit)")
# use_condaenv("C:/Users/22367/AppData/Local/Continuum/anaconda3")
# 
# py_config()#安装的python版本环境查看，显示anaconda和numpy的详细信息。
# py_available()#[1] TRUE   #检查您的系统是否安装过Python
py_module_available("pandas")#检查“pandas”是否安装
# Sys.setenv(RETICULATE_PYTHON = "C:/Users/22367/AppData/Local/Microsoft/WindowsApps/python.exe")
# Sys.setenv(RETICULATE_PYTHON = "C:/Users/22367/AppData/Local/Continuum/anaconda3/python.exe")

py_install('pandas')
py_install('selenium')

repl_python()
# 加载“panda”数据集

#import pandas as pd
# 载入数据集

SNPid='rs78378222'
xpath='/html/body/div[1]/div[2]/div/div[4]/div[2]/div/div/div[1]/a[2]'

#安装对应浏览器驱动
# 我们可以通过下面的方式先看看Selenium.Webdriver支持哪些浏览器
from selenium import webdriver
import time
print(help(webdriver))
driver = webdriver.Chrome()
# #访问页面
driver.get("https://atlas.ctglab.nl/PheWAS")
print(driver.page_source)
print(driver.current_url)
print(driver.title)
search=driver.find_element_by_id("searchText")
flag=search.is_displayed()
if flag:
  print("Display is normal");
else:
  print("Display is abnormal");
time.sleep(2)
search.clear()
search.send_keys(SNPid)
# 获取文本框内容并打印
print(search.get_attribute("value"))

#driver.find_element_by_id("searchText").click()
#search.click()
#driver.forward()
#driver.back()
#模拟回车操作
#search.submit()
#try:
#    search.submit()
#    print('success')
#except Exception as e:
#    print('fail')

from selenium.webdriver.common.keys import Keys
# 回车按键
driver.find_element_by_id("searchText").send_keys(Keys.ENTER)
#search.submit()
download=driver.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[4]/div[2]/div/div/div[1]/a[2]')
download.click()
# 回到R
exit
