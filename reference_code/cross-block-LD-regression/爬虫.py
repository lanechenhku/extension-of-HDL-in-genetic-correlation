# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 06:18:18 2021

@author: 22367
"""

#安装selenium库
SNPid='rs78378222'
pip install selenium

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
xpath='/html/body/div[1]/div[2]/div/div[4]/div[2]/div/div/div[1]/a[2]'
download=driver.find_element_by_xpath(xpath)
download.click()
