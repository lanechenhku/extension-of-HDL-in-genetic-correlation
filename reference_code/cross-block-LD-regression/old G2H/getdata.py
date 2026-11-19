import pandas as pd
import time
from selenium import webdriver
from selenium.webdriver.common.keys import Keys

with open("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/old G2H/r.txt") as file:
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]
    URL=lines[2]
  #from webdriver_manager.chrome import ChromeDriverManager
  driver = webdriver.Chrome('C:/Users/22367/Downloads/selenium-java-4.0.0/chromedriver.exe')
  #driver = webdriver.Chrome(URL)
  driver.get(URL)
  search=driver.find_element_by_id("searchText")
  # flag=search.is_displayed()
  # if flag:
  #   print("Display is normal");
  # else:
  #   print("Display is abnormal");
  time.sleep(2)
  search.clear()
  search.send_keys(lines[0])
  # 获取文本框内容并打印
  print(search.get_attribute("value"))
  driver.find_element_by_id("searchText").send_keys(Keys.ENTER)
  time.sleep(5)
  download=driver.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[4]/div[2]/div/div/div[1]/a[2]')
  download.click()
  
  time.sleep(2)
  search.clear()
  search.send_keys(lines[1])
  # 获取文本框内容并打印
  print(search.get_attribute("value"))
  driver.find_element_by_id("searchText").send_keys(Keys.ENTER)
  time.sleep(5)
  download=driver.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[4]/div[2]/div/div/div[1]/a[2]')
  download.click()
  time.sleep(8)
