#!/opt/bin/python3

from selenium import webdriver
from selenium.webdriver.common.by import By

browser = webdriver.Chrome()

from selenium import webdriver
from selenium.webdriver.common.keys import Keys

browser = webdriver.Chrome()
browser.get('https://umgear.org/')

search_box = browser.find_element(By.ID, 'search_gene_symbol_intro')
search_box.send_keys('Sox2')

exact_match = browser.find_element(By.ID, 'exact_match_icon')
exact_match.click()

