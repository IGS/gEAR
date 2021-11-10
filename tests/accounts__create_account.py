#!/opt/bin/python3

"""

This testing script performs the following actions:

1. Open main page
2. Clicks link to create account
3. Fills out account creation form (data in gear.ini)
4. Attempts to submit without required field (email)
5. Submits completed form

Future things to test

- Test error handling if passwords don't match

"""

from selenium import webdriver
browser = webdriver.Chrome()
import configparser
import time

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By

config = configparser.ConfigParser()
config.read('../gear.ini')

browser = webdriver.Chrome()
browser.get(config['test']['host'])

create_acct_link = browser.find_element(By.ID, 'create_acct_link')
create_acct_link.click()

name_box = browser.find_element(By.ID, 'inputName')
name_box.send_keys(config['test']['user_name'])

inst_box = browser.find_element(By.ID, 'inputInstitution')
inst_box.send_keys(config['test']['user_institution'])

pass_box = browser.find_element(By.ID, 'inputPassword')
pass_box.send_keys(config['test']['password'])

repass_box = browser.find_element(By.ID, 'retypePassword')
repass_box.send_keys(config['test']['password'])

# skip required e-mail and test for error

submit_box = browser.find_element(By.ID, 'btn_account_creation_submit')
submit_box.click()

email_warning =  browser.find_element(By.ID, 'email_invalid')

if email_warning.is_displayed():
    print("Initial E-mail wasn't provided")
else:
    print("Initial E-mail was provided")

email_box = browser.find_element(By.ID, 'inputEmail')
email_box.send_keys(config['test']['user_email'])

submit_box.click()

if email_warning.is_displayed():
    print("Secondary E-mail wasn't provided")
else:
    print("Secondary E-mail was provided")



time.sleep(2)
browser.close()
