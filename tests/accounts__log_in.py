#!/opt/bin/python3

"""

This testing script performs the following actions:

1. Open main page
2. Fills in the user/password credentials (incorrectly)
3. Tests that failed login is reported
4. Fills in the user/password credentials (correctly)
5. Checks that logged in state is shown


Future things to test

- None yet

"""

from selenium import webdriver
import configparser
import time

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By

def main():
    browser = webdriver.Chrome()
    results = list()
    
    config = configparser.ConfigParser()
    config.read('../gear.ini')

    browser = webdriver.Chrome()
    browser.get(config['test']['host'])

    user_box = browser.find_element(By.ID, 'user_email')
    pass_box = browser.find_element(By.ID, 'user_pass')
    submit_box = browser.find_element(By.ID, 'btn_sign_in')

    # Send the wrong ones first
    user_box.send_keys('slartibartfast@magrathea.com')
    pass_box.send_keys('hangthesenseofit')
    submit_box.click()

    if email_warning.is_displayed():
        results.append({"success": 1, "label": "Missing input validation"})
    else:
        results.append({"success": 0, "label": "Missing input validation"})

    # Now send the corrent ones
    user_box.send_keys(config['test']['user_name'])
    pass_box.send_keys(config['test']['password'])



    time.sleep(2)
    browser.close()

    return results

if __name__ == '__main__':
    main()
