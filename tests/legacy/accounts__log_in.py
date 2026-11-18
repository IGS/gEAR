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
import sys


from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By

# This is the color the UI makes text within text boxes to show an error
ERROR_COLOR = 'rgba(255, 0, 0, 1)'

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
    user_logged_in = browser.find_element(By.ID, 'user_logged_in')

    # Send the wrong ones first
    user_box.send_keys('slartibartfast@magrathea.com')
    pass_box.send_keys('hangthesenseofit')
    submit_box.click()
    time.sleep(2)

    if user_box.value_of_css_property('color') == ERROR_COLOR:
        results.append({"success": 1, "label": "Incorrect user name validation"})
    else:
        results.append({"success": 0, "label": "Incorrect user name validation"})

    # Now send the corrent user name, password should still be wrong
    user_box.clear()
    user_box.send_keys(config['test']['user_email'])
    
    submit_box.click()
    time.sleep(1)

    if pass_box.value_of_css_property('color') == ERROR_COLOR:
        results.append({"success": 1, "label": "Incorrect password validation"})
    else:
        results.append({"success": 0, "label": "Incorrect password validation"})

    pass_box.clear()
    pass_box.send_keys(config['test']['password'])
    time.sleep(1)

    # have to reconnect to the submit box since it was pulled from the DOM
    submit_box = browser.find_element(By.ID, 'btn_sign_in')
    submit_box.click()
    time.sleep(1)

    if user_logged_in.is_displayed():
        results.append({"success": 1, "label": "Complete login"})
    else:
        results.append({"success": 0, "label": "Complete login"})

    time.sleep(2)
    browser.close()

    return results

if __name__ == '__main__':
    main()
