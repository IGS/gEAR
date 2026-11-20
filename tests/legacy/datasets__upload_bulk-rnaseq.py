#!/opt/bin/python3

"""

This testing script performs the following actions:

1. Open main page
2. Navigate to upload page
3. Check that login warning was displayed
4. Log in
5. Make sure upload page is displayed
6. Copy metadata template
7. Edit metadata template with values from gear.ini
8. Submit template for upload
9. Verify template submission worked
10. Submit demo datasets for upload
11. Verify upload completed

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

    logged_in = log_in(config, browser, results)

    if not logged_in:
        return results
        sys.exit(1)

    # Go to the upload page
    browser.get("{0}/upload_dataset.html".format(config['test']['host']))
    metadata_upload_c = browser.find_element(By.ID, 'metadata_upload_c')

    if metadata_upload_c.is_displayed():
        results.append({"success": 1, "label": "Navigate to upload page"})
    else:
        results.append({"success": 0, "label": "Navigate to upload page"})
    
        
    time.sleep(2)
    browser.close()

    import json
    print(json.dumps(results))

    return results

def log_in(config, browser, results):
    user_box = browser.find_element(By.ID, 'user_email')
    pass_box = browser.find_element(By.ID, 'user_pass')
    submit_box = browser.find_element(By.ID, 'btn_sign_in')
    user_logged_in = browser.find_element(By.ID, 'user_logged_in')

    # Now send the corrent user name, password 
    user_box.send_keys(config['test']['user_email'])
    pass_box.send_keys(config['test']['password'])
    submit_box.click()
    time.sleep(1)

    if user_logged_in.is_displayed():
        results.append({"success": 1, "label": "Complete login"})
    else:
        results.append({"success": 0, "label": "Complete login"})

    return results[-1]["success"]

if __name__ == '__main__':
    main()
