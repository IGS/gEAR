= Testing =

Manual UI testing is time consuming and is one of the fastest ways to make people look for another job. To keep our great people, we've created this automated testing framework to exercise as many features of the site as possible. This should free us to do more development and identify issues faster.

Generally testing progresses in a few phases

== API testing ==

Starting with API testing allows us to make sure the API is working properly to do things like create accounts, insert datasets, etc.  This also sets the stage for UI testing with the data created during these first API steps.  The UI is then tested, and then we remove the testing data.

=== Current and pending API tests ===



== UI testing ==

We use Selenium for this and these steps can take a while. Automated UI testing isn't necessarily quick, and there are a lot of pages/features to check.


=== Current and pending UI tests ===

- Create account
- Log in




=== Clean up ===

- Delete account


== Selenium cheat sheet ==

Common imports

```
from selenium import webdriver
browser = webdriver.Chrome()
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
```

Getting an element by ID

```
name_box = browser.find_element(By.ID, 'inputName')
```

Typing things

```
name_box.send_keys('Foo')
```

Clicking something (link, button, etc.)

```
submit_box = browser.find_element(By.ID, 'btn_submit')
submit_box.click()
```

Checking if an element is visible

```
email_warning =  browser.find_element(By.ID, 'email_invalid')

if email_warning.is_displayed():
    print("Initial E-mail wasn't provided")
```
