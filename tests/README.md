# Testing

Manual UI testing is time consuming and is one of the fastest ways to make people look for another job. To keep our great people, we've created this automated testing framework to exercise as many features of the site as possible. This should free us to do more development and identify issues faster.

Generally testing progresses in a few phases

## API testing

Starting with API testing allows us to make sure the API is working properly to do things like create accounts, insert datasets, etc.  This also sets the stage for UI testing with the data created during these first API steps.  The UI is then tested, and then we remove the testing data.

### Completed API tests

None yet listed

### Current and pending API tests

* Deleting an account

## UI testing ##

We use Selenium for this and these steps can take a while. Automated UI testing isn't necessarily quick, and there are a lot of pages/features to check.

### Completed UI tests

#### Account creation

* Navigate from home page to account creation
* Test error handling
  * Leaving out required field (e-mail)
  * Submitting with duplicate account info
* Saving of new account info

#### Main Page (index)

* Logging in
  * With incorrect credentials
  * With correct credentials

### Current and pending UI tests

#### Main Page (index)

#### Main Page - Display Panel mode

#### Contact Us

#### Analysis (single-cell) workbench

#### Comparision tool

#### Dataset (single-gene) curator

* Create a bar plot
* Create a scatter plot
  * Facet rows and cols
  * color
  * marker size
* Create a line plot
* Create a violin plot
* Create a tSNE/uMAP plot
  * split and color by a category
* Load an existing plot

#### Multigene curator

* Create a heatmap
  * Click some genes to bring up secondary violin plot
  * Cluster observations too
  * Second heatmap with axes flipped
* Create a violin
* Create a volcano
* Create a dotplot
* Create a quadrant plot
* Load an existing plot
* Save a plot
* Save a new gene cart

#### Manual Documentation

#### Dataset Uploader - Expression Data

#### Dataset Uploader - Epigenetic Data

#### Dataset Explorer

#### Gene Cart Explorer

#### Epiviz Panel Designer


## Selenium cheat sheet ##

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
