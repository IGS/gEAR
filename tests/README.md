# Testing

Manual UI testing is time consuming and is one of the fastest ways to make people look for another job. To keep our great people, we've created this automated testing framework to exercise as many features of the site as possible. This should free us to do more development and identify issues faster.

Generally testing progresses in a few phases

## API testing

Starting with API testing allows us to make sure the API is working properly to do things like create accounts, insert datasets, etc.  This also sets the stage for UI testing with the data created during these first API steps.  The UI is then tested, and then we remove the testing data.

### Completed API tests

None yet listed

### Current and pending API tests

* Deleting an account

## UI testing

We use Selenium for this and these steps can take a while. Automated UI testing isn't necessarily quick, and there are a lot of pages/features to check.

## Visual regression testing

These kind of tests take a screenshot of a particular HTML element and compare it to a baseline screenshot to ensure the images have not changed.  This is useful to ensure plots have not changed over time (via algorithm or parameters, etc.).

### Testing using SeleniumBase

The preferred way to do visual regression testing is to use SeleniumBase to take the screenshot and HTML data.  Data is stored in the "visual_baseline" directory.  For each test name, a directory will be created for that test name.  The first time `seleniumbase.BaseCase.check_window()` is run, it will save a screenshot as "baseline.png" and a series of text files corresponding to various levels of HTML granularity.

There are three levels of HTML logs:

1. Nested HTML elements
2. Number 1, plus the attribute properties present
3. Number 2, plus the attribute values present

Calling `check_window(<test_name>, <log_level>)` after a baseline image is present will compare the current browser window's HTML tags to those in the specified level chosen, and will also screenshot a "latest" image for viewing. If the current HTML tags do not match those from the baseline run, the test fails.

Read https://seleniumbase.io/examples/visual_testing/ReadMe/ for more infomation about
other processes that happen when this is called.

### Manual screenshots

To save a screenshot from Chrome:

1. Right-click page -> Inspect
2. Find div/element/HTML subset that contains the image you want to save
3. Right-click element -> Capture Node Screenshot.
4. The PNG image should download to a specified default location.  You can then rename it and move it to a "visual_regression_screenshots" directory for future use.

This could be used for image regression purposes, but currently we have no implementation of this.

## Writing SeleniumBase tests

The paradigm that seems easiest to work with is to use two classes. The first class represents the page being tested, including properties to use in testing. This class will also contain the code that deals with page navigation and manipulation. The second class represents the tests and assertions, and is an extension of the SeleniumBase.BaseCase class.  Do note that an instance from the second class will be passed to methods in the first class, as that object is the SeleniumBase driver itself.

For SeleniumBase tests, we will use pytest to run the tests. To run these tests, run `pytest <script>`.  It will run all tests with "test_" as the function name. To run in a localhost environment (to test on Docker images), pass in `--data=localhost` as a option after the script name, which gets stored in `SeleniumBase.BaseCase.data`.

Reference: https://seleniumbase.io/help_docs/syntax_formats/ (see #5)

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

#### Main Page - Display Panel mode

#### Contact Us

#### Analysis (single-cell) workbench

#### Comparision tool

* Select dataset
* Select conditions
* Ensure condition labels are reflected in the plot
* Ensure plot can be generated
  * Default options
  * Significance test
    * Filter
    * Color
* Gene highlighting
  * Found genes show in plot
  * Not found genes show in that div
* Select genes from plot
  * Ensure they show in table
  * Ensure highlighted genes are colored
* Name and save gene cart
* Download table
* Use visual tool like needle (python package) for visual regression testing (plot doesn't differ)
  * TODO: https://github.com/python-needle/needle
* colorblind mode (check colorized filter colorscheme)

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
* colorblind mode

#### Multigene curator

* Ensure plot is loaded when dataset is chosen (loaded from saved displays or default volcano plot)
* Volcanoes should be disabled when no categories have 2+ groups
* Quadrants should be disabled when no categories have 3+ groups
* Create a heatmap
  * Must have 2+ genes
  * Alt heatmap with cluster observations checkbox
  * Alt heatmap with cluster genes checkbox
  * Alt heatmap with axes flipped
  * Distance metric for clustering observations/genes
* Create a violin
  * Stacked violin plot
  * With jitter
* Create a volcano
  * Categories should be not available from select if they have <2 groups
  * REQUIRED - query/ref conditions
    * Window alert if not chosen or category is different
  * DE Algorithm
  * Annotate non-signficant p-values
  * Use adjusted p-vals
* Create a dotplot
* Create a quadrant plot
  * Categories should be not available from select if they have <3 groups
  * REQUIRED - query1/query2/ref conditions
    * Window alert if not chosen or category is different
  * DE Algorithm
  * Foldchange cutuff
  * FDR cutoff
* Misc.
  * Primary category
  * Secondary category (may need to choose a new dataset)
  * Sort category (either primary or secondary)
* Load an existing plot
* Save a plot
* Save a new gene cart
  * From volcano or quadrant
* Use visual tool like needle (python package) for visual regression testing (plot doesn't differ)
* Ensure heatmap and matrixplot expression value for a single gene and observation/celltype is correct
  * This tests that sorting was fine
  * For heatmaps/violins/dotplots
* colorblind mode (check all plots

#### Manual Documentation

#### Dataset Uploader - Expression Data

#### Dataset Uploader - Epigenetic Data

#### Dataset Explorer

#### Gene Cart Explorer

#### Epiviz Panel Designer

## Selenium cheat sheet

Common imports

```python3
from selenium import webdriver
browser = webdriver.Chrome()
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
```

Getting an element by ID

```python3
name_box = browser.find_element(By.ID, 'inputName')
```

Typing things

```python3
name_box.send_keys('Foo')
```

Clicking something (link, button, etc.)

```python3
submit_box = browser.find_element(By.ID, 'btn_submit')
submit_box.click()
```

Checking if an element is visible

```python3
email_warning =  browser.find_element(By.ID, 'email_invalid')

if email_warning.is_displayed():
    print("Initial E-mail wasn't provided")
```

Clear a form element

```python3
name_box = browser.find_element(By.ID, 'inputName')
name_box.clear()
```
