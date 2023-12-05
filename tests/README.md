# Testing

Manual UI testing is time consuming and is one of the fastest ways to make people look for another job. To keep our great people, we've created this automated testing framework to exercise as many features of the site as possible. This should free us to do more development and identify issues faster.

Generally testing progresses in a few phases

## API testing

Starting with API testing allows us to make sure the API is working properly to do things like create accounts, insert datasets, etc.  This also sets the stage for UI testing with the data created during these first API steps.  The UI is then tested, and then we remove the testing data.

Reference:  https://towardsdatascience.com/unit-testing-python-data-visualizations-18e0250430

### Completed API tests

None yet listed

### Current and pending API tests

* Deleting an account

## UI (acceptance) testing

NOTE: This will be outdated as soon as the new tests are written in `<root>/www/js/test`

We use Selenium for this and these steps can take a while. Automated UI testing isn't necessarily quick, and there are a lot of pages/features to check.

Reference: https://medium.com/empathyco/the-front-end-testing-of-data-visualizations-29a5644b9e0e

TODO: Should we monkeypatch API and CGI calls to speed things along? - https://docs.pytest.org/en/7.1.x/how-to/monkeypatch.html

### Visual regression testing

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

Reference: https://seleniumbase.io/help_docs/syntax_formats/ (see #5)

For SeleniumBase tests, we will use pytest to run the tests. To run these tests, run `pytest <script>`. It will run all tests with "test_" as the function name.

To run all scripts in the directory, omit the `<script>` in the pytest arguments, or pass in a directory
. It will run all tests from files with "test_*.py" or "*_test.py" as the filename

To run in a localhost environment (to test on Docker images), pass in `--data=localhost` as a option after the script name, which gets stored in `SeleniumBase.BaseCase.data`.

If a test fails, the default tracebacks can be pretty long, so you can also pass in `tb=short`, `tb=line`, or `tb=no` to shorten the traceback or remove it entirely. Adding the option `-rA` will print a summary table of passes and fails by test.  You can also pass in `--demo` which runs the test at a slower rate, and gives visual indicators on which elements are being interacted with (i.e. clicks).

**RECOMMENDED** Adding `--headless` runs the tests with a headless browser, which can be good if you are running tests on a server, or while you are doing work (when the browser will pop up in front of the window you are working in).

Reference: https://seleniumbase.io/help_docs/customizing_test_runs/#seleniumbase-methods-api-reference

You can also add "browser capabilities" to the Selenium Webdriver by passing them in a JSON string using the `--cap-string` argument, or in a file using `--cap-file`.

References: https://www.selenium.dev/documentation/webdriver/capabilities/shared/
https://seleniumbase.io/help_docs/desired_capabilities/#seleniumbase-methods-api-reference

### Note about commonalities with Selenium

SeleniumBase runs Selenium methods under the hood, and provides default timeout settings, and by default uses CSS selectors... both of which are normally written out explicitly in a Selenium test.  They can each be modified in pretty much every function with a "timeout" argument or a "by" argument, respectively.

You can also use method like `BaseCase.find_element(selector)` to return a WebElement object that can be used with Selenium methods.  A method like `BaseCase.find_elements(selector)` returns a list of WebElement objects to iterate over.  There are also variations of these commands, which handle specfic situations such as if the element is visible, clickable, present, or not.

References: https://seleniumbase.io/help_docs/method_summary/
https://github.com/seleniumbase/SeleniumBase/blob/master/seleniumbase/fixtures/base_case.py (for breakdowns of the methods themselves)

### Note about assert statements

SeleniumBase has it's own "assert" statements, like "assert_element_visible". Unlike the traditional pytest assert statements, these do not take an optional message argument in case of assertion failure.  Don't be like me and spend hours trying to figure out what the test wasn't passing when the element clearly existed.

SeleniumBase also has a way to defer the assertion failure on a test case.  This is very useful if you want to find all the bugs in the test before failing the test.  To do this, you use one of the following functions:

* `BaseCase.deferred_assert_element`
* `BaseCase.deferred_assert_element_present`
* `BaseCase.deferred_assert_text`
* `BaseCase.deferred_check_window` (for visual regression testing)

To process these at the end of the test, call `BaseCase.process_deferred_asserts()`

It is also worth noting that the "deferred_assert_element*" methods have a conditional where if the URL has not changed, the "timeout" parameter is set to 1, which can be problematic if one is waiting for a plot to load. In this case, it is better to set a sleep event after the plot creation.

### Note about matching colors

When trying to check a property for the correct color (to indicate success, error, etc.), match the color using RGB values. From my experience, hex-codes, CSS shorthand names ("red", "darkgrey", etc.), and RGBA values do not match.

### Note about alert boxes

Seems that the Selenium Webdriver capabilities for choosing how to handle "alert" modals is ignored, and the default of dismissing and notifying is forced in SeleniumBase.  For this reason, I would advise not trying to handle alerts yourself and just assumed they were already clicked

https://github.com/seleniumbase/SeleniumBase/discussions/1284

### Note on retrying failed tests

Currently I have had some issues with a test succeeding or failing on an inconsistent basis. Seleniumbase is designed to auto-wait for an element to render or appear visible (with an adjustable timeout).  If you wish to have a test rerun, you can add `from seleniumbase import decorators` and add the `@retry_on_exception` fixture to any pytest function.  Alternatively in the `pytest` command, pass in `--reruns=<NUM> --reruns-delay=<SECONDS>`

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

#### Comparison tool

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
  * Ensure highlighted genes are colored (Pou4f3)
* Visual regression testing
* Download gene selection table
* Name and save gene cart

#### Multigene Curator

* Create a heatmap
  * Must have 2+ genes
  * Alt heatmap with cluster observations checkbox
  * Alt heatmap with cluster genes checkbox
  * Alt heatmap with axes flipped
  * Distance metric for clustering observations/genes
  * Matrix plot
  * Sort by primary category
* Create a violin
  * Stacked violin plot
  * With jitter
* Create a volcano
  * REQUIRED - query/ref conditions
    * Window alert if not chosen or category is different
  * DE Algorithm
  * Annotate non-signficant p-values
  * Use adjusted p-vals
* Create a dotplot
* Create a quadrant plot
  * REQUIRED - query1/query2/ref conditions
    * Window alert if not chosen or category is different
  * DE Algorithm
  * Foldchange cutuff
  * FDR cutoff

### Current and pending UI tests

#### Main Page - Display Panel mode

#### Contact Us

#### Analysis (single-cell) workbench

#### Colorblind mode

* Comparision tool
* Dataset curator
* Multigene curator
* Main page - display panel

#### Tests that use a private dataset instead of a public one (requires login)

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

* Ensure plot is loaded when dataset is chosen (loaded from saved displays or default volcano plot)
* Volcanoes should be disabled when no categories have 2+ groups
* Quadrants should be disabled when no categories have 3+ groups
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

#### Analysis (single-cell) Workbench

* Go through all steps in new analysis
* Resume unsaved analysis
* Resume saved analysis
* Load primary analysis

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
