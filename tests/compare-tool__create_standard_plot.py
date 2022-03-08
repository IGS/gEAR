#!/opt/bin/python3

"""

This testing script performs the following actions:

1. Open compare datasets page
2. Logs user in (if necessary)
3. Select a dataset
4. Choose X and Y conditions
6. Select no options
7. Verify plot was generated
- NOTE: Does not currently verify accuracy of generated plot

"""

import argparse, sys, time

from selenium import webdriver

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support.select import Select

import common.compare_datasets as compare


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--localhost', action="store_true", help="Run tests on localhost instead of umgear.org")
    args = parser.parse_args()

    # Determine location to test
    url = 'https://umgear.org/'
    if args.localhost:
        url = "http://localhost:8080/"

    url += "compare_datasets.html"

    results = []

    browser = webdriver.Chrome()
    compare_test = compare.CompareTest(browser)

    try:
        compare_test.browser.get(url)

        # Check if logged in, and do so
        # Dataset selection
        if compare_test.test_dataset_selection():
            results.append({"success": 1, "label": "Dataset selected from tree"})
        else:
            results.append({"success": 0, "label": "Dataset selected from tree"})

        time.sleep(compare_test.timeout)

        if compare_test.test_condition_selection():
            results.append({"success": 1, "label": "X and Y conditions selected"})
        else:
            results.append({"success": 0, "label": "X and Y conditions selected"})


        # Create Plot
        if compare_test.test_plot_creation():
            results.append({"success": 1, "label": "Plot successfully made"})
        else:
            results.append({"success": 0, "label": "Plot successfully made"})

    finally:
        compare_test.browser.quit()
        return results

if __name__ == "__main__":
    results = main()
    print(results)