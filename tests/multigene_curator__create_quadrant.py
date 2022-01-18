#!/opt/bin/python3

"""

This testing script performs the following actions:

1. Open multigene curator page
2. Logs user in (if necessary)
3. Select a dataset
4. Choose plot type
5. Choose some genes
6. Select some options
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

import common.multigene_curator as mg

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--localhost', action="store_true", help="Run tests on localhost instead of umgear.org")
    args = parser.parse_args()

    # Determine location to test
    url = 'https://umgear.org/'
    if args.localhost:
        url = "http://localhost:8080/"

    url += "multigene_curator.html"

    results = []

    browser = webdriver.Chrome()
    mg_test = mg.MGTest("Quadrant", browser)

    try:
        mg_test.browser.get(url)

        # Check if logged in, and do so
        # Dataset selection
        if mg_test.test_dataset_selection():
            results.append({"success": 1, "label": "Dataset selected from tree"})
        else:
            results.append({"success": 0, "label": "Dataset selected from tree"})

        time.sleep(mg_test.timeout)

        if mg_test.test_plot_load_after_dataset_selection():
            results.append({"success": 1, "label": "Default plot loaded after dataset selection"})
        else:
            results.append({"success": 0, "label": "Default plot loaded after dataset selection"})

        # Choose some genes
        if mg_test.test_gene_entry():
            results.append({"success": 1, "label": "Genes typed in manually"})
        else:
            results.append({"success": 0, "label": "Genes typed in manually"})


        # Select plot type
        if mg_test.test_plot_type_selection():
            results.append({"success": 1, "label": "Plot type selected from select2 dropdown"})
        else:
            results.append({"success": 0, "label": "Plot type selected from select2 dropdown"})

        # Choose some options

        # Create Plot
        if mg_test.test_plot_creation():
            results.append({"success": 1, "label": "Quadrant successfully made"})
        else:
            results.append({"success": 0, "label": "Quadrant successfully made"})

    finally:
        mg_test.browser.quit()
        return results

if __name__ == "__main__":
    results = main()
    print(results)