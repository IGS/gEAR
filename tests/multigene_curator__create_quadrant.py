#!/opt/bin/python3

"""

This testing script performs the following actions:

1. Open multigene curator page
2. Logs user in (if necessary)
3. Select a dataset
4. Choose plot type
5. Choose some genes
6. Select some options

"""

import argparse, sys, time

from selenium import webdriver

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support.select import Select

TIMEOUT_PERIOD = 5
DATASET_TITLE = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
PLOT_TYPE_TEXT = "Quadrant"
GENECART_NAME = "sadkins_savetest"

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

    try:
        browser.get(url)

        # Check if logged in, and do so
        # Dataset selection
        print("-- DATASET SELECTION")
        # NOTE: Uses JSTree
        dataset_box = WebDriverWait(browser, timeout=TIMEOUT_PERIOD).until(lambda d: d.find_element(By.ID, 'dataset'))
        dataset_box.click()
        dataset_search_box = browser.find_element(By.ID, "dataset_tree_q")
        dataset_search_box.send_keys("kelley")
        dataset_tree_items = WebDriverWait(browser, timeout=TIMEOUT_PERIOD).until(lambda d: d.find_elements(By.CLASS_NAME, "jstree-search"))
        for elt in dataset_tree_items:
            if elt.text == DATASET_TITLE:
                elt.click()
                break

        if dataset_box.text == DATASET_TITLE:
            results.append({"success": 1, "label": "Dataset selected from tree"})
        else:
            results.append({"success": 0, "label": "Dataset selected from tree"})

        time.sleep(TIMEOUT_PERIOD)

        # Select plot type
        print("-- PLOT TYPE SELECTION")
        # NOTE: Select2 is actually used in the page, and uses a different set of HTML tags to abstract the select element
        select2_plot = browser.find_element(By.ID, 'select2-plot_type_select-container')
        select2_plot.click()
        select2_plot_list = browser.find_element(By.ID, 'select2-plot_type_select-results')
        select2_plot_list_elts = select2_plot_list.find_elements(By.TAG_NAME, "li")
        for elt in select2_plot_list_elts:
            if elt.text == PLOT_TYPE_TEXT:
                elt.click()
                break

        # For some reason the correct plot is selected, but not displayed in the select2 closed dropdown
        if select2_plot.text == PLOT_TYPE_TEXT:
            results.append({"success": 1, "label": "Plot type selected from select2 dropdown"})
        else:
            results.append({"success": 0, "label": "Plot type selected from select2 dropdown"})

        # Choose some genes
        print("-- GENE SELECTION - VIA GENE CART")
        # NOTE: Uses JSTree
        genecart_box = WebDriverWait(browser, timeout=TIMEOUT_PERIOD).until(lambda d: d.find_element(By.ID, 'gene_cart'))
        genecart_box.click()
        genecart_search_box = browser.find_element(By.ID, "gene_cart_tree_q")
        genecart_search_box.send_keys(GENECART_NAME)
        genecart_tree_items = WebDriverWait(browser, timeout=TIMEOUT_PERIOD).until(lambda d: d.find_elements(By.CLASS_NAME, "jstree-search"))
        for elt in genecart_tree_items:
            if elt.text == GENECART_NAME:
                elt.click()
                break
        select2_gene_list = browser.find_element(By.ID, "select2-gene-dropdown-container")
        select2_gene_list_elts = select2_gene_list.find_elements(By.TAG_NAME, "li")
        if len(select2_gene_list_elts):
            results.append({"success": 1, "label": "Gene cart selected from tree"})
        else:
            results.append({"success": 0, "label": "Gene cart selected from tree"})

        # Choose some options

        #except Exception as e:
        #print(str(e), file=sys.stderr)
    finally:
        browser.quit()
        return results

if __name__ == "__main__":
    results = main()
    print(results)