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

from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support.select import Select

import common.multigene_curator as mg

TIMEOUT_PERIOD = 5
DATASET_TITLE = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
PLOT_TYPE_TEXT = "Heatmap"
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
        if mg.test_dataset_selection(browser, DATASET_TITLE, TIMEOUT_PERIOD):
            results.append({"success": 1, "label": "Dataset selected from tree"})
        else:
            results.append({"success": 0, "label": "Dataset selected from tree"})

        time.sleep(TIMEOUT_PERIOD)

        # Select plot type
        if mg.test_plot_type_selection(browser, PLOT_TYPE_TEXT):
            results.append({"success": 1, "label": "Plot type selected from select2 dropdown"})
        else:
            results.append({"success": 0, "label": "Plot type selected from select2 dropdown"})

        # Choose some genes
        if mg.test_gene_entry(browser, TIMEOUT_PERIOD):
            results.append({"success": 1, "label": "Genes typed in manually"})
        else:
            results.append({"success": 0, "label": "Genes typed in manually"})

        """
        print("-- GENE SELECTION - VIA GENE CART")
        # NOTE: Uses JSTree and requires login
        genecart_box = WebDriverWait(browser, timeout=TIMEOUT_PERIOD).until(lambda d: d.find_element(By.ID, 'gene_cart'))
        genecart_box.click()
        genecart_search_box = browser.find_element(By.ID, "gene_cart_tree_q")
        genecart_search_box.send_keys(GENECART_NAME)
        genecart_tree_items = WebDriverWait(browser, timeout=3).until(lambda d: d.find_elements(By.CLASS_NAME, "jstree-search"))
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
        """

        # Choose some options
        """
        print("-- FILTER_BY SELECTION")
        # In this case, all groups in all observations are included.  Need to click 'close' on some groups
        select2_cluster_filter_by_box = WebDriverWait(browser, timeout=TIMEOUT_PERIOD).until(lambda d: d.find_element(By.ID,'select2-cluster_dropdown-container'))
        select2_cluster_filter_by_textarea = select2_cluster_filter_by_box.find_element(By.XPATH,"//span/textarea")
        select2_cluster_filter_by_textarea.click()
        select2_cluster_filter_by_textarea.send_keys("HC (i)" + Keys.ENTER)
        select2_cluster_filter_by_textarea.send_keys("SC (i)" + Keys.ENTER)
        select2_cluster_filter_by_textarea.send_keys("TEC" + Keys.ENTER)
        select2_cluster_filter_by_box_elts = select2_cluster_filter_by_box.find_elements(By.TAG_NAME, "li")
        if len(select2_cluster_filter_by_box_elts):
            results.append({"success": 1, "label": "Filter by 'cluster' category set"})
        else:
            results.append({"success": 0, "label": "Filter by 'cluster' category set"})
        """

        print("-- GROUP_BY SELECTION")
        cluster_group_by_radio = browser.find_element(By.ID, "cluster_groupby")
        cluster_group_by_radio.click()

        # Not worrying about distance metric - Euclidean is default

        # Create Plot
        if mg.test_plot_creation(browser, TIMEOUT_PERIOD):
            results.append({"success": 1, "label": "Heatmap successfully made"})
        else:
            results.append({"success": 0, "label": "Heatmap successfully made"})

    except Exception as e:
        print(str(e), file=sys.stderr)
    finally:
        browser.quit()
        return results

if __name__ == "__main__":
    results = main()
    print(results)