from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait

from dataclasses import dataclass
@dataclass
class MGTest:
    plot_type: str
    dataset: str = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
    genecart_to_load: str = "sadkins_savetest"
    genes: list = ["Pou4f3", "Rfx7", "Sox2"]
    timeout: int = 5

    def test_dataset_selection(self, browser) -> bool:
        print("-- DATASET SELECTION")
        try:
            # NOTE: Uses JSTree
            dataset_box = WebDriverWait(browser, timeout=self.timeout).until(lambda d: d.find_element(By.ID, 'dataset'))
            dataset_box.click()
            dataset_search_box = browser.find_element(By.ID, "dataset_tree_q")
            dataset_search_box.send_keys("kelley")
            dataset_tree_items = WebDriverWait(browser, timeout=self.timeout).until(lambda d: d.find_elements(By.CLASS_NAME, "jstree-search"))
            for elt in dataset_tree_items:
                if elt.text == self.dataset:
                    elt.click()
                    break
            return True if dataset_box.text == self.dataset else False
        except:
            return False

    def test_gene_entry(self, browser) -> bool:
        print("-- GENE SELECTION - MANUAL INPUT")
        try:
            select2_gene_box = WebDriverWait(browser, timeout=self.timeout).until(lambda d: d.find_element(By.ID,'select2-gene_dropdown-container'))
            select2_gene_textarea = select2_gene_box.find_element(By.XPATH,"//span/textarea")
            select2_gene_textarea.click()
            for gene in self.genes:
                select2_gene_textarea.send_keys(gene + Keys.ENTER)
            select2_gene_box_elts = select2_gene_box.find_elements(By.TAG_NAME, "li")
            return True if len(select2_gene_box_elts) else False
        except:
            return False

    def test_plot_creation(self, browser) -> bool:
        print("-- PLOT CREATION")
        try:
            create_plot_btn = browser.find_element(By.ID, "create_plot")
            create_plot_btn.click()
            plot_container = WebDriverWait(browser, timeout=self.timeout).until(lambda d: d.find_element(By.CLASS_NAME,'plotly-container'))
            return True if plot_container else False
        except:
            return False

    def test_plot_type_selection(self, browser) -> bool:
        print("-- PLOT TYPE SELECTION")
        try:
            # NOTE: Select2 is actually used in the page, and uses a different set of HTML tags to abstract the select element
            select2_plot = browser.find_element(By.ID, 'select2-plot_type_select-container')
            select2_plot.click()
            select2_plot_list = browser.find_element(By.ID, 'select2-plot_type_select-results')
            select2_plot_list_elts = select2_plot_list.find_elements(By.TAG_NAME, "li")
            for elt in select2_plot_list_elts:
                if elt.text == self.plot_type:
                    elt.click()
                    break

            # For some reason the correct plot is selected, but not displayed in the select2 closed dropdown
            return True if select2_plot.text == self.plot_type else False
        except:
            return False
