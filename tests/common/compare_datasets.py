from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait

from dataclasses import dataclass, field
@dataclass(frozen=True)
class CompareTest:
    browser: webdriver = webdriver.Chrome()
    dataset: str = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
    genecart_to_save: str = "sadkins_selenium"
    genes: list = field(default_factory=lambda: ["Pou4f3", "Rfx7", "Sox2"])
    invalid_gene: str = "abc123"
    condition_cat: str = "cluster"
    x_group: str = "HC (i)"
    y_group: str = "TEC"
    timeout: int = 5

    def test_dataset_selection(self) -> bool:
        print("-- DATASET SELECTION")
        try:
            # NOTE: Uses JSTree
            dataset_box = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_element(By.ID, 'dataset'))
            dataset_box.click()
            dataset_search_box = self.browser.find_element(By.ID, "dataset_tree_q")
            dataset_search_box.send_keys("kelley")
            dataset_tree_items = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_elements(By.CLASS_NAME, "jstree-search"))
            for elt in dataset_tree_items:
                if elt.text == self.dataset:
                    elt.click()
                    break
            return True if dataset_box.text == self.dataset else False
        except:
            return False

    def test_condition_selection(self) -> bool:
        return False

    def test_condition_selection_same_conditions(self) -> bool:
        # Should throw error if both X and Y conditions are the same
        return False

    def test_condition_label_adjustment(self) -> bool:
        # Should reflect in the plot's axes
        return False

    def test_significance_test_colorize(self) -> bool:
        return False

    def test_significance_test_filter(self) -> bool:
        return False

    def test_found_gene_highlighting(self) -> bool:
        print("-- GENE HIGHLIGHTING - FOUND GENES")
        try:
            gene_box = self.browser.find_element(By.ID,'highlighted_genes')
            gene_box.sendKeys(', '.join(self.genes))


            if not self.test_plot_creation():
                return False

            # Check if genes are a) in plot, b) in "not found" text, and/or c) in selection table
            # All three of these genes are in this plot
            gene_annotation = self.browser.find_element(By.CSS_SELECTOR, '[data-unformatted="{}"'.format(self.genes[0]))
            return True if gene_annotation else False
        except:
            return False

    def test_notfound_gene_highlighting(self) -> bool:
        print("-- GENE HIGHLIGHTING - NOT FOUND GENES")
        try:
            gene_box = self.browser.find_element(By.ID,'highlighted_genes')
            gene_box.sendKeys(self.invalid_gene)


            if not self.test_plot_creation():
                return False

            # Check if genes are a) in plot, b) in "not found" text, and/or c) in selection table
            # This gene does not exist in the plot
            gene_not_found = self.browser.find_element(By.ID, 'genes_not_found')
            return True if gene_not_found and self.invalid_gene in gene_not_found.getText() else False
        except:
            return False

    def test_plot_creation(self) -> bool:
        print("-- PLOT CREATION")
        try:
            create_plot_btn = self.browser.find_element(By.ID, "btn_apply_dataset_changes")
            create_plot_btn.click()
            plot_container = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_element(By.CLASS_NAME,'plotly'))
            return True if plot_container else False
        except:
            return False

    def test_selected_genes_in_table(self) -> bool:
        return False

    def test_highlighted_genes_in_table_highlighted(self) -> bool:
        return False
