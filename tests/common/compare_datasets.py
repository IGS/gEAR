from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support.select import Select
from selenium.webdriver.support.color import Color

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
        print("-- TEST CONDITION SELECTION - VALID SELECTIONS")
        try:
            condition_x_tab = self.browser.find_element(By.ID, "condition_x_tab")
            condition_y_tab = self.browser.find_element(By.ID, "condition_y_tab")
            condition_collapse = self.browser.find_element(By.ID, "{}_collaps".format(self.condition_cat))
            # Select condition x's group
            condition_x_tab.click()
            condition_collapse.click()
            x_cond_checkbox = self.browser.find_element(By.CSS_SELECTOR, 'data-group="{};-;{}"'.format(self.condition_cat, self.x_group))
            x_cond_checkbox.click()

            # Select condition y's group
            condition_y_tab.click()
            condition_collapse.click()
            y_cond_checkbox = self.browser.find_element(By.CSS_SELECTOR, 'data-group="{};-;{}"'.format(self.condition_cat, self.y_group))
            y_cond_checkbox.click()

            # Click back to x-tab to see if condition is still enabled
            condition_x_tab.click()
            if not x_cond_checkbox.is_selected():
                print("ERROR: X condition was not selected after clicking back")
                return False

            # Repeat for y-tab
            condition_y_tab.click()
            if not y_cond_checkbox.is_selected():
                print("ERROR: Y condition was not selected after clicking back")
                return False

            return True
        except:
            return False

    def test_condition_selection_same_conditions(self) -> bool:
        # Should throw error if both X and Y conditions are the same
        print("-- TEST CONDITION SELECTION - EQUAL SELECTIONS")
        try:
            condition_x_tab = self.browser.find_element(By.ID, "condition_x_tab")
            condition_y_tab = self.browser.find_element(By.ID, "condition_y_tab")
            condition_collapse = self.browser.find_element(By.ID, "{}_collaps".format(self.condition_cat))
            # Select condition x's group
            condition_x_tab.click()
            condition_collapse.click()
            x_cond_checkbox = self.browser.find_element(By.CSS_SELECTOR, 'data-group="{};-;{}"'.format(self.condition_cat, self.x_group))
            x_cond_checkbox.click()

            # Select condition y's group
            condition_y_tab.click()
            condition_collapse.click()
            y_cond_checkbox = self.browser.find_element(By.CSS_SELECTOR, 'data-group="{};-;{}"'.format(self.condition_cat, self.x_group))
            y_cond_checkbox.click()

            # plot should fail
            if not self.test_plot_creation():
                error_container = self.browser.find_element(By.ID, "error_loading_c")
                return True if error_container else False
            return False
        except:
            return False

    def test_condition_label_adjustment(self) -> bool:
        # Should reflect in the plot's axes
        print("-- VERIFY CONDITION LABELS SHOW ON AXES")
        x_label = self.browser.find_element(By.ID, "x_label").getText()
        y_label = self.browser.find_element(By.ID, "y_label").getText()

        if not self.test_plot_creation():
            return False

        # Assumes the plot axes titles are the first/only elements with this class
        xtitle = self.browser.find_element(By.CLASS_NAME, "xtitle").getText()
        ytitle = self.browser.find_element(By.CLASS_NAME, "ytitle").getText()

        if xtitle == x_label and ytitle == y_label:
            return True
        return False

    def test_significance_test_colorize(self) -> bool:
        # Perform a colorize significance test
        print("-- SIGNIFICANCE TEST - COLORIZE MODE")
        try:
            sig_select = Select(self.browser.find_element(By.ID, "statistical_test"))
            sig_select.select_by_value('t-test')
            if not self.test_plot_creation():
                return False
            points = self.browser.find_elements(By.CLASS_NAME, "points")
            # Looking for at least one red-colored point
            RED = Color.from_string('red')
            for p in points:
                if p.value_of_css_property("fill") == RED:
                    return True
            return False
        except:
            return False

    def test_significance_test_filter(self) -> bool:
        # Perform a colorize significance test
        print("-- SIGNIFICANCE TEST - FILTER MODE")
        try:
            sig_select = Select(self.browser.find_element(By.ID, "statistical_test"))
            sig_select.select_by_value('t-test')
            if not self.test_plot_creation():
                return False
            points = self.browser.find_elements(By.CLASS_NAME, "points")
            num_colorized_points = len(points)

            # click filter button
            filter_radio = self.browser.find_element(By.ID, "stat_action_filter")
            filter_radio.click()
            if not self.test_plot_creation():
                return False
            points = self.browser.find_elements(By.CLASS_NAME, "points")
            num_filtered_points = len(points)

            # Colored points should be filtered out so the counts should differ
            return True if not num_colorized_points == num_filtered_points else False
        except:
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
