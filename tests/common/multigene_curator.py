from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait

from dataclasses import dataclass, field
@dataclass(frozen=True)
class MGTest:
    plot_type: str
    browser: webdriver = webdriver.Chrome()
    dataset: str = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
    genecart_to_load: str = "sadkins_savetest"
    genes: list = field(default_factory=lambda: ["Pou4f3", "Rfx7", "Sox2"])
    filter_by: list = field(default_factory=lambda: ["HC (i)", "SC (i)", "TEC"])
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

    def test_filter_by_selection(self) -> bool:
        print("-- FILTER_BY SELECTION")
        try:
            # In this case, all groups in all observations are included.  Need to click 'close' on some groups
            select2_cluster_filter_by_box = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_element(By.ID,'select2-cluster_dropdown-container'))
            select2_cluster_filter_by_textarea = select2_cluster_filter_by_box.find_element(By.XPATH,"//span/textarea")
            select2_cluster_filter_by_textarea.click()
            for cat in self.filter_by:
                select2_cluster_filter_by_textarea.send_keys(cat + Keys.ENTER)
                select2_cluster_filter_by_textarea.send_keys(cat + Keys.ENTER)
                select2_cluster_filter_by_textarea.send_keys(cat + Keys.ENTER)
            select2_cluster_filter_by_box_elts = select2_cluster_filter_by_box.find_elements(By.TAG_NAME, "li")
            return True if len(select2_cluster_filter_by_box_elts) else False
        except:
            return False

    def test_gene_entry(self) -> bool:
        print("-- GENE SELECTION - MANUAL INPUT")
        try:
            select2_gene_box = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_element(By.ID,'select2-gene_dropdown-container'))
            select2_gene_textarea = select2_gene_box.find_element(By.XPATH,"//span/textarea")
            select2_gene_textarea.click()
            for gene in self.genes:
                select2_gene_textarea.send_keys(gene + Keys.ENTER)
            select2_gene_box_elts = select2_gene_box.find_elements(By.TAG_NAME, "li")
            return True if len(select2_gene_box_elts) else False
        except:
            return False

    def test_gene_cart_entry(self) -> bool:
        print("-- GENE SELECTION - VIA GENE CART")
        try:
            # NOTE: Uses JSTree and requires login
            genecart_box = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_element(By.ID, 'gene_cart'))
            genecart_box.click()
            genecart_search_box = self.browser.find_element(By.ID, "gene_cart_tree_q")
            genecart_search_box.send_keys(self.genecart_to_load)
            genecart_tree_items = WebDriverWait(self.browser, timeout=3).until(lambda d: d.find_elements(By.CLASS_NAME, "jstree-search"))
            for elt in genecart_tree_items:
                if elt.text == self.genecart_to_load:
                    elt.click()
                    break
            select2_gene_list = self.browser.find_element(By.ID, "select2-gene-dropdown-container")
            select2_gene_list_elts = select2_gene_list.find_elements(By.TAG_NAME, "li")
            return True if len(select2_gene_list_elts) else False
        except:
            return False

    def test_plot_creation(self) -> bool:
        print("-- PLOT CREATION")
        try:
            create_plot_btn = self.browser.find_element(By.ID, "create_plot")
            create_plot_btn.click()
            plot_container = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_element(By.CLASS_NAME,'plotly-container'))
            return True if plot_container else False
        except:
            return False

    def test_plot_load_after_dataset_selection(self) -> bool:
        print("-- PLOT LOADING AFTER SELECTION OF DATASET")
        try:
            plot_container = WebDriverWait(self.browser, timeout=self.timeout).until(lambda d: d.find_element(By.CLASS_NAME,'plotly-container'))
            return True if plot_container else False
        except:
            return False

    def test_plot_type_selection(self) -> bool:
        print("-- PLOT TYPE SELECTION")
        select_id = "plot_type_select"
        try:
            # NOTE: Select2 is actually used in the page, and uses a different set of HTML tags to abstract the select element
            select2 = self.browser.find_element(By.ID, 'select2-{}-container'.format(select_id))
            select2.click()
            select2_list = self.browser.find_element(By.ID, 'select2-{}-results'.format(select_id))
            select2_list_elts = select2_list.find_elements(By.TAG_NAME, "li")
            for elt in select2_list_elts:
                if elt.text == self.plot_type:
                    elt.click()
                    break

            # For some reason the correct plot is selected, but not displayed in the select2 closed dropdown
            return True if select2.text == self.plot_type else False
        except:
            return False

    def test_quadrant_group_selection(self) -> bool:
        pass

    def test_volcano_group_selection(self) -> bool:
        pass

    ### OPTIONS TESTING

    def test_annotate_nonsignficant(self) -> bool:
        print("-- ANNOTATE NONSIGNIFICANT POINTS CHECKBOX SELECTION")
        try:
            annot_nonsig_checkbox = self.browser.find_element(By.ID, "annot_nonsig")
            annot_nonsig_checkbox.click()
            return True
        except:
            return False

    def test_clusterbar(self) -> bool:
        print(" -- CLUSTERBAR CATEGORY SELECTION")

    def test_cluster_observations(self) -> bool:
        print("-- CLUSTER OBSERVATIONS CHECKBOX SELECTION")
        try:
            cluster_obs_checkbox = self.browser.find_element(By.ID, "cluster_obs")
            cluster_obs_checkbox.click()
            return True
        except:
            return False

    def test_cluster_genesw(self) -> bool:
        print("-- CLUSTER GENES CHECKBOX SELECTION")
        try:
            cluster_genes_checkbox = self.browser.find_element(By.ID, "cluster_genes")
            cluster_genes_checkbox.click()
            return True
        except:
            return False

    def test_de_algo(self) -> bool:
        print("-- DE ALGORITHM SELECTION")
        select_id = "de_test_select"
        try:
            # NOTE: Select2 is actually used in the page, and uses a different set of HTML tags to abstract the select element
            select2 = self.browser.find_element(By.ID, 'select2-{}-container'.format(select_id))
            select2.click()
            select2_list = self.browser.find_element(By.ID, 'select2-{}-results'.format(select_id))
            select2_list_elts = select2_list.find_elements(By.TAG_NAME, "li")
            for elt in select2_list_elts:
                if elt.text == "Welch's t-test":
                    elt.click()
                    break

            return True if select2.text == "Welch's t-test" else False
        except:
            return False

    def test_distance_metric(self) -> bool:
        print("-- DISTANCE METRIC SELECTION")
        select_id = "distance_select"
        try:
            # NOTE: Select2 is actually used in the page, and uses a different set of HTML tags to abstract the select element
            select2 = self.browser.find_element(By.ID, 'select2-{}-container'.format(select_id))
            select2.click()
            select2_list = self.browser.find_element(By.ID, 'select2-{}-results'.format(select_id))
            select2_list_elts = select2_list.find_elements(By.TAG_NAME, "li")
            for elt in select2_list_elts:
                if elt.text == "Euclidean":
                    elt.click()
                    break

            return True if select2.text == "Euclidean" else False
        except:
            return False

    def test_fdr_cutoff(self) -> bool:
        print("-- FDR CUTOFF INPUT")
        try:
            fdr_input = self.browser.find_element(By.ID, "quadrant_foldchange_cutoff")
            fdr_input.send_keys("0.05")   # The default
            fdr_input.click()
            return True
        except:
            return False

    def test_flip_axes(self) -> bool:
        print("-- FLIP AXES CHECKBOX SELECTION")
        try:
            flip_axes_checkbox = self.browser.find_element(By.ID, "flip_axes")
            flip_axes_checkbox.click()
            return True
        except:
            return False

    def test_fold_change_cutoff(self) -> bool:
        print("-- FOLDCHANGE CUTOFF INPUT")
        try:
            foldchange_input = self.browser.find_element(By.ID, "quadrant_foldchange_cutoff")
            foldchange_input.send_keys("2")   # The default
            foldchange_input.click()
            return True
        except:
            return False

    def test_include_zero_foldchange(self) -> bool:
        print("-- INCLUDE ZERO FOLDCHANGE CHECKBOX SELECTION")
        try:
            include_zero_foldchange_checkbox = self.browser.find_element(By.ID, "include_zero_foldchange")
            include_zero_foldchange_checkbox.click()
            return True
        except:
            return False

    def test_matrixplot(self) -> bool:
        print("-- MATRIXPLOT CHECKBOX SELECTION")
        try:
            matrixplot_checkbox = self.browser.find_element(By.ID, "matrixplot")
            matrixplot_checkbox.click()
            return True
        except:
            return False

    def test_primary_category(self) -> bool:
        pass

    def test_secondary_category(self) -> bool:
        pass

    def test_sortable_primary_order(self) -> bool:
        pass

    def test_stacked_violin_plot(self) -> bool:
        print("-- STACKED VIOLIN CHECKBOX SELECTION")
        try:
            stacked_violin_checkbox = self.browser.find_element(By.ID, "stacked_violin")
            stacked_violin_checkbox.click()
            return True
        except:
            return False

    def test_use_adjusted_pvals(self) -> bool:
        print("-- USE ADJUSTED P-VALUES CHECKBOX SELECTION")
        try:
            adj_pvals_checkbox = self.browser.find_element(By.ID, "adj_pvals")
            adj_pvals_checkbox.click()
            return True
        except:
            return False

    def test_violin_jitter(self) -> bool:
        print("-- JITTER POINTS CHECKBOX SELECTION")
        try:
            violin_add_points_checkbox = self.browser.find_element(By.ID, "violin_add_points")
            violin_add_points_checkbox.click()
            return True
        except:
            return False

    ### POST-PLOT OPTIONS

    def test_save_new_display(self) -> bool:
        pass

    def test_save_new_genecart(self) -> bool:
        pass