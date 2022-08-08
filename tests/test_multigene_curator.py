#!/opt/bin/python3

"""
To run these tests, run `pytest <script>`.  It will run all tests with "test_" as the function name.

To run as localhost (to test on Docker images), pass in --data=localhost as a option after the script name.
"""

import configparser, random

from seleniumbase import BaseCase
from selenium.webdriver.common.by import By

from dataclasses import dataclass, field

config = configparser.ConfigParser()
config.read('../gear.ini')

RED = 'rgb(255, 0, 0)'
PALE_GREEN = 'rgb(195, 230, 203)'
GREY = 'rgb(128, 128, 128)'

@dataclass(frozen=True)
class MultigenePage:
    dataset: str = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
    genecart_to_load: str = "sadkins_savetest"
    genecart_to_save: str = "multigene_genecart_selenium_{}".format(random.randint(0, 999999))
    display_to_save: str = "multigene_display_selenium_{}".format(random.randint(0, 999999))
    genes: list = field(default_factory=lambda: ["Pou4f3", "Rfx7", "Sox2"])
    invalid_gene: str = "abc123"
    condition_cat: str = "cluster"  # Currently this dataset has only "cluster" and "louvain" which are the same set of groups
    filter_by: list = field(default_factory=lambda: ["HC (i)", "SC (i)", "TEC"])

    select2_container = "#select2-{}-container"
    select2_results = "#select2-{}-results"

    # TODO: Add test for checking indiviudal adata[obs, ensmbl_id] value and df[ensembl_id, obs] value to check for sort preservation
    # TODO: Add a similar test for checking the mean of a cell type for the same ensembl_id
    # TODO: Use a different dataset to check for secondary categories

    def nav_to_url(self, sb):
        sb.get("http://localhost:8080/multigene_curator.html" if sb.data == "localhost" else "https://umgear.org/multigene_curator.html")

    def login(self, sb, user=config['test']['user_email'], pw=config['test']['password']):
        sb.type("#user_email", user)
        sb.type("#user_pass", pw)
        sb.click("#btn_sign_in")
        sb.sleep(1)

    def click_by_select2_text(self, sb, select_id, text, optgroup=None):
        select2_results_id = self.select2_results.format(select_id)
        # This HTML fragment renders under the <scripts> tags when the Select2 dropdown is opened
        select2_list_elts = sb.find_visible_elements("{} .select2-results__option--selectable".format((select2_results_id)))
        for elt in select2_list_elts:
            # In cases where the text can appear multiple times in the dropdown, we need to check for the optgroup
            # The optgroup name is stored in the parents sibling element.
            if optgroup:
                parent = elt.find_element(by=By.XPATH, value="..")    # Should be <ul>
                sibling = parent.find_element(by=By.XPATH, value="../strong[@class='select2-results__group']")    # Should be <strong>
                if not sibling.text == optgroup:
                    continue
            if elt.text == text:
                elt.click()
                return

    def enter_genes(self, sb):
        gene_dropdown_select = "gene_dropdown"
        select2_gene_box = self.select2_container.format(gene_dropdown_select)
        select2_gene_textarea = "{} + span textarea".format(select2_gene_box)
        sb.click(select2_gene_textarea)
        for gene in self.genes:
            sb.add_text(select2_gene_textarea, gene)
            # Clicking "Enter" would add the first returned gene to the list, which may be incorrect (ex: search Sox2, get Qsox2)
            # Find the gene in the list and click on it
            # (I think there is a shorter way to do this with XPath but XPaths get confusinng in my opinion)
            self.click_by_select2_text(sb, gene_dropdown_select, gene)
        return select2_gene_box # return element "id"

    def enter_gene_carts(self, sb):
        sb.click("#gene_cart")
        sb.type("#gene_cart_tree_q", self.genecart_to_load)
        sb.click(".jstree-search:contains('{}')".format(self.genecart_to_load))
        return self.select2_container.format("gene_dropdown")   # return element "id"

    def filter_by_selection(self, sb):
        # First clear the select2 dropdown of all selected groups
        sb.click("#{}_clear".format(self.condition_cat))
        # Now start adding some.
        cluster_dropdown = "{}_dropdown".format(self.condition_cat)
        select2_filter_condition = self.select2_container.format(cluster_dropdown)
        select2_filter_textarea = "{} + span textarea".format(select2_filter_condition)
        sb.click(select2_filter_textarea)
        for cat in self.filter_by:
            sb.add_text(select2_filter_textarea, cat)
            self.click_by_select2_text(sb, cluster_dropdown, cat)
        return select2_filter_condition # return element "id"

    def plot_creation(self, sb):
        sb.click("#create_plot")
        sb.sleep(10)

    def plot_visual_regression(self, sb, img_name, level=3, deferred=False):
        """
        Test that a screenshot of a plot element is the same as the screenshot on disk

        :param img_name:
            Unique name parameter to establish or compare a baseline to

        Read https://seleniumbase.io/examples/visual_testing/ReadMe/ for more infomation about
        other processes that happen when this is called.
        """
        if deferred:
            sb.deferred_check_window(name=img_name, level=level)
            return
        sb.check_window(name=img_name, level=level)

    def select_dataset(self, sb):
        sb.wait_for_element_not_visible("#pre_dataset_spinner")   # Give time for the API to load datasets
        sb.sleep(2)  # A little more buffer room
        sb.click("#dataset")    # Worth noting this is "dataset_id" on other pages
        sb.type("#dataset_tree_q", "kelley")
        sb.click(".jstree-search:contains('{}')".format(self.dataset))

    def select_de_algo(self, sb, algo="Welch's t-test"):
        de_algo_select  = "de_test_select"
        select_id = self.select2_container.format(de_algo_select)
        sb.click(select_id)
        self.click_by_select2_text(sb, de_algo_select, algo)
        return sb.find_element(select_id)   # return Selenium WebElement

    def select_distance_metric(self, sb, metric="Euclidean"):
        distance_select  = "distance_select"
        select_id = self.select2_container.format(distance_select)
        sb.click(select_id)
        self.click_by_select2_text(sb, distance_select, metric)
        return sb.find_element(select_id)   # return Selenium WebElement

    def select_plot_type(self, sb, plot_type):
        plot_type_select = "plot_type_select"
        select_id = self.select2_container.format(plot_type_select)
        sb.click(select_id)
        self.click_by_select2_text(sb, plot_type_select, plot_type)
        return sb.find_element(select_id)   # return Selenium WebElement

    def select_quadrant_conditions(self, sb):
        compare1 = "quadrant_compare1_condition"
        compare1_container = self.select2_container.format(compare1)
        sb.click(compare1_container)
        self.click_by_select2_text(sb, compare1, "HC (i)", optgroup=self.condition_cat)

        compare2 = "quadrant_compare2_condition"
        compare2_container = self.select2_container.format(compare2)
        sb.click(compare2_container)
        self.click_by_select2_text(sb, compare2, "SC (i)", optgroup=self.condition_cat)

        ref = "quadrant_ref_condition"
        ref_container = self.select2_container.format(ref)
        sb.click(ref_container)
        self.click_by_select2_text(sb, ref, "TEC", optgroup=self.condition_cat)

    def select_volcano_conditions(self, sb):
        query = "volcano_query_condition"
        query_container = self.select2_container.format(query)
        sb.click(query_container)
        self.click_by_select2_text(sb, query, "HC (i)", optgroup=self.condition_cat)

        ref = "volcano_ref_condition"
        ref_container = self.select2_container.format(ref)
        sb.click(ref_container)
        self.click_by_select2_text(sb, ref, "SC (i)", optgroup=self.condition_cat)


    ### OPTIONS TESTING

    def add_clusterbar(self) -> bool:
        print(" -- CLUSTERBAR CATEGORY SELECTION")


    def sortable_primary_order(self):
        pass

    ### POST-PLOT OPTIONS

    def save_new_display(self):
        pass

    def save_new_genecart(self):
        pass


class MultigeneTests(BaseCase):

    def test_valid_login(self):
        """Test that user can successfully log in."""
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.login(self)
        self.assert_element_visible("#user_logged_in")

    def test_invalid_user_email(self):
        """Test that user can't log in with invalid user email."""
        mg = MultigenePage()
        mg.nav_to_url(self)
        # Incorrect user and password
        mg.login(self, "slartibartfast@magrathea.com", "hangthesenseofit")
        # Cannot get color from CSS property since it wasn't computed.  Instead use style attribute string
        style = self.get_attribute("#user_email", "style")
        self.assert_true("color: {}".format(RED) in style, "User should not be logged in")

    def test_invalid_password(self):
        """Test that user can't log in with invalid password."""
        mg = MultigenePage()
        mg.nav_to_url(self)
        # Correct user, incorrect password
        mg.login(self, config['test']['user_email'], "hangthesenseofit")
        style = self.get_attribute("#user_pass", "style")
        self.assert_true("color: {}".format(RED) in style, "User should have incorrect password")

    # TODO: Test colorpalette, and colorblindness
    # TODO: Test if no genes are selected for dotplot (bad), or volcano plot (good).  Test if 1 gene is selected for heatmap (bad)

    #! Currently does not work as genecart JSTree is not loaded after login
    def test_genecart_entry(self):
        self.skip(reason="Not working as genecart JSTree is not loaded after login")
        # ---
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.login(self)  # Need to login to retrieve gene cart
        mg.select_dataset(self)
        select2_gene_box = mg.enter_gene_carts(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))

    def test_filter_selection(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        # Don't need to worry about selecting genes
        select2_elt = mg.select_plot_type(self, "Heatmap")
        self.assert_true(select2_elt.text == "Heatmap", "Heatmap plot type should be selected")
        select2_filter_box = mg.filter_by_selection(self)
        select2_filter_box_elts = self.find_elements("{} .select2-selection__choice__display ".format(select2_filter_box))
        self.assert_true(len(select2_filter_box_elts), "Filter selection box should have {} conditions".format(len(mg.filter_by)))
        # ? Not sure where to go next with this

    def test_plot_heatmap(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        select2_gene_box = mg.enter_genes(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))
        select2_elt = mg.select_plot_type(self, "Heatmap")
        self.assert_true(select2_elt.text == "Heatmap", "Heatmap plot type should be selected")
        self.click("#matrixplot")   # Enabled by default so turn off
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        # Flip axes
        self.click("#flip_axes")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_flip_axes", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#flip_axes")
        # Change distance metric
        select2_elt = mg.select_distance_metric(self, "Manhattan (Cityblock)")
        self.assert_true(select2_elt.text == "Manhattan (Cityblock)", "Manhattan distance metric should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_distance_metric", 2, True)   # id values in plot get randomized, so stick with level 2
        select2_elt = mg.select_distance_metric(self, "Euclidean")
        self.assert_true(select2_elt.text == "Euclidean", "Euclidean distance metric should be selected")
        # Cluster genes
        self.click("#cluster_genes")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_cluster_genes", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#cluster_genes")
        # Cluster observations
        self.click("#cluster_obs")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_cluster_obs", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#cluster_obs")
        # Cluster both
        self.click("#cluster_genes")
        self.click("#cluster_obs")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_cluster_both", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#cluster_genes")
        self.click("#cluster_obs")
        # Sort by Primary Category
        self.click("#{}_primary".format(mg.condition_cat))
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_primary", 2, True)   # id values in plot get randomized, so stick with level 2
        # ---
        self.process_deferred_asserts()

    def test_plot_matrixplot(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        select2_gene_box = mg.enter_genes(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))
        select2_elt = mg.select_plot_type(self, "Heatmap")
        self.assert_true(select2_elt.text == "Heatmap", "Heatmap plot type should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_matrixplot", 2, True)   # id values in plot get randomized, so stick with level 2
        # Flip axes
        self.click("#flip_axes")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "matrixplot_flip_axes", 2, True)   # id values in plot get randomized, so stick with level 2
        # Change distance metric
        select2_elt = mg.select_distance_metric(self, "Manhattan (Cityblock)")
        self.assert_true(select2_elt.text == "Manhattan (Cityblock)", "Manhattan distance metric should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "matrixplot_distance_metric", 2, True)   # id values in plot get randomized, so stick with level 2
        select2_elt = mg.select_distance_metric(self, "Euclidean")
        self.assert_true(select2_elt.text == "Euclidean", "Euclidean distance metric should be selected")
        # Cluster genes
        self.click("#cluster_genes")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "matrixplot_cluster_genes", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#cluster_genes")
        # Cluster observations
        self.click("#cluster_obs")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "matrixplot_cluster_obs", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#cluster_obs")
        # Cluster both
        self.click("#cluster_genes")
        self.click("#cluster_obs")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "matrixplot_cluster_both", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#cluster_genes")
        self.click("#cluster_obs")
        # Sort by Primary Category
        self.click("#{}_primary".format(mg.condition_cat))
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "matrixplot_primary", 2, True)   # id values in plot get randomized, so stick with level 2
        # ---
        self.process_deferred_asserts()

    def test_plot_dotplot(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        select2_gene_box = mg.enter_genes(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))
        select2_elt = mg.select_plot_type(self, "Dotplot")
        self.assert_true(select2_elt.text == "Dotplot", "Dotplot plot type should be selected")
        # Primary Category is mandatory
        self.click("#{}_primary".format(mg.condition_cat))
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "dotplot_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        # ---
        self.process_deferred_asserts()

    def test_plot_quadrant(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        select2_gene_box = mg.enter_genes(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))
        select2_elt = mg.select_plot_type(self, "Quadrant")
        self.assert_true(select2_elt.text == "Quadrant", "Quadrant plot type should be selected")
        mg.select_quadrant_conditions(self)
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "quadrant_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        # Instead of checking for the presence of grey points, let's check the legend instead.  This should save time
        traces = self.find_visible_elements(".traces text")
        grey_found = False
        for t in traces:
            if "NONE/NONE" in t.text:
                grey_found = True
                break
        self.assert_true(grey_found, msg="No points were colored grey in standard quadrant plot.")
        # Test fold change cutoff
        self.click("#quadrant_foldchange_cutoff")
        self.type("#quadrant_foldchange_cutoff", "1.5")   # default value is 2
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "quadrant_foldchange_cutoff_adjusted", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#quadrant_foldchange_cutoff")
        self.type("#quadrant_foldchange_cutoff", "2")
        # Test excluding zero-foldchange genes (normally in gray)
        self.click("#include_zero_foldchange")  # disable
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "quadrant_zero_foldchange_off", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#include_zero_foldchange")  # enable
        traces = self.find_visible_elements(".traces text")
        grey_found = False
        for t in traces:
            if "NONE/NONE" in t.text:
                grey_found = True
                break
        self.assert_false(grey_found, msg="Found points that were colored grey when zero fold-change was disabled.")
        # Test FDR cutoff
        self.click("#quadrant_fdr_cutoff")
        self.type("#quadrant_fdr_cutoff", "0.1")   # default value is 0.05
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "quadrant_fdr_cutoff_adjusted", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#quadrant_fdr_cutoff")
        self.type("#quadrant_fdr_cutoff", "0.05")
        # Change DE algorithm
        select2_elt = mg.select_de_algo(self, "Wilcoxon rank-sum test")
        self.assert_true(select2_elt.text == "Wilcoxon rank-sum test", "Wilcoxon rank-sum test should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "quadrant_de_algo", 2, True)   # id values in plot get randomized, so stick with level 2
        select2_elt = mg.select_de_algo(self, "Welch's t-test")
        self.assert_true(select2_elt.text == "Welch's t-test", "Welch's t-test should be selected")
        # ---
        self.process_deferred_asserts()

    def test_plot_violin(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        select2_gene_box = mg.enter_genes(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))
        select2_elt = mg.select_plot_type(self, "Violin")
        self.assert_true(select2_elt.text == "Violin", "Violin plot type should be selected")
        # Primary Category is mandatory
        self.click("#{}_primary".format(mg.condition_cat))
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "violin_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        # Add jitter
        self.click("#violin_add_points")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "violin_jitter", 2, True)   # id values in plot get randomized, so stick with level 2
        # Sort by Primary Category
        self.click("#{}_primary".format(mg.condition_cat))
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "violin_primary", 2, True)   # id values in plot get randomized, so stick with level 2
        # ---
        self.process_deferred_asserts()

    def test_plot_stacked_violin(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        select2_gene_box = mg.enter_genes(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))
        select2_elt = mg.select_plot_type(self, "Violin")
        self.assert_true(select2_elt.text == "Violin", "Violin plot type should be selected")
        self.click("#stacked_violin")
        # Primary Category is mandatory
        self.click("#{}_primary".format(mg.condition_cat))
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "stacked_violin", 2, True)   # id values in plot get randomized, so stick with level 2
        # Add jitter
        self.click("#violin_add_points")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "stacked_violin_jitter", 2, True)   # id values in plot get randomized, so stick with level 2
        # Sort by Primary Category
        self.click("#{}_primary".format(mg.condition_cat))
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "stacked_violin_primary", 2, True)   # id values in plot get randomized, so stick with level 2
        # ---
        self.process_deferred_asserts()

    def test_plot_volcano(self):
        mg = MultigenePage()
        mg.nav_to_url(self)
        mg.select_dataset(self)
        select2_gene_box = mg.enter_genes(self)
        select2_gene_box_elts = self.find_elements("{} .select2-selection__choice__display".format(select2_gene_box))
        self.assert_true(len(select2_gene_box_elts), "Gene selection box should have {} genes".format(len(mg.genes)))
        select2_elt = mg.select_plot_type(self, "Volcano")
        self.assert_true(select2_elt.text == "Volcano", "Volcano plot type should be selected")
        mg.select_volcano_conditions(self)
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "volcano_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        # Disable annotation of nonsignificant genes
        self.click("#annot_nonsig")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "volcano_annot_nonsig_off", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#annot_nonsig")
        # Now disable adjusted p-values
        self.click("#adj_pvals")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "volcano_adj_pvals_off", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#adj_pvals")
        # Change DE algorithm
        select2_elt = mg.select_de_algo(self, "Wilcoxon rank-sum test")
        self.assert_true(select2_elt.text == "Wilcoxon rank-sum test", "Wilcoxon rank-sum test should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "volcano_de_algo", 2, True)   # id values in plot get randomized, so stick with level 2
        select2_elt = mg.select_de_algo(self, "Welch's t-test")
        self.assert_true(select2_elt.text == "Welch's t-test", "Welch's t-test should be selected")
        # ---
        self.process_deferred_asserts()