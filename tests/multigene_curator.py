#!/opt/bin/python3

"""
To run these tests, run `pytest <script>`.  It will run all tests with "test_" as the function name.

To run as localhost (to test on Docker images), pass in --data=localhost as a option after the script name.
"""

import configparser, random

from seleniumbase import BaseCase
from selenium.webdriver.common.keys import Keys

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
    condition_cat: str = "cluster"
    filter_by: list = field(default_factory=lambda: ["HC (i)", "SC (i)", "TEC"])

    select2_container = "#select2-{}-container"
    select2_results = "#select2-{}-results"

    # TODO: Add test for checking indiviudal adata[obs, ensmbl_id] value and df[ensembl_id, obs] value to check for sort preservation
    # TODO: Add a similar test for checking the mean of a cell type for the same ensembl_id

    def nav_to_url(self, sb):
        sb.get("http://localhost:8080/multigene_curator.html" if sb.data == "localhost" else "https://umgear.org/multigene_curator.html")

    def login(self, sb, user=config['test']['user_email'], pw=config['test']['password']):
        sb.type("#user_email", user)
        sb.type("#user_pass", pw)
        sb.click("#btn_sign_in")
        sb.sleep(1)

    def enter_genes(self, sb):
        select2_gene_box = self.select2_container.format("gene_dropdown")
        select2_gene_textarea = "{} + span textarea".format(select2_gene_box)
        sb.click(select2_gene_textarea)
        for gene in self.genes:
            sb.add_text(select2_gene_textarea, gene)
            sb.add_text(select2_gene_textarea, Keys.ENTER)
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
        select2_filter_condition = self.select2_container.format(self.condition_cat)
        select2_filter_textarea = "{} + span textarea".format(select2_filter_condition)
        sb.click(select2_filter_textarea)
        for cat in self.filter_by:
            sb.add_text(select2_filter_textarea, cat)
        return sb.find_element(select2_filter_condition)

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
            sb.deferred_check_window("#{}".format(img_name))
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
        select2_list = self.select2_results.format(de_algo_select)
        select2_list_elts = sb.find_elements("{} .select2-results__option--selectable".format(select2_list))
        for elt in select2_list_elts:
            if elt.text == algo:
                elt.click()
                break
        return sb.find_element(select_id)   # return Selenium WebElement

    def select_distance_metric(self, sb, metric="Euclidean"):
        distance_select  = "distance_select"
        select_id = self.select2_container.format(distance_select)
        sb.click(select_id)
        select2_list = self.select2_results.format(distance_select)
        select2_list_elts = sb.find_elements("{} .select2-results__option--selectable".format(select2_list))
        for elt in select2_list_elts:
            if elt.text == metric:
                elt.click()
                break
        return sb.find_element(select_id)   # return Selenium WebElement

    def select_plot_type(self, sb, plot_type):
        plot_type_select = "plot_type_select"
        select_id = self.select2_container.format(plot_type_select)
        sb.click(select_id)
        select2_list = self.select2_results.format(plot_type_select)
        select2_list_elts = sb.find_elements("{} .select2-results__option--selectable".format(select2_list))
        for elt in select2_list_elts:
            if elt.text == plot_type:
                elt.click()
                break
        return sb.find_element(select_id)   # return Selenium WebElement

    def add_quadrant_group_selection(self) -> bool:
        pass

    def add_volcano_group_selection(self) -> bool:
        pass

    ### OPTIONS TESTING

    def add_annotate_nonsignficant(self) -> bool:
        print("-- ANNOTATE NONSIGNIFICANT POINTS CHECKBOX SELECTION")
        try:
            annot_nonsig_checkbox = self.browser.find_element(By.ID, "annot_nonsig")
            annot_nonsig_checkbox.click()
            return True
        except:
            return False

    def add_clusterbar(self) -> bool:
        print(" -- CLUSTERBAR CATEGORY SELECTION")

    def add_cluster_observations(self) -> bool:
        print("-- CLUSTER OBSERVATIONS CHECKBOX SELECTION")
        try:
            cluster_obs_checkbox = self.browser.find_element(By.ID, "cluster_obs")
            cluster_obs_checkbox.click()
            return True
        except:
            return False

    def primary_category(self) -> bool:
        pass

    def secondary_category(self) -> bool:
        pass

    def sortable_primary_order(self) -> bool:
        pass

    ### POST-PLOT OPTIONS

    def save_new_display(self) -> bool:
        pass

    def save_new_genecart(self) -> bool:
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

    def test_genecart_entry(self):
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
        self.assert_true(select2_elt.text == "Manhattan", "Manhattan distance metric should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "heatmap_distance_metric", 2, True)   # id values in plot get randomized, so stick with level 2
        select2_elt = mg.select_distance_metric(self, "Euclidean")
        self.assert_text(select2_elt.text == "Euclidean", "Euclidean distance metric should be selected")
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
        self.click("#matrixplot")
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
        self.assert_true(select2_elt.text == "Manhattan", "Manhattan distance metric should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "matrixplot_distance_metric", 2, True)   # id values in plot get randomized, so stick with level 2
        select2_elt = mg.select_distance_metric(self, "Euclidean")
        self.assert_text(select2_elt.text == "Euclidean", "Euclidean distance metric should be selected")
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
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        # Test fold change cutoff
        mg.plot_visual_regression(self, "quadrant_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        self.type("#quadrant_foldchange_cutoff", "1.5")   # default value is 2
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "quadrant_foldchange_cutoff", 2, True)   # id values in plot get randomized, so stick with level 2
        self.type("#quadrant_foldchange_cutoff", "2")
        # Test including zero-foldchange genes (in gray)
        self.click("#include_zero_foldchange")  # enable
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        points = self.find_visible_elements(".points path")
        # The "fill" property is not computed but hardcoded by Plotly into the "style" attribute
        # Unfortunately this will take a while if none of the points are colored grey
        grey_found = False
        for p in points:
            p_style = p.get_attribute("style")
            if "fill: {}".format(GREY) in p_style:
                grey_found = True
                break
        self.assert_true(grey_found, msg="No points were colored grey when zero fold-change was allowed.")
        mg.plot_visual_regression(self, "quadrant_zero_foldchange", 2, True)   # id values in plot get randomized, so stick with level 2
        self.click("#include_zero_foldchange")  # enable
        # Test FDR cutoff
        self.type("#quadrant_fdr_cutoff", "0.1")   # default value is 0.05
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "quadrant_fdr_cutoff", 2, True)   # id values in plot get randomized, so stick with level 2
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
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "violin_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        # Add jitter
        self.click("#violin_add_points")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "violin_jitter", 2, True)   # id values in plot get randomized, so stick with level 2
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
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "stacked_violin", 2, True)   # id values in plot get randomized, so stick with level 2
        # Add jitter
        self.click("#violin_add_points")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "stacked_violin_jitter", 2, True)   # id values in plot get randomized, so stick with level 2
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
        print(select2_elt.text)
        self.assert_true(select2_elt.text == "Volcano", "Volcano plot type should be selected")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "volcano_no_opts", 2, True)   # id values in plot get randomized, so stick with level 2
        # Now use adjusted p-values
        self.click("#adj_pvals")
        mg.plot_creation(self)
        self.deferred_assert_element(".plotly")
        mg.plot_visual_regression(self, "volcano_adj_pvals", 2, True)   # id values in plot get randomized, so stick with level 2
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