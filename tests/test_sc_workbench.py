#!/opt/bin/python3

"""
To run these tests, run `pytest <script>`.  It will run all tests with "test_" as the function name.

To run as localhost (to test on Docker images), pass in --data=localhost as a option after the script name.
"""

import configparser, random

from seleniumbase import BaseCase
import seleniumbase.common.exceptions
from selenium.webdriver.common.by import By

from dataclasses import dataclass, field

config = configparser.ConfigParser()
config.read('../gear.ini')

RED = 'rgb(255, 0, 0)'

# This just verifies the HTML tags are identical.  With Plotly IDs are randomized
VISUAL_LEVEL = 2

@dataclass(frozen=True)
class SCAnalysisPage:
    dataset: str = "P2, mouse, scRNA-seq, cochlea (Hertzano/Ament)"
    genes: list = field(default_factory=lambda: ["Pou4f3", "Pou4f1", "Pou4f2"])
    condition_cat: str = "cluster"  # Currently this dataset has only "cluster" and "louvain" which are the same set of groups
    deferred_errors: list = field(default_factory=list[str])

    query_compare_gene="Test1"  # Changed clustering labels to something I can reliably compare
    ref_compare_gene="Test2"

    def nav_to_url(self, sb):
        sb.get("http://localhost:8080/analyze_dataset.html" if sb.data == "localhost" else "https://umgear.org/analyze_dataset.html")

    def login(self, sb, user=config['test']['user_email'], pw=config['test']['password']):
        sb.type("#user_email", user)
        sb.type("#user_pass", pw)
        sb.click("#btn_sign_in")
        sb.sleep(1)

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
        sb.click("#dataset_id")    # Worth noting this is "dataset_id" on other pages
        sb.type("#dataset_tree_q", "P2, mouse, scRNA-seq, cochlea")
        sb.click(".jstree-search:contains('{}')".format(self.dataset))

class MultigeneTests(BaseCase):

    def test_valid_login(self):
        """Test that user can successfully log in."""
        sc = SCAnalysisPage()
        sc.nav_to_url(self)
        sc.login(self)
        self.assert_element_visible("#user_logged_in")

    def test_invalid_user_email(self):
        """Test that user can't log in with invalid user email."""
        sc = SCAnalysisPage()
        sc.nav_to_url(self)
        # Incorrect user and password
        sc.login(self, "slartibartfast@magrathea.com", "hangthesenseofit")
        # Cannot get color from CSS property since it wasn't computed.  Instead use style attribute string
        style = self.get_attribute("#user_email", "style")
        self.assert_true("color: {}".format(RED) in style, "User should not be logged in")

    def test_invalid_password(self):
        """Test that user can't log in with invalid password."""
        sc = SCAnalysisPage()
        sc.nav_to_url(self)
        # Correct user, incorrect password
        sc.login(self, config['test']['user_email'], "hangthesenseofit")
        style = self.get_attribute("#user_pass", "style")
        self.assert_true("color: {}".format(RED) in style, "User should have incorrect password")

    def test_from_scratch(self):
        """Test that the entire analysis can be completed from a "New" analysis."""
        sc = SCAnalysisPage()
        sc.nav_to_url(self)
        sc.login(self)
        sc.select_dataset(self)

        ### Prelim steps
        self.click_chain([
            "#filter_cells_lt_n_genes_selected"
            , "#filter_genes_lt_n_cells_selected"
            , "#btn_apply_primary_filter"
            ], timeout=15)
        self.assert_element("#primary_top_genes_c")

        ### QC by mitochondrial content
        self.click_chain([
            "#toggle_qc_by_mito + .toggle-group label"
            , "#btn_do_analysis_qc_by_mito"
            ])
        self.assert_elements([
            "#qbm_violin_c"
            , "#qbm_scatter_percent_mito_c"
            , "#qbm_scatter_n_genes_c"
            ])
        self.click("#btn_qbm_save")
        self.sleep(5)   # Give enough time to save

        ### Identify highly-variable genes
        self.click_chain([
            "#toggle_select_variable_genes + .toggle-group label"
            , "#btn_do_analysis_select_variable_genes"
            ])
        self.assert_element("#top_genes strong", timeout=10) # ? Should I check for the actual top genes? List may be random if there are ties.
        self.assert_element("#asvg_plot_c")
        self.click("#btn_asvg_save")
        self.sleep(5)   # Give enough time to save

        ### Principal Component Analysis (PCA)
        self.click_chain([
            "#toggle_pca + .toggle-group label"
            , "#btn_pca_run"
            ], timeout=10)
        self.assert_elements([
            "#pca_scatter_c"
            , "#pca_variance_c"
            ], timeout=10) # grayscale PCA since no genes searched
        self.type("#pca_genes_to_color", ",".join(sc.genes))
        self.click("#btn_pca_run")
        self.assert_elements([
            "#pca_scatter_c"
            , "#pca_variance_c"
            ])
        self.type("#top_pca_genes", "1,3,5")
        self.click("#btn_pca_top_genes")
        self.assert_element("#pca_top_genes_c")

        ### tSNE/UMAP
        self.click("#toggle_tsne + .toggle-group label")
        self.type("#tsne_genes_to_color", ",".join(sc.genes))
        self.type("#dredux_n_neighbors", "5")   # Keep these numbers low to speed up plot generation
        self.type("#tsne_n_pcs", "5")
        self.click_chain([
            "#dimensionality_reduction_method_tsne" # UMAP already clicked
            , "#btn_tsne_run"
            ])
        self.sleep(10)
        self.assert_elements([
            "#tsne_plot_c"
            , "#umap_plot_c"
            ], timeout=15)

        ### Clustering (Louvain)
        self.click("#toggle_louvain + .toggle-group label")
        self.type("#louvain_resolution", "0.5") # Setting this to a value where a smallish number of clusters arise
        self.click("#btn_louvain_run")
        self.assert_elements([
            "#louvain_tsne_plot_c"
            , "#louvain_umap_plot_c"
            ], timeout=10)

        ### Find Marker Genes
        self.click("#toggle_marker_genes + .toggle-group label")
        self.type("#marker_genes_n_genes", "5")
        self.click("#btn_marker_genes_run")
        self.assert_elements([
            "#marker_genes_plot_c"
            , "#marker_genes_table"
            ], timeout=10)
        self.click("#btn_download_marker_genes")
        # Currently there is no deferred assertion for downloaded file.  Make one.
        try:
            self.assert_downloaded_file("marker_genes.xls")
        except seleniumbase.common.exceptions.NoSuchFileException as e:
            log = "DEFERRED ASSERT\n\tAssert: {}\n\tError: {}\n".format('assert_downloaded_file("marker_genes.xls")', str(e))
            sc.deferred_errors.append(log)
        self.type("#marker_genes_manually_entered", ",".join(sc.genes))
        self.assert_text("3", "#marker_genes_entered_count")
        self.assert_text("3", "#marker_genes_unique_count")
        self.click("#marker_genes_table tr:nth-of-type(2) .js-row-idx")    # Select row 1 of table which selects all genes
        self.assert_text("11", "#marker_genes_selected_count")    # This assumes genes in row 1 do not change
        self.assert_text("14", "#marker_genes_unique_count")
        self.click("#btn_visualize_marker_genes")
        self.assert_elements([
            "#marker_genes_dotplot_c"
            , "#marker_genes_violin_c"
        ])

        ### Louvain Part 2
        # The labels tend to be inconsistent across runs (score ties and so on), so let's give them artificial names.
        for idx, elt in enumerate(self.find_elements("#marker_genes_group_labels tbody tr")):
            selector = "#marker_genes_group_labels tbody tr:nth-of-type({}) .group_user_label input".format(idx+1)
            self.type(selector, "Test{}".format(idx+1))
        self.click("#btn_louvain_rerun_with_groups")
        self.assert_elements([
            "#louvain_tsne_plot_c"
            , "#louvain_umap_plot_c"
            ])

        ### Compare genes and clusters
        self.click("#toggle_compare_genes + .toggle-group label")
        # NOTE: Since "query_cluster" and "ref_cluster" have options under an optgroup,
        # using self.select_option_by_value will not work.  Need to find element by xpath instead

        # 1 vs rest compare
        query_elt = self.find_element("#query_cluster")
        query_elt.find_element(by=By.XPATH, value="//optgroup/option[@value='{}']".format(sc.query_compare_gene)).click()
        self.select_option_by_value("#reference_cluster", "all-reference-clusters")
        self.click("#btn_compare_genes_run")
        self.assert_elements([
            "#compare_genes_ranked_c"
            , "#compare_genes_violin_c"
            ])
        self.click("#cg_download_table_f")
        # Currently there is no deferred assertion for downloaded file.  Make one.
        try:
            self.assert_downloaded_file("cluster_comparison_{}_vs_all-reference-clusters.xls".format(sc.query_compare_gene))
        except seleniumbase.common.exceptions.NoSuchFileException as e:
            log = "DEFERRED ASSERT\n\tAssert: {}\n\tError: {}\n".format('assert_downloaded_file("marker_genes.xls")', str(e))
            sc.deferred_errors.append(log)
        self.click("#cg_show_table_f")
        self.assert_element("#compare_genes_table_f")
        sc.plot_visual_regression(self, "cluster_1_vs_rest", VISUAL_LEVEL, True)
        # 1 vs 1 compare
        ref_elt = self.find_element("#reference_cluster")
        ref_elt.find_element(by=By.XPATH, value="//optgroup/option[@value='{}']".format(sc.ref_compare_gene)).click()
        self.click("#btn_compare_genes_run")
        self.assert_elements([
            "#compare_genes_ranked_c"
            , "#compare_genes_violin_c"
            ])
        sc.plot_visual_regression(self, "cluster_1_vs_1", VISUAL_LEVEL, True)
        # Changing method
        for val in ["t-test", "wilcoxon"]:
            self.select_option_by_value("#compare_genes_method", val)
            self.click("#btn_compare_genes_run")
            self.deferred_assert_element("#compare_genes_ranked_c")
            self.deferred_assert_element("#compare_genes_violin_c")
        self.select_option_by_value("#compare_genes_method", "t-test_overestim_var")    # Reset
        # Changing correction method
        self.select_option_by_value("#compare_genes_corr_method", "bonferroni")
        self.click("#btn_compare_genes_run")
        self.deferred_assert_element("#compare_genes_ranked_c")
        self.deferred_assert_element("#compare_genes_violin_c")

        #---
        # Custom deferred asserts we could not use built-in methods with.
        # Need to be added before "process_deferred_asserts since that raises the assertions"
        if sc.deferred_errors:
            print("\n".join(sc.deferred_errors))
            self.deferred_assert_element("#custom_failures_present")    # dummy variable to force failure
        self.process_deferred_asserts()

    def test_from_primary_analysis(self):
        """Test that a saved primary analysis can be resumed."""
        pass

    def test_from_saved_analysis(self):
        """Test that a saved analysis can be resumed."""
        pass

    def test_from_unsaved_analysis(self):
        """Test that an unsaved analysis can be resumed."""
        pass