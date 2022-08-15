#!/opt/bin/python3

"""
To run these tests, run `pytest <script>`.  It will run all tests with "test_" as the function name.

To run as localhost (to test on Docker images), pass in --data=localhost as a option after the script name.
"""

import configparser, random

from seleniumbase import BaseCase

from dataclasses import dataclass, field

config = configparser.ConfigParser()
config.read('../gear.ini')

RED = 'rgb(255, 0, 0)'

# This just verifies the HTML tags are identical.  With Plotly IDs are randomized
VISUAL_LEVEL = 2

@dataclass(frozen=True)
class FrontPage:
    genes: list = field(default_factory=lambda: ["Pou4f3", "Pou4f1", "Pou4f2"])
    profile: str = "Ear (diverse variety)"
    permalink_url = "{}/index.html?multigene_plots=0&gene_symbol_exact_match=1&gene_symbol=Sox2"
    pattern: str = "P2_cochlea_PCA"

    def nav_to_url(self, sb):
        sb.get("http://localhost:8080/" if sb.data == "localhost" else "https://umgear.org/")

    def nav_to_permalink(self, sb):
        sb.get(self.permalink_url.format("http://localhost:8080/" if sb.data == "localhost" else "https://umgear.org/"))

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

    def select_profile(self, sb):
        sb.sleep(2)  # A little more buffer room
        sb.click("#search_param_profile")    # Worth noting this is "dataset_id" on other pages
        sb.type("#profile_tree_q", self.profile)
        sb.click(".jstree-search:contains('{}')".format(self.profile))

    def select_results_profile(self, sb):
        sb.sleep(2)  # A little more buffer room
        sb.click("#selected_profile")    # Worth noting this is "dataset_id" on other pages
        sb.type("#selected_profile_tree_q", self.profile)
        sb.click(".jstree-search:contains('{}')".format(self.profile))

    def select_results_pattern(self, sb):
        sb.sleep(2)  # A little more buffer room
        sb.click("#projection_source")    # Worth noting this is "dataset_id" on other pages
        sb.type("#projection_source_tree_q", self.pattern)
        sb.click(".jstree-search:contains('{}')".format(self.pattern))

class FrontPageSearchTests(BaseCase):
    """Tests where the search occurs on the front page."""
    def test_valid_login(self):
        """Test that user can successfully log in."""
        fp = FrontPage()
        fp.nav_to_url(self)
        fp.login(self)
        self.assert_element_visible("#user_logged_in")

    def test_invalid_user_email(self):
        """Test that user can't log in with invalid user email."""
        fp = FrontPage()
        fp.nav_to_url(self)
        # Incorrect user and password
        fp.login(self, "slartibartfast@magrathea.com", "hangthesenseofit")
        # Cannot get color from CSS property since it wasn't computed.  Instead use style attribute string
        style = self.get_attribute("#user_email", "style")
        self.assert_true("color: {}".format(RED) in style, "User should not be logged in")

    def test_invalid_password(self):
        """Test that user can't log in with invalid password."""
        fp = FrontPage()
        fp.nav_to_url(self)
        # Correct user, incorrect password
        fp.login(self, config['test']['user_email'], "hangthesenseofit")
        style = self.get_attribute("#user_pass", "style")
        self.assert_true("color: {}".format(RED) in style, "User should have incorrect password")

    def test_index_search_exact_genes(self):
        fp = FrontPage()
        fp.nav_to_url(self)
        self.type("#search_gene_symbol_intro", ",".join(fp.genes))
        fp.select_profile(self)
        # Exact match is already checked by default
        self.click("#intro_search_icon")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_element("#links_out_c")
        self.assert_element_not_visible("#functional_not_supported_alert")
        self.assert_text("3", "#search_result_count")
        # TODO: Test link-outs, test links in dataset panel box, test selecting a different gene

    def test_index_search_inexact_genes(self):
        fp = FrontPage()
        fp.nav_to_url(self)
        self.type("#search_gene_symbol_intro", "Pou4f")
        fp.select_profile(self)
        # Turn off exact match
        self.click("#exact_match_icon")
        self.click("#intro_search_icon")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_element("#links_out_c")
        self.assert_element_not_visible("#functional_not_supported_alert")
        self.assert_text("3", "#search_result_count")

    def test_index_search_multigenes(self):
        fp = FrontPage()
        fp.nav_to_url(self)
        self.type("#search_gene_symbol_intro", ",".join(fp.genes))
        fp.select_profile(self)
        # Exact match is already checked by default
        self.click("#multigene_search_icon")
        self.click("#intro_search_icon")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_element_not_visible("#links_out_c")
        self.assert_element("#functional_not_supported_alert")
        self.assert_element_not_visible("#search_result_count")

    def test_projection(self):
        self.skip(reason="ProjectR not merged into 'main' branch")
        fp = FrontPage()
        fp.nav_to_url(self)
        self.click("[data-tool-name='projection']")
        fp.select_results_pattern(self)
        fp.select_results_profile(self)
        # Single pattern is already checked by default
        # PCA is already checked by default
        self.click("submit_search_projection")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_element_not_visible("#links_out_c")
        self.assert_element("#functional_not_supported_alert")
        self.assert_text("50", "#search_result_count")
        # NMF
        self.click("#nmf_algo")
        self.click("submit_search_projection")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        # Multi pattern
        self.click("#pca_algo")
        self.click("#multi_pattern")
        self.click("#projection_pattern_deselect_all")
        self.click(".js-projection-pattern-elts-check:nth-of-type(-n+10)")  # Select 10 patterns (n starts at 0)
        self.click("submit_search_projection")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible

class ResultsPageSearchTests(BaseCase):
    """Tests where the search occurs on the results page."""

    def test_index_search_exact_genes(self):
        fp = FrontPage()
        fp.nav_to_permalink(self)
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_text("1", "#search_result_count")
        self.type("#search_gene_symbol", ",".join(fp.genes))
        fp.select_results_profile(self)
        # Exact match is already checked by default
        self.click("#submit_search")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_element("#links_out_c")
        self.assert_element_not_visible("#functional_not_supported_alert")
        self.assert_text("3", "#search_result_count")

    def test_index_search_exact_genes(self):
        fp = FrontPage()
        fp.nav_to_permalink(self)
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_text("1", "#search_result_count")
        self.type("#search_gene_symbol", "Pou4f")
        fp.select_results_profile(self)
        # Turn off exact match
        self.click("#exact_match_input")
        self.click("#submit_search")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_element("#links_out_c")
        self.assert_element_not_visible("#functional_not_supported_alert")
        self.assert_text("3", "#search_result_count")

    def test_index_search_multigenes(self):
        fp = FrontPage()
        fp.nav_to_permalink(self)
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_text("1", "#search_result_count")
        self.type("#search_gene_symbol", ",".join(fp.genes))
        fp.select_results_profile(self)
        # Exact match is already checked by default
        self.click("#multigene_plots_input")
        self.click("#submit_search")
        self.assert_element(".js-plotly-plot")  # At least one plot is visible
        self.assert_element_not_visible("#links_out_c")
        self.assert_element("#functional_not_supported_alert")
        self.assert_element_not_visible("#search_result_count")
