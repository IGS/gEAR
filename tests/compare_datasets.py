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
PALE_GREEN = 'rgb(195, 230, 203)'

@dataclass(frozen=True)
class ComparePage:
    dataset: str = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
    genecart_to_save: str = "compare_selenium_{}".format(random.randint(0, 999999))
    genes: list = field(default_factory=lambda: ["Pou4f3", "Rfx7", "Sox2"])
    invalid_gene: str = "abc123"
    condition_cat: str = "cluster"
    x_group: str = "HC (i)"
    y_group: str = "TEC"
    x_selection: str = "{};-;{}".format(condition_cat, x_group)
    y_selection: str = "{};-;{}".format(condition_cat, y_group)

    # sb => SeleniumBase "test" class object

    def nav_to_url(self, sb):
        sb.get("http://localhost:8080/compare_datasets.html" if sb.data == "localhost" else "https://umgear.org/compare_datasets.html")

    def login(self, sb, user=config['test']['user_email'], pw=config['test']['password']):
        sb.type("#user_email", user)
        sb.type("#user_pass", pw)
        sb.click("#btn_sign_in")
        sb.sleep(1)

    def plot_creation(self, sb):
        sb.click("#btn_apply_dataset_changes")  # I believe the default timeout will be enough to generate the plot
        sb.sleep(10)

    def select_conditions(self, sb, condition_category, x_selection, y_selection):
        sb.click("#condition_x_tab")
        sb.click("#" + "{}_collapse".format(condition_category))
        sb.click("input[data-group='{}']".format(x_selection))
        sb.click("#condition_y_tab")
        #sb.click("#" + "{}_collapse".format(condition_category))   # Should already be expanded
        sb.click("input[data-group='{}']".format(y_selection))
        sb.sleep(1)

    def select_dataset(self, sb):
        sb.wait_for_element_not_visible("#pre_dataset_spinner")   # Give time for the API to load datasets
        sb.sleep(2)  # A little more buffer room
        sb.click("#dataset_id")
        sb.type("#dataset_tree_q", "kelley")
        sb.click(".jstree-search:contains('{}')".format(self.dataset))

    def select_genes_in_plot(self, sb):
        # In order to make the modebar options visible, we need to hover over the plot
        sb.hover_and_click(".plotly", ".modebar-btn[data-title='Box Select']")
        # Use drag and drop to select the genes (x to right, y up, start in middle of plot)
        sb.drag_and_drop_with_offset(".nsewdrag", 400, -400)
        # For now, I think it's easier to select all genes with a double click (but the table will be larger)
        #sb.double_click(".nsewdrag")

    def select_significance_test(self, sb, select_value):
        sb.select_option_by_value("#statistical_test", select_value)

    def plot_visual_regression(self, sb, img_name, level=3):
        """
        Test that a screenshot of a plot element is the same as the screenshot on disk

        :param img_name:
            Unique name parameter to establish or compare a baseline to

        Read https://seleniumbase.io/examples/visual_testing/ReadMe/ for more infomation about
        other processes that happen when this is called.
        """
        sb.check_window(name=img_name, level=level)


class CompareTests(BaseCase):

    def test_valid_login(self):
        """Test that user can successfully log in."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.login(self)
        self.assert_element_visible("#user_logged_in")

    def test_invalid_user_email(self):
        """Test that user can't log in with invalid user email."""
        cp = ComparePage()
        cp.nav_to_url(self)
        # Incorrect user and password
        cp.login(self, "slartibartfast@magrathea.com", "hangthesenseofit")
        # Cannot get color from CSS property since it wasn't computed.  Instead use style attribute string
        style = self.get_attribute("#user_email", "style")
        self.assert_true("color: {}".format(RED) in style, "User should not be logged in")

    def test_invalid_password(self):
        """Test that user can't log in with invalid password."""
        cp = ComparePage()
        cp.nav_to_url(self)
        # Correct user, incorrect password
        cp.login(self, config['test']['user_email'], "hangthesenseofit")
        style = self.get_attribute("#user_pass", "style")
        self.assert_true("color: {}".format(RED) in style, "User should have incorrect password")

    def test_create_standard_plot(self):
        """Tests that a basic plot can be created."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        self.assert_true(self.get_text("#dataset_id") == cp.dataset, "Dataset should be {}".format(cp.dataset))
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        cp.plot_creation(self)
        self.assert_element(".plotly")
        cp.plot_visual_regression(self, "normal_plot", 2)   # id values in plot get randomized, so stick with level 2

    def test_invalid_plot_duplicate_conditions(self):
        """Tests that a basic plot fails if X- and Y- conditions are the same."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.x_selection)
        cp.plot_creation(self)
        self.assert_element("#error_loading_c")

    def test_conditions_are_preserved_upon_change(self):
        """Test that selected conditions are saved even when choosing the opposite axis' set of conditions."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        self.click("#condition_x_tab")
        x_selected = self.is_selected("input[data-group='{}']".format(cp.x_selection))
        self.assert_true(x_selected, msg="'x' conditions were not saved after the switch.")
        self.click("#condition_y_tab")
        y_selected = self.is_selected("input[data-group='{}']".format(cp.y_selection))
        self.assert_true(y_selected, msg="'y' conditions were not saved after the switch.")

    def test_axis_label_adjustment(self):
        """Test that the axes labels are adjusted when conditions are changed."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        # "input" tags have values, not text.
        x_label = self.get_value("#x_label")
        y_label = self.get_value("#y_label")
        self.assert_equal(x_label, cp.x_group, msg="X-axis label was not adjusted correctly.")
        self.assert_equal(y_label, cp.y_group, msg="Y-axis label was not adjusted correctly.")
        cp.plot_creation(self)
        # Labels should reflect in the plot's axes
        # Assumes the plot axes titles are the first/only elements with this class
        x_title = self.get_text(".xtitle")
        y_title = self.get_text(".ytitle")
        self.assert_equal(x_label, x_title, msg="X-axis title was not adjusted in plot.")
        self.assert_equal(y_label, y_title, msg="Y-axis title was not adjusted in plot.")

    def test_significance_test_colorize(self):
        """Test that a colorized plot is created when a significance test is performed with the "colorize" radio checked."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        cp.select_significance_test(self, "t-test") # colorize should already be checked by default
        cp.plot_creation(self)
        points = self.find_visible_elements(".points path")
        # The "fill" property is not computed but hardcoded by Plotly into the "style" attribute
        # Unfortunately this will take a while if none of the points are colored
        red_found = False
        for p in points:
            p_style = p.get_attribute("style")
            if "fill: {}".format(RED) in p_style:
                red_found = True
                break
        self.assert_true(red_found, msg="No points were colored red when a colorized significance test was performed.")

    def test_significance_test_filter(self):
        """Test that a colorized plot is created when a significance test is performed with the "filter" radio checked."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        cp.plot_creation(self)
        points = self.find_visible_elements(".points path")
        num_unfiltered_points = len(points)
        # Get filtered points
        cp.select_significance_test(self, "t-test")
        # NOTE: The "custom-control-input" class on the radio button seems to make it not visible (z-index: -1)
        # Initially I was going to use the normal "form-check-input" class, but I found that clicking the label works
        self.click("#stat_action_filter_label")
        cp.plot_creation(self)
        points = self.find_visible_elements(".points path")
        num_filtered_points = len(points)
        self.assert_not_equal(num_unfiltered_points, num_filtered_points, msg="Points were not filtered when a filtered significance test was performed.")

    def test_found_gene_highlighting(self):
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        self.type('#highlighted_genes', ', '.join(cp.genes))
        cp.plot_creation(self)
        # All three genes are in this plot
        gene_annotation = self.find_visible_elements('[data-unformatted="{}"'.format(cp.genes[0]))
        self.assert_true(len(gene_annotation), msg="Genes that should have been in plot were not found in the plot.")

    def test_notfound_gene_highlighting(self):
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        self.type('#highlighted_genes', cp.invalid_gene)
        cp.plot_creation(self)
        gene_not_found = self.get_text('#genes_not_found')
        self.assert_true((gene_not_found and cp.invalid_gene in gene_not_found), msg="Genes that should not have been in plot were found in the plot.")

    def test_genes_in_selection_table(self):
        """Test that highlighted genes are in the selection table."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        self.type('#highlighted_genes', ', '.join(cp.genes))
        cp.plot_creation(self)
        cp.select_genes_in_plot(self)
        self.assert_element_visible('#tbl_selected_genes')
        # Was highlighted gene in selection table?
        table_elts = self.find_visible_elements('#tbl_selected_genes td.text-truncate')
        highlighted_gene_found = False
        for e in table_elts:
            if e.text in cp.genes:
                highlighted_gene_found = True
                break
        self.assert_true(highlighted_gene_found, msg="Genes that should have been in table were not found in the table.")
        # Was green text found in table for a highlighted gene? Applied to .table-success class
        success_elts = self.find_visible_elements('#tbl_selected_genes tr.table-success')
        self.assert_true(len(success_elts), msg="Highlighted genes in table were not colored appropriately.")

    def test_download_table(self):
        """Test that gene selection table can be downloaded."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        cp.plot_creation(self)
        cp.select_genes_in_plot(self)
        self.assert_element_visible('#tbl_selected_genes')
        self.click("#download_gene_table")
        self.assert_downloaded_file("selected_genes.tsv")
        self.delete_downloaded_file_if_present("selected_genes.tsv")

    def test_save_gene_cart(self):
        """Test that a gene cart can be named and saved."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.login(self)  # Need to login to do this
        cp.select_dataset(self)
        cp.select_conditions(self, cp.condition_cat, cp.x_selection, cp.y_selection)
        cp.plot_creation(self)
        cp.select_genes_in_plot(self)
        self.assert_element_visible('#tbl_selected_genes')
        self.click("#create_gene_cart")
        # Name cart
        self.type("#gene_cart_name", cp.genecart_to_save)
        # Save cart
        self.click("#save_gene_cart")
        self.assert_element_visible('#gene_cart_member_count')
        # TODO: Ensure cart is in the geardb database
