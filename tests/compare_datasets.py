"""
To run these tests, run `pytest <script>`.  It will run all tests with "test_" as the function name.

To run as localhost (to test on Docker images), pass in --data=localhost as a option after the script name.
"""

import configparser

from seleniumbase import BaseCase

from dataclasses import dataclass, field

config = configparser.ConfigParser()
config.read('../gear.ini')

@dataclass(frozen=True)
class ComparePage:
    dataset: str = "P1, mouse, scRNA-seq, utricle, hair cells, supporting cells, and transitional epithelial cells (Kelley)"
    genecart_to_save: str = "sadkins_selenium"
    genes: list = field(default_factory=lambda: ["Pou4f3", "Rfx7", "Sox2"])
    invalid_gene: str = "abc123"
    condition_cat: str = "cluster"
    x_group: str = "HC (i)"
    y_group: str = "TEC"
    x_selection: str = "{};-;{}".format(condition_cat, x_group)
    y_selection: str = "{};-;{}".format(condition_cat, y_group)

    # sb => SeleniumBase "test" class object

    def get_class_text(self, sb, class_name):
        return sb.get_text("." + class_name)

    def get_id_text(self, sb, label_id):
        return sb.get_text("#" + label_id)

    def nav_to_url(self, sb):
        sb.get("http://localhost:8080/compare_datasets.html" if sb.data == "localhost" else "https://umgear.org/compare_datasets.html")

    def highlighted_genes_in_table_highlighted(self, sb):
        pass

    def plot_creation(self, sb):
        sb.click("#btn_apply_dataset_changes")

    def select_conditions(self, sb):
        sb.click("#condition_x_tab")
        sb.click("#" + "{}_collapse".format(self.condition_cat))
        sb.click("input[data-group='{}']".format(self.x_selection))
        sb.click("#condition_y_tab")
        sb.click("#" + "{}_collapse".format(self.condition_cat))
        sb.click("input[data-group='{}']".format(self.y_selection))

    def select_same_conditions(self, sb):
        sb.click("#condition_x_tab")
        sb.click("#" + "{}_collapse".format(self.condition_cat))
        sb.click("input[data-group='{}']".format(self.x_selection))
        sb.click("#condition_y_tab")
        sb.click("#" + "{}_collapse".format(self.condition_cat))
        sb.click("input[data-group='{}']".format(self.x_selection))

    def select_dataset(self, sb):
        #self.wait_for_element_not_visible("#pre_dataset_spinner")   # Give time for the API to load datasets
        sb.click("#dataset_id")
        sb.type("#dataset_tree_q", "kelley")
        sb.click(".jstree-search:contains('{}')".format(self.dataset))

    def select_significance_test(self, sb, select_value):
        sb.select_option_by_value("#statistical_test", select_value)

    def selected_genes_in_table(self, sb):
        pass

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
        self.send_keys("#user_email", config['test']['user_email'])
        self.send_keys("#user_pass", config['test']['password'])
        self.click("#btn_sign_in")
        self.assert_element_visible("#user_logged_in", "User should be logged in")

    def no_test_invalid_login(self):
        """Test that user can't log in with invalid credentials."""
        cp = ComparePage()
        cp.nav_to_url(self)
        # Incorrect user and password
        self.send_keys("#user_email", "slartibartfast@magrathea.com")
        self.send_keys("#user_pass", "hangthesenseofit")
        self.click("#btn_sign_in")
        RED = 'rgba(255, 0, 0, 1)'
        color = self.get_property_value()("#user_email", "color")
        self.assert_true(color == RED, "User should not be logged in")
        # Correct user, incorrect password
        self.send_keys("#user_email", config['test']['user_email'])
        self.click("#btn_sign_in")
        color = self.get_property_value()("#user_pass", "color")
        self.assert_true(color == RED, "User should have incorrect password")

    def no_test_create_standard_plot(self):
        """Tests that a basic plot can be created."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        cp.plot_creation(self)
        self.assert_element(".plotly")
        cp.plot_visual_regression(self, "normal_plot", 2)   # id values in plot get randomized, so stick with level 2

    def no_test_invalid_plot_duplicate_conditions(self):
        """Tests that a basic plot fails if X- and Y- conditions are the same."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_same_conditions(self)
        cp.plot_creation(self)
        self.assert_element("#error_loading_c")

    def no_test_conditions_are_preserved_upon_change(self):
        """Test that selected conditions are saved even when choosing the opposite axis' set of conditions."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        self.click("#condition_x_tab")
        self.assert_true(self.is_selected("input[data-group='{}']".format(cp.x_selection)), msg="'x' conditions were not saved after the switch.")
        self.click("#condition_y_tab")
        self.assert_true(self.is_selected("input[data-group='{}']".format(cp.y_selection)), msg="'y' conditions were not saved after the switch.")

    def no_test_axis_label_adjustment(self):
        """Test that the axes labels are adjusted when conditions are changed."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        x_label = cp.get_id_text(self, "x_label")
        y_label = cp.get_id_text(self, "y_label")
        cp.plot_creation(self)
        # Labels should reflect in the plot's axes
        # Assumes the plot axes titles are the first/only elements with this class
        x_title = cp.get_class_text(self, "xtitle")
        y_title = cp.get_class_text(self, "ytitle")
        self.assert_equal(x_label, x_title, msg="X-axis label was not adjusted when conditions were changed.")
        self.assert_equal(y_label, y_title, msg="Y-axis label was not adjusted when conditions were changed.")

    def no_test_significance_test_colorize(self):
        """Test that a colorized plot is created when a significance test is performed with the "colorize" radio checked."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        cp.select_significance_test(self, "t-test") # colorize should already be checked by default
        cp.plot_creation(self)
        points_fill = self.get_property_value(".points", "fill")
        red_points = list(filter(lambda f: f == "#FF0000", points_fill))
        self.assert_true(len(red_points), msg="No points were colored red when a colorized significance test was performed.")


    def no_test_significance_test_filter(self):
        """Test that a colorized plot is created when a significance test is performed with the "filter" radio checked."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        cp.select_significance_test(self, "t-test") # colorize should already be checked by default
        cp.plot_creation(self)
        points = self.get_element((".points"))
        num_colorized_points = len(points)
        self.click("#stat_action_filter")
        cp.plot_creation(self)
        points = self.get_element(".points")
        num_filtered_points = len(points)
        self.assert_not_equal(num_colorized_points, num_filtered_points, msg="Points were not filtered when a filtered significance test was performed.")

    def no_test_found_gene_highlighting(self):
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        self.send_keys('#highlighted_genes', ', '.join(self.genes))
        cp.plot_creation(self)
        # All three genes are in this plot
        gene_annotation = self.get_element('[data-unformatted="{}"'.format(self.genes[0]))
        self.assert_true(len(gene_annotation) > 0, msg="Genes that should have been in plot were not found in the plot.")

    def no_test_notfound_gene_highlighting(self):
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        self.send_keys('#highlighted_genes', ', '.join(self.invalid_gene))
        cp.plot_creation(self)
        gene_not_found = cp.get_id_text(self, "genes_not_found")
        self.assertTrue((gene_not_found and self.invalid_gene in gene_not_found), msg="Genes that should not have been in plot were found in the plot.")

    def no_test_genes_in_selection_table(self):
        """Test that highlighted genes are in the selection table."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        self.send_keys('#highlighted_genes', ', '.join(self.genes))
        cp.plot_creation(self)
        # Select some points
        self.click("[data-title='Box Select']")
        # Drag across all points (start in upper left of chart and grab all points 300px right and down)
        # This should hopefully grab the Pou4f3 gene.
        self.drag_and_drop_with_offset(".points", 300, 300)
        self.assert_element_visible('#tbl_selected_genes', msg="Selected genes table was not visible.")
        table_elts = self.get_elements('#tbl_selected_genes td.text_truncate')
        table_genes = self.get_text(table_elts)
        self.assertTrue(self.genes[0] in table_genes, msg="Genes that should have been in table were not found in the table.")
        bg_color = self.get_property_value(table_elts, "background-color")
        green_bg = list(filter(lambda c: c == "#c3e6cb", bg_color))
        self.assert_true(len(green_bg), msg="Highlighted genes in table were not colored appropriately.")

    def no_test_download_table(self):
        """Test that gene selection table can be downloaded."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        cp.plot_creation(self)
        # Select some points
        self.click("[data-title='Box Select']")
        # Drag across all points (start in upper left of chart and grab all points 300px right and down)
        self.drag_and_drop_with_offset(".points", 300, 300)
        self.assert_element_visible('#tbl_selected_genes', msg="Selected genes table was not visible.")
        self.click("#download_gene_table")
        self.assert_downloaded_file("selected_genes.tsv")
        self.delete_downloaded_file_if_present("selected_genes.tsv")

    def no_test_save_gene_cart(self):
        """Test that a gene cart can be named and saved."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        cp.plot_creation(self)
        # Select some points
        self.click("[data-title='Box Select']")
        # Drag across all points (start in upper left of chart and grab all points 300px right and down)
        self.drag_and_drop_with_offset(".points", 300, 300)
        self.assert_element_visible('#tbl_selected_genes', msg="Selected genes table was not visible.")
        self.click("#create_gene_cart")
        # Name cart
        import random
        random_int = random.randint(0, 999999)
        cart_name = "Test Cart {}".format(random_int)
        self.send_keys("#gene_cart_name", cart_name)
        # Save cart
        self.click("#save_gene_cart")
        self.assert_element_visible('#gene_cart_member_count', msg="Gene cart was not saved.")
        # TODO: Ensure cart is in the geardb database
