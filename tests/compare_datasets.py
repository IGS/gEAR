"""
To run these tests, run `pytest <script>`.  It will run all tests with "test_" as the function name.

To run as localhost (to test on Docker images), pass in --data=localhost as a option after the script name.
"""

from seleniumbase import BaseCase

from dataclasses import dataclass, field

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

    """

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

    """

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

    def test_create_standard_plot(self):
        """Tests that a basic plot can be created."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        cp.plot_creation(self)
        self.assert_element(".plotly")
        cp.plot_visual_regression(self, "normal_plot", 2)   # id values in plot get randomized, so stick with level 2

    def test_invalid_plot_duplicate_conditions(self):
        """Tests that a basic plot fails if X- and Y- conditions are the same."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_same_conditions(self)
        cp.plot_creation(self)
        self.assert_element("#error_loading_c")

    def test_conditions_are_preserved_upon_change(self):
        """Test that selected conditions are saved even when choosing the opposite axis' set of conditions."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_conditions(self)
        self.click("#condition_x_tab")
        self.assert_true(self.is_selected("input[data-group='{}']".format(cp.x_selection)), msg="'x' conditions were not saved after the switch.")
        self.click("#condition_y_tab")
        self.assert_true(self.is_selected("input[data-group='{}']".format(cp.y_selection)), msg="'y' conditions were not saved after the switch.")

    def test_axis_label_adjustment(self) -> None:
        """Test that the axes labels are adjusted when conditions are changed."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        x_label = cp.get_id_text(self, "x_label")
        y_label = cp.get_id_text(self, "y_label")
        cp.plot_creation(self)
        # Labels should reflect in the plot's axes
        # Assumes the plot axes titles are the first/only elements with this class
        x_title = cp.get_class_text(self, "xtitle")
        y_title = cp.get_class_text(self, "ytitle")
        self.assert_equal(x_label, x_title, msg="X-axis label was not adjusted when conditions were changed.")
        self.assert_equal(y_label, y_title, msg="Y-axis label was not adjusted when conditions were changed.")

    def test_significance_test_colorize(self) -> bool:
        """Test that a colorized plot is created when a significance test is performed with the "colorize" radio checked."""
        cp = ComparePage()
        cp.nav_to_url(self)
        cp.select_dataset(self)
        cp.select_significance_test(self, "t-test") # colorize should already be checked by default
        cp.plot_creation(self)
        points_fill = self.get_property_value(".points", "fill")
        red_points = list(filter(self.is_attribute_present("fill", "#FF0000"), points_fill))
        self.assert_true(len(red_points) > 0, msg="No points were colored red when a colorized significance test was performed.")
