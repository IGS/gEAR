"""
These are tests that rely on Needle (https://needle.readthedocs.io/)
to assist with visual regression testing.
"""

from needle.cases import NeedleTestCase
from needle.driver import NeedleChrome

from dataclasses import dataclass, field
@dataclass(frozen=True)
class CompareNeedleTest(NeedleTestCase):
    baseline_directory: str = "visual_regression_screenshots"

    @classmethod
    def get_web_driver(cls):
        return NeedleChrome

    def test_plot_visual_regression(self, plot_img) -> bool:
        """
        Test that a screenshot of a plot element is the same as the screenshot on disk

        :param plot_img:
            Basename of plot image saved on disk, minus the ".png" part which will be added
        """
        print("-- PLOT VISUAL REGRESSION")
        try:
            self.compareScreenshot("plotly", plot_img)
            return True
        except AssertionError as e:
            print(str(e))
            return False
        except Exception:
            return False