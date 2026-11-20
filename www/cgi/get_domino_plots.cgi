#!/opt/bin/python3

"""
Generates cell-cell communication plots (dot plots and river/Sankey diagrams) from a
DominoSignal output .parquet file.  Ported from notebooks/domino_plotting.ipynb; the
plotting itself lives in lib/gear/dominosignal.py.

The parquet file is looked up by bare filename within www/cache/domino_signal.

Request parameters (all optional except filename):
    filename            - Name of the .parquet file within the domino_signal directory
    selected_tf         - TF for the dot_tf and river_tf plots
    selected_ligand     - Ligand for the river_ligand plot
    selected_receptor   - Receptor for the river_receptor plot
    sending_clusters, receiving_clusters
                        - Lists for dot_clusters and dot_clusters_filtered
    selected_tfs, selected_ligands, selected_receptors,
    selected_sending, selected_receiving
                        - Lists for river_selected (and river_selected_notf, minus TFs)
    exp_range_threshold - Filter threshold for dot_clusters_filtered (default 0.5)
    exp_max_threshold   - Filter threshold for dot_clusters_filtered (default 1.5)
    river_apply_thresholds
                        - If "1" or "true", river_selected and river_selected_notf only include ligand/receptor
                          pairs passing the above thresholds, computed over the river plot's own selection
    plots               - List of plots to generate (default all)
    colorblind_mode     - If "1" or "true", use the color-blind safe palette (ds.PALETTE_2)
    dot_size_scale      - Multiplier for the dot plots' default dot sizes (default 1, clamped to 0.25-3)
    dot_white_background
                        - If "1" or "true", dot plots use a white background instead of Plotly's default
    image_format        - If "pdf", return the single plot in "plots" as a PDF file download instead of JSON
    width, height       - Image size in pixels when image_format is set (defaults to the plot's own size)

Lists may be passed as JSON arrays or comma-separated strings.

Data structure returned:

{
    success: 1,
    plots: { dot_tf: {data, layout} | null, dot_clusters: ..., ... },
    messages: { dot_tf: "Select a TF to generate this plot", ... }
}

With image_format set, the image file is returned, or on error {success: -1, message: ...}.

"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

from gear import dominosignal as ds

ALL_PLOTS = [
    "dot_tf", "dot_clusters", "dot_clusters_filtered",
    "river_tf", "river_ligand", "river_receptor",
    "river_selected", "river_selected_notf",
]

# Image formats that can be requested via image_format, with their content types
IMAGE_FORMATS = {
    "pdf": "application/pdf",
}
IMAGE_MIN_SIZE = 100
IMAGE_MAX_SIZE = 20000
KALEIDO_DEFAULT_WIDTH = 700
KALEIDO_DEFAULT_HEIGHT = 500

# (width, height) pixel changes to try, in order, until a PDF fits on a single page (see render_single_page_pdf)
PDF_SIZE_ADJUSTMENTS = [(0, 0), (0, -1), (-1, 0), (-1, -1), (0, -2), (-2, 0), (-2, -1), (-1, -2), (-2, -2)]


class InvalidRequest(Exception):
    """Raised when the request itself is unusable."""


class PlotImage:
    """A single plot rendered as an image file, to be returned as a download."""

    def __init__(self, filename, image_format, content):
        self.filename = filename
        self.content_type = IMAGE_FORMATS[image_format]
        self.content = content


def main():
    form = cgi.FieldStorage()

    domino_df = ds.load_domino_df(form.getvalue('filename'))

    selected_tf = form.getvalue('selected_tf')
    selected_ligand = form.getvalue('selected_ligand')
    selected_receptor = form.getvalue('selected_receptor')

    selected_tfs = get_list_param(form, 'selected_tfs')
    selected_ligands = get_list_param(form, 'selected_ligands')
    selected_receptors = get_list_param(form, 'selected_receptors')
    selected_sending = get_list_param(form, 'selected_sending')
    selected_receiving = get_list_param(form, 'selected_receiving')

    sending_clusters = get_list_param(form, 'sending_clusters')
    receiving_clusters = get_list_param(form, 'receiving_clusters')

    exp_range_threshold = float(form.getvalue('exp_range_threshold') or ds.DEFAULT_EXP_RANGE_THRESHOLD)
    exp_max_threshold = float(form.getvalue('exp_max_threshold') or ds.DEFAULT_EXP_MAX_THRESHOLD)

    river_apply_thresholds = get_flag_param(form, 'river_apply_thresholds')

    requested_plots = get_list_param(form, 'plots') or ALL_PLOTS

    palette = ds.PALETTE_2 if get_flag_param(form, 'colorblind_mode') else ds.PALETTE_1

    # Advanced plotting options for the dot plots
    dot_size_scale = min(ds.DOT_SIZE_SCALE_MAX, max(ds.DOT_SIZE_SCALE_MIN, float(form.getvalue('dot_size_scale') or 1)))
    dot_options = dict(palette=palette, dot_size_scale=dot_size_scale,
                       white_background=get_flag_param(form, 'dot_white_background'))

    builders = {
        "dot_tf": lambda: ds.domino_dot_plot_tf(domino_df, selected_tf, **dot_options),
        "dot_clusters": lambda: ds.domino_dot_plot(domino_df, sending_clusters, receiving_clusters, **dot_options),
        "dot_clusters_filtered": lambda: ds.domino_dot_plot(domino_df, sending_clusters, receiving_clusters, filter=True,
                                                            exp_range_threshold=exp_range_threshold, exp_max_threshold=exp_max_threshold,
                                                            **dot_options),
        "river_tf": lambda: ds.domino_river_plot_tf(domino_df, selected_tf, palette=palette),
        "river_ligand": lambda: ds.domino_river_plot_ligand(domino_df, selected_ligand, palette=palette),
        "river_receptor": lambda: ds.domino_river_plot_receptor(domino_df, selected_receptor, palette=palette),
        "river_selected": lambda: ds.domino_river_plot_selected(domino_df, selected_ligands, selected_receptors,
                                                                selected_tfs, selected_sending, selected_receiving,
                                                                apply_thresholds=river_apply_thresholds,
                                                                exp_range_threshold=exp_range_threshold, exp_max_threshold=exp_max_threshold,
                                                                palette=palette),
        "river_selected_notf": lambda: ds.domino_river_plot_selected_notf(domino_df, selected_ligands, selected_receptors,
                                                                          selected_sending, selected_receiving,
                                                                          apply_thresholds=river_apply_thresholds,
                                                                          exp_range_threshold=exp_range_threshold, exp_max_threshold=exp_max_threshold,
                                                                          palette=palette),
    }

    image_format = form.getvalue('image_format')
    if image_format:
        return build_plot_image(form, builders, requested_plots, image_format)

    result = {'success': 1, 'plots': {}, 'messages': {}}

    for plot_name in requested_plots:
        if plot_name not in builders:
            result['messages'][plot_name] = "Unknown plot type"
            continue

        try:
            fig = builders[plot_name]()
            result['plots'][plot_name] = json.loads(fig.to_json())
        except ds.PlotSkipped as e:
            result['plots'][plot_name] = None
            result['messages'][plot_name] = str(e)

    return result


def build_plot_image(form, builders, requested_plots, image_format):
    """Renders the single requested plot as an image file (via Kaleido)."""
    if image_format not in IMAGE_FORMATS:
        raise InvalidRequest("Unsupported image format: {0}".format(image_format))

    if len(requested_plots) != 1 or requested_plots[0] not in builders:
        raise InvalidRequest("Request exactly one plot to download as an image")

    plot_name = requested_plots[0]
    fig = builders[plot_name]()

    # Kaleido can only draw the river plots' label text shadows as embedded bitmaps (one per label),
    #  so drop them to keep the file fully vector and small.  They stay on screen for readability.
    fig.update_traces(textfont_shadow="none", selector=dict(type="sankey"))

    width = get_size_param(form, 'width')
    height = get_size_param(form, 'height')
    if image_format == "pdf":
        content = render_single_page_pdf(fig, width, height)
    else:
        content = fig.to_image(format=image_format, width=width, height=height)

    # Name the file after the DominoSignal output and plot, e.g. domino_network_p21_dot_clusters.pdf
    file_stem = os.path.splitext(os.path.basename(form.getvalue('filename')))[0]
    file_stem = "".join(c if c.isalnum() or c in "-_." else "_" for c in file_stem)
    return PlotImage("{0}_{1}.{2}".format(file_stem, plot_name, image_format), image_format, content)


def render_single_page_pdf(fig, width, height):
    """
    Renders the figure as a one-page PDF.

    Kaleido rounds the PDF page size, so at about half of all image sizes the plot overflows the page
    by a fraction of a pixel, adding a blank second page.  When that happens, re-render at a size a
    pixel or two smaller (an invisible difference), which always gives a single page.
    """
    # Kaleido's defaults when no size is given
    width = width or fig.layout.width or KALEIDO_DEFAULT_WIDTH
    height = height or fig.layout.height or KALEIDO_DEFAULT_HEIGHT

    for width_change, height_change in PDF_SIZE_ADJUSTMENTS:
        content = fig.to_image(format="pdf", width=width + width_change, height=height + height_change)
        if count_pdf_pages(content) == 1:
            break

    return content


def count_pdf_pages(content):
    # Page objects are "/Type /Page"; the page tree itself is "/Type /Pages"
    return content.count(b"/Type /Page") - content.count(b"/Type /Pages")


def get_size_param(form, name):
    """Image width or height in pixels, clamped to a sensible range, or None to use the plot's own size."""
    value = form.getvalue(name)
    if not value:
        return None

    return int(min(IMAGE_MAX_SIZE, max(IMAGE_MIN_SIZE, float(value))))


def get_flag_param(form, name):
    # Flags are sent from the client as form data, so arrive as strings
    return str(form.getvalue(name, '')).lower() in ('1', 'true')


def get_list_param(form, name):
    value = form.getvalue(name)
    if not value:
        return []

    if value.startswith('['):
        return [str(v) for v in json.loads(value)]

    return [v.strip() for v in value.split(',') if v.strip()]


if __name__ == '__main__':
    try:
        result = main()
    except (ds.DominoFileError, ds.PlotSkipped, InvalidRequest) as e:
        result = {'success': -1, 'message': str(e)}
    except Exception as e:
        result = {'success': -1, 'message': "Error generating DominoSignal plots: {0}".format(e)}

    sys.stdout = original_stdout

    if isinstance(result, PlotImage):
        sys.stdout.write("Content-Type: {0}\r\n".format(result.content_type))
        sys.stdout.write('Content-Disposition: attachment; filename="{0}"\r\n\r\n'.format(result.filename))
        sys.stdout.flush()
        sys.stdout.buffer.write(result.content)
        sys.stdout.buffer.flush()
    else:
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
