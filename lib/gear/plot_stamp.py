"""
plot_stamp.py - Site watermark and dataset citation stamped onto downloaded plot images.

The stamp is a few lines of small text placed below the plot: the dataset title, then its
PubMed ID (or permalink) with a "Made with gEAR" line. The JS counterpart for plots
downloaded on the client is Citation.plotStampLines in www/js/classes/citation.js.
"""

import geardb
import plotly.graph_objects as go

MAX_TITLE_LENGTH = 100
STAMP_COLOR = "#555555"

# Plotly legend ID reserved for the stamp, high enough not to collide with legends a plot already uses
STAMP_LEGEND = "legend99"


def build_stamp_lines(dataset) -> list[str]:
    """
    Build the stamp text lines for a dataset.

    Args:
        dataset (geardb.Dataset): Dataset being plotted.

    Returns:
        list[str]: Lines of plain text, top to bottom.
    """
    site_url = (geardb.domain_url or "https://umgear.org").rstrip("/")
    site_label = geardb.domain_short_label or "gEAR"

    title = (dataset.title or "").strip()
    if len(title) > MAX_TITLE_LENGTH:
        title = title[: MAX_TITLE_LENGTH - 1].rstrip() + "…"

    # Point back to the dataset permalink when there is no publication to cite
    pubmed_id = str(dataset.pubmed_id or "").strip()
    if pubmed_id.isdigit():
        site_line = "PMID {} · Made with {} · {}".format(
            pubmed_id, site_label, site_url.removeprefix("https://").removeprefix("http://")
        )
    elif dataset.share_id:
        site_line = f"Made with {site_label} · {site_url}/p?id=d.{dataset.share_id}"
    else:
        site_line = f"Made with {site_label} · {site_url}"

    return [line for line in (title, site_line) if line]


def stamp_plotly_figure(fig: go.Figure, lines: list[str]) -> None:
    """
    Add the stamp lines below a Plotly figure. Updates 'fig' inplace.

    Plotly annotations cannot push the margins, so they would overlap long or rotated tick
    labels. Instead the text is the title of an otherwise empty legend anchored to the bottom
    of the container, which makes Plotly grow the bottom margin to fit it.
    """
    if not lines:
        return

    # An explicit showlegend=False would hide the stamp legend too, so hide the plot's legend
    # entries one trace at a time instead.
    if fig.layout.showlegend is False:
        fig.for_each_trace(lambda trace: trace.update(showlegend=False))
        fig.update_layout(showlegend=True)

    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker={"color": "rgba(0,0,0,0)", "size": 1},
            name="",
            showlegend=True,
            legend=STAMP_LEGEND,
            hoverinfo="skip",
        )
    )
    fig.update_layout(
        {
            STAMP_LEGEND: {
                "orientation": "h",
                "xref": "container",
                "yref": "container",
                "x": 0.01,
                "y": 0,
                "xanchor": "left",
                "yanchor": "bottom",
                "bgcolor": "rgba(0,0,0,0)",
                "entrywidth": 1,
                "title": {"text": "<br>".join(lines), "font": {"size": 10, "color": STAMP_COLOR}},
            }
        }
    )


def stamp_matplotlib_figure(fig, lines: list[str], fontsize: int = 8) -> None:
    """
    Add the stamp lines below everything already drawn on a matplotlib figure. Updates 'fig' inplace.

    The text can land outside the figure canvas, so save the figure with bbox_inches="tight".
    """
    if not lines:
        return

    renderer = fig.canvas.get_renderer()
    content = fig.get_tightbbox(renderer)  # in inches
    pad = 0.05
    fig.text(
        content.x0 / fig.get_figwidth(),
        (content.y0 - pad) / fig.get_figheight(),
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=fontsize,
        color=STAMP_COLOR,
    )
