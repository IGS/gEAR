#!/opt/bin/python3

"""
generate_static_layouts.py - Generate a static plot based on the dataset's layout configurations stored in database

Shaun Adkins - sadkins@som.umaryland.edu
"""

DATASET_PREVIEWS_DIR = "/var/www/img/dataset_previews"

import os, json, requests, sys

import plotly.graph_objects as go

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
import geardb

def get_all_displays(cursor):
    """Get all dataset displays out of the database."""

    query  = "SELECT dataset_id, id, plot_type, plotly_config from dataset_display"

    cursor.execute(query)

    displays = dict()

    for (dataset_id, display_id, plot_type, plotly_config) in cursor:
        displays.setdefault(dataset_id, dict())
        displays[dataset_id].setdefault(display_id, {"plot_type":plot_type, "config":plotly_config, "default":False})

    # Determine default display images (to be symlinked later)
    defaults_query  = "SELECT d.id AS dataset_id, dd.id AS display_id " \
        "FROM dataset d " \
        "JOIN dataset_display dd ON dd.dataset_id=d.id " \
        "JOIN dataset_preference dp ON dp.display_id=dd.id " \
        "WHERE dp.user_id=d.owner_id"

    cursor.execute(defaults_query)

    for (dataset_id, display_id) in cursor:
        displays[dataset_id][display_id]["default"] = True

    return displays

def make_static_plotly_graph(filename, config, url):
    """Create a static plotly PNG image using the existing config."""
    # WARNING: Disabling SSL verification in the POST call
    # Throw error if things went awry (check apache ssl_error logs)
    result = requests.post(url, json=config, verify=False)
    result.raise_for_status()

    decoded_result = result.json()

    # If plotly API threw an error, report as failed
    if "success" in decoded_result and decoded_result["success"] < 0:
        return False

    #plot_config = result["plot_config"]
    plot_json = decoded_result["plot_json"]

    fig = go.Figure(data=plot_json["data"], layout=plot_json["layout"])

    fig.write_image(filename)
    return True

def make_static_svg(filename, config, url):
    return False
    requests.post("http://127.0.0.1/api/plot/{}/svg?gene={}".format(dataset_id, config["gene_symbol"]))

def make_static_tsne_graph(filename, config, url):
    """Create (or overwrite) a static tsne PNG image using the existing config."""

    # WARNING: Disabling SSL verification in the POST call
    result = requests.post(url, json=config, verify=False)
    result.raise_for_status()

    decoded_result = result.json()

    # If plotly API threw an error, report as failed
    if not "success" in decoded_result:
        return False
    if "success" in decoded_result and decoded_result["success"] < 0:
        return False

    import base64
    img_data = base64.urlsafe_b64decode(decoded_result["image"])
    with open(filename, "wb") as fh:
        fh.write(img_data)

    return True

def main():

    if not os.path.isdir(DATASET_PREVIEWS_DIR):
        sys.exit("Directory to store static dataset images does not exist.")

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    # Return all layouts
    displays = get_all_displays(cursor)

    # Plot and save as static image
    for dataset_id in displays:
        for display_id, props in displays[dataset_id].items():
            filename = os.path.join(DATASET_PREVIEWS_DIR, "{}.{}.png".format(dataset_id, display_id))

            # Add to config so it can be passed in the POST cmd as well
            config = json.loads(props["config"])
            config['plot_type'] = props['plot_type']

            gene = "single"

            if os.path.isfile(filename):
                print("Overwriting file {}".format(filename))

            success = False
            url = "https://localhost/api/plot/{}".format(dataset_id)
            try:
                # Plotly
                if props["plot_type"] in ['bar', 'scatter', 'violin', 'line', 'contour', 'tsne_dynamic', 'tsne/umap_dynamic']:
                    url += "/"
                    success = make_static_plotly_graph(filename, config, url)
                elif props["plot_type"] in ["mg_violin", "dotplot", "volcano", "heatmap", "quadrant"]:
                    url += "/mg_dash"
                    gene = "multi"
                    success = make_static_plotly_graph(filename, config, url)
                # tSNE (todo later)
                elif props["plot_type"] in ["tsne_static", "umap_static", "pca_static", "tsne"]:
                    url += "/tsne"
                    success = make_static_tsne_graph(filename, config, url)
                # SVG (todo later)
                elif props["plot_type"] in ["svg"]:
                    url += "/svg"
                    success = make_static_svg(filename, config, url)
                # Epiviz (todo later)
                elif props["plot_type"] in ["epiviz"]:
                    url += "/epiviz"
                    pass
                else:
                    print("Plot type {} for display id {} is not recognizable".format(props["plot_type"], display_id))
            except Exception as e:
                print("Error with plotting dataset {}".format(dataset_id), file=sys.stderr)
                print(str(e), file=sys.stderr)

            if not success:
                print("Could not create static image for display id {}".format(display_id))
                continue

            # Changing to group write permissions, so that edited curator displays (by www-data or jorvis user) can be used to generate new images
            try:
                os.chmod(filename, 0o666)
            except Exception as e:
                print("Could not chmod file {}. Reason: {}".format(filename, str(e)))

            # Create symlink to filename for all the designated 'default' displays
            if props["default"]:
                symlink_path = os.path.join(DATASET_PREVIEWS_DIR, "{}.{}.default.png".format(dataset_id, gene))
                # If running multiple times, we shouldn't need to recreate the symlink
                try:
                    os.symlink(filename, symlink_path)
                except:
                    print("Symlink for {} to default already exists. Skipping".format(filename))

    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()
