#!/opt/bin/python3

import cgi
import json
import os
import sys

import requests

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

DATASET_PREVIEWS_DIR = os.path.abspath(os.path.join('..', 'img', 'dataset_previews'))

def make_static_plotly_graph(filename, config, url):
    """Create (or overwrite) a static plotly PNG image using the existing config."""

    # WARNING: Disabling SSL verification in the POST call
    result = requests.post(url, json=config, verify=False)
    result.raise_for_status()

    decoded_result = result.json()

    # If plotly API threw an error, report as failed
    if "success" not in decoded_result:
        return False
    if "success" in decoded_result:
        # If success is a list, flatten and take the first element
        if isinstance(decoded_result["success"], list):
            decoded_result["success"] = decoded_result["success"][0]
        if decoded_result["success"] < 0:
            return False

    plot_json = decoded_result["plot_json"]

    import plotly.graph_objects as go
    # We know the figure is valid, so skip potential illegal property issues.
    fig = go.Figure(data=plot_json["data"], layout=plot_json["layout"], skip_invalid=True)
    fig.write_image(filename)
    try:
        os.chmod(filename, 0o666)
    except Exception as e:
        print("Could not chmod {} for reason: {}".format(filename, str(e)), file=sys.stderr)
    return True

def make_static_tsne_graph(filename, config, url):
    """Create (or overwrite) a static tsne PNG image using the existing config."""

    # WARNING: Disabling SSL verification in the POST call
    result = requests.post(url, json=config, verify=False)
    result.raise_for_status()

    decoded_result = result.json()

    # If plotly API threw an error, report as failed
    if "success" not in decoded_result:
        return False
    if "success" in decoded_result:
        # If success is a list, flatten and take the first element
        if isinstance(decoded_result["success"], list):
            decoded_result["success"] = decoded_result["success"][0]
        if decoded_result["success"] < 0:
            return False

    import base64
    img_data = base64.urlsafe_b64decode(decoded_result["image"])
    try:
        with open(filename, "wb") as fh:
            fh.write(img_data)
    except Exception as e:
        print("Could not write {} for reason: {}".format(filename, str(e)), file=sys.stderr)
        return False

    return True

def make_static_svg(filename, dataset_id):
    """Create (or overwrite) a static svg PNG image using the existing config."""

    svg_filepath = "../datasets_uploaded/{}.svg".format(dataset_id)

    import cairosvg
    cairosvg.svg2png(url=svg_filepath, write_to=filename)

    try:
        os.chmod(filename, 0o666)
    except Exception as e:
        print("Could not chmod {} for reason: {}".format(filename, str(e)), file=sys.stderr)
    return True

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    display_id = form.getvalue('id')
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')
    label = form.getvalue('label')
    plot_type = form.getvalue('plot_type')
    plotly_config = form.getvalue('plotly_config')
    is_local = form.getvalue('is_local', False)

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    user = geardb.get_user_from_session_id(session_id=session_id)

    if not user:
        print('User not found for session_id {}'.format(session_id), file=sys.stderr)
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(dict(display_id=None, success=False)))
        return

    user_id = user.id

    if user_id is None:
        print('User ID is None for session_id {}'.format(session_id), file=sys.stderr)
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(dict(display_id=None, success=False)))
        return

    if display_id:
        # display_id exists, so update

        query = "SELECT user_id FROM dataset_display WHERE id = %s"
        cursor.execute(query, (display_id,))
        row = cursor.fetchone()
        if row is not None:
            (display_owner,) = row

            if int(user_id) == display_owner:
                # A user must be the owner of the particular
                # display in order to save/update. We check here
                # in case a user maliciously changes the HTML id of the display
                # and then saves
                query = """
                    UPDATE dataset_display
                    SET label = %s, plot_type = %s, plotly_config = %s
                    WHERE id = %s;
                """
                cursor.execute(query, (label, plot_type, plotly_config, display_id))
                result = dict(display_id=display_id, success=True)
            else:
                print('UPDATE DIDNT HAPPEN?', file=sys.stderr)
                result = dict(display_id=display_id, success=False)
        else:
            print('Display ID {} not found in database.'.format(display_id), file=sys.stderr)
            result = dict(display_id=display_id, success=False)
    else:
        # Display doesn't exist yet, insert new

        query = """
            INSERT INTO dataset_display
            (dataset_id, user_id, label, plot_type, plotly_config)
            VALUES (%s, %s, %s, %s, %s)
        """
        cursor.execute(query,
            (dataset_id, user_id, label, plot_type, plotly_config))
        result = dict(success=True)

        # Retrieve display ID so we can generate the static image
        query = """
            SELECT id FROM dataset_display
            WHERE dataset_id = %s
                AND user_id = %s
                AND label = %s
                AND plot_type = %s
                AND plotly_config = %s
            ORDER BY id DESC LIMIT 1
        """
        cursor.execute(query,
            (dataset_id, user_id, label, plot_type, plotly_config))

        row = cursor.fetchone()
        if row:
            (display_id,) = row
            result["display_id"] = display_id

        if not display_id:
            print('Display ID not found after insert.', file=sys.stderr)
            result = dict(success=False)

    filename = os.path.join(DATASET_PREVIEWS_DIR, "{}.{}.png".format(dataset_id, display_id))

    # Normalize path to avoid directory traversal attacks (e.g. ../../../etc/passwd) and validate
    filename = os.path.normpath(filename)
    if not filename.startswith(DATASET_PREVIEWS_DIR):
        print("Invalid filename: {}".format(filename), file=sys.stderr)
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        result["success"] = False
        print(json.dumps(result))
        return

    # Add plot_type to config so it can be passed in the POST cmd as well
    config = json.loads(plotly_config)
    config["plot_type"] = plot_type

    alt_url = "http://localhost/api/plot/{}".format(dataset_id)

    image_success = False
    url = "https://localhost/api/plot/{}".format(dataset_id)

    if is_local:
        url = alt_url

    try:
        if plot_type.lower() in ['bar', 'scatter', 'violin', 'line', 'contour', 'tsne_dynamic', 'tsne/umap_dynamic']:
            image_success = make_static_plotly_graph(filename, config, url)
        elif plot_type.lower() in ["mg_violin", "dotplot", "volcano", "heatmap", "quadrant"]:
            url += "/mg_plotly"
            image_success = make_static_plotly_graph(filename, config, url)
        elif plot_type.lower() in ["tsne_static", "umap_static", "pca_static", "tsne"]:
            url += "/tsne"
            image_success = make_static_tsne_graph(filename, config, url)
        elif plot_type.lower() in ["svg"]:
            url += "/svg"
            image_success = make_static_svg(filename, dataset_id)
        elif plot_type.lower() in ["epiviz"]:
            url += "/epiviz"
            pass
        else:
            print("Plot type {} for display id {} is not recognizable".format(plot_type, display_id))
    except Exception as e:
        print("Error with plotting dataset {}".format(dataset_id), file=sys.stderr)
        print(str(e), file=sys.stderr)

    if not image_success:
        try:
            print("Attempting to use alt_url: {}".format(alt_url), file=sys.stderr)
            url = alt_url
            if plot_type.lower() in ['bar', 'scatter', 'violin', 'line', 'contour', 'tsne_dynamic', 'tsne/umap_dynamic']:
                image_success = make_static_plotly_graph(filename, config, url)
            elif plot_type.lower() in ["mg_violin", "dotplot", "volcano", "heatmap", "quadrant"]:
                url += "/mg_plotly"
                image_success = make_static_plotly_graph(filename, config, url)
            elif plot_type.lower() in ["tsne_static", "umap_static", "pca_static", "tsne"]:
                url += "/tsne"
                image_success = make_static_tsne_graph(filename, config, url)
            elif plot_type.lower() in ["mg_tsne_static", "mg_umap_static", "mg_pca_static"]:
                url += "/mg_tsne"
                image_success = make_static_tsne_graph(filename, config, url)
            elif plot_type.lower() in ["svg"]:
                url += "/svg"
                image_success = make_static_svg(filename, dataset_id)
        except Exception as e:
            print("Error with plotting dataset {} using alt_url".format(dataset_id), file=sys.stderr)
            print(str(e), file=sys.stderr)

    if not image_success:
        print("Could not create static image file for display id {}".format(display_id), file=sys.stderr)

    cnx.commit()
    cursor.close()
    cnx.close()


    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
