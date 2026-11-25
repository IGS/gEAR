#!/opt/bin/python3

"""
save_dataset_display.cgi - Create or update a saved dataset display and regenerate its static preview PNG.

Input: id (display_id; omit to insert new), session_id (required), dataset_id, label, plot_type,
       plotly_config (JSON string), is_local.
Output: JSON {display_id, success}; writes img/dataset_previews/<dataset_id>.<display_id>.png
        (only when the display was saved).
"""

import cgi
import json
import os
import sys

import requests

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear import plot_types

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

def has_preview_image(plot_type):
    """Return True if this plot type gets a preview image (gosling and spatial displays do not)."""
    return plot_types.family(plot_type) == "svg" or plot_types.preview_endpoint(plot_type) is not None

def make_preview_image(filename, plot_type, config, dataset_id, base_url):
    """Draw a display's preview PNG. Returns True if the image was written.

    The plot-type registry (www/plot_types.json) decides which API endpoint under base_url draws it.
    SVG previews are drawn from the dataset's uploaded SVG instead.
    """
    family = plot_types.family(plot_type)
    if family == "svg":
        return make_static_svg(filename, dataset_id)

    url = base_url + plot_types.preview_endpoint(plot_type)
    if family in ("scanpy", "mg_scanpy"):
        return make_static_tsne_graph(filename, config, url)
    return make_static_plotly_graph(filename, config, url)

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    display_id = form.getfirst('id')
    session_id = form.getfirst('session_id')
    dataset_id = form.getfirst('dataset_id')
    label = form.getfirst('label')
    plot_type = form.getfirst('plot_type')
    plotly_config = form.getfirst('plotly_config')
    is_local = form.getfirst('is_local', False)

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

        # Use the new row's ID directly. Looking it up again by matching the saved columns missed
        #  rows with a blank label (stored as NULL, and "label = NULL" is never true), which named
        #  the preview image "<dataset_id>.None.png"
        display_id = cursor.lastrowid
        if display_id:
            result = dict(display_id=display_id, success=True)
        else:
            print('Display ID not found after insert.', file=sys.stderr)
            result = dict(display_id=None, success=False)

    # Don't generate a preview for a display that wasn't saved (e.g. one this user doesn't own)
    if not result["success"]:
        cnx.commit()
        cursor.close()
        cnx.close()
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

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

    url = "https://localhost/api/plot/{}".format(dataset_id)

    if is_local:
        url = alt_url

    if plot_types.normalize(plot_type) is None:
        print("Plot type {} for display id {} is not recognizable".format(plot_type, display_id), file=sys.stderr)

    # Draw the preview, retrying once on the plain-http URL if the first attempt fails
    image_success = False
    if has_preview_image(plot_type):
        for attempt_url in dict.fromkeys([url, alt_url]):   # a single attempt when both URLs are the same
            if attempt_url == alt_url and attempt_url != url:
                print("Attempting to use alt_url: {}".format(alt_url), file=sys.stderr)
            try:
                image_success = make_preview_image(filename, plot_type, config, dataset_id, attempt_url)
            except Exception as e:
                print("Error with plotting dataset {} using {}".format(dataset_id, attempt_url), file=sys.stderr)
                print(str(e), file=sys.stderr)
            if image_success:
                break

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
