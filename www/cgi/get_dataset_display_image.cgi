#!/opt/bin/python3

""" Given a display ID, return the static image """

import cgi, json
import os, sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

WEB_IMAGE_ROOT = '../img/dataset_previews'  # origin is cgi script dir

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    dataset_id = form.getfirst("dataset_id")
    display_id = form.getfirst('display_id')

    image_root = os.path.normpath(WEB_IMAGE_ROOT)
    image_preview_url = os.path.normpath(os.path.join(image_root, "{}.{}.png".format(dataset_id, display_id)))

    # The IDs must not lead outside the preview directory, into a sibling like dataset_previews_old,
    #  or into a subdirectory. Kept inline (normpath, then startswith) so CodeQL sees the path guard.
    if image_preview_url.startswith(image_root + os.sep) and os.path.dirname(image_preview_url) == image_root:
        display = geardb.get_display_by_id(display_id)
        if display and display.get("plot_type", None) == "gosling":
            image_preview_url = "{0}/gosling.png".format(WEB_IMAGE_ROOT)
        elif not os.path.exists(image_preview_url):
            image_preview_url = "{0}/missing.png".format(WEB_IMAGE_ROOT)
    else:
        image_preview_url = "{0}/missing.png".format(WEB_IMAGE_ROOT)


    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(image_preview_url))

if __name__ == '__main__':
    main()