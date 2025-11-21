#!/opt/bin/python3

""" Given a display ID, return the static image """

import cgi, json
import os, sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

WEB_IMAGE_ROOT = '../img/dataset_previews'  # origin is cgi script dir

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    dataset_id = form.getvalue("dataset_id")
    display_id = form.getvalue('display_id')

    image_preview_url = os.path.join(WEB_IMAGE_ROOT, "{}.{}.png".format(dataset_id, display_id))

    image_preview_url = os.path.normpath(image_preview_url)
    if not image_preview_url.startswith(WEB_IMAGE_ROOT):
        raise Exception("Invalid filename: {}".format(image_preview_url))

    #print(image_preview_url, file=sys.stderr)

    display = geardb.get_display_by_id(display_id)
    if display and display.get("plot_type", None) == "gosling":
        image_preview_url = "{0}/gosling.png".format(WEB_IMAGE_ROOT)

    if not os.path.exists(image_preview_url):
        image_preview_url = "{0}/missing.png".format(WEB_IMAGE_ROOT)


    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(image_preview_url))

if __name__ == '__main__':
    main()