#!/opt/bin/python3

"""
Creates a Github issue from a submitted comment on the gEAR page.

"""

import cgi, json, os, requests, socket, sys
from dotenv import load_dotenv
from requests.exceptions import HTTPError
from pathlib import Path
from uuid import uuid4

env_path = Path('..') / '.env'  # .env file is in "www" directory
load_dotenv(dotenv_path=env_path)


GITHUB_ACCESS_TOKEN=os.getenv("GITHUB_ACCESS_TOKEN")
GEAR_GIT_URL="https://api.github.com/repos/jorvis/gEAR/issues"
ASSIGNEES=["songeric1107"]

SITE_COMMENTS_PROJ_URL="https://api.github.com/projects/columns/8150789/cards" # Corresponds to jorvis/gEAR

SCREENSHOT_DIR = "contact_screenshots"
SCREENSHOT_URL = 'https://umgear.org/{}'.format(SCREENSHOT_DIR)

def main():

    print('Content-Type: application/json\n\n')
    result = {'error': [], 'success': 0 }

    form = cgi.FieldStorage()
    firstname = form.getvalue('submitter_firstname')
    lastname = form.getvalue('submitter_lastname')
    email = form.getvalue('submitter_email')
    title = form.getvalue('comment_title')
    comment = form.getvalue('comment')
    tag = form.getvalue('comment_tag')
    screenshot = form.getvalue('screenshot', None)
    private = form.getvalue('private_check')

    if not tag:
        tag = ''

    # Get IP address of the server hosting this CGI script (to determine gEAR flavor)
    hostname = socket.gethostname()
    ip_address = socket.gethostbyname(hostname)

    # TODO: this should probably go in the config file.
    # Also this assumes all screenshots comes from production
    if "nemo" in hostname:
        SCREENSHOT_URL = "https://nemoanalytics.org/{}".format(SCREENSHOT_DIR)
    elif "gcid" in hostname:
        SCREENSHOT_URL = "https://gcid.umgear.org/{}".format(SCREENSHOT_DIR)

    # If screenshot was provided, get URL and eventually assign to body
    screenshot_url = "None"
    if screenshot and not screenshot == "null":
        ext = os.path.splitext(screenshot)[1]
        new_basename = str(uuid4()) + ext
        src = f"../{SCREENSHOT_DIR}/files/{screenshot}"
        dst = f"../{SCREENSHOT_DIR}/{new_basename}" # Synlink is up a directory
        os.symlink(src, dst)
        screenshot_url = "{}/{}".format(SCREENSHOT_URL, new_basename)

    # In an effort to not blow up the "tags" field in github, I will just indicate the tags in the body of the Github issue
    body = (f"**From:** {firstname} {lastname}\n\n"
           f"**Email:** {email}\n\n"
           f"**Server IP:** {ip_address}\n\n"
           f"**Msg:** {comment}\n\n"
           f"**Tags:** {tag.split(', ')}\n\n"
           f"**Screenshot:** {screenshot_url}"
           )

    # Headers data (i.e. authentication)
    headers = { "Authorization": "token {}".format(GITHUB_ACCESS_TOKEN) }

    # Issue metadata
    data = {
        "title":title,
        "body":body,
        "labels":["site_comment"],
        "assignees":ASSIGNEES
    }

    # If user clicked "private" checkbox, send to private git repo
    git_url = GEAR_GIT_URL
    site_comments_url = SITE_COMMENTS_PROJ_URL
    if private == "false":
        # "false" is javascript string "false"
        git_url = GEAR_GIT_URL.replace("jorvis", "IGS")
        site_comments_url = SITE_COMMENTS_PROJ_URL.replace("8150789", "15595055")

    # Code from https://realpython.com/python-requests/
    try:
        response = requests.post(git_url, data = json.dumps(data), headers = headers)
        # If the response was successful (200- and 300- level status codes), no Exception will be raised
        response.raise_for_status()
    except HTTPError as http_err:
        result["error"] = f'HTTP error occurred: {http_err}'
    except Exception as err:
        result["error"] = f'Other error occurred: {err}'
    else:
        result["success"] = 1
    print(json.dumps(result))

    # Submission failed :-(
    if result["success"] == 0:
        sys.exit(0)

    # Add to 'Site comments' projects board
    headers["Accept"] = 'application/vnd.github.v3+json'
    content_id = response.json()["id"]   # ID of ticket just created
    data = {
        "content_id": content_id,
        "content_type": "Issue"
    }
    try:
        response = requests.post(site_comments_url, data = json.dumps(data), headers = headers )
        # If the response was successful (200- and 300- level status codes), no Exception will be raised
        response.raise_for_status()
    except HTTPError as http_err:
        print("HTTP error occurred{}.\n Count not assign id {} to 'Site Comments' project card".format(http_err, content_id), file=sys.stderr)
    except Exception as err:
        print("Other error occurred{}.\n Count not assign id {} to 'Site Comments' project card".format(err, content_id), file=sys.stderr)



if __name__ == '__main__':
    main()
