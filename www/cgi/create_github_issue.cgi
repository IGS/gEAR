#!/opt/bin/python3

"""
Creates a Github issue from a submitted comment on the gEAR page.

"""

import cgi
import json
import os
import re
import sys
from pathlib import Path
from uuid import uuid4

import requests
from requests.exceptions import HTTPError

lib_path = Path(__file__).resolve().parents[2] / 'lib'
sys.path.append(str(lib_path))
import geardb

# Repository information for creating issues in Github
GEAR_GIT_URL="https://api.github.com/repos/jorvis/gEAR/issues"
PRIVATE_OWNER = "jorvis"
PUBLIC_OWNER = "IGS"
ASSIGNEES=["jorvis"]

SCREENSHOT_DIR = "contact_screenshots"
# www/contact_screenshots. The upload handler (index.php) saves files into its "files" subdirectory.
SCREENSHOT_ROOT = Path(__file__).resolve().parents[1] / SCREENSHOT_DIR

GITHUB_TIMEOUT = 15     # seconds to wait for the Github API before giving up
SECURITY_CHECK_ANSWER = "28"    # 7 x 4, must match the contact form

# Github limits issue titles to 256 characters and bodies to 65536 characters.
# The comment limit leaves room for the rest of the issue body.
MAX_LENGTHS = {
    "submitter_fullname": 256,
    "submitter_email": 256,
    "comment_title": 256,
    "comment": 60000,
    "comment_tag": 1000,
}
EMAIL_PATTERN = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")

def validate_form(fields, security_check):
    """
    Server-side copy of the contact form's checks, so direct POSTs to this script
    cannot bypass the browser validation.

    Args:
        fields: dict of form field name to stripped value.
        security_check: the submitted answer to the "7 x 4" question.

    Returns:
        A list of error messages (empty if the form is valid).
    """
    errors = []

    required = {
        "submitter_fullname": "Name",
        "submitter_email": "Email address",
        "comment_title": "Title",
        "comment": "Question / Comment",
    }
    for name, label in required.items():
        if not fields[name]:
            errors.append(f"{label} is required")

    if fields["submitter_email"] and not EMAIL_PATTERN.match(fields["submitter_email"]):
        errors.append("Email address is not valid")

    for name, max_length in MAX_LENGTHS.items():
        if len(fields[name]) > max_length:
            label = required.get(name, "Tags")
            errors.append(f"{label} is longer than {max_length} characters")

    if security_check != SECURITY_CHECK_ANSWER:
        errors.append("Security check answer is incorrect")

    return errors

def link_screenshot(screenshot, screenshot_url):
    """
    Creates a public, randomly named symlink to an uploaded screenshot.

    Args:
        screenshot: file name returned by the upload handler (a name inside contact_screenshots/files).
        screenshot_url: public base URL of the contact_screenshots directory.

    Returns:
        The public URL of the symlink, or a short note if the screenshot could not be linked.
    """
    files_dir = (SCREENSHOT_ROOT / "files").resolve()
    src = (files_dir / screenshot).resolve()

    # Only accept a plain file name inside files/ (no "../" or other path tricks)
    if src.parent != files_dir or not src.is_file():
        print(f"Screenshot '{screenshot}' not found in {files_dir}", file=sys.stderr)
        return "None (uploaded screenshot not found)"

    new_basename = f"{uuid4()}{src.suffix}"
    try:
        # Relative target so the link still works if the web root moves
        os.symlink(Path("files") / src.name, SCREENSHOT_ROOT / new_basename)
    except OSError as err:
        print(f"Could not link screenshot '{screenshot}': {err}", file=sys.stderr)
        return "None (could not link screenshot)"

    return f"{screenshot_url}/{new_basename}"

def main():

    print('Content-Type: application/json\n\n')
    result = {'error': [], 'success': 0 }

    form = cgi.FieldStorage()
    fields = {name: (form.getfirst(name) or "").strip() for name in MAX_LENGTHS}
    security_check = (form.getfirst('super_impressive_security_check') or "").strip()
    screenshot = form.getfirst('screenshot', None)
    private = form.getfirst('private_check')

    errors = validate_form(fields, security_check)
    if errors:
        result["error"] = "; ".join(errors)
        print(json.dumps(result))
        return

    # Verify that the access token is set in gear.ini (the templates ship with an empty value)
    access_token = geardb.servercfg.get("github", "access_token", fallback="").strip()
    if not access_token:
        result["error"] = "Github access token not set in gear.ini"
        print(json.dumps(result))
        return

    # Get domain URL of this server, so we can ensure this is about to be read locally and for the UCSC exporting
    if os.getenv("ENVIRONMENT", "production").lower() == "development":
        domain_url = "http://localhost:8080"
    else:
        domain_url = geardb._read_domain_url()

    SCREENSHOT_URL = f'{domain_url}/{SCREENSHOT_DIR}'

    # If screenshot was provided, get URL and eventually assign to body
    screenshot_url = "None"
    if screenshot and not screenshot == "null":
        screenshot_url = link_screenshot(screenshot, SCREENSHOT_URL)

    # In an effort to not blow up the "tags" field in github, I will just indicate the tags in the body of the Github issue
    # (contact.js already sends the tags as one comma-separated string)
    body = (f"**From:** {fields['submitter_fullname']}\n\n"
           f"**Email:** {fields['submitter_email']}\n\n"
           f"**URL:** {domain_url}\n\n"
           f"**Msg:** {fields['comment']}\n\n"
           f"**Tags:** [{fields['comment_tag'] or 'None'}]\n\n"
           f"**Screenshot:** {screenshot_url}"
           )

    # Headers data (i.e. authentication)
    headers = { "Authorization": f"token {access_token}" }

    # Issue metadata
    data = {
        "title": fields["comment_title"],
        "body":body,
        "labels":["site_comment"],
        "assignees":ASSIGNEES
    }

    # If user clicked "private" checkbox, send to private git repo.
    # Anything other than "false" (including a missing value) stays private, to be safe.
    git_url = GEAR_GIT_URL
    if private == "false":
        # "false" is javascript string "false"
        git_url = GEAR_GIT_URL.replace(PRIVATE_OWNER, PUBLIC_OWNER)

    # Code from https://realpython.com/python-requests/
    try:
        response = requests.post(git_url, json=data, headers=headers, timeout=GITHUB_TIMEOUT)
        # If the response was successful (200- and 300- level status codes), no Exception will be raised
        response.raise_for_status()
    except HTTPError as http_err:
        # Include Github's own message (e.g. "Resource not accessible by personal access token"),
        # which the status line alone does not show
        try:
            github_message = response.json().get("message", "")
        except ValueError:
            github_message = response.text
        result["error"] = f'HTTP error occurred: {http_err}. Github says: {github_message}'
        print(f"Github issue creation failed for {git_url}: {result['error']}", file=sys.stderr)
    except Exception as err:
        result["error"] = f'Other error occurred: {err}'
        print(f"Github issue creation failed for {git_url}: {result['error']}", file=sys.stderr)
    else:
        result["success"] = 1
    print(json.dumps(result))


if __name__ == '__main__':
    main()
