#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
import mysql.connector


def attempt_symlink(cursor, user_id, dataset_id, display_id, is_multigene):
    """Attempt to create a symlink if the user of the saved display is also the dataset owner."""

    DATASET_PREVIEWS_DIR = "/var/www/img/dataset_previews"

    gene = "single"
    if is_multigene > 0:
        gene = "multi"

    filename = os.path.join(DATASET_PREVIEWS_DIR, "{}.{}.png".format(dataset_id, display_id))

    if not os.path.isfile(filename):
        print("File to static image does not exist.  Skipping attempted symlink.", file=sys.stderr)
        return

    # Symlink as a default static image if dataset owner curated the image
    defaults_query = """
        SELECT d.id AS dataset_id, dp.display_id AS display_id
        FROM dataset d, dataset_preference dp
        WHERE d.id = dp.dataset_id
        AND dp.user_id = d.owner_id
        AND dp.user_id = %s AND d.id = %s AND dp.display_id = %s
    """
    cursor.execute(defaults_query, (user_id, dataset_id, display_id))

    # If a row was returned, dataset owner saved the display.  Symlink to newest display
    row = cursor.fetchone()
    if row:
        symlink_path = os.path.join(DATASET_PREVIEWS_DIR, "{}.{}.default.png".format(dataset_id, gene))

        # If symlink exists, we want to attempt to remove, so that it can be remade again (to point to a new display ID if applicable)
        try:
            os.remove(symlink_path)
        except:
            pass

        os.symlink(filename, symlink_path)
    else:
        print("Will not create 'default static image' symlink for display id {}. " \
            "Most likely the user requesting to make the display default is not the dataset owner".format(display_id), file=sys.stderr)


def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    user_id = form.getvalue('user_id')
    dataset_id = form.getvalue('dataset_id')
    display_id = form.getvalue('display_id')
    is_multigene = int(form.getvalue('is_multigene', 0))

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    try:
        query = """
            INSERT INTO dataset_preference (user_id, dataset_id, display_id, is_multigene)
            VALUES (%s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE display_id = VALUES(display_id);
        """
        cursor.execute(query, (user_id, dataset_id, display_id, is_multigene))
        result = dict(success=True)
    except mysql.connector.Error as err:
        result = dict(success=False)

    attempt_symlink(cursor, user_id, dataset_id, display_id, is_multigene)

    cnx.commit()
    cursor.close()

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
