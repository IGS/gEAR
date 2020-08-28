#!/opt/bin/python3

import cgi, json
import configparser
import os, sys
import mysql.connector
from mysql.connector import errorcode

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    display_id = form.getvalue('id')
    user_id = form.getvalue('user_id')
    dataset_id = form.getvalue('dataset_id')
    label = form.getvalue('label')
    plot_type = form.getvalue('plot_type')
    plotly_config = form.getvalue('plotly_config')

    config = configparser.ConfigParser()
    config.read('../../gear.ini')

    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'],
                                      buffered=True)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)

    cursor = cnx.cursor()

    if display_id:
        # display_id exists, so update

        query = "SELECT user_id FROM dataset_display WHERE id = %s"
        cursor.execute(query, (display_id,))
        (display_owner,) = cursor.fetchone()

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
            result = dict(success=True)
        else:
            print('UPDATE DIDNT HAPPEN?', file=sys.stderr)
            result = dict(success=False)
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

    cnx.commit()
    cursor.close()
    cnx.close()

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
