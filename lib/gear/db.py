import mysql.connector
from mysql.connector import errorcode
"""

"""

import os
import sys
import configparser

class MySQLDB:
    """
    Provides a connection to the MySQL instance.

    Returns: a mysql.connector connection object
    """
    def __init__(self):
        connection = None

    def connect(self):
        ini_path = "{0}/../../gear.ini".format(os.path.dirname(__file__))
        config = configparser.ConfigParser()
        config.read(ini_path)
        
        try:
            cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                          host=config['database']['host'], database=config['database']['name'], buffered=True)
            
            return cnx
        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password", file=sys.stderr)
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist", file=sys.stderr)
            else:
                print(err, file=sys.stderr)

