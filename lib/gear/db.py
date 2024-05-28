import mysql.connector
from mysql.connector import errorcode

# This resolves some "no localization support" error
from mysql.connector.locales.eng import client_error

"""

"""

import os
import sys
import configparser

from gear.serverconfig import ServerConfig

class MySQLDB:
    """
    Provides a connection to the MySQL instance.

    Returns: a mysql.connector connection object
    """
    def __init__(self):
        self.connection = None

    def connect(self):
        """
        TODO: Investigate effect of reusing this connection when populated
        vs. creating new ones each time.
        """
        config = ServerConfig().parse()

        try:
            cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                          host=config['database']['host'], database=config['database']['name'],
                                          buffered=True, use_pure=True)
            self.connection = cnx
            return cnx
        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password", file=sys.stderr)
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist", file=sys.stderr)
            else:
                print(err, file=sys.stderr)

