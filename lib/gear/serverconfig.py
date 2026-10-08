"""
serverconfig.py - Access to the gEAR server configuration.

Provides ServerConfig, which parses the gear.ini file at the repository root.
"""

import os
import sys
import configparser

class ServerConfig:
    """
    Provides the parsed version of the gear.ini file

    Returns: a configparser.ConfigParser() object
    """
    def __init__(self):
        self.config = None

    def parse(self) -> configparser.ConfigParser:
        """
        Read gear.ini and return the parsed configuration.
        """
        ini_path = "{0}/../../gear.ini".format(os.path.dirname(__file__))
        self.config = configparser.ConfigParser()
        self.config.read(ini_path)
        return self.config
