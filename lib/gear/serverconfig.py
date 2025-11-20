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
        ini_path = "{0}/../../gear.ini".format(os.path.dirname(__file__))
        self.config = configparser.ConfigParser()
        self.config.read(ini_path)
        return self.config
