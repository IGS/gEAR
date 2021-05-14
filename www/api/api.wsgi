#!env/bin/python
import os
import sys

abs_dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, abs_dir_path)

from api import app as application
