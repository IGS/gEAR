#!env/bin/python
import os
import sys

# WSGI + Virtual Environment
# Working with Virtual Environments
# http://flask.pocoo.org/docs/1.0/deploying/mod_wsgi/#configuring-apache
abs_dir_path = os.path.dirname(os.path.realpath(__file__))
#rel_path_to_activate_this = '/env/bin/activate_this.py'
#activate_this = ''.join((abs_dir_path, rel_path_to_activate_this))
#with open(activate_this) as file_:
#    exec(file_.read(), dict(__file__=activate_this))
sys.path.insert(0, abs_dir_path)

from api import app as application
