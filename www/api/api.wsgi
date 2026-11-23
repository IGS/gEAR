#!/opt/bin/python3
import os
import sys

# mod_wsgi does not always set sys.executable to the interpreter running this app (on one server
# it was the system /usr/bin/python3). Libraries that start a Python helper process with
# sys.executable need this interpreter and its packages. For example, kaleido (Plotly image export)
# launches its browser through a Python wrapper script that imports kaleido's dependencies.
# See docs/developer/setup/plotly_image_export.md
python_exe = os.path.join(sys.prefix, "bin", f"python{sys.version_info.major}.{sys.version_info.minor}")
if os.access(python_exe, os.X_OK) and os.path.realpath(sys.executable) != os.path.realpath(python_exe):
    sys.executable = python_exe

abs_dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, abs_dir_path)

from api import app as application
