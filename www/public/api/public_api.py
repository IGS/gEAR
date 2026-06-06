import os

from flask import Flask
from routes.datasets import datasets_bp
from routes.status import status_bp

# from routes.users import users_bp  # Example of another noun entity

debug = os.environ.get('DEBUG', False)
debug = debug in ('1', 'true', 'True', 'TRUE')

app = Flask(__name__)

# Register the dataset blueprint with the public API prefix
app.register_blueprint(datasets_bp, url_prefix='/')
app.register_blueprint(status_bp, url_prefix='/')


if __name__ == '__main__':
    app.run(debug=debug, threaded=True)