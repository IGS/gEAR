from flask import Flask
from routes.datasets import datasets_bp
from routes.status import status_bp

# from routes.users import users_bp  # Example of another noun entity

app = Flask(__name__)

# Register the dataset blueprint with the public API prefix
app.register_blueprint(datasets_bp, url_prefix='/api/public')
app.register_blueprint(status_bp, url_prefix='/api/public')

if __name__ == '__main__':
    app.run()