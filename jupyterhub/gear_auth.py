import os
from typing import Any, Dict, Optional

from itsdangerous import URLSafeTimedSerializer, BadSignature, SignatureExpired
from jupyterhub.auth import Authenticator
from jupyterhub.handlers import BaseHandler
from tornado import web


class GearLoginHandler(BaseHandler):
    async def get(self):
        self.log.warning("GearLoginHandler GET called")
        self.log.warning("launch_token present: %s", "launch_token" in self.request.arguments)

        user = await self.login_user()
        self.log.warning("login_user returned: %r", user)

        if user is None:
            raise web.HTTPError(403)

        next_url = self.get_argument("next", default=None)
        if next_url:
            self.redirect(next_url)
        else:
            self.redirect(self.hub.base_url)


class GearLaunchTokenAuthenticator(Authenticator):
    token_param_name = "launch_token"
    max_age_seconds = 3600
    auto_login = True

    def _serializer(self) -> URLSafeTimedSerializer:
        secret = os.environ.get("GEAR_LAUNCH_SECRET")
        if not secret:
            raise RuntimeError("GEAR_LAUNCH_SECRET is not set")
        return URLSafeTimedSerializer(secret, salt="gear-launch")

    def login_url(self, base_url):
        return base_url + "gear-login"

    def get_handlers(self, app):
        return [
            (r"/gear-login", GearLoginHandler),
        ]

    async def check_allowed(self, username, authentication=None):
        self.log.warning("check_allowed called for username=%r", username)
        return True

    async def authenticate(self, handler, data=None) -> Optional[Dict[str, Any]]:
        self.log.warning("authenticate called")
        token = handler.get_argument(self.token_param_name, default=None)
        self.log.warning("token present: %s", bool(token))

        if not token:
            self.log.warning("No token found in request")
            return None

        serializer = self._serializer()

        try:
            payload = serializer.loads(token, max_age=self.max_age_seconds)
            self.log.warning("Token decoded successfully: %r", payload)
        except SignatureExpired:
            self.log.warning("Token expired")
            return None
        except BadSignature:
            self.log.warning("Bad token signature")
            return None
        except Exception as e:
            self.log.warning("Unexpected token decode error: %r", e)
            return None

        username = payload.get("username")
        datasets = payload.get("datasets", [])
        notebook_env = payload.get("notebook_env", "python")
        selected_dataset = payload.get("selected_dataset")

        if not isinstance(username, str) or not username:
            self.log.warning("Invalid username in payload: %r", username)
            return None

        if not isinstance(datasets, list):
            self.log.warning("Invalid datasets payload: %r", datasets)
            return None

        for path in datasets:
            if not isinstance(path, str) or not path.startswith("/"):
                self.log.warning("Invalid dataset path in payload: %r", path)
                return None

        if selected_dataset is not None:
            if not isinstance(selected_dataset, str) or not selected_dataset.startswith("/"):
                self.log.warning("Invalid selected_dataset in payload: %r", selected_dataset)
                return None

        auth_model = {
            "name": username,
            "auth_state": {
                "gear_datasets": datasets,
                "gear_selected_dataset": selected_dataset,
                "gear_payload": payload,
                "gear_notebook_env": notebook_env,
            },
        }

        self.log.warning("Returning auth_model: %r", auth_model)
        return auth_model
