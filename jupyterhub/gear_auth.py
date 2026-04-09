import os
from typing import Any, Dict, Optional

from itsdangerous import URLSafeTimedSerializer, BadSignature, SignatureExpired
from jupyterhub.auth import Authenticator


class GearLaunchTokenAuthenticator(Authenticator):
    """
    Authenticate users via a short-lived launch token minted by gEAR.

    Expected token payload:
      {
        "username": "jorvis",
        "datasets": ["/var/www/datasets/DS123.h5ad"],
        "selected_dataset": "/var/www/datasets/DS123.h5ad"   # optional
      }
    """

    token_param_name = "launch_token"
    max_age_seconds = 300

    def _serializer(self) -> URLSafeTimedSerializer:
        secret = os.environ.get("GEAR_LAUNCH_SECRET")
        if not secret:
            raise RuntimeError("GEAR_LAUNCH_SECRET is not set")
        return URLSafeTimedSerializer(secret, salt="gear-launch")

    async def authenticate(self, handler, data=None) -> Optional[Dict[str, Any]]:
        token = handler.get_argument(self.token_param_name, default=None)
        if not token:
            return None

        serializer = self._serializer()

        try:
            payload = serializer.loads(token, max_age=self.max_age_seconds)
        except SignatureExpired:
            return None
        except BadSignature:
            return None

        username = payload.get("username")
        datasets = payload.get("datasets", [])
        selected_dataset = payload.get("selected_dataset")

        if not isinstance(username, str) or not username:
            return None

        if not isinstance(datasets, list):
            return None

        for path in datasets:
            if not isinstance(path, str) or not path.startswith("/"):
                return None

        if selected_dataset is not None:
            if not isinstance(selected_dataset, str) or not selected_dataset.startswith("/"):
                return None

        return {
            "name": username,
            "auth_state": {
                "gear_datasets": datasets,
                "gear_selected_dataset": selected_dataset,
                "gear_payload": payload,
            },
        }
