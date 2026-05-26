= Notes on setting up Jupyter Notebook support

== For pre-existing servers where this feature is being added:

If install gEAR for the first time you can skip this section, as a new install will have these done already.

=== Directory creation:

```bash
export GEARROOT=/home/jorvis/git/gEAR

mkdir $GEARROOT/jupyterhub
mkdir $GEARROOT/jupyterhub/state
mkdir $GEARROOT/jupyterhub/images
mkdir $GEARROOT/jupyterhub/notebooks
mkdir $GEARROOT/jupyterhub/userhomes


```

== For new gEAR installs:

=== Apache configuration

You'll need to create a config file to support the vhost:

/etc/apache2/sites-available/hub-ssl.conf

```
<VirtualHost *:443>
  ServerName hub.umgear.org

  SSLEngine on
  SSLCertificateFile /etc/apache2/ssl_umgear/umgear_ca_bundle.crt
  SSLCertificateKeyFile /etc/apache2/ssl_umgear/server.key

  ProxyPreserveHost On
  ProxyRequests Off

  RequestHeader set X-Forwarded-Proto "https"
  RequestHeader set X-Forwarded-Port "443"

  RewriteEngine On
  RewriteCond %{HTTP:Upgrade} =websocket [NC]
  RewriteRule /(.*) ws://127.0.0.1:8000/$1 [P,L]

  ProxyPass / http://127.0.0.1:8000/
  ProxyPassReverse / http://127.0.0.1:8000/

  ErrorLog /var/log/apache2/jupyterhub_error.log
  CustomLog /var/log/apache2/jupyterhub_access.log combined
</VirtualHost>

```

Then enable needed modules and make sure that vhost is enabled:

```bash
sudo a2enmod proxy proxy_http proxy_wstunnel headers rewrite
audo a2ensite hub-ssl.conf
sudo systemctl reload apache2
```

=== DNS configuration

You need to create a DNS entry for the hub which points to same IP as the main server. For example:

hub.umgear.org  ->  same public IP as www.umgear.org

=== Generate secrets

In the /jupyterhub/.env file change the JUPYTERHUB_CRYPT_KEY and GEAR_LAUNCH_SECRET values using:

```bash
openssl rand -hex 32
```

=== Build Docker images

```bash
cd jupyterhub/
docker build -t gear-notebook:py -f images/python/Dockerfile .
docker build -t gear-notebook:r -f images/r/Dockerfile .

docker compose up -d
docker logs -f gear-jupyterhub
```
