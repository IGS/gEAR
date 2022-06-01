## Apache configuration notes for gEAR instances

This document is intended to capture all the customizations to the apache2 config files needed for a gEAR instance to operate properly.  The paths given are those for a standard Ubuntu installation, but adjust as needed for another platform.

### Enabling mod-rewrite, CGI, PHP and WSGI

    $ sudo a2enmod rewrite
    $ sudo a2enmod cgi

    $ sudo apt install libapache2-mod-wsgi apache2-dev
    $ sudo a2enmod wsgi
    $ sudo a2enmod proxy
    $ sudo a2dismod mpm_event && sudo a2enmod mpm_prefork && sudo a2enmod php8.1
    
### PHP (used by uploader)

    $ sudo apt install php

The php.ini file also needs updating on some systems to get around file upload limitations. Visit
the test PHP page in gEAR, look for the 'Configuration file (php.ini) Path' setting, then modify
that one.

     http://localhost/cgi/test.php

Modifications to make:

    post_max_size = 3000M
    upload_max_filesize = 3000M



### /etc/apache2/apache2.conf

Aside from the default things, these sections are important

    # This is needed starting in scanpy 1.4.x because the Numba module attempts to write
    #  cacha data where it doesn't have permission due to brazen developers
    #  assuming everyone will use virtualenv all the time.
    #  https://github.com/ska-sa/montblanc/issues/253
    SetEnv NUMBA_CACHE_DIR /tmp

    <Directory /var/www/>
            Options FollowSymLinks
            AllowOverride None
            Require all granted

            # Need to do this first: sudo a2enmod rewrite
            RewriteEngine on
            #LogLevel alert rewrite:trace2
            RewriteCond %{DOCUMENT_ROOT}%{REQUEST_FILENAME} !-f
            RewriteRule ^(.+)\.(\d+)\.(js|css)$ $1.$3 [L]
    </Directory>

    #Enable redirector script 'p' in www dir
    <FilesMatch ^p$>
            Options ExecCGI FollowSymLinks
            SetHandler cgi-script
    </FilesMatch>

    WSGIPythonHome "/opt/Python-3.10.4"
    LoadModule wsgi_module "/opt/Python-3.10.4/lib/python3.10/site-packages/mod_wsgi/server/mod_wsgi-py310.cpython-310-x86_64-linux-gnu.so"

### /etc/apache2/sites-available/000-default.conf

Assuming you have SSL setup (below) you want to make sure apache redirects http traffic
to https, so you need to add this line to the <VirtualHost *:80> block, adjusted for your
domain, of course:

    Redirect permanent / https://umgear.org/

### Disable apache's PrivateTmp

gEAR lets users write datafiles such as analyses in an area under /tmp until they want to
save, when they are moved to a directory within the application.  If Apache has PrivateTmp
on this causes an error like this:

   OSError: [Errno 18] Invalid cross-device link

For now, the solution is to disable this.  If you run this:

    $ systemctl show apache2 | grep PrivateTmp

And you get this:

    PrivateTmp=yes

You need to turn it off.  On Ubuntu 20, you can find the setting in this file:

    /usr/lib/systemd/system/apache2.service

Set it to false, then restart apache.

### Victor's SSL instructions

    $ sudo a2enmod ssl
    $ sudo service apache2 restart

Copied /etc/apache2/ssl_umgear from production server to this server
Created /etc/apache2/sites-available/umgear-ssl.conf file with contents from production server
In /etc/apache2/sites-available/umgear-ssl.conf, edited servername, removed server aliases

Used a2ensite to enable the umgear-ssl site.

     $ sudo a2ensite
     $ sudo service apache2 restart

#### Custom config needed for Flask API

Resources:
- https://pypi.org/project/mod_wsgi/
- http://flask.pocoo.org/docs/1.0/deploying/mod_wsgi/
- https://www.digitalocean.com/community/tutorials/how-to-deploy-a-flask-application-on-an-ubuntu-vps

### /etc/apache2/sites-available/umgear-ssl.conf

This needs to be tailored for each machine's resources to match the processors (cores) present
and number of threads within each process.  If not running under SSL, this goes in 000-default.conf

There's a lot in here, but the CGI-related addition is:

    <Directory /var/www/cgi>
         Options ExecCGI FollowSymLinks
         SetHandler cgi-script
    </Directory>

    # FLASK API
    WSGIDaemonProcess api user=www-data group=www-data processes=8 threads=2
    WSGIScriptAlias /api /var/www/api/api.wsgi
    <Directory /var/www/api>
            WSGIProcessGroup api
            WSGIApplicationGroup %{GLOBAL}
            WSGIScriptReloading On
            <IfVersion < 2.4>
                    Order allow,deny
                    Allow from all
            </IfVersion>
            <IfVersion >= 2.4>
                     Require all granted
            </IfVersion>
    </Directory>

### /etc/apache2/mods-available/wsgi.load

## The version numbers here need to coincide with the python version installed

LoadModule wsgi_module "/opt/Python-3.10.4/lib/python3.10/site-packages/mod_wsgi/server/mod_wsgi-py310.cpython-310-x86_64-linux-gnu.so"

### /etc/apache2/mods-enabled/wsgi.conf

Add the line `WSGIPythonHome "/opt/Python-3.10.4"` into the IfModule block.

Then, finally restart apache again.

      $ sudo service apache2 restart

## Optional configurations

I don't like the new private tmp directories employed lately, as they don't servive service restarts
and we make use of the tmp area while the users are building pipelines.  Turn this off by editing
the following file:

### /etc/systemd/system/multi-user.target.wants/apache2.service

And then setting this line:

    PrivateTmp=false


## Common errors

If you get errors like this:

```
RuntimeError: cannot cache function '_downsample_array': no locator available for
file '/home/jorvis/git/gEAR/www/api/env/lib/python3.7/site-packages/scanpy/preprocessing/_simple.py'
```

The solution is to chmod 777 the __pycache__ directory in that same folder.

If you get errors like module can't be loaded:

```
ModuleNotFoundError: No module named 'encodings'
```

WSGI is probably not loading the correct Python.  Make sure the following are pointing to the right env:

/etc/apache2/mods-available/wsgi*
