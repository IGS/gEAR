# Setting up HiGlass VM for epigenomics

There is currently one VM used for HiGlass services, and it is being used across all gEAR platforms.  Currently it stores annotation files and Hi-C chromatin files.  Please see that [documentation](../misc/higlass_file_upload.md) for how to prep and upload those files.

## Setting up HiGlass VM

* E2-small VM (2 CPUs, 2 Gb memory)
* 20Gb storage
* For higlass-manage-prod I made a persistent image
  * URL [http://34.75.224.1:8989](http://34.75.224.1:8989/admin/)
  * This is the static external URL
  * domain will be [https://higlass.umgear.org:8989](https://higlass.umgear.org:8989)
* Chose Debian Bookwork (just went with defaults here)
* For the network, export TCP port 8989 across IPv4 0.0.0.0/0
  * The URL will be http (i.e.[https://higlass.umgear.org/app](https://higlass.umgear.org/app))
  * Gemini suggested I remove this after creating the load balancer below but it seems this prevents access within the VM (blocked the load balancer from talking to the server
  * TODO: Only allow from certain Google IPs
    * Load Balancer Health-check probes
    * `130.211.0.0/22` and `35.191.0.0/16` (VERIFY)
* Used a Google Cloud External Application Global Load Balancer to redirect HTTPS (443) to HTTP under port 8989\.
  * Made a "higlass" backend instance group as the instance group was what the load balancer connects to for its backend.
  * Google managed certificates were created under the domain [higlass.umgear.org](http://higlass.umgear.org) (which Joshua had to add as a DNS CNAME A-type)

* Debian would not let me use "pip" to install packages system-wide so need "venv"
* Install python3-pip, python3-venv git
* `python3 -m venv higlass`
* `cd ~/higlass/bin; source ./activate`
* The next steps are not traditional due to various errors from version differences
* `cd ~; git clone https://github.com/higlass/higlass-manage.git`
  * We want to use the scripts from this repo as they are more recently updated (v0.8.2) but not in the PyPI version (v0.8.0) which seems to not work for me.
  * Make sure the "higlass" python package is not installed as it is a predecessor to higlass-python
* `cd higlass-manage; pip install -r requirements.txt`
* `pip install -U higlass-python==1.3.4`
* `python setup.py install`
  * This installs 0.8.2 of higlass-manage
  * This step is optional I think.
* `pip install docker==7.1.0`
  * This fixes a bug in the requirements.txt version of docker related to the requests module not supporting the http+docker URL scheme
* Install docker ([https://docs.docker.com/engine/install/debian/](https://docs.docker.com/engine/install/debian/))
* To get test files cd into \~/higlass-manage and run "./get\_test\_data.sh".  This will pull test files into `~/hg-data/media/uploads/` which you can put wherever, like in the home directory.
* Run ./[test.sh](http://test.sh) to test ingesting (to ensure things are set up correctly)
* To start the server `higlass-manage start –use-redis –site-url=higlass.umgear.org –version=local`
  * Test with [https://higlass.umgear.org/api/v1/tilesets/](https://higlass.umgear.org/api/v1/tilesets/) (replace the IPv4)
  * The instructions for API calls do not have a slash at the end, but I found adding the slash is the only way I can view results
  * If performing some debugging, create a DEBUG=True environment variable so that error lines are not cached
  * Add `--use-redis` to use the Redis caching service for performance.  By default it stores in \~/redis-data.
  * By default temp data is stored in /tmp/higlass-docker and files are stored in \~/hg-data
  * The documentation suggests using the "site-url" option but I have not seen a change with or without it
  * I added a Dockerfile.extended file to the home directory that updates uWSGI timeout parameters to allow for larger HiC datasets to be ingested without a 504 Gateway Timeout. When building the image must be named "image-default" and the `--version local` must be passed.
  * Sometimes I get a message about no ReadTimout (which is a typo in their code). Everything still seems to work fine despite that, but I fixed the type in the higlass-manage/[start.py](http://start.py) code anyways, which gives a successful start
* Now that Higlass is started, the docker container is running.  Unfortunately the docker image has some issues with clodius.
  * `higlass-manage shell`
  * replicates a docker exec /bin/bash
  * If you write print statements, check the /tmp/uwsgi-stderr-supervisor-\*.log. These prints will not show when you do `higlass-manage logs`
* Easy way to ingest a dataset
  * `higlass-manage ingest <file>`
  * This will print the "ingest tileset" CLI command with the UUID inside the "-uid" argument
* Create superuser
  * `higlass-manage create superuser`
    * User \- gear
    * Email \- [sadkins@som.umaryland.edu](mailto:sadkins@som.umaryland.edu)
    * Pass \- gearadmin
  * We can now use this to upload data via POST using curl.  Pass in `-u gear:gearadmin` as user credentials in the request.
  * This also lets us use the admin panel at [https://higlass.umgear.org/aadmin](https://higlass.umgear.org/admin)

## API endpoint URLs

* /tilesets/
  * Multiple datasets
  * You can find the UUID in each entry
* /tilesets/\<uuid\>
  * Individual dataset properties
* /tileset\_info/?d=\<uuid\>
  * Individual dataset
  * [https://docs.higlass.io/jupyter.html?highlight=tileset\_info\#serving-custom-data](https://docs.higlass.io/jupyter.html?highlight=tileset_info#serving-custom-data)
  * I believe this is what is used for Gosling
* /tiles/?d\<uuid\>.z.x
  * Tile information for an individual dataset at zoom level "z" and position "x"
  * [https://docs.higlass.io/jupyter.html?highlight=tileset\_info\#serving-custom-data](https://docs.higlass.io/jupyter.html?highlight=tileset_info#serving-custom-data)
  * There are optional "dot" parameters too
  * For a 2D matrix you need to add a position "y" as well. So z.x.y