import sys

#from flask import request
#from flask_restful import Resource, reqparse

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
#this = sys.modules[__name__]
#from gear.serverconfig import ServerConfig
#this.servercfg = ServerConfig().parse()

def uploader_callback(session_id, share_id, fh):
    success = 1
    message = ""

    if not fh:
        fh=sys.stderr

    ## Do the work here
   

    return {
        "success": success
        , "message": message
    }
