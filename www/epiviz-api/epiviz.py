from epivizfileserver import setup_app, create_fileHandler, MeasurementManager
from epivizfileserver.trackhub import TrackHub
import os
import numpy
import pickle
import math

# create measurements to load multiple trackhubs or configuration files
mMgr = MeasurementManager()

# create file handler, enables parallel processing of multiple requests
mHandler = create_fileHandler()

# add genome. - for supported genomes 
# check https://obj.umiacs.umd.edu/genomes/index.html
genome = mMgr.add_genome("mm10", url = os.getcwd() + "/genomes/")
genome = mMgr.add_genome("hg19", url = os.getcwd() + "/genomes/")
genome = mMgr.add_genome("hg38", url = os.getcwd() + "/genomes/")
genome = mMgr.add_genome("marmoset", url = os.getcwd() + "/genomes/")

# setup the app from the measurements manager 
# and run the app
app = setup_app(mMgr)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000)