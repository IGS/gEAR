"""Stand-in for gear.spatialhandler (which imports heavy spatial libraries) in upload CGI tests."""

SPATIALTYPE2CLASS = {"visium": object, "visiumhd": object, "xenium": object, "curio": object}
