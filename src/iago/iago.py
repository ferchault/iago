# standard modules
import json
import os
import re


# custom modules
import LocationProvider


def get_location_group(conffile=None):
	lg = LocationProvider.LocationGroup()
	lg.from_file(conffile)
	return lg
