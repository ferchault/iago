import pandas as pd
import os
import re
import json
import glob

import cp2k
import DatabaseProvider
import LocationProvider
import StorageProvider

class Config(dict):
	def __init__(self):
		self['basepath'] = '.'

	def __getattr__(self, item):
		return self[item]


class Store(object):
	def __init__(self):
		hierarchy = os.getcwd().split(os.sep)
		self._bucket = [_ for _ in hierarchy if re.match('^.*-[0-9a-f]{32}$', _)][0]
		self._basepath = os.sep + os.path.join(*hierarchy[:hierarchy.index(self._bucket) + 1])

	def load_database(self, clear=False):
		try:
			fh = open(os.path.join(self._basepath, 'db.json'))
			data = json.load(fh)
		except:
			clear = True
		if clear:
			data = dict()
		df = pd.DataFrame(data.pop('ts', {}))
		return data, df

	def save_database(self, db, df):
		db['ts'] = df.to_dict()
		json.dump(db, open(os.path.join(self._basepath, 'db.json'), 'w'))

	def get_basepath(self):
		return self._basepath




def get_location_group():
	lg = LocationProvider.LocationGroup()
	lg.from_file()
	return lg

# def _load_cp2k(self, data, df):
# 	cp = cp2k.LogFileReader()
# 	ci = cp2k.InputFileReader()
# 	data['input'] = {}
# 	dfs = [df]
# 	for run in _find_runs():
# 		r = cp.read_run(self, run)
# 		i = ci.read_run(self, run)
# 		data['input'][run] = i
# 		dfs.append(r)
# 	df = pd.concat(dfs).reset_index(drop=True)
# 	return data, df

#def _find_runs():
#	runs = glob.glob(os.path.join(self._basepath, 'run-*'))
#	return [os.path.basename(_) for _ in runs if os.path.isdir(_)]

def parse(code_type):
	if code_type == 'cp2k':
		pass
	else:
		raise NotImplementedError()