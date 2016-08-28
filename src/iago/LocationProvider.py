import os
import re
import pandas as pd

class LocationProvider(object):
	@staticmethod
	def create(url):
		protocol = url.split(':')[0]
		if protocol == 'file':
			return FileLocationProvider(url[len('file://'):])
		else:
			raise NotImplementedError()

	def get_bucket_list(self):
		raise NotImplementedError()


class FileLocationProvider(LocationProvider):
	def __init__(self, path):
		self._basepath = path

	def get_bucket_list(self):
		directories = [os.path.join(self._basepath, o) for o in os.listdir(self._basepath) if os.path.isdir(os.path.join(self._basepath, o))]
		buckets = [os.path.basename(_) for _ in directories if re.match('^.*-[0-9a-f]{32}$', _)]
		names = [_[:-33] for _ in buckets]
		ids = [_[-32:] for _ in buckets]
		return pd.DataFrame({'name': names, 'id': ids, 'rawname': buckets}).sort_values(by='name')
