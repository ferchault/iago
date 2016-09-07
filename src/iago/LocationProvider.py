import os
import re
import stat
try:
	import ConfigParser as cp
except:
	import configparser as cp
import pandas as pd
try:
	import paramiko
	HAS_PARAMIKO = True
except:
	HAS_PARAMIKO = False
	pass


class LocationGroup(object):
	""" Bundles a set of locations to abstract access to raw data no matter where the data is stored. Behaves like a set of buckets.
	"""
	def __init__(self):
		self._hosts = dict()
		self._buckets = None

	def from_file(self):
		config = cp.ConfigParser()
		config.read([os.path.expanduser('~/.iago.conf')])

		self._hosts = dict()
		for hostalias in config.sections():
			self._hosts[hostalias] = config.get(hostalias, 'url')

	def get_bucket_list(self, force_update=False):
		for hostalias, host in self._hosts.iteritems():
			if not isinstance(host, LocationProvider):
				host = LocationProvider.create(host)
				self._hosts[hostalias] = host
			else:
				if not force_update:
					continue

			df = host.get_bucket_list()
			df['hostalias'] = hostalias
			if self._buckets is None:
				self._buckets = df
			else:
				self._buckets = self._buckets[self._buckets.hostalias != hostalias]
				self._buckets = self._buckets.append(df)
		self._buckets.reset_index()
		return self._buckets


class LocationProvider(object):
	""" Base class for all supported protocols.

	Implements the factory pattern for new locations.
	"""
	@staticmethod
	def create(url):
		try:
			protocol, url = re.match(r'^([^:]*)://(.*)$', url).groups()
		except:
			raise ValueError('Invalid location URL.')

		links = {'file': FileLocationProvider, 'ssh': SSHLocationProvider}
		if protocol not in links:
			raise NotImplementedError()

		return links[protocol](url)

	def get_bucket_list(self):
		raise NotImplementedError()

	def _buckets_to_df(self):
		names = [_[:-33] for _ in self._buckets]
		ids = [_[-32:] for _ in self._buckets]
		return pd.DataFrame({'name': names, 'id': ids, 'rawname': self._buckets}).sort_values(by='name')

	def __init__(self):
		self._buckets = None


class FileLocationProvider(LocationProvider):
	""" Location class for raw file access on the local machine.
	"""
	def __init__(self, path):
		self._basepath = path

	def get_bucket_list(self):
		directories = [os.path.join(self._basepath, o) for o in os.listdir(self._basepath) if os.path.isdir(os.path.join(self._basepath, o))]
		self._buckets = [os.path.basename(_) for _ in directories if re.match('^.*-[0-9a-f]{32}$', _)]
		return self._buckets_to_df()


class SSHLocationProvider(LocationProvider):
	""" Location class for remote file access using SSH.
	"""
	def __init__(self, path):
		if not HAS_PARAMIKO:
			raise RuntimeError('The paramiko python module is required for this location provider.')

		match = re.match(r'^(.+)@(.+):(.*)$', path)
		if match is None:
			raise ValueError('Incorrect SSH URL.')

		self._username, self._hostname, self._basepath = match.groups()

		client = paramiko.SSHClient()
		client.load_system_host_keys()
		client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		client.connect(self._hostname, username=self._username)

		self._client = client
		self._sftp = None

	def _connect_sftp(self):
		if self._sftp is not None:
			return
		self._sftp = self._client.open_sftp()

	def get_bucket_list(self):
		self._connect_sftp()
		try:
			self._sftp.chdir(self._basepath)
		except IOError:
			raise ValueError('Path does not exist on SFTP server. Check your SSH URL.')

		directories = []
		for inode in self._sftp.listdir_attr():
			if stat.S_ISDIR(inode.st_mode):
				directories.append(inode.filename)
		self._buckets = [_ for _ in directories if re.match('^.*-[0-9a-f]{32}$', _)]
		return self._buckets_to_df()
