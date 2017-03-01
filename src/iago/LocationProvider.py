# standard modules
import os
import getpass
import re
import stat
import shlex
import shutil
import subprocess
import tempfile
import warnings
try:
	import ConfigParser as cp
except:
	import configparser as cp

# third-party modules
import pandas as pd
try:
	import paramiko
	HAS_PARAMIKO = True
except:
	HAS_PARAMIKO = False
	pass

# custom modules
import DatabaseProvider


class LocationGroup(object):
	""" Bundles a set of locations to abstract access to raw data no matter where the data is stored.

	Behaves like a set of buckets.
	"""
	def __init__(self):
		self._hosts = dict()
		self._buckets = None

	def from_file(self, filename=None):
		""" Reads the user-defined configuration of storage location.

		:param filename: Filename to read the config from. If not present, `.iago.conf` in the home directory will be used instead.
		"""
		config = cp.ConfigParser()
		if filename is None:
			filenames = [os.path.expanduser('~/.iago.conf')]
		else:
			filenames = [filename]
		config.read(filenames)

		self._hosts = dict()
		for hostalias in config.sections():
			skip = False
			try:
				skip = config.getboolean(hostalias, 'skip')
			except cp.NoOptionError:
				pass
			if not skip:
				self._hosts[hostalias] = config.get(hostalias, 'url')

	def get_bucket_list(self, force_update=False, with_status=False):
		for hostalias, host in self._hosts.iteritems():
			if not isinstance(host, LocationProvider):
				host = LocationProvider.create(host)
				self._hosts[hostalias] = host
			else:
				if not force_update:
					continue
			try:
				df = host.get_bucket_list()
			except ValueError as e:
				warnings.warn('Unable to load host %s: %s' % (hostalias, str(e)))
				continue

			df['hostalias'] = hostalias
			if self._buckets is None:
				self._buckets = df
			else:
				self._buckets = self._buckets[self._buckets.hostalias != hostalias]
				self._buckets = self._buckets.append(df)
		if self._buckets is None:
			return None
		self._buckets = self._buckets.reset_index(drop=True)

		if with_status:
			return self._get_bucket_status(self._buckets)
		else:
			return self._buckets

	def _get_bucket_status(self, bucket_list=None):
		""" Checks for bucket database status.

		:param bucket_list: Pandas dataframe with buckets to query.
		:return: Pandas dataframe.
		"""
		if bucket_list is None:
			bucket_list = self.get_bucket_list()

		bucket_list['hasMeta'] = False
		bucket_list['hasDB'] = False
		for idx, bucket in bucket_list.iterrows():
			lp = self._hosts[bucket.hostalias]
			if lp.has_file(bucket.rawname, '.iago.conf'):
				bucket_list.loc[idx, 'hasMeta'] = True
			if lp.has_file(bucket.rawname, 'iagodb.json'):
				bucket_list.loc[idx, 'hasDB'] = True

		return bucket_list

	def _get_bucket_from_name_or_id(self, bucket):
		if bucket in self._buckets.id.values:
			rows = self._buckets[self._buckets.id == bucket]
		elif bucket in self._buckets.name.values:
			rows = self._buckets[self._buckets.name == bucket]

		if len(rows) > 1:
			raise ValueError('Bucket name or id %s is not unique' % bucket)
		if len(rows) == 0:
			raise ValueError('No bucket with this id or name: %s' % bucket)

		return rows.iloc[0]

	def fetch_database(self, bucket):
		row = self._get_bucket_from_name_or_id(bucket)
		try:
			fh = self._hosts[row.hostalias].open_file(row.rawname, 'iagodb.json')
		except:
			fh = None
		try:
			fn, ephemeral = self._hosts[row.hostalias].local_file(row.rawname, 'iagodb.h5')
			if fh is not None:
				fh.close()
				fh = None
		except:
			fn = None

		if fh is None and fn is None:
			raise RuntimeError('Unable to open remote database.')

		db = DatabaseProvider.DB()
		db.read(handle=fh, name=fn)
		if ephemeral:
			os.remove(fn)
		return db

	def build_database(self, bucket):
		row = self._get_bucket_from_name_or_id(bucket)
		self._hosts[row.hostalias].build_database(row.rawname)


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

		links = {'file': FileLocationProvider, 'ssh': SSHLocationProvider, 'cloud': CloudLocationProvider}
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

	def has_file(self, bucket, filename):
		raise NotImplementedError()

	def open_file(self, bucket, filename):
		raise NotImplementedError()

	def local_file(self, bucket, filename):
		""" If necessary, copies a remote file to local file system.

		:return: Tuple. Absolute filename and flag whether file is ephemeral."""
		raise NotImplementedError()

	def build_database(self, bucket):
		raise NotImplementedError()

	def _get_temp_filename(self):
		fh = tempfile.NamedTemporaryFile()
		name = fh.name
		fh.close()
		return name

class FileLocationProvider(LocationProvider):
	""" Location class for raw file access on the local machine.
	"""
	def __init__(self, path):
		self._basepath = path

	def get_bucket_list(self):
		try:
			directories = [os.path.join(self._basepath, o) for o in os.listdir(self._basepath) if os.path.isdir(os.path.join(self._basepath, o))]
		except OSError:
			raise ValueError('Folder does not exist')
		self._buckets = [os.path.basename(_) for _ in directories if re.match('^.*-[0-9a-f]{32}$', _)]
		return self._buckets_to_df()

	def has_file(self, bucket, filename):
		return os.path.isfile(os.path.join(self._basepath, bucket, filename))

	def open_file(self, bucket, filename):
		return open(os.path.join(self._basepath, bucket, filename))

	def local_file(self, bucket, filename):
		return os.path.join(self._basepath, bucket, filename), False


class SSHLocationProvider(LocationProvider):
	""" Location class for remote file access using SSH.
	"""
	def __init__(self, path):
		if not HAS_PARAMIKO:
			raise RuntimeError('The paramiko python module is required for this location provider.')

		client, self._port, self._sock, self._username, self._hostname, self._basepath = self._prepare_connect(path)
		client.connect(self._hostname, username=self._username, port=self._port, sock=self._sock)

		self._client = client
		self._sftp = None

	@staticmethod
	def _prepare_connect(path):
		""" Reads user configurations related to SSH."""
		match = re.match(r'^(.+)@(.+):(.*)$', path)
		_port = 22
		_sock = None
		if match is None:
			# OpenSSH config
			match = re.match(r'^(.+):(.*)$', path)
			if match is None:
				raise ValueError('Incorrect SSH URL.')
			else:
				_username = getpass.getuser()
				_hostname, _basepath = match.groups()
				ssh_config = paramiko.SSHConfig()
				user_config_file = os.path.expanduser("~/.ssh/config")
				if os.path.exists(user_config_file):
					with open(user_config_file) as fh:
						ssh_config.parse(fh)
				user_config = ssh_config.lookup(_hostname)
				if len(user_config) != 0:
					if 'hostname' in user_config:
						_hostname = user_config['hostname']
					if 'user' in user_config:
						_username = user_config['user']
					if 'port' in user_config:
						_port = user_config['port']
					if 'proxycommand' in user_config:
						_sock = paramiko.ProxyCommand(user_config['proxycommand'])
		else:
			_username, _hostname, _basepath = match.groups()

		client = paramiko.SSHClient()
		client.load_system_host_keys()
		client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		return client, _port, _sock, _username, _hostname, _basepath

	def _connect_sftp(self):
		if self._sftp is not None:
			self._sftp.chdir(self._basepath)
			return
		self._sftp = self._client.open_sftp()
		try:
			self._sftp.chdir(self._basepath)
		except IOError:
			raise ValueError('Path does not exist on SFTP server. Check your SSH URL.')

	def get_bucket_list(self):
		self._connect_sftp()

		directories = []
		for inode in self._sftp.listdir_attr():
			if stat.S_ISDIR(inode.st_mode):
				directories.append(inode.filename)
		self._buckets = [_ for _ in directories if re.match(r'^.*-[0-9a-f]{32}$', _)]
		return self._buckets_to_df()

	def has_file(self, bucket, filename):
		self._connect_sftp()
		try:
			s = self._sftp.stat(bucket + '/' + filename)
			return True
		except:
			return False

	def open_file(self, bucket, filename):
		self._connect_sftp()
		self._sftp.chdir(bucket)
		return self._sftp.file(filename)

	def local_file(self, bucket, filename):
		fn = self._get_temp_filename()
		self._connect_sftp()
		self._sftp.chdir(bucket)
		self._sftp.get(filename, fn)
		return fn, True

	def build_database(self, bucket):
		if not self.has_file(bucket, 'iago-analysis.py'):
			raise RuntimeError('No analysis script specified for this bucket.')

		stdin, stdout, stderr = self._client.exec_command('cd %s; cd %s; python iago-analysis.py' % (self._basepath, bucket))
		if stdout.channel.recv_exit_status() != 0:
			print stderr.read()
			raise RuntimeError('Remote analysis failed.')


class CloudDelegatedFile(object):
	def __init__(self, fh, cloudtemppath):
		self.__dict__['file'] = fh
		self._cloudtemppath = cloudtemppath

	def __close__(self):
		self.file.close(self)
		shutil.rmtree(self._cloudtemppath)

	def __getattr__(self, attr):
		return getattr(self.file, attr)

	def __setattr__(self, attr, value):
		return setattr(self.file, attr, value)


class CloudLocationProvider(LocationProvider):
	""" Location class for remote file access using rclone.
	"""
	def __init__(self, path):
		# cloud://remote:folder
		match = re.match(r'^(.+):(.*)$', path)
		if match is None:
			raise ValueError('Incorrect cloud URL.')

		self._remote, self._basepath = match.groups()

	def get_bucket_list(self):
		directories = []

		command = 'rclone -q --max-depth 1 lsd %s:%s' % (self._remote, self._basepath)
		command = shlex.split(command)
		try:
			output = subprocess.check_output(command)
		except subprocess.CalledProcessError:
			raise RuntimeError('rclone not installed or not configured properly.')

		lines = output.splitlines()
		for line in lines:
			directories.append(line.split()[-1])
		self._buckets = [_ for _ in directories if re.match(r'^.*-[0-9a-f]{32}$', _)]
		return self._buckets_to_df()

	def has_file(self, bucket, filename):
		command = 'rclone -q --max-depth 1 ls %s:%s/%s' % (self._remote, self._basepath, bucket)
		command = shlex.split(command)
		try:
			output = subprocess.check_output(command)
		except subprocess.CalledProcessError:
			raise RuntimeError('rclone not installed or not configured properly.')

		lines = output.splitlines()
		for line in lines:
			if line.split()[-1] == filename:
				return True
		return False

	def open_file(self, bucket, filename):
		tdir = tempfile.mkdtemp()

		raise NotImplementedError()

		return CloudDelegatedFile(fh, tdir)