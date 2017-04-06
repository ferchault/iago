# standard modules
import gzip
import os
import os.path
import utils
import warnings
import re

# third-party modules
try:
	import MDAnalysis as mda
	HAS_MDA = True
except ImportError:
	HAS_MDA = False
	pass
try:
	import mdtraj
	HAS_MDTRAJ = True
except ImportError:
	HAS_MDTRAJ = False
import pandas as pd

# custom modules
import cp2k

class UniverseFactory(object):
	def __init__(self, topology, trajectory):
		self._topology = topology
		self._trajectory = trajectory
		if not HAS_MDA and not HAS_MDTRAJ:
			raise RuntimeError('Neither mdtraj nor MDAnalysis available. Unable to read trajectory data.')


class EmptyUniverse(object):
	""" Class for calculations that do not have trajectory data.

	Mimics the 'MDAnalysis.Universe' interface as far as necessary to be a drop-in replacement for an empty Universe.
	"""
	def __len__(self):
		return 0

	def __init__(self):
		self.trajectory = []


class Reader(object):
	""" Base class for any reader. Defines the API the rest of iago is expecting."""
	def __init__(self, path):
		""" Prepares the data structure.

		:param path: Absolute path to the single run.
		"""
		#: Absolute path to the run directory.
		self._path = path
		#: Holds options of the subclass. No relevance to the base class.
		self._options = dict()
		#: MDAnalysis.Universe instance for trajectory data.
		self._universe = EmptyUniverse()

	def get_options(self):
		""" Lists all available options supported by this reader subclass.

		:return: List of strings.
		"""
		raise NotImplementedError()

	def read(self):
		""" Parse the input and output as well as load the trajectory.

		Any data accessible in the run directory should be stored in attributes of the subclass.
		"""
		raise NotImplementedError()

	def get_input(self):
		""" Fetch the parsed input file data.

		:return: :class:`dict` like or :class:`iago.utils.Map` instance
		"""
		raise NotImplementedError()

	def get_output(self):
		""" Fetch the parsed output file data.

		:return: :class:`pandas.DataFrame`
		"""
		raise NotImplementedError()

	def get_universe(self):
		""" Fetch the trajectory data.

		When no trajectory data is present, this function must return an :class:`iago.Reader.EmptyUniverse` instance.

		:return: :class:`MDAnalysis.core.AtomGroup.Universe` instance or :class:`iago.Reader.EmptyUniverse` instance.
		"""
		return self._universe

	def get_trajectory_frames(self):
		""" Fetch a list of frame numbers for the trajectory.

		Particularly helpful if the log file and the trajectory data are written in different intervals. Frame numbers
		should be consistently one-based or zero-based. Frame numbers coincide with time steps from molecular dynamics
		calculations.

		:return: Iterable.
		"""
		return range(len(self._universe.trajectory))

	def _get_first_file_matching(self, filenames):
		""" Helper function checking for the first existing file in a list.

		Paths can be absolute or relative to the run path. If no file is found, function raises a 'KeyError'.

		:param filenames: Iterable of strings. Filenames to check.
		:return: First entry that points to an existing file as absolute path.
		"""
		for filename in filenames:
			if os.path.isabs(filename):
				mergepath = filename
			else:
				mergepath = os.path.join(self._path, filename)
			if os.path.isfile(mergepath):
				return os.path.realpath(mergepath)
		raise KeyError('No file matches.')

	def claims(self):
		""" Checks whether the run path specified when creating the Reader instance looks like a valid run for this code.

		:return: True if the reader is sure it can handle this run given the filename options."""
		raise NotImplementedError()

class NAMDReader(Reader):
	""" Parses NAMD runs."""
	def __init__(self, *args, **kwargs):
		""" Prepare internal data structure and set the default file paths.
                """
		super(NAMDReader, self).__init__(*args, **kwargs)

		#: List of supported options
		self._options = 'inputnames logs'

		#: List of filenames to test for input files. Precedence in order of the list.
		self.inputnames = ['namd.conf', '*.nmd']

		#: List of logfiles to test for input files. Precedence in order of the list.
		self.logs = ['run.log', 'run.log.gz']

		#: Holds the atom-frame table.
		self._output = pd.DataFrame()

		#: Holds the relevant sections of the input file
		self._config = dict()

	@staticmethod
	def _recognize_value(value):
		try:
			return float(re.match(r'([.0-9]*).*', value).groups()[0])
		#to remove the , after some number and convert the number into a float.
		except:
			return value

	@staticmethod
	def _parse_input_file(lines):
		dict = {'ensemble': 'NVE'}
		for line in lines:
			value=[]
			line = line.strip()
			# check line not remark and not empty line
			if (not line.startswith('#')) and (line):
				#regular expression to match lines do not have ending remark';#'
				#return a list with one string element
				line=re.match(r'([^;#]*);?#?.*',line).groups()
				field = line[0].split()
				# parse all setted parameters
				if field[0] == 'set':
					dict.update({field[1]: NAMDReader._recognize_value(field[2])})
					continue
				# stop reading after first 'minimize' or 'run' command
				if (field[0] == 'minimize') or (field[0] == 'run'):
					if field[1].startswith('$'):
						var = field[1][1:]
						dict.update({field[0]: dict[var]})
					else:
						dict.update({field[0]: NAMDReader._recognize_value(field[1])})
					break
				# detect ensemble
				if field[0] == 'langevin':
					if (field[1] == 'on' and dict['ensemble'] == 'NVE'):
						dict.update({'ensemble': 'NVT'})
				if field[0] == 'langevinPiston':
					if (field[1] == 'on' and dict['ensemble'] == 'NVT'):
						dict.update({'ensemble': 'NPT'})
					# default parse
					# substitute variable with their real value
					# in case have multi entry, value can be a list
				for number in field[1:]:
					if number.startswith('$'):
						var = number[1:]
						value.append(NAMDReader._recognize_value(dict[var]))
					else:
						value.append(NAMDReader._recognize_value(number))
				# if single value entry, convert value list to a number
				if len(value) == 1:
					dict.update({field[0]: value[0]})
				else:
					dict.update({field[0]: value})
		return dict

	@staticmethod
	def _parse_output_file(lines):
		namd_to_iago = dict((
			('TS', 'frame'),
			('BOND', 'bond'),
			('ANGLE', 'angle'),
			('DIHED', 'dihedral'),
			('IMPRP', 'improper'),
			('ELECT', 'electrostatic'),
			('VDW', 'vdw'),
			('BOUNDARY', 'boundary'),
			('MISC', 'external'),
			('KINETIC', 'kinetic'),
			('TOTAL', 'total'),
			('TOTAL2', 'conserved'),
			('TOTAL3', 'NAMD_total3'),
			('POTENTIAL', 'potential'),
			('PRESSURE', 'pressure'),
			('GPRESSURE', 'NAMD_GPRESSURE'),
			('VOLUME', 'volume'),
			('TEMP', 'temperature'),
			('PRESSAVG', 'NAMD_PRESSAVG'),
			('GPRESSAVG', 'NAMD_GPRESSAVG'),
			('TEMPAVG', 'NAMD_TEMPAVG')
		))
		frames = []
		properties = []
		kcal_to_eV = 0.0433634
		energy_columns = 'BOND ANGLE DIHED IMPRP ELECT VDW BOUNDARY MISC KINETIC TOTAL TOTAL2 TOTAL3 POTENTIAL'.split()
		other_columns = 'PRESSURE GPRESSURE VOLUME PRESSAVG GPRESSAVG TEMP TEMPAVG'.split()
		for line in lines:
			line = line.strip()
			if line.startswith('ETITLE:'):
				#to update the order of properties of interest in one line
				properties = line.split()
				continue
			if line.startswith('ENERGY'):
				data = {}
				values = line.split()
				#would assume the first column is always TS
				data['frame'] = int(values[1])
				lookup = dict(zip(properties[2:], values[2:]))
				for energy in energy_columns:
					try:
						data[namd_to_iago[energy]] = float(lookup[energy])*kcal_to_eV
					except KeyError:
						pass
				for other in other_columns:
					try:
						data[namd_to_iago[other]] = float(lookup[other])
					except KeyError:
						pass
				frames.append(data)
		return frames

	def get_options(self):
		return self._options.split()

	def get_input(self):
		return self._config

	def get_output(self):
		return self._output

	def read(self):
		configname = None
		try:
			configname = self._get_first_file_matching(self.inputnames)
		except KeyError:
			warnings.warn('No input file for run in path %s' % self._path)

		logname = None
		try:
			logname = self._get_first_file_matching(self.logs)
		except KeyError:
			warnings.warn('No log file for run in path %s' % self._path)

		# parse input
		if configname is not None:
			fh = open(configname)
			self._config = self._parse_input_file(fh.readlines())
			fh.close()

		# load universe
		if configname is not None:
			try:
				topologyfile = self._config['structure']
				self._get_first_file_matching([topologyfile, ])
			except KeyError:
				warnings.warn('No topology file for run in path %s' % self._path)
				topologyfile = None

			trajectoryfile = self._config['outputname'] + '.dcd'
			if 'DCDfile' in self._config:
				trajectoryfile = self._config['DCDfile']

			try:
				self._universe = mda.Universe(topologyfile, trajectoryfile)
			except OSError:
				self._universe = EmptyUniverse()

		# parse logfile
		if logname is not None:
			if logname[-3:] == '.gz':
				fh = gzip.GzipFile(logname)
			else:
				fh = open(logname)

			if fh is not None:
				configlines = fh.readlines()
				fh.close()
				self._output = pd.DataFrame(self._parse_output_file(configlines))

	def claims(self):
		""" Checks whether this run is indeed a NAMD run.

		:return: True if this certainly is a NAMD run."""
		try:
			logname = self._get_first_file_matching(self.logs)
		except KeyError:
			return False

		line_limit = 50
		if logname is not None:
			if logname[-3:] == '.gz':
				fh = gzip.GzipFile(logname)
			else:
				fh = open(logname)
		for lineno in xrange(line_limit):
			line = fh.readline()
			if line.startswith('Info:') and 'NAMD' in line:
				return True
		return False

class CP2KReader(Reader):
	""" Parses CP2K Quickstep runs.
	"""
	def __init__(self, *args, **kwargs):
		""" Prepare internal data structure and set the default file paths.
		"""
		super(CP2KReader, self).__init__(*args, **kwargs)

		#: List of supported options
		self._options = 'inputnames topologies logs'

		#: List of filenames to test for input files. Precedence in order of the list.
		self.inputnames = ['run.inp', ]

		#: List of topologies to test for input files. Precedence in order of the list.
		self.topologies = ['input.psf', '../input/input.psf', 'input-1.psf']

		#: List of logfiles to test for input files. Precedence in order of the list.
		self.logs = ['run.log', 'run.log.gz']

		#: Holds the atom-frame table.
		self._output = pd.DataFrame()

		#: Holds the tree-based input file
		self._config = utils.Map()

	def get_options(self):
		return self._options.split()

	def get_input(self):
		return self._config

	def get_output(self):
		return self._output

	@staticmethod
	def _recognize_line(line):
		parts = line.split()
		if len(parts) == 1:
			return parts[0], None
		if len(parts) == 2:
			try:
				return parts[0], float(parts[1])
			except ValueError:
				return tuple(parts)
		try:
			return parts[0], tuple(map(float, parts[1:]))
		except ValueError:
			return parts[0], tuple(parts[1:])

	@staticmethod
	def _parse_input_file(lines):
		sections = []
		result = utils.Map()
		for no, line in enumerate(lines):
			line = line.strip()
			if len(line) == 0:
				continue
			if line.startswith('&END'):
				ending = line[4:].strip()
				if ending != '' and (len(sections) == 0 or ending != sections[-1]):
					raise ValueError('Trying to end the %s section with %s' % (sections[-1], ending))
				sections = sections[:-1]
				continue
			if line.startswith('&'):
				section = (line.split()[0][1:]).strip()
				if section in result.traverse(sections):
					if isinstance(result.traverse(sections)[section], utils.Map):
						result.traverse(sections)[section] = [result.traverse(sections)[section], utils.Map()]
					else:
						result.traverse(sections)[section].append(utils.Map())
				else:
					result.traverse(sections)[section] = utils.Map()
				sections.append(section)
				continue
			if not (line.startswith('#') or line.startswith('!')):
				k, v = CP2KReader._recognize_line(line)
				q = result.traverse(sections)
				if k in q:
					if isinstance(q[k], list):
						q[k].append(v)
					else:
						q[k] = [q[k], v]
				else:
					q[k] = v
				continue
		return result

	def read(self):
		configname = None
		try:
			configname = self._get_first_file_matching(self.inputnames)
		except KeyError:
			warnings.warn('No input file for run in path %s' % self._path)

		topname = None
		try:
			topname = self._get_first_file_matching(self.topologies)
		except KeyError:
			warnings.warn('No topology file for run in path %s' % self._path)

		logname = None
		try:
			logname = self._get_first_file_matching(self.logs)
		except KeyError:
			warnings.warn('No log file for run in path %s' % self._path)

		# parse input
		if configname is not None:
			fh = open(configname)
			self._config = self._parse_input_file(fh.readlines())

		# load universe
		if topname is not None and configname is not None:
			try:
				format = self._config.MOTION.PRINT.TRAJECTORY.FORMAT.upper()
			except:
				format = 'XMOL'  # default from CP2K manual
			fileextensions = {'DCD': 'dcd', 'DCD_ALIGNED_CELL': 'dcd', 'PDB': 'pdb', 'XMOL': 'xyz', 'XYZ': 'xyz'}
			if self._config.MOTION.PRINT.TRAJECTORY.FILENAME is not None:
				outputname = self._config.MOTION.PRINT.TRAJECTORY.FILENAME
			else:
				outputname = self._config.GLOBAL.PROJECT_NAME
			try:
				outputfile = os.path.join(self._path, outputname + '-pos-1.' + fileextensions[format])
				self._universe = mda.Universe(topname, outputfile)

			except OSError:
				self._universe = EmptyUniverse()

		# parse logfile
		if logname is not None:
			openmethod = open
			if logname[-3:] == '.gz':
				openmethod = gzip.GzipFile
			fh = openmethod(logname)
			if fh is not None:
				c = cp2k.LogFile(fh)
				self._output = c.parse()

	def claims(self):
		""" Checks whether this run is indeed a CP2K run.

		:return: True if this certainly is a CP2K run."""
		try:
			logname = self._get_first_file_matching(self.logs)
		except KeyError:
			return False

		line_limit = 50
		if logname is not None:
			if logname[-3:] == '.gz':
				fh = gzip.GzipFile(logname)
			else:
				fh = open(logname)
		for lineno in xrange(line_limit):
			line = fh.readline()
			if line.startswith(' CP2K|'):
				return True
		return False