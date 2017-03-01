# standard modules
import gzip
import os
import os.path
import utils
import warnings

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
		return self._options

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
			format = self._config.MOTION.PRINT.TRAJECTORY.FORMAT
			format_dict = {'XYZ': 'xyz', 'DCD': 'dcd'}
			if self._config.MOTION.PRINT.TRAJECTORY.FILENAME is not None:
				outputname = self._config.MOTION.PRINT.TRAJECTORY.FILENAME
			else:
				outputname = self._config.GLOBAL.PROJECT_NAME
			outputfile = os.path.join(self._path, outputname + '-pos-1.' + format_dict.get(format))
			try:
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
