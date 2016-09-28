# standard modules
import os
import utils
import gzip

# third-party modules
import MDAnalysis as mda
import pandas as pd

# custom modules
import cp2k


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
		self._universe = None

	def set_options(self, **kwargs):
		""" Updates the internally used options for the reader subclass.
		:param kwargs: Any options the user may want to configure for the reader subclass.
		:return:
		"""
		self._options.update(kwargs)

	def read(self):
		""" Parse the input and output as well as load the trajectory.

		Any data accessible in the run directory should be stored in attributes of the subclass.
		"""
		raise NotImplementedError()

	def get_input(self):
		""" Fetch the parsed input file data.
		:return: Dict-like or 'utils.Map' instance.
		"""
		raise NotImplementedError()

	def get_output(self):
		""" Fetch the parsed output file data.
		:return: 'pandas.DataFrame'
		"""
		raise NotImplementedError()

	def get_universe(self):
		""" Fetch the trajectory data.

		When no trajectory data is present, this function must return an EmptyUniverse instance.
		:return: 'MDAnalysis.Universe' instance or 'EmptyUniverse' instance.
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


class CP2KReader(Reader):
	""" Parses CP2K Quickstep runs.
	"""
	def __init__(self, *args, **kwargs):
		""" Prepare internal data structure.
		"""
		super(CP2KReader, self).__init__(*args, **kwargs)
		self._config = None
		self._output = pd.DataFrame()

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
		# TODO make configurable
		configname = 'run.inp'
		psffile = os.path.join(self._path, '../input/input.psf')
		logfile = 'run.log.gz'

		# parse input
		fh = open(os.path.join(self._path, configname))
		self._config = self._parse_input_file(fh.readlines())

		# load universe
		outputfile = os.path.join(self._path, self._config.MOTION.PRINT.TRAJECTORY.FILENAME + '-pos-1.dcd')
		try:
			self._universe = mda.Universe(psffile, outputfile)
		except OSError:
			self._universe = EmptyUniverse()

		# parse logfile
		logpath = os.path.join(self._path, logfile)
		if logfile[-3:] == '.gz':
			fh = gzip.GzipFile(logpath)
		else:
			fh = open(logpath)
		c = cp2k.LogFile(fh)
		self._output = c.parse()
