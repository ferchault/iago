import utils
import os
import MDAnalysis as mda


class EmptyUniverse(object):
	def __len__(self):
		return 0

	def __init__(self):
		self.trajectory = []


class Reader(object):
	def __init__(self, path):
		self._path = path
		self._options = dict()
		self._universe = None

	def set_options(self, **kwargs):
		self._options.update(kwargs)

	def read(self):
		raise NotImplementedError()

	def get_universe(self):
		return self._universe

	def get_trajectory_frames(self):
		return range(len(self._universe.trajectory))


class CP2KReader(Reader):
	def __init__(self, *args, **kwargs):
		super(CP2KReader, self).__init__(*args, **kwargs)
		self._config = None

	def _recognize_line(self, line):
		parts = line.split()
		if len(parts) == 1:
			return parts, None
		if len(parts) == 2:
			try:
				return parts[0], float(parts[1])
			except:
				return parts
		try:
			return parts[0], map(float, parts[1:])
		except:
			return parts[0], parts[1:]

	def _parse_input_file(self, lines):
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
				result.traverse(sections)[section] = utils.Map()
				sections.append(section)
				continue
			if not (line.startswith('#') or line.startswith('!')):
				k, v = self._recognize_line(line)
				result.traverse(sections)[k] = v
				continue
		return result

	def read(self):
		# TODO make configurable
		configname = 'run.inp'
		psffile = os.path.join(self._path, '../input/input.psf')

		# parse input
		fh = open(os.path.join(self._path, configname))
		self._config = self._parse_input_file(fh.readlines())

		# load universe
		outputfile = os.path.join(self._path, self._config.MOTION.PRINT.TRAJECTORY.FILENAME + '-pos-1.dcd')
		try:
			self._universe = mda.Universe(psffile, outputfile)
		except OSError:
			self._universe = EmptyUniverse()

