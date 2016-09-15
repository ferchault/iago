import pandas as pd
import utils
import json
import pint


class DB(object):
	def __init__(self):
		ureg = pint.UnitRegistry()

		self._groups = utils.SafeDict()
		self.planes = utils.AnnotatedDataFrame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Plane name', None),
			'normal_x': ('Normal vector: x component', None),
			'normal_y': ('Normal vector: y component', None),
			'normal_z': ('Normal vector: z component', None),
			'support_x': ('Support point: x component', ureg.angstrom),
			'support_y': ('Support point: y component', ureg.angstrom),
			'support_z': ('Support point: z component', ureg.angstrom),
		})

	@property
	def groups(self):
		""" Known static atom groups.
		:return: Dictionary with group names as key, list of 0-based atom indices as value.
		"""
		return self._groups

	def write(self, fh):
		""" Writes the database to disk or stream.
		:param fh: File handle or filename.
		:return:
		"""
		if not hasattr(fh, 'write'):
			fh = open(fh, 'w')

		data = dict()

		# groups
		data['groups'] = self.groups
		# planes
		data['planes'] = self.planes.to_dict()
		data['planes-meta'] = self.planes.annotations_to_dict()

		# finalise
		fh.write(json.dumps(data, separators=(',',':')))

	def read(self, fh):
		""" Reads the database from disk or stream.
		:param fh: File handle or filename.
		:return:
		"""
		if not hasattr(fh, 'write'):
			fh = open(fh, 'w')

		data = json.load(fh)

		# groups
		self._groups = utils.SafeDict(data['groups'])
		# planes
		self.planes = pd.DataFrame(data['planes-meta'], data['planes'])

class Map(dict):
	def __init__(self, *args, **kwargs):
		super(Map, self).__init__(*args, **kwargs)
		for arg in args:
			if isinstance(arg, dict):
				for k, v in arg.iteritems():
					self[k] = v

		if kwargs:
			for k, v in kwargs.iteritems():
				self[k] = v

	def __getattr__(self, attr):
		return self.get(attr)

	def __setattr__(self, key, value):
		self.__setitem__(key, value)

	def __setitem__(self, key, value):
		super(Map, self).__setitem__(key, value)
		self.__dict__.update({key: value})

	def __delattr__(self, item):
		self.__delitem__(item)

	def __delitem__(self, key):
		super(Map, self).__delitem__(key)
		del self.__dict__[key]

	def traverse(self, keys):
		return reduce(lambda d, k: d[k], keys, self)


class DatabaseProvider(object):
	@staticmethod
	def get_annotated_collections():
		return 'atommeta', 'atomframemeta', 'framemeta'


class MemoryDatabaseProvider(DatabaseProvider):
	def __init__(self, run, atom, atommeta, atomframe, atomframemeta, frame, framemeta):
		self.run = Map(run)
		self.atom = pd.DataFrame(atom)
		self.atommeta = Map(atommeta)
		self.atomframe = pd.DataFrame(atomframe)
		self.atomframemeta = Map(atomframemeta)
		self.frame = pd.DataFrame(frame)
		self.framemeta = Map(framemeta)
