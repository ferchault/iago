import pandas as pd
import utils
import json
import pint


class DB(object):
	def __init__(self):
		self._groups = utils.SafeDict()
		self.planes = utils.AnnotatedDataFrame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Plane name', None),
			'normal_x': ('Normal vector: x component', None),
			'normal_y': ('Normal vector: y component', None),
			'normal_z': ('Normal vector: z component', None),
			'support_x': ('Support point: x component', 'angstrom'),
			'support_y': ('Support point: y component', 'angstrom'),
			'support_z': ('Support point: z component', 'angstrom'),
		})
		self.distances = utils.AnnotatedDataFrame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Plane name', None),
			'atom1': ('First atom index', None),
			'atom2': ('Second atom index', None),
			'dist': ('Distance', 'angstrom')
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
		fh.write(json.dumps(data, separators=(',', ':')))

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
		self.planes = utils.AnnotatedDataFrame(data['planes-meta'], data['planes'])

		# cleanup
		fh.close()


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
