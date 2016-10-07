# standard modules
import utils
import json

# third-party modules
import pandas as pd


class DB(object):
	def __init__(self):
		self._groups = utils.SafeDict()

		#: Plane data as calculated by :func:`iago.Analyser.Analyser.dynamic_plane`
		self.planes = utils.annotated_data_frame({
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
		self.distances = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Distance set name', None),
			'atom1': ('First atom index', None),
			'atom2': ('Second atom index', None),
			'dist': ('Distance', 'angstrom')
		})
		self.planedistances = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Distance set name', None),
			'plane': ('Plane name', None),
			'atom1': ('First atom index', None),
			'dist': ('Distance', 'angstrom')
		})
		self.input = utils.Map()
		self.output = pd.DataFrame()

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
		data['planes'] = self.planes.to_dict('list')
		data['planes-meta'] = self.planes.annotations_to_dict()
		# distances
		data['distances'] = self.distances.to_dict('list')
		data['distances-meta'] = self.distances.annotations_to_dict()
		data['planedistances'] = self.planedistances.to_dict('list')
		data['planedistances-meta'] = self.planedistances.annotations_to_dict()
		# input / output
		data['input'] = self.input
		data['output'] = self.output.to_dict('list')
		# data['output-meta'] = self.output.annotations_to_dict()

		# finalise
		fh.write(json.dumps(data, separators=(',', ':')))

	def read(self, fh):
		""" Reads the database from disk or stream.
		:param fh: File handle or filename.
		:return:
		"""
		if not hasattr(fh, 'read'):
			fh = open(fh)

		data = json.load(fh)

		# groups
		try:
			self._groups = utils.SafeDict(data['groups'])
		except KeyError:
			pass
		# planes
		try:
			self.planes = utils.annotated_data_frame(data['planes-meta'], data['planes'])
		except KeyError:
			pass
		# distances
		try:
			self.distances = utils.annotated_data_frame(data['distances-meta'], data['distances'])
		except KeyError:
			pass
		try:
			self.planedistances = utils.annotated_data_frame(data['planedistances-meta'], data['planedistances'])
		except KeyError:
			pass
		# input / output
		try:
			self.input = utils.Map(data['input'])
		except KeyError:
			pass
		# self.output = utils.AnnotatedDataFrame(data['output-meta'], data['output'])
		try:
			self.output = pd.DataFrame.from_dict(data['output'])
		except KeyError:
			pass

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
