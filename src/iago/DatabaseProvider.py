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

		#: Atom-atom distance data as calculated by :func:`iago.Analyser.Analyser.dynamic_distance`
		self.distances = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Distance set name', None),
			'atom1': ('First atom index', None),
			'atom2': ('Second atom index', None),
			'dist': ('Distance', 'angstrom')
		})

		#: Atom-plane distance data as calculated by :func:`iago.Analyser.Analyser.dynamic_distance`
		self.planedistances = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Distance set name', None),
			'plane': ('Plane name', None),
			'atom1': ('First atom index', None),
			'dist': ('Distance', 'angstrom')
		})

		self.energies = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'total': ('Total energy', 'hartree'),
			'conserved': ('Conserved quantity', 'hartree'),
			'coreself': ('Core-Self energy', 'hartree'),
			'corehamiltonian': ('Core Hamiltonian', 'hartree'),
			'hartree': ('Hartree energy', 'hartree'),
			'xc': ('Exchange-Correlation energy', 'hartree'),
			'hfx': ('Hartree-Fock Exchange energy', 'hartree'),
			'dispersion': ('Dispersion energy', 'hartree'),
			'potential': ('Potential energy', 'hartree'),
			'kinetic': ('Kinetic energy', 'hartree'),
			'drift': ('Energy drift per atom', 'kelvin')
		})

		self.cells = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'a': ('First cell length', 'angstrom'),
			'b': ('Second cell length', 'angstrom'),
			'c': ('Third cell length', 'angstrom'),
			'alpha': ('First cell angle', 'degrees'),
			'beta': ('Second cell angle', 'degrees'),
			'gamma': ('Third cell angle', 'degrees'),
			'volume': ('Cell volume', 'angstrom**3'),
		})

		self.ensembles = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'temperature': ('Temperature', 'kelvin'),
			'pressure': ('Pressure', 'bar'),
		})

		self.meta = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'iasd': ('Integrated absolute spin density', None),
			's2': ('Determinant S**2', None),
			'scfcycles': ('Number of SCF cycles', None),
			'otcycles': ('Number of outer SCF cycles', None),
			'globaleri': ('Number of ERI evaluated', None)
		})

		self.points = utils.annotated_data_frame({
			'run': ('Run', None),
			'frame': ('Frame number', None),
			'name': ('Point set name', None),
			'x': ('X position', 'angstrom'),
			'y': ('Y position', 'angstrom'),
			'z': ('Z position', 'angstrom'),
			'fractional_x': ('Fractional x position', None),
			'fractional_y': ('Fractional y position', None),
			'fractional_z': ('Fractional z position', None),
		})

		self.input = utils.Map()
		self.output = pd.DataFrame()

		self._data_tables = 'points meta ensembles cells energies'.split()

	@staticmethod
	def _transfer_missing_units(origin, destination, columns):
		for column in columns:
			destination._iago_units[column] = origin._iago_units[column]
			destination._iago_comments[column] = origin._iago_comments[column]

	def assign_output_columns(self):
		for table in self._data_tables:
			overlap = (set(self.output.columns) & set(self.__dict__[table].columns)) - set(['run', 'frame'])
			if len(overlap) == 0:
				continue

			newtable = self.__dict__[table].append(self.output[['run', 'frame'] + list(overlap)])
			DB._transfer_missing_units(self.__dict__[table], newtable, ['run', 'frame'] + list(overlap))
			self.__dict__[table] = newtable
			self.output.drop(list(overlap), inplace=True, axis=1)

	@property
	def groups(self):
		""" Known static atom groups.

		:return: Dictionary with group names as key, list of 0-based atom indices as value.
		"""
		return self._groups

	def write(self, fh):
		""" Writes the database to disk or stream.

		:param fh: File handle or filename.
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
		# points
		data['points'] = self.points.to_dict('list')
		data['points-meta'] = self.points.annotations_to_dict()
		# input / output
		data['input'] = self.input
		data['output'] = self.output.to_dict('list')
		# data['output-meta'] = self.output.annotations_to_dict()
		# data tables
		for table in self._data_tables:
			data[table] = self.__dict__[table].to_dict('list')
			data[table + '-meta'] = self.__dict__[table].annotations_to_dict()

		# finalise
		fh.write(json.dumps(data, separators=(',', ':')))

	def read(self, fh):
		""" Reads the database from disk or stream.

		:param fh: File handle or filename.
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
		# points
		try:
			self.points = utils.annotated_data_frame(data['points-meta'], data['points'])
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
		# data tables
		for table in self._data_tables:
			try:
				self.__dict__[table] = utils.annotated_data_frame(data[table + '-meta'], data[table])
			except KeyError:
				pass
		# cleanup
		fh.close()

	def monitor(self, quantity):
		for table in self._data_tables:
			tobj = self.__dict__[table]
			if quantity in tobj.columns:
				# keep dependency local to this function
				import matplotlib.pyplot as plt
				plt.plot(tobj['frame'], tobj[quantity])
				unit = tobj.explain(quantity).Unit.values[0]
				if unit == 'No unit available.':
					unit = 'n/a'
				label = '%s in %s' % (tobj.explain(quantity).Comment.values[0], unit)
				plt.ylabel(label)
				plt.xlabel('Frame')
				return
		raise ValueError('No such quantity.')