import StorageProvider
import pandas as pd
import types
import itertools as it
import utils
import DatabaseProvider

class Analyser(object):
	def __init__(self):
		self._db = DatabaseProvider.DB()

	def setup(self):
		""" Points out the required settings and parsers. Has to be user-defined."""
		raise NotImplementedError()

	def calculated_columns(self):
		""" Holds the definitions of calculated columns as derived from the parsed result. May be user-defined."""
		pass

	def define_groups(self):
		""" Holds the definitions of static atom groups. May be user-defined."""
		pass

	def static_load_groups(self, filename):
		try:
			lines = open(filename).readlines()
		except:
			raise ValueError('Unable to read %s', filename)

		self._db.groups.update(utils.parse_ndx(lines))

	def static_group(self, *args):
		if len(args) < 1:
			raise RuntimeError('No empty groups allowed.')

		first, rest = args[:1], args[1:]

		# only dict given as argument
		if len(args) == 1:
			self._db.groups.update(first)

		# selector string or list of elements
		if len(args) == 2:
			if isinstance(args[1], types.StringTypes):
				self._parser.get_atom_indices(args[1])
			else:
				self._db.groups.update({first: args[1]})

		# list of elements
		self._db.groups.update({first: rest})

	def dynamic_plane(self, name, selector, normal=None, comment=None):
		""" Not suitable for wrapped trajectories.
		:param name: Name of the plane.
		:param selector: Which atoms are to be approximated by a plane.
		:param normal: Preferred orientation of the normal vector. 3-tuple.
		:param comment: Description of the plane.
		:return:
		"""
		try:
			maxidx = max(self._db.planes.index)
		except ValueError:
			maxidx = 0
		for run in self.parser.get_runs():
			try:
				u = self.parser.get_universe(run)
			except:
				continue

			steps = self.parser.get_trajectory_frames(run)

			for step, frame in it.izip(steps, u.trajectory):
				ag = u.selectAtoms(selector)
				normal_vector, center_of_geometry = utils.fit_plane(ag.positions, normal=normal)
				self._db.planes.loc[maxidx] = {
					'run': run,
					'frame': step,
					'name': name,
					'normal_x': normal_vector[0],
					'normal_y': normal_vector[1],
					'normal_z': normal_vector[2],
					'support_x': center_of_geometry[0],
					'support_y': center_of_geometry[1],
					'support_z': center_of_geometry[2],
				}
				maxidx += 1





	def run(self):
		self.setup()
		self.parser.run()
		self.define_groups()
		self.calculated_columns()
		return self.storage_provider.save(self._db)

	def _compare_predicate(self, a, p, b):
		return (p == 'eq' and a == b) or (p == 'gt' and a > b) or (p == 'lt' and a < b) or (p == 'ge' and a >= b) or (p == 'le' and a <= b)

	def dynamic_atomlist(self, columnname, filter, where_count=None, comment=None):
		run_col = []
		frame_col = []
		result_col = []
		for run in self.parser.get_runs():
			universe, framenumbers = self.parser.get_coordinates(run)
			if universe is None:
				continue

			for frame in universe.trajectory:
				ag = universe.select_atoms(filter)
				frameidx = next(framenumbers)

				if where_count is not None:
					filtered = []
					condition, mode, value = where_count
					for atom in ag:
						subset = universe.select_atoms(condition % atom.index)
						if self._compare_predicate(len(subset), mode, value):
							filtered.append(atom.index)
					ag = universe.atoms[filtered]

				run_col.append(run)
				frame_col.append(frameidx)
				result_col.append([_.index for _ in ag])

		df = pd.DataFrame({'run': run_col, 'frame': frame_col, 'result': result_col})
		self._add_per_frame(df, columnname)
