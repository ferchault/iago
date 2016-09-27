# standard modules
import itertools as it
import os
import types
import utils

# third-party modules
import MDAnalysis.analysis.distances as mdaad
import pandas as pd

# custom modules
import DatabaseProvider
import Parser


class Analyser(object):
	def __init__(self):
		self._db = DatabaseProvider.DB()
		self.parser = Parser.Parser()
		self.path = None

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
		if not os.path.isabs(filename):
			filename = os.path.abspath(os.path.join(self.path, filename))
		try:
			lines = open(filename).readlines()
		except:
			raise ValueError('Unable to read %s', filename)

		self._db.groups.update(utils.parse_ndx(lines))

	def static_group(self, *args):
		if len(args) < 1:
			raise RuntimeError('No empty groups allowed.')

		first, rest = args[0], args[1:]

		# only dict given as argument
		if len(args) == 1:
			self._db.groups.update(first)

		# selector string or list of elements
		if len(args) == 2:
			if isinstance(args[1], types.StringTypes):
				idx = self.parser.get_atom_indices(args[1])
			else:
				idx = args[1:][0]
			self._db.groups.update({first: idx})

		# list of elements
		if len(args) > 2:
			self._db.groups.update({first: list(rest)})

	def get_groups(self):
		return self._db.groups

	def dynamic_plane(self, name, selector, normal=None, comment=None, framesel=None):
		""" Not suitable for wrapped trajectories.
		:param name: Name of the plane.
		:param selector: Which atoms are to be approximated by a plane.
		:param normal: Preferred orientation of the normal vector. 3-tuple.
		:param comment: Description of the plane.
		:return:
		"""
		try:
			maxidx = max(self._db.planes.index) + 1
		except ValueError:
			maxidx = 0
		for run in self.parser.get_runs():
			try:
				u = self.parser.get_universe(run)
			except:
				continue

			steps = self.parser.get_trajectory_frames(run)
			tgroups = self.parser.get_groups(run, self.get_groups())

			if framesel is None:
				framesel = slice(None)

			for step, frame in it.izip(steps[framesel], u.trajectory[framesel]):
				ag = u.select_atoms(selector, **tgroups)
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

	def dynamic_distance(self, name, selectorA, selectorB, cutoff=None, comment=None, framesel=None, signed=False):
		"""
		:param name: Name of the distance set between groups A and B
		:param selectorA: Selector for the group A
		:param selectorB: Selector for the grouop B
		:param cutoff: Maximum included distance in Angstrom. None to disable.
		:param comment: Description of the distance set
		:return:
		"""
		try:
			maxidx = max(self._db.distances.index) + 1
		except ValueError:
			maxidx = 0
		for run in self.parser.get_runs():
			try:
				u = self.parser.get_universe(run)
			except:
				continue

			steps = self.parser.get_trajectory_frames(run)
			tgroups = self.parser.get_groups(run, self.get_groups())

			if framesel is None:
				framesel = slice(None)

			if not (utils.is_plane_selector(selectorA) or utils.is_plane_selector(selectorB)):
				# simple atom-atom distance
				for step, frame in it.izip(steps[framesel], u.trajectory[framesel]):
					agA = u.select_atoms(selectorA, **tgroups)
					agB = u.select_atoms(selectorB, **tgroups)
					distances = mdaad.distance_array(agA.positions, agB.positions, box=frame.dimensions)
					for iidx, i in enumerate([_.id for _ in agA]):
						for jidx, j in enumerate([_.id for _ in agB]):
							if cutoff is None or distances[iidx, jidx] < cutoff:
								self._db.distances.loc[maxidx] = {
									'run': run,
									'frame': step,
									'name': name,
									'atom1': iidx,
									'atom2': jidx,
									'dist': distances[iidx, jidx],
								}
								maxidx += 1
			elif utils.is_plane_selector(selectorA) and utils.is_plane_selector(selectorB):
				raise ValueError('Distance between planes not defined.')
			else:
				# atom-plane distances
				if utils.is_plane_selector(selectorB):
					planeselector, atomselector = selectorB, selectorA
				else:
					planeselector, atomselector = selectorA, selectorB

				planename = planeselector.split()[1]
				df = self._db.planes[(self._db.planes.run == run) & (self._db.planes.name == planename)]
				for step, frame in it.izip(steps[framesel], u.trajectory[framesel]):
					ag = u.select_atoms(atomselector, **tgroups)

					# fetch plane
					try:
						s = df[df.frame == step].iloc[0]
					except IndexError:
						# no entry, no data
						continue
					ds = utils.plane_point_distance([s.normal_x, s.normal_y, s.normal_z], (s.support_x, s.support_y, s.support_z), ag.positions)
					if signed is False:
						ds = abs(ds)
					for iidx, i in enumerate([_.id for _ in ag]):
						if cutoff is None or ds[iidx] < cutoff:
							self._db.distances.loc[maxidx] = {
								'run': run,
								'frame': step,
								'name': name,
								'plane': planename,
								'atom1': iidx,
								'dist': ds[iidx],
							}
							maxidx += 1

	def run(self):
		self.setup()
		self.parser.run(self.path)
		self.define_groups()
		self.calculated_columns()
		self._db.write(os.path.join(self.path, 'iagodb.json'))

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
