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
		self.runmatch = {}

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
		""" Reads groups of atoms from a GROMACS ndx file.

		The 'ndx' file format looks like this:

		.. code-block:: none

			[ group_A ]
			1 2 3
			[ group_B ]
			4 5 6

		with indices being one-based. This group has the same composition for the whole simulation.

		:param filename: Path or filename of the 'ndx' file to load.
		"""
		if not os.path.isabs(filename):
			filename = os.path.abspath(os.path.join(self.path, filename))
		try:
			lines = open(filename).readlines()
		except:
			raise ValueError('Unable to read %s', filename)

		self._db.groups.update(utils.parse_ndx(lines))

	def static_group(self, *args):
		""" Creates one or multiple static groups from an atom selector or indices list.

		If a dictionary is the single argument, its keys are taken as group names while its values are considered to be
		the zero-based atom index list. The following example creates the group 'cas' with atoms 1-3.

		.. code-block:: python
			:emphasize-lines: 3

			import iago.Analyser as Analyser
			c = Analyser()
			c.static_group({'cas': (0, 1, 2)})

		If the group is to be created by a atom selector string, then two arguments are required: a group name first
		followed by the selection string. The following example groups all carbon atoms of atom type CA into one group
		with the name cas:

		.. code-block:: python
			:emphasize-lines: 3

			import iago.Analyser as Analyser
			c = Analyser()
			c.static_group('cas', 'type CA')

		Instead of the selector, a list of zero-based atom indices is accepted, as well.

		.. code-block:: python
			:emphasize-lines: 3

			import iago.Analyser as Analyser
			c = Analyser()
			c.static_group('cas', (1, 2, 3))

		Finally, a list of atoms directly as parameters is accepted:

		.. code-block:: python
			:emphasize-lines: 3

			import iago.Analyser as Analyser
			c = Analyser()
			c.static_group('cas', 1, 2, 3)

		:param args: Dict or name and selector string or name and list of indices.
		"""
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

	def dynamic_plane(self, name, selector, normal=None, framesel=None, comment=None):
		""" Fits a plane through at least three atoms in each frame.

.. warning::

  Be careful **not to use this function with wrapped trajectories**, i.e. coordinates that have been wrapped into the
  simulation box. Imagine points on a plane that is not parallel to any of the simulation box surfaces. Any point
  outside the box will work fine for unwrapped trajectories but will give results other than the expected.

The *name* of the plane is informational only and should be kept concise since it is the identifier of the plane.

The *selector* can be any :ref:`atom selector <selection-atom>` but must select at least three atoms in every frame.
The selector will be evaluated for any frame which means that selected atoms may change over the course
of the trajectory. If you do not want this behaviour, define a static group in your :func:`define_groups`
implementation in your :ref:`analyser script <whatis-analyser>`.

The *normal* defines the surface orientation. Since the plane is stored in the database by recording a support vector
and the surface normal vector, the sign of the normal vector can be important, e.g. if you want to calculate the atom
distance from a certain plane with :func:`dynamic_distance`. The *normal* vector defines the orientation of the normal
vectors calculated by this function such that they point in the same direction. More precisely, the dot product of the
*normal* parameter and the calculated normal vectors is guaranteed to be positive (or zero).

The *framesel* selects the frames that are to be analysed. If follows the regular python slice syntax. Its application
in the context of *iago* is explained in greater depth i

The *comment* can be as verbose as needed to describe the plane in a human-readable format. It will be stored in the
database alongside the results for documentation purposes. Together with the plane *name*, the user should be able to
understand the meaning and origin of this plane.

The resulting data is stored in the database as :attr:`iago.DatabaseProvider.DB.planes`.

Technical detail: the plane is calculated using the SVD approach. There is no mass-weighting of the atoms in question.

:param name: Name of the plane.
:param selector: Which atoms are to be approximated by a plane.
:param normal: Preferred orientation of the normal vector. 3-tuple.
:param framesel: A python slice selecting only certain frames.
:param comment: Description of the plane.
"""
		resultdata = []
		for run, alias in self.parser.get_runs().iteritems():
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
				resultdata.append({
					'run': alias,
					'frame': step,
					'name': name,
					'normal_x': normal_vector[0],
					'normal_y': normal_vector[1],
					'normal_z': normal_vector[2],
					'support_x': center_of_geometry[0],
					'support_y': center_of_geometry[1],
					'support_z': center_of_geometry[2],
				})
		self._db.planes = self._db.planes.append(pd.DataFrame(resultdata), ignore_index=True)

	def dynamic_distance(self, name, selectorA, selectorB, cutoff=None, comment=None, framesel=None, signed=False):
		""" Calculates a distance between two objects in each frame.

.. warning::

  Currently, the minimum image convention is only respected for rectangular unit cells. This is likely to be fixed soon.
  This warning only affects atom-atom distances.

This routine supports atom-atom distances, and atom-plane distances. For atom-atom distances, the cartesian product of
the two selections is built and distances are calculated for each tuple of atoms. For atom-plane calculations, the
shortest distance between each atom and the plane is measured. No wrapping of atoms back into the unit cell is happening
in this case. In the case of atom-plane distances, signed distances are stored in the database. Positive values denote
atom positions in the direction of the plane normal vector (i.e. above the plane).

The *name* of the plane is informational only and should be kept concise since it is the identifier of the distance set.

The two selectors *selectorA* and *selectorB* can either be both :ref:`atom selectors <selection-atom>` or at most one
can be a :ref:`plane selector <selection-plane>`. The selectors are evaluated dynamically, i.e. for each :term:`frame`
separately.

The *cutoff* parameter limits entries to be included in the database. While all combinations from the cartesian product
of the two selections are calculated, only those are stored in the database where the calculated absolute distance is
below the cutoff.

The *comment* can be as verbose as needed to describe the distance set in a human-readable format. It will be stored
in the database alongside the results for documentation purposes. Together with the *name*, the user should be able to
understand the meaning and origin of this distance set.

:param name: Name of the distance set between groups A and B
:param selectorA: Selector for the group A
:param selectorB: Selector for the group B
:param cutoff: Maximum included distance in Angstrom. None to disable.
:param comment: Description of the distance set
"""
		resultdata = []
		for run, alias in self.parser.get_runs().iteritems():
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
								resultdata.append({
									'run': alias,
									'frame': step,
									'name': name,
									'atom1': iidx,
									'atom2': jidx,
									'dist': distances[iidx, jidx],
								})
				df = pd.DataFrame(resultdata)
				self._db.distances = self._db.distances.append(df, ignore_index=True)
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
							resultdata.append({
								'run': alias,
								'frame': step,
								'name': name,
								'plane': planename,
								'atom1': iidx,
								'dist': ds[iidx],
							})
				df = pd.DataFrame(resultdata)
				self._db.planedistances = self._db.planedistances.append(df, ignore_index=True)

	def collect_input_output(self):
		input = utils.Map()
		output = []
		for run, alias in self.parser.get_runs().iteritems():
			try:
				input[alias] = self.parser.get_input(run)
			except NotImplementedError:
				continue
			try:
				output.append(self.parser.get_output(run, alias))
			except NotImplementedError:
				continue

		try:
			output = pd.concat(output, ignore_index=True)
		except ValueError:
			output = pd.DataFrame()
		self._db.input = input
		self._db.output = output

	def run(self):
		self.setup()
		self.parser.run(self.path, self.runmatch)
		self.define_groups()
		self.calculated_columns()
		self.collect_input_output()
		self._db.write(os.path.join(self.path, 'iagodb.json'))

	def _compare_predicate(self, a, p, b):
		return (p == 'eq' and a == b) or (p == 'gt' and a > b) or (p == 'lt' and a < b) or (p == 'ge' and a >= b) or (p == 'le' and a <= b)

	def dynamic_atomlist(self, columnname, filter, where_count=None, comment=None):
		run_col = []
		frame_col = []
		result_col = []
		for run, alias in self.parser.get_runs().iteritems():
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

				run_col.append(alias)
				frame_col.append(frameidx)
				result_col.append([_.index for _ in ag])

		df = pd.DataFrame({'run': run_col, 'frame': frame_col, 'result': result_col})
		self._add_per_frame(df, columnname)
