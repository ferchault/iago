# standard modules
import os

# custom modules
import Reader


class Parser(object):
	def __init__(self):
		self._readers = dict()
		self.path = None

	def get_atom_indices(self, selector):
		"""
		:param selector: Valid selection string.
		:return: List of 0-based atom indices.
		"""
		u = self._readers.itervalues().next().get_universe()
		if isinstance(u, Reader.EmptyUniverse):
			return []

		ag = u.select_atoms(selector)
		return [_.index for _ in ag]

	def get_runs(self):
		"""
		:return: List of run names available in this bucket.
		"""
		inodes = os.listdir(self.path)
		directories = [_ for _ in inodes if os.path.isdir(os.path.join(self.path, _))]
		runs = [_ for _ in directories if _.startswith('run-')]
		return runs

	def get_universe(self, run):
		return self._readers[run].get_universe()

	def get_groups(self, run, groups):
		u = self.get_universe(run)
		if isinstance(u, Reader.EmptyUniverse):
			return {key: [] for (key, value) in groups.iteritems()}
		return {key: u.atoms[value] for (key, value) in groups.iteritems()}

	def get_trajectory_frames(self, run):
		return self._readers[run].get_trajectory_frames()

	def get_run_code(self, run):
		return 'cp2k'

	def run(self, path):
		""" Parses all runs of a certain bucket.
		:return:
		"""
		self.path = path

		for run in self.get_runs():
			code = self.get_run_code(run)
			if code == 'cp2k':
				self._readers[run] = Reader.CP2KReader(os.path.join(path, run))
			else:
				raise NotImplementedError()
			self._readers[run].read()
