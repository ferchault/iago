import Reader
import os

class Parser(object):
	def __init__(self):
		self._readers = dict()
		self.path = None

	def get_atom_indices(self, selector):
		"""
		:param selector: Valid selection string.
		:return: List of 0-based atom indices.
		"""
		raise NotImplementedError()

	def get_runs(self):
		"""
		:return: List of run names available in this bucket.
		"""
		inodes = os.listdir(self.path)
		directories = [_ for _ in inodes if os.path.isdir(os.path.join(self.path, _))]
		runs = [_ for _ in directories if _.startswith('run-')]
		return runs

	def get_universe(self, run):
		raise NotImplementedError()

	def get_trajectory_frames(self, run):
		raise NotImplementedError()

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
