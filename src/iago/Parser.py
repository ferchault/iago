import Reader
import os

class Parser(object):
	def __init__(self):
		self._readers = dict()

	def get_atom_indices(self, selector):
		"""
		:param selector: Valid selection string.
		:return: List of 0-based atom indices.
		"""
		raise NotImplementedError()

	def get_runs(self):
		raise NotImplementedError()

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
		for run in self.get_runs():
			code = self.get_run_code(run)
			if code == 'cp2k':
				self._readers[run] = Reader.CP2KReader(os.path.join(path, run))
			else:
				raise NotImplementedError()
			self._readers[run].read()
