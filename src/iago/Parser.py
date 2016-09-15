class Parser(object):
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