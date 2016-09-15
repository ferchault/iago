import StorageProvider
import pandas as pd


class Analyser(object):
	def __init__(self):
		self.storage_provider = StorageProvider.JSONProvider()

	def setup(self):
		""" Points out the required settings and parsers. Has to be user-defined."""
		raise NotImplementedError()

	def calculated_columns(self):
		""" Holds the definitions of calculated columns as derived from the parsed result. Has to be user-defined."""
		raise NotImplementedError()

	def run(self):
		self.setup()
		self.parser.run()
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
