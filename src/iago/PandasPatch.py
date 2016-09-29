""" The purpose of this file is to monkey-patch the pandas.DataFrame class to include metadata about columns like
units.
"""

# system modules
import types

# third-party modules
import pandas as pd


def explain(self, columnnames=None):
	if isinstance(columnnames, types.StringTypes):
		columnnames = [columnnames]

	if columnnames is None:
		columnnames = self.columns

	comments, units = [], []
	for columnname in columnnames:
		try:
			comment = self._iago_comments[columnname]
		except KeyError:
			comment = 'No description available.'

		try:
			unit = self._iago_units[columnname]
			if unit is None:
				unit = 'Dimensionless.'
		except KeyError:
			unit = 'No unit available.'

		comments.append(comment)
		units.append(unit)

	return pd.DataFrame({'Name': columnnames, 'Comment': comments, 'Unit': units})


def modified_init(self, *args, **kwargs):
	self._iago_old_init(*args, **kwargs)
	self._iago_comments = dict()
	self._iago_units = dict()

def annotations_to_dict(self):
	columns = self._iago_units.keys()
	return {column: (self._iago_comments[column], str(self._iago_units[column])) for column in columns}


pd.DataFrame._iago_old_init = pd.DataFrame.__init__
pd.DataFrame.__init__ = modified_init
pd.DataFrame.explain = explain
pd.DataFrame.annotations_to_dict = annotations_to_dict
