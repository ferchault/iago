import re
import types
import pandas as pd
import numpy as np


def parse_ndx(lines):
	""" Parses a GROMACS ndx file.
	:param lines: Iterable of strings.
	:return: Dictionary with group names as key, list of 0-based atom indices as value.
	"""
	groups = dict()
	lastgroup = None
	thisvalues = []

	for line in lines:
		line = line.strip()

		# skip comments and empty lines
		if line.startswith('#') or line.startswith(';') or len(line) == 0:
			continue

		# group definition
		if line.startswith('['):
			if lastgroup is not None:
				groups[lastgroup] = thisvalues
				thisvalues = []
			lastgroup = re.match(r'\[(.*)\]', line).groups()[0].strip()
			if len(lastgroup) == 0:
				raise ValueError('Empty group names not supported.')
		else:
			if lastgroup is None:
				raise ValueError('First non-comment line has to contain a group definition.')
			try:
				new = map(lambda _: int(_)-1, line.split())
				if len([_ for _ in new if _ < 0]) > 0:
					raise ValueError('Only non-negative and one-based indices supported.')
				thisvalues += new
			except:
				raise ValueError('Non-numeric atom numbers.')
	if lastgroup is not None:
		groups[lastgroup] = thisvalues

	return groups


def fit_plane(positions, normal=None):
	""" Not suited for wrapped trajectories.
	:param positions:
	:param normal:
	:return:
	"""
	rows, cols = positions.shape
	p = np.ones((rows, 1))
	ab = np.hstack([positions, p])
	u, d, v = np.linalg.svd(ab, True)
	b = v[3, :]
	nn = np.linalg.norm(b[0:3])
	b /= nn

	normal_vector = b[0:3]
	center_of_geometry = np.average(positions, axis=0)

	if normal is not None:
		normal = np.array(normal)
		if np.dot(normal, normal_vector) < 0:
			normal_vector *= -1

	return normal_vector, center_of_geometry


class SafeDict(dict):
	def __init__(self, other=None, **kwargs):
		super(SafeDict, self).__init__()
		self.update(other, **kwargs)

	def __setitem__(self, key, value):
		if key in self:
			raise ValueError('Overwriting an existing entry is prohibited.')
		super(SafeDict, self).__setitem__(key, value)

	def update(self, other=None, **kwargs):
		if other is not None:
			for k, v in other:
				self[k] = v
		for k, v in kwargs.items():
			self[k] = v


class AnnotatedDataFrame(pd.DataFrame):
	def __init__(self, labels, data=None):
		if data is None:
			super(AnnotatedDataFrame, self).__init__(columns=labels.keys())
		else:
			super(AnnotatedDataFrame, self).__init__(data)

		self._comments = dict()
		self._units = dict()
		for k, v in labels.iteritems():
			comment, unit = v
			self._comments[k] = comment
			self._units[k] = unit

	def explain(self, columnnames):
		if isinstance(columnnames, types.StringTypes):
			columnnames = [columnnames]

		comments, units = [], []
		for columnname in columnnames:
			try:
				comment = self._comments[columnname]
			except KeyError:
				comment = 'No description available.'

			try:
				unit = self._units[columnname]
				if unit is None:
					unit = 'Dimensionless.'
			except KeyError:
				unit = 'No unit available.'

			comments.append(comment)
			units.append(units)

		return pd.DataFrame({'Name': columnnames, 'Comment': comments, 'Unit': units})