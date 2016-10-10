# standard modules
import re
import types

# third-party modules
import numpy as np
import pandas as pd
import pint


class Map(dict):
	def __init__(self, *args, **kwargs):
		super(Map, self).__init__(*args, **kwargs)
		for arg in args:
			if isinstance(arg, dict):
				for k, v in arg.iteritems():
					if isinstance(v, (list, tuple)):
						v = [Map(_) if isinstance(_, dict) else _ for _ in v]
					if isinstance(v, dict):
						v = Map(v)
					self[k] = v

		if kwargs:
			for k, v in kwargs.iteritems():
				self[k] = v

	def __getattr__(self, attr):
		return self.get(attr)

	def __setattr__(self, key, value):
		self.__setitem__(key, value)

	def __setitem__(self, key, value):
		super(Map, self).__setitem__(key, value)
		if type(value) is dict:
			value = Map(value)
		self.__dict__.update({key: value})

	def __delattr__(self, item):
		self.__delitem__(item)

	def __delitem__(self, key):
		super(Map, self).__delitem__(key)
		del self.__dict__[key]

	def traverse(self, keys):
		r = self
		for k in keys:
			r = r[k]
			if isinstance(r, list):
				r = r[-1]
		return r


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
				new = map(lambda _: int(_) - 1, line.split())
				if len([_ for _ in new if _ < 0]) > 0:
					raise ValueError('Only non-negative and one-based indices supported.')
				thisvalues += new
			except:
				raise ValueError('Non-numeric atom numbers.')
	if lastgroup is not None:
		groups[lastgroup] = thisvalues

	return groups


def fit_plane(positions, normal=None):
	""" Calculates a fitted plane given atom positions.

	.. warning::
	 	Not suited for wrapped trajectories. Imagine points on a plane that is not parallel to any of the simulation
	 	box surfaces. Any point outside the box will work fine for unwrapped trajectories but will give results other
	 	than the expected.

	:param positions: :class:`numpy.array` with atom positions.
	:param normal: Tuple of the preferred normal vector orientation
	:return: Tuple of normal vector and center of geometry as :class:`numpy.array`
	"""
	if np.linalg.matrix_rank(positions) <= 1:
		raise ValueError('Collinear positions.')

	rows, cols = positions.shape
	if rows < 3:
		raise ValueError('Plane undetermined - need at least three points.')

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


def is_plane_selector(selector):
	return re.match('^plane [^ ]+$', selector)


def plane_point_distance(normal_vector, plane_point, data_point):
	nv = np.array(normal_vector)
	pp = np.array(plane_point)
	d = -np.dot(nv, pp)
	signed_distance = (np.dot(normal_vector, np.array(data_point).T) + d)/np.linalg.norm(nv)
	return signed_distance


class SafeDict(dict):
	def __init__(self, other=None, **kwargs):
		super(SafeDict, self).__init__()
		self.update(other, **kwargs)

	def __setitem__(self, key, value):
		if key in self:
			raise ValueError('Overwriting an existing entry (%s) is prohibited.' % key)
		super(SafeDict, self).__setitem__(key, value)

	def update(self, other=None, **kwargs):
		if other is not None:
			for k, v in other.iteritems():
				self[k] = v
		for k, v in kwargs.items():
			self[k] = v


def annotated_data_frame(definition, datadict=None):
	columns = definition.keys()
	if datadict is None:
		df = pd.DataFrame(columns=columns)
	else:
		df = pd.DataFrame.from_dict(datadict)
	df._iago_comments = {k: definition[k][0] for k in definition.iterkeys()}
	df._iago_units = {k: definition[k][1] for k in definition.iterkeys()}
	return df


def compare(reference, other=None, labels=None):
	""" Compares two or more dictionaries in the ipython / jupyter context."""
	import deepdiff
	from IPython.core.display import display, HTML

	if not isinstance(other, (list, tuple)):
		if other is None:
			other = reference[1:]
			reference = reference[0]
		else:
			other = [other]

	if not isinstance(labels, (list, tuple)):
		if labels is None:
			labels = map(str, range(len(other) + 1))
		else:
			labels = [labels(reference)] + map(labels, other)

	diffs = [deepdiff.DeepDiff(reference, _) for _ in other]

	output = '''<style>
        table.iago-diff > tr {
            border: 0px
        }
        td.iago-diff-changed {
            background-color: #f8cbcb;
        }
        td.iago-diff-unchanged {
            background-color: #a6f3a6;
            font-size: 60%;
        }
        td.iago-diff-reference {
            background-color: lightgrey;
        }
        td.iago-diff-key {
            text-align: center;
        }
    </style>'''
	output += '<table class="iago-diff"><tr>'
	output += ''.join(map(lambda _: '<th>%s</th>' % _, labels))
	output += '</tr>'

	headline = '<tr><td colspan="%d" class="iago-diff-key">%%s</td></tr>' % len(labels)

	# changed values
	all_changed_values = []
	for diff in diffs:
		try:
			all_changed_values += diff['values_changed'].keys()
		except KeyError:
			pass
		pass

	for value in set(all_changed_values):
		output += headline % path_to_html(value)

		for diff in diffs:
			try:
				output += '<td class="iago-diff-reference">%s</td>' % diff['values_changed'][value]['old_value']
			except KeyError:
				pass
			break

		for didx, diff in enumerate(diffs):
			try:
				output += '<td class="iago-diff-changed">%s</td>' % diff['values_changed'][value]['new_value']
			except KeyError:
				output += '<td class="iago-diff-unchanged">(unchanged)</td>'

	output += headline % 'other differences'
	output += '<tr><td colspan="%d" >' % len(labels)
	for didx, diff in enumerate(diffs):
		missing = []
		for k in diff.keys():
			if k != 'values_changed':
				missing.append(k)
		if len(missing) == 0:
			continue
		output += '<b>%s</b><br/>' % labels[didx + 1]
		for m in missing:
			output += '<i>%s</i>: %s </br>' % (m, diff[m])

	output += '</td></tr></table>'
	return display(HTML(output))


def path_to_html(path):
	if path.startswith('root'):
		path = path[4:]

	output_parts = []
	parts = path.split('][')
	for idx, part in enumerate(parts):
		if (idx < len(parts) - 1 and part[-1] == "'") or (idx == len(parts) - 1 and part[-2] == "'"):
			if idx == 0:
				start = 2
			else:
				start = 1
			if idx == len(parts) - 1:
				end = -2
			else:
				end = -1
			output_parts.append(part[start:end])
		else:
			if idx == len(parts) - 1:
				output_parts.append('[' + part)
			else:
				output_parts.append('[%s]' % part)

	output = ''
	for part in output_parts:
		if len(output) > 0 and output[-1] != ']':
			output += '.' + part
		else:
			if part.startswith('[') or len(output) == 0:
				output += part
			else:
				output += '.' + part
	return output
