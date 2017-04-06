# standard modules
import os
import re

# custom modules
import Reader


class Parser(object):
	def __init__(self):
		self._readers = dict()
		self.path = None
		self.runmatch = dict()
		self._runcache = None

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
		""" Discovers all available runs in this bucket.

		:return: Dict of run names available in this bucket. Keys: paths, values: names.
		"""
		if self._runcache is not None:
			return self._runcache

		# regular runs
		inodes = os.listdir(self.path)
		directories = [_ for _ in inodes if os.path.isdir(os.path.join(self.path, _))]
		runs = {_: _ for _ in directories if _.startswith('run-')}

		# alternative run directories
		for root, dirs, files in os.walk(self.path):
			relpath = os.path.relpath(root, self.path)
			for regex, replace in self.runmatch.iteritems():
				g = re.match(regex, relpath)
				if g is not None:
					runs[relpath] = replace.format(**g.groupdict())

		self._runcache = runs
		return self._runcache

	def get_universe(self, run):
		return self._readers[run].get_universe()

	def get_input(self, run):
		return self._readers[run].get_input()

	def get_output(self, run, alias):
		o = self._readers[run].get_output()
		o['run'] = alias
		return o

	def get_groups(self, run, groups):
		u = self.get_universe(run)
		if isinstance(u, Reader.EmptyUniverse):
			return {key: [] for (key, value) in groups.iteritems()}
		return {key: u.atoms[value] for (key, value) in groups.iteritems()}

	def get_trajectory_frames(self, run):
		return self._readers[run].get_trajectory_frames()

	def get_run_code(self, runpath, topologyfiles, configfiles, logfiles):
		readers = {'cp2k': Reader.CP2KReader, 'namd': Reader.NAMDReader}
		for label, reader in readers.iteritems():
			r = reader(runpath)
			if 'inputnames' in r.get_options():
				r.inputnames = configfiles + r.inputnames
			if 'topologies' in r.get_options():
				r.topologies = topologyfiles + r.topologies
			if 'logs' in r.get_options():
				r.logs = logfiles + r.logs

			if r.claims():
				return label

	def run(self, path, runmatch=dict(), topologyfiles=[], configfiles=[], logfiles=[]):
		""" Parses all runs of a certain bucket.

		:param path: Basepath of all runs in this bucket.
		:param runmatch: For run autodiscovery: dict of regular expressions matching relative paths from bucket root as
			keys and named group replacements as values.
		"""
		self.path = path
		self.runmatch = runmatch

		for run in self.get_runs():
			code = self.get_run_code(os.path.join(path, run), topologyfiles, configfiles, logfiles)
			if code == 'cp2k':
				self._readers[run] = Reader.CP2KReader(os.path.join(path, run))
			elif code == 'namd':
				self._readers[run] = Reader.NAMDReader(os.path.join(path, run))
			else:
				raise NotImplementedError()
			self._readers[run].read()
