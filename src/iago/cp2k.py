# system modules
import os

# third-party modules
import pandas as pd
import numpy as np


class InputFileReader(object):
	def read_run(self, *args, **kwargs):
		i = InputFile()
		i.read_run(*args, **kwargs)
		return i


class LogFileReader2(object):
	def read_run(self, store, run='run-0'):
		lines = open(os.path.join(store.get_basepath(), run, 'run.log')).readlines()
		steps = []
		times = []
		cons = []
		global_eri = []
		for line in lines:
			if "STEP NUMBER" in line:
				steps.append(int(line.strip().split()[-1]))
			if 'TIME [fs]' in line:
				times.append(float(line.strip().split()[-1]))
			if 'CONSERVED QUANTITY' in line:
				cons.append(float(line.strip().split()[-1]))
			if 'Global ERI' in line:
				global_eri.append(float(line.strip().split()[-1]))
		return pd.DataFrame({'stepno': steps, 'time_fs': times, 'cons_h': cons, 'global_eri': global_eri[1:]})


class LogFileReader(object):
	def read_run(self, *args, **kwargs):
		i = LogFile()
		return i.read_run(*args, **kwargs)


class LogFile(object):
	def __init__(self, *args, **kwargs):

		self._skip = []
		self._framedata = None
		self._noatoms = None
		self._pos = 0
		self._keywords = dict((
			('stepnumber', 'STEP NUMBER'),
			('steptime', 'TIME [fs]'),
			('cell_a', 'CELL LNTHS[bohr]'),
			('cell_b', 'CELL LNTHS[bohr]'),
			('cell_c', 'CELL LNTHS[bohr]'),
			('cell_alpha', 'CELL ANGLS[deg]'),
			('cell_beta', 'CELL ANGLS[deg]'),
			('cell_gamma', 'CELL ANGLS[deg]'),
			('temperature', 'TEMPERATURE [K]'),
			('pressure', 'PRESSURE [bar]'),
			('volume', 'VOLUME[bohr^3]'),
			('conserved', 'CONSERVED QUANTITY [hartree]'),
			('coreselfenergy', 'Self energy of the core charge distribution:'),
			('corehamiltonian', 'Core Hamiltonian energy:'),
			('hartree', 'Hartree energy:'),
			('xc', 'Exchange-correlation energy:'),
			('hfx', 'Hartree-Fock Exchange energy:'),
			('dispersion', 'Dispersion energy:'),
			('etot', 'ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):'),
			('epot', 'POTENTIAL ENERGY[hartree]'),
			('ekin', 'KINETIC ENERGY [hartree]'),
			('drift', 'ENERGY DRIFT PER ATOM [K]'),
			('iasd', 'Integrated absolute spin density  :'),
			('s2', 'Ideal and single determinant S**2 :'),
			('scfcycles', 'outer SCF loop converged in'),
			('otcycles', 'outer SCF loop converged in'),
			('globaleri', 'Global ERI counter')
		))

	def read_run(self, store, run='run-0'):
		self._basepath = store.get_basepath()
		self._run = run
		self._loglines = open(os.path.join(store.get_basepath(), run, 'run.log')).readlines()

		self.skip_header()
		results = []
		while self.find_next_frame():
			rowdata = {'run': self._run}
			for keyword in self._keywords:
				retval = self.get_frame_desired(keyword)
				if retval is not None:
					rowdata[keyword] = retval
			rowdata.update(self.get_frame_multi())
			results.append(rowdata)
		return pd.DataFrame(results)

	def skip_header(self):
		""" Set the internal cursor to the beginning of the MD run. """
		if self._pos != 0:
			return
		for idx in range(len(self._loglines)):
			if 'GO CP2K GO!' in self._loglines[idx]:
				self._pos = idx
				return

	def find_next_frame(self):
		""" Sets the next complete frame for analysis.

		:return: Boolean. True if a frame could be loaded.
		"""
		first = False
		for idx in range(self._pos, len(self._loglines)):
			if '*' * 79 in self._loglines[idx]:
				if first:
					self._framedata = self._loglines[self._pos:idx]
					self._pos = idx + 1
					return True
				else:
					first = True
		self._framedata = None
		return False

	def get_frame_multi(self):
		lc, lh, lr, fs = [], [], [], []
		for line in self._framedata:
			line = line.strip()
			if line.startswith('Local Correction array'):
				parts = map(float, line.split()[6:])
				lc.append(parts)
			if line.startswith('Local history array'):
				parts = map(float, line.split()[6:])
				lh.append(parts)
			if line.startswith('a and b'):
				parts = map(float, line.split()[6:])
				lr.append(parts)
			if ' fock_4c ' in line:
				parts = map(float, line.split()[-3:])
				fs.append(parts)

		def _reformat(ln):
			a = np.array(ln)
			numhist = a.shape[-1]
			if numhist == 0:
				return a
			return a.reshape((-1, 3, numhist)).tolist()

		lc, lh, lr = map(_reformat, (lc, lh, lr))

		retvals = {}
		if len(lc) > 0:
			retvals['hfx_corr_a'] = lc
		if len(lh) > 0:
			retvals['hfx_hist_a'] = lh
		if len(lr) > 0:
			retvals['hfx_hist_reg'] = lr
		if len(fs) > 0:
			retvals['fock_forces'] = fs
		return retvals

	def get_frame_desired(self, kind):
		""" Looks up specified properties from the input file tree.

		:param kind: The key to look for. Check the source for known keys.
		:return: Converted value in Angstrom/eV/K/mu_B
		"""

		def keep_float(var):
			try:
				val = float(var)
			except:
				return None
			return val

		keyword = self._keywords[kind]
		converter = dict((
			('cell_a', lambda _: _ * 0.52917721967),
			('cell_b', lambda _: _ * 0.52917721967),
			('cell_c', lambda _: _ * 0.52917721967),
			('volume', lambda _: _ * 0.52917721967 ** 3),
			('conserved', lambda _: _ * 27.21138602),
			('coreselfenergy', lambda _: _ * 27.21138602),
			('corehamiltonian', lambda _: _ * 27.21138602),
			('hartree', lambda _: _ * 27.21138602),
			('xc', lambda _: _ * 27.21138602),
			('hfx', lambda _: _ * 27.21138602),
			('etot', lambda _: _ * 27.21138602),
			('epot', lambda _: _ * 27.21138602),
			('ekin', lambda _: _ * 27.21138602),
			('drift', lambda _: _ * 27.21138602),
		))
		indexs = dict((
			('cell_a', 0),
			('cell_b', 1),
			('cell_c', 2),
			('cell_alpha', 0),
			('cell_beta', 1),
			('cell_gamma', 2),
			('scfcycles', 1),
			('otcycles', 0),
			('temperature', 0),
			('drift', 0),
			('ekin', 0),
			('epot', 0),
			('pressure', 0),
			('s2', 0)
		))
		if kind in converter:
			converter = converter[kind]
		else:
			converter = lambda _: _
		try:
			pos = indexs[kind]
		except:
			pos = None
		if self._framedata is None:
			return None
		for idx in range(len(self._framedata))[::-1]:
			if keyword in self._framedata[idx]:
				parts = self._framedata[idx].split()
				parts = [converter(_) for _ in map(keep_float, parts) if _ is not None]
				if pos is not None:
					parts = [parts[pos]]
				if len(parts) == 1:
					return parts[0]
				return parts
		return None

	def get_frame_mulliken_charges(self):
		if self._framedata is None:
			raise ValueError('No frame loaded.')

		chargelines = []
		begin = None
		for idx in range(len(self._framedata)):
			if 'MULLIKEN POPULATION ANALYSIS' in self._framedata[idx]:
				begin = idx + 3
			if '# Total charge and spin' in self._framedata[idx] and begin is not None:
				chargelines = self._framedata[begin:idx]
		if len(chargelines) == 0:
			return None

		objs = list()
		for line in chargelines:
			parts = line.strip().split()
			objs.append({
				'alpha': float(parts[3]),
				'beta': float(parts[4]),
				'charge': float(parts[5]),
				'spin': float(parts[6]),
			})
		return objs

	def get_frame_dimensions(self):
		if self._framedata is None:
			raise ValueError('No frame loaded.')

		a, b, c = self.get_frame_desired('celllengths')
		alpha, beta, gamma = self.get_frame_desired('cellangles')
		return a, b, c, alpha, beta, gamma
