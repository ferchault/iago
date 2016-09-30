# third-party modules
import pandas as pd


class LogFile(object):
	def _skip_header(self):
		""" Set the internal cursor to the beginning of the MD run. """
		for line in self._fh:
			if 'GO CP2K GO!' in line:
				return

	def __init__(self, fh):
		if not hasattr(fh, 'read'):
			fh = open(fh)

		self._fh = fh
		self._skip_header()
		self._results = []
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
		self._converters = dict((
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
		self._parse_positions = dict((
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

	def _extract_value(self, keyword, line):
		def keep_float(var):
			try:
				val = float(var)
			except:
				return None
			return val

		keyword = self._keywords[keyword]

		try:
			converter = self._converters[keyword]
		except KeyError:
			def converter(x): x
		try:
			pos = self._parse_positions[keyword]
		except KeyError:
			pos = None

		parts = line.split()
		parts = [converter(_) for _ in map(keep_float, parts) if _ is not None]
		if pos is not None:
			parts = [parts[pos]]
		if len(parts) == 1:
			return parts[0]
		return parts

	def parse(self):
		parselines = []
		sep = '*' * 79
		for line in self._fh:
			if sep in line:
				self._readframe(parselines)
				parselines = []
			else:
				parselines.append(line)
		self._readframe(parselines)

		return pd.DataFrame(self._results)

	def _readframe(self, parselines):
		rowdata = {}
		for line in parselines:
			for keyword in self._keywords:
				if keyword in line:
					retval = self._extract_value(keyword, line)
					if retval is not None:
						rowdata[keyword] = retval
					break
		self._results.append(rowdata)