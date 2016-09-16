class Reader(object):
	def __init__(self, path):
		self._path = path
		self._options = dict()

	def set_options(self, **kwargs):
		self._options.update(kwargs)

	def read(self):
		raise NotImplementedError()


class CP2KReader(Reader):
	def read(self):
		# TODO make configurable
		filename = 'run.inp'

