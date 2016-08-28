import pandas as pd


class Map(dict):
	def __init__(self, *args, **kwargs):
		super(Map, self).__init__(*args, **kwargs)
		for arg in args:
			if isinstance(arg, dict):
				for k, v in arg.iteritems():
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
		self.__dict__.update({key: value})

	def __delattr__(self, item):
		self.__delitem__(item)

	def __delitem__(self, key):
		super(Map, self).__delitem__(key)
		del self.__dict__[key]

	def traverse(self, keys):
		return reduce(lambda d, k: d[k], keys, self)


class DatabaseProvider(object):
	@staticmethod
	def get_annotated_collections():
		return 'atommeta', 'atomframemeta', 'framemeta'


class MemoryDatabaseProvider(DatabaseProvider):
	def __init__(self, run, atom, atommeta, atomframe, atomframemeta, frame, framemeta):
		self.run = Map(run)
		self.atom = pd.DataFrame(atom)
		self.atommeta = Map(atommeta)
		self.atomframe = pd.DataFrame(atomframe)
		self.atomframemeta = Map(atomframemeta)
		self.frame = pd.DataFrame(frame)
		self.framemeta = Map(framemeta)
