# system modules
import json
import os
import pickle

# third-party modules
import pint

# custom modules
import DatabaseProvider as dp


class StorageProvider(object):
	def load(self, bucket_id, config, stream=None):
		raise NotImplementedError()

	def save(self, db, stream=None):
		raise NotImplementedError()


class JSONProvider(StorageProvider):
	@staticmethod
	def _get_default_filename(self, config):
		return os.path.join(config.basepath, 'db.json')

	def load(self, bucket_id, config, stream=None):
		if stream is None:
			fh = open(self._get_default_filename())
		else:
			fh = stream

		data = json.load(fh)
		ureg = pint.UnitRegistry()
		annotated = dict()
		for collection in dp.DatabaseProvider.get_annotated_collections():
			annotated[collection] = {k: ureg.Quantity.from_tuple(pickle.loads(v)) for k, v in data[collection].items()}

		db = dp.MemoryDatabaseProvider(data['run'], data['atom'], annotated['atommeta'], data['atomframe'],
										annotated['atomframemeta'], data['frame'], annotated['framemeta'])
		return db

	def save(self, db, stream=None):
		fn = None
		if stream is None:
			fh = open(self._get_default_filename(), 'w')
			fn = self._get_default_filename()
		else:
			fn = stream.filename
			fh = stream

		data = {'run': db.run, 'atom': db.atom, 'atomframe': db.atomframe, 'frame': db.frame}
		for collection in dp.DatabaseProvider.get_annotated_collections():
			data[collection] = {k: pickle.dumps(v.to_tuple()) for k, v in getattr(db, collection)}

		return fn
