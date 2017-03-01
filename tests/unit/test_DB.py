# standard modules
import StringIO
from unittest import TestCase

# custom modules
from iago.DatabaseProvider import DB


class TestDB(TestCase):
	def test_read_empty(self):
		s = StringIO.StringIO('{}')
		d = DB()
		try:
			d.read(s, format='json')
		except KeyError:
			self.fail('DB cannot handle empty JSON files.')
