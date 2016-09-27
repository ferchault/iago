# system modules
from unittest import TestCase

# custom modules
import iago


class TestAnalyser(TestCase):
	def test_static_group_args(self):
		a = iago.Analyser()
		a.static_group('debug', 1, 3, 4, 5)
		self.assertDictEqual(a.get_groups(), {'debug': [1, 3, 4, 5]})

	def test_static_group_listarg(self):
		a = iago.Analyser()
		a.static_group('debug', [1, 3, 4, 5])
		self.assertDictEqual(a.get_groups(), {'debug':  [1, 3, 4, 5]})