from unittest import TestCase
import iago


class TestAnalyser(TestCase):
	def test_static_group(self):
		a = iago.Analyser()
		a.static_group('debug', 1, 3, 4, 5)
		self.assertDictEqual(a.get_groups(), {'debug': (1, 3, 4, 5)})
