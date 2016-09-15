from unittest import TestCase
import iago.utils as u


class TestUtils(TestCase):
	def test_parse_ndx_basic_load(self):
		lines = ['[demo]', '1 2 3']
		groups = u.parse_ndx(lines)

		self.assertDictEqual(groups, {'demo': [0, 1, 2]})

	def test_parse_ndx_skip_comments(self):
		lines = ['#test', '; test', '[demo]', '1 2 3']
		groups = u.parse_ndx(lines)

		self.assertDictEqual(groups, {'demo': [0, 1, 2]})

	def test_parse_ndx_skip_initial_whitespace(self):
		lines = ['  [demo]', '  1 2 3']
		groups = u.parse_ndx(lines)

		self.assertDictEqual(groups, {'demo': [0, 1, 2]})

	def test_parse_ndx_skip_name_whitespace(self):
		lines = ['[ demo	]', '1 2 3']
		groups = u.parse_ndx(lines)

		self.assertDictEqual(groups, {'demo': [0, 1, 2]})
