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

	def test_parse_ndx_raise_empty_group(self):
		lines = ['[ ]', '1 2 3']

		self.assertRaises(ValueError, u.parse_ndx, lines)

	def test_parse_ndx_raise_missing_group(self):
		lines = ['1 2 3']

		self.assertRaises(ValueError, u.parse_ndx, lines)

	def test_parse_ndx_raise_zero_based(self):
		lines = ['[ demo	]', '0 2 3']

		self.assertRaises(ValueError, u.parse_ndx, lines)

	def test_parse_ndx_raise_noninteger(self):
		lines = ['[ demo	]', 'a b c']

		self.assertRaises(ValueError, u.parse_ndx, lines)
