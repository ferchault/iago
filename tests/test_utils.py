from unittest import TestCase
import iago.utils as u
import math
import numpy as np

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

	def test_is_plane_selector(self):
		self.assertTrue(u.is_plane_selector('plane foobar'))
		self.assertFalse(u.is_plane_selector('plane foo bar'))
		self.assertFalse(u.is_plane_selector('group plane'))

	def test_plane_point_distance(self):
		nv, pp, dp = (0, 0, 2), (0, 0, 0), (0, 0, 1)
		self.assertEqual(u.plane_point_distance(nv, pp, dp), 1)
		nv, pp, dp = (0, 0, 2), (0, 0, 0), (0, 0, -1)
		self.assertEqual(u.plane_point_distance(nv, pp, dp), -1)
		nv, pp, dp = (7, 11, 13), (2, 3, 5), (17, 19, 23)
		self.assertAlmostEqual(u.plane_point_distance(nv, pp, dp), 515*math.sqrt(339)/339)

	def test_plane_point_distance_multiple(self):
		nv, pp, dp = (0, 0, 2), (0, 0, 0), ((0, 0, 1), (0, 0, -1))
		self.assertTrue(np.allclose(u.plane_point_distance(nv, pp, dp), np.array((1., -1.))))
