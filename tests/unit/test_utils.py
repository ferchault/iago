# standard modules
from unittest import TestCase
import math

# third-party modules
import numpy as np

# custom modules
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

	def test_map(self):
		m = u.Map()
		m['foo'] = 3
		self.assertEqual(m['foo'], 3)
		self.assertEqual(m.foo, 3)
		m.foo = 4
		self.assertEqual(m['foo'], 4)
		self.assertEqual(m.foo, 4)
		m = u.Map({'foo': 2})
		self.assertEqual(m, {'foo': 2})
		m = u.Map(foo=2)
		self.assertEqual(m, {'foo': 2})
		m = u.Map(foo=2, bar=3)
		del m.bar
		self.assertEqual(m, {'foo': 2})
		m = u.Map(foo=2, bar=3)
		del m['bar']
		self.assertEqual(m, {'foo': 2})
		self.assertRaises(KeyError, m.__getattr__('bar'))

	def test_map_traverse(self):
		m = u.Map(foo=u.Map({'bar': 3}))
		self.assertEqual(m.foo.bar, 3)
		self.assertEqual(m.traverse('foo bar'.split()), 3)

		m = u.Map()
		m.foo = [u.Map(bar=3), u.Map(bar=4)]
		self.assertEqual(m.traverse('foo bar'.split()), 4)

	def test_map_traverse_assignment(self):
		m = u.Map()
		m.foo = [u.Map(bar=3), u.Map(bar=4)]
		q = m.traverse('foo'.split())
		q.bar = 5
		self.assertEqual(m.foo[1].bar, 5)

	def test_map_dict_assignment(self):
		m = u.Map()
		m.foo = 3
		m.bar = {'foo': 4}
		self.assertTrue(isinstance(m.bar, u.Map))

	def test_map_list_dicts(self):
		data = {'foo': [{'a': 1}, {'a': 2}]}
		m = u.Map(data)
		self.assertIsInstance(m.foo[0], u.Map)

		data = {'input': {'run-a': {'a': 1}, 'run-b': {'a': 2}}}
		m = u.Map(data)
		self.assertIsInstance(m.input['run-a'], u.Map)

	def test_path_to_html(self):
		self.assertEqual(u.path_to_html("root['FORCE_EVAL']['SUBSYS']['COORD']['He'][1][0]['foo']"),
						 "FORCE_EVAL.SUBSYS.COORD.He.[1][0].foo")
		self.assertEqual(u.path_to_html("root['FORCE_EVAL']['SUBSYS']['COORD']['He'][1][0]"),
						 "FORCE_EVAL.SUBSYS.COORD.He.[1][0]")
		self.assertEqual(u.path_to_html("root['FORCE_EVAL']['SUBSYS']['COORD']['He']"),
						 "FORCE_EVAL.SUBSYS.COORD.He")

	def test_safe_dict_raises(self):
		s = u.SafeDict()
		s['test'] = 1
		try:
			s['test'] = 2
		except ValueError:
			self.assertTrue(True)
			return
		self.fail('SafeDict allows overwriting entries.')

	def test_safe_dict_update(self):
		s = u.SafeDict()
		s.update(test=2)
		self.assertEqual(s['test'], 2)

	def test_fit_plane(self):
		nv, cog = u.fit_plane(np.array([[-1, 0, 0], [0, -1, 0], [1, 0, 0], [0, 1, 0]]), normal=(0, 0, 1))
		self.assertTrue(np.allclose(nv, (0, 0, 1)))
		self.assertTrue(np.allclose(cog, (0, 0, 0)))

	def test_fit_plane_undef(self):
		positions_in_line = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, 3]])
		positions_too_few = np.array([[0, 0, 0], [0, 0, 1]])
		self.assertRaises(ValueError, u.fit_plane, positions_in_line)
		self.assertRaises(ValueError, u.fit_plane, positions_too_few)

	def test_minimum_image_distance(self):
		d = 2.446792
