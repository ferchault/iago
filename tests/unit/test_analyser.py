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

	def test_empty_input_output(self):
		# ensure empty buckets do not break iago
		def _mock_get_runs():
			return {}
		a = iago.Analyser()
		a.parser.get_runs = _mock_get_runs
		a.collect_input_output()

	def test_mandatory_and_optional_overrides(self):
		a = iago.Analyser()

		self.assertRaises(NotImplementedError, a.setup)
		try:
			a.define_groups()
			a.calculated_columns()
		except:
			self.fail('Mistakenly mandatory routines.')

	def test_fail_missing_ndx(self):
		a = iago.Analyser()
		self.assertRaises(ValueError, a.static_load_groups, 'nonexisting.ndx')
		a.path = '/'
		self.assertRaises(ValueError, a.static_load_groups, 'nonexisting.ndx')

	def test_no_empty_groups(self):
		a = iago.Analyser()
		self.assertRaises(ValueError, a.static_group, 'test')
		self.assertRaises(ValueError, a.static_group)

	def test__compare_predicate(self):
		func = iago.Analyser._compare_predicate

		self.assertEqual(func(2, 'eq', 2), True)
		self.assertEqual(func(2, 'eq', 3), False)
		self.assertEqual(func(2, 'gt', 2), False)
		self.assertEqual(func(2, 'gt', 1), True)
		self.assertEqual(func(2, 'lt', 2), False)
		self.assertEqual(func(2, 'lt', 3), True)
		self.assertEqual(func(2, 'ge', 2), True)
		self.assertEqual(func(2, 'ge', 1), True)
		self.assertEqual(func(2, 'le', 2), True)
		self.assertEqual(func(2, 'le', 3), True)
