# system modules
from unittest import TestCase

# third-party modules
import pandas as pd

# custom modules
import iago


class TestPandasPatch(TestCase):
	def test_explain(self):
		df = pd.DataFrame(columns='a b'.split())
		self.assertEqual(df.explain('a').to_dict(), {
			'Comment': {0: 'No description available.'},
			'Name': {0: 'a'},
			'Unit': {0: 'No unit available.'}})
		self.assertEqual(df.explain('c').to_dict(), {
			'Comment': {0: 'No description available.'},
			'Name': {0: 'c'},
			'Unit': {0: 'No unit available.'}})
		self.assertEqual(df.explain('a b'.split()).to_dict(), {
			'Comment': {0: 'No description available.', 1: 'No description available.'},
			'Name': {0: 'a', 1: 'b'},
			'Unit': {0: 'No unit available.', 1: 'No unit available.'}})
		self.assertEqual(df.explain().to_dict(), {
			'Comment': {0: 'No description available.', 1: 'No description available.'},
			'Name': {0: 'a', 1: 'b'},
			'Unit': {0: 'No unit available.', 1: 'No unit available.'}})
		df._iago_units['a'] = 'bar'
		df._iago_comments['a'] = 'foo'
		self.assertEqual(df.explain('a').to_dict(), {
			'Comment': {0: 'foo'},
			'Name': {0: 'a'},
			'Unit': {0: 'bar'}})
		df._iago_units['a'] = None
		self.assertEqual(df.explain('a').to_dict(), {
			'Comment': {0: 'foo'},
			'Name': {0: 'a'},
			'Unit': {0: 'Dimensionless.'}})

	def test_annotations_to_dict(self):
		df = pd.DataFrame(columns='a b'.split())
		df._iago_units['a'] = 'bar'
		df._iago_comments['a'] = 'foo'
		self.assertEqual(df.annotations_to_dict(), {'a': ('foo', 'bar')})

	def test_instance_attributes(self):
		df1 = pd.DataFrame(columns='a b'.split())
		df1._iago_units['a'] = 'test'
		df2 = pd.DataFrame(columns='a b'.split())
		self.assertTrue(df1.explain().to_dict() != df2.explain().to_dict())
