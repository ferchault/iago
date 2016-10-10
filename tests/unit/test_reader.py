# system modules
import os
from unittest import TestCase

# custom module
from iago import Reader

class TestReader(TestCase):
	def test_emtpy_universe(self):
		e = Reader.EmptyUniverse()
		self.assertEqual(len(e), 0)
		self.assertEqual(len(e.trajectory), 0)

	def test_reader(self):
		r = Reader.Reader(os.getcwd())
		self.assertRaises(NotImplementedError, r.get_options)
		self.assertRaises(NotImplementedError, r.read)
		self.assertRaises(NotImplementedError, r.get_input)
		self.assertRaises(NotImplementedError, r.get_output)
		self.assertIsInstance(r.get_universe(), Reader.EmptyUniverse)
		self.assertEqual(r.get_trajectory_frames(), [])

	def test_reader_file_matching(self):
		r = Reader.Reader(os.getcwd())
		self.assertRaises(KeyError, r._get_first_file_matching, ['foo', ])
		file_absolute = os.path.abspath(os.path.realpath(__file__))
		file_relative = os.path.relpath(file_absolute)
		self.assertEqual(r._get_first_file_matching([file_absolute, ]), file_absolute)
		self.assertEqual(r._get_first_file_matching([file_relative, ]), file_absolute)
