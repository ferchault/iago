# system modules
import os
from unittest import TestCase

# custom module
from iago import Reader

BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


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

	def test_reader_run_discovery(self):
		# CP2K run
		path = os.path.join(BASE_DIR, 'fixtures', 'debug-88ac57b3e437fa1d5d26d00d6c768324', 'run-1')
		r_cp2k = Reader.CP2KReader(path)
		r_namd = Reader.NAMDReader(path)
		self.assertTrue(r_cp2k.claims())
		self.assertFalse(r_namd.claims())

		# NAMD run
		path = os.path.join(BASE_DIR, 'fixtures', 'ubqTUTORIAL-1aa3d370076a9db6c69c587ba561ecd0', 'run-1')
		r_cp2k = Reader.CP2KReader(path)
		r_namd = Reader.NAMDReader(path)
		self.assertFalse(r_cp2k.claims())
		self.assertTrue(r_namd.claims())
