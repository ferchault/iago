# system modules
from unittest import TestCase

# custom modules
import iago.Reader


class TestCP2KReader(TestCase):
	def test__parse_input_file(self):
		# simple case
		lines = '''
		&MOTION
			&MD
				STEPS 42
			&END MD
		&END'''.split('\n')
		r = iago.Reader.CP2KReader._parse_input_file(lines)
		self.assertEqual(r, {'MOTION': {'MD': {'STEPS': 42.0}}})

		# duplicate entries
		lines = '''
		&FORCE_EVAL
			&DFT
				BASIS_SET_FILE_NAME a
				BASIS_SET_FILE_NAME b
				BASIS_SET_FILE_NAME c
			&END
		&END'''.split('\n')
		r = iago.Reader.CP2KReader._parse_input_file(lines)
		self.assertEqual(r, {'FORCE_EVAL': {'DFT': {'BASIS_SET_FILE_NAME': ['a', 'b', 'c']}}})

		# duplicate tree
		lines = '''
		&FORCE_EVAL
			METHOD QuickstepA
		&END
		&FORCE_EVAL
			METHOD QuickstepB
		&END
		&FORCE_EVAL
			METHOD QuickstepC
		&END'''.split('\n')
		r = iago.Reader.CP2KReader._parse_input_file(lines)
		self.assertEqual(r, {'FORCE_EVAL': [{'METHOD': 'QuickstepA'}, {'METHOD': 'QuickstepB'}, {'METHOD': 'QuickstepC'}]})

		# mismatching end clause
		lines = '''
		&MOTION
			&MD
				STEPS 42
			&END FOOBAR
		&END'''.split('\n')
		self.assertRaises(ValueError, iago.Reader.CP2KReader._parse_input_file, lines)

	def test_recognize_line(self):
		self.assertEqual(iago.Reader.CP2KReader._recognize_line('FOO bar'), ('FOO', 'bar'))
		self.assertEqual(iago.Reader.CP2KReader._recognize_line('FOO 3'), ('FOO', 3.))
		self.assertEqual(iago.Reader.CP2KReader._recognize_line('FOO 3.1'), ('FOO', 3.1))
		self.assertEqual(iago.Reader.CP2KReader._recognize_line('FOO 3e7'), ('FOO', 3e7))
		self.assertEqual(iago.Reader.CP2KReader._recognize_line('FOO 1 2'), ('FOO', (1., 2.)))
		self.assertEqual(iago.Reader.CP2KReader._recognize_line('FOO'), ('FOO', None))
		self.assertEqual(iago.Reader.CP2KReader._recognize_line('FOO bar bar'), ('FOO', ('bar', 'bar')))
