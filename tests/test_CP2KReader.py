# system modules
from unittest import TestCase

# custom modules
import iago.Reader


class TestCP2KReader(TestCase):
	def test__parse_input_file(self):
		c = iago.Reader.CP2KReader(None)

		# simple case
		lines = '''
		&MOTION
			&MD
				STEPS 42
			&END MD
		&END'''.split('\n')
		r = c._parse_input_file(lines)
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
		r = c._parse_input_file(lines)
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
		r = c._parse_input_file(lines)
		self.assertEqual(r, {'FORCE_EVAL': [{'METHOD': 'QuickstepA'}, {'METHOD': 'QuickstepB'}, {'METHOD': 'QuickstepC'}]})

		# mismatching end clause
		lines = '''
		&MOTION
			&MD
				STEPS 42
			&END FOOBAR
		&END'''.split('\n')
		self.assertRaises(ValueError, c._parse_input_file, lines)
