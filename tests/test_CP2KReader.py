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
			&END
		&END'''.split('\n')
		r = c._parse_input_file(lines)
		self.assertEqual(r, {'FORCE_EVAL': {'DFT': {'BASIS_SET_FILE_NAME': ['a', 'b']}}})

		# duplicate tree
		lines = '''
		&FORCE_EVAL
			METHOD QuickstepA
		&END
		&FORCE_EVAL
			METHOD QuickstepB
		&END'''.split('\n')
		r = c._parse_input_file(lines)
		self.assertEqual(r, {'FORCE_EVAL': [{'METHOD': 'QuickstepA'}, {'METHOD': 'QuickstepB'}]})
