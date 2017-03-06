from unittest import TestCase

# custom modules
import iago.Reader

class TestNAMDReader(TestCase):
    def test__parse_input_file_last_one_wins(self):
        lines = '''
        outputPressure 200
        outputPressure 100
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['outputPressure'], 100)

    def test__parse_input_file_float_variable(self):
        lines = '''
        set pressure_repeat 200
        outputPressure $pressure_repeat
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['outputPressure'], 200)