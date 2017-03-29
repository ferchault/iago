from unittest import TestCase

# custom modules
import iago.Reader

class TestNAMDReader(TestCase):
    def test__parse_input_file_last_one_wins(self):
        lines = '''
        outputPressure off
        outputPressure on
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['outputPressure'], 'on')

    def test__parse_input_file_float_variable(self):
        lines = '''
        set pressure_repeat 200
        outputPressure $pressure_repeat
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['outputPressure'], 200)

    def test__parse_input_file_break_after_minimize(self):
        lines = '''
        outputPressure off
        minimize 100
        outputPressure on
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['outputPressure'], 'off')
        self.assertEqual(result['minimize'], 100)

    def test__parse_input_file_run_with_variable(self):
        lines = '''
        set steps 100
        run $steps
        run 200
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['run'], 100)

    def test__parse_input_file_test_NVT(self):
        lines = '''
        langevin on
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['ensemble'], 'NVT')

    def test__parse_input_file_test_NPT(self):
        lines = '''
        langevin on
        langevinPiston on
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['ensemble'], 'NPT')

    def test__parse_input_file_test_remove_ending_comment(self):
        lines = '''
        restartfreq   500;#500ps
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['restartfreq'], 500)

    def test__parse_input_file_test_multiple_entries(self):
        lines = '''
        sphericalBCcenter   30.3081743413, 28.8049907121, 15.353994423
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['sphericalBCcenter'], [30.3081743413, 28.8049907121, 15.353994423])

    def test__parse_logfile_test_readin_one_step(self):
        lines = '''
        ETITLE:      TS           BOND          ANGLE
        ENERGY: 0 1.5 2.5
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_output_file(lines)
        self.assertEqual(result[0]['TS'], 0)
        self.assertEqual(result[0]['ANGLE'], 2.5)

    def test__parse_logfile_test_change_ETITLE(self):
        lines = '''
        ETITLE:  TS  BOND  ANGLE
        ENERGY: 0 1.5 2.5
        ENERGY: 1 11.5 12.5
        ETITLE: TS DIHEDRAL
        ENERGY: 2 21.5
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_output_file(lines)
        self.assertEqual(result[1]['ANGLE'], 12.5)
        self.assertEqual(result[2]['TS'], 2)
        self.assertEqual(result[2]['DIHEDRAL'], 21.5)
