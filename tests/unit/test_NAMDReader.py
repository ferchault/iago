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
        set temp 50
        reinitvels ${temp}
        run $steps
        run 200
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['run'], 100)
        self.assertEqual(result['reinitvels'], 50)

    def test__parse_input_file_run_with_curve_bracket_variable(self):
        lines = '''
        set steps 100
        set temp 50
        reinitvels ${temp}
        run ${steps}
        run 200
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['run'], 100)
        self.assertEqual(result['reinitvels'], 50)

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
        sphericalBCcenter   30.3081743413, 28.8049907121, 15.353994423+3
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['sphericalBCcenter'], [30.3081743413, 28.8049907121, 15.353994423])

    def test__parse_input_file_test_multiple_varialbes(self):
        lines = '''
        set path testpath
        set directory testdirectory
        restartpath I/am/${directory}/${path}/${path}/tests
        anotherpath I/am/$directory/$path$path${path}/tests'''.split('\n')
        result = iago.Reader.NAMDReader._parse_input_file(lines)
        self.assertEqual(result['restartpath'], 'I/am/testdirectory/testpath/testpath/tests')
        self.assertEqual(result['anotherpath'], 'I/am/testdirectory/testpathtestpathtestpath/tests')

    def test__parse_logfile_test_readin_one_step(self):
        lines = '''
        ETITLE:      TS           BOND          ANGLE  NOTRELEVANT PRESSAVG
        ENERGY: 0 1.5 1.0 42 3.5
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_output_file(lines)
        self.assertEqual(result[0]['frame'], 0)
        self.assertEqual(result[0]['angle'], 0.0433634)
        self.assertNotIn('NOTRELEVANT', result[0])
        self.assertEqual(result[0]['NAMD_PRESSAVG'], 3.5)

    def test__parse_logfile_test_change_ETITLE(self):
        lines = '''
        ETITLE:  TS  BOND  ANGLE
        ENERGY: 0 1.5 2.5
        ENERGY: 1 11.5 12.5
        ETITLE: TS DIHED
        ENERGY: 2 21.5
        '''.split('\n')
        result = iago.Reader.NAMDReader._parse_output_file(lines)
        self.assertEqual(result[1]['angle'], 12.5*0.0433634)
        self.assertEqual(result[2]['frame'], 2)
        self.assertEqual(result[2]['dihedral'], 21.5*0.0433634)
