# standard modules
import os
from unittest import TestCase

# custom modules
import iago

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class Analyser(iago.Analyser):
	def setup(self):
		self.path = os.path.join(BASE_DIR, 'fixtures', 'debug-88ac57b3e437fa1d5d26d00d6c768324')

	def define_groups(self):
		self.static_load_groups('index.ndx')
		self.static_group('debug', 1, 3, 4, 5)
		self.static_group('test2', [1, 2, 3])
		self.static_group('iron', 'type FE1 or type FE2')
		self.static_group({'test3': [1, 2, 3]})

	def calculated_columns(self):
		self.dynamic_plane(
			'O3A',
			'group O3A',
			normal=(0, 0, 1),
			framesel=slice(2),
			comment='Plane defined by the last triply coordinated oxygens on the A side.')
		self.dynamic_distance(
			'OH',
			'type O',
			'type H',
			cutoff=1.25,
			framesel=slice(2),
			comment='OH-distances up to transition state bond length'
		)
		#self.dynamic_distance(
		#	'O3A',
		#	'plane O3A',
		#	where='element O',
		#	comment='Height above dynamic surface.')
		#self.dynamic_hbonds(
		#	'all',
		#	'element O',
		#	'element O',
		#	comment='All hydrogen bonds.')


class TestUtils(TestCase):
	def test_create_analyser(self):
		a = Analyser()
		a.run()
