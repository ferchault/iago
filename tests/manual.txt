import iago

class Analyser(iago.Analyser):
	def setup(self):
		self.parser = iago.CP2KParser(self)

	def define_groups(self):
		self.static_load_groups('index.ndx')
		self.static_group('test', 1, 3, 4, 5)
		self.static_group('test2', [1, 2, 3])
		self.static_group('iron', 'type FE1 or type FE2')
		self.static_group({'test3': [1, 2, 3]})

	def calculated_columns(self):
		self.dynamic_plane(
			'O3A',
			'group O3A',
			normal=(0, 0, 1),
			comment='Plane defined by the last triply coordinated oxygens on the A side.')
		self.dynamic_distance(
			'O3A',
			'plane O3A',
			where='element O',
			comment='Height above dynamic surface.')
		self.dynamic_hbonds(
			'all',
			'element O',
			'element O',
			comment='All hydrogen bonds.')

if __name__ == '__main__':
	a = Analyser()
	a.run()


# calculate results
import iago
db = iago.fetch_bucket('name or id')

# hbonds donated per atom per time
hbts = db.hbonds.group_by('run frame donor').count()
dts = db.distances.O3A
combo = pd.merge(hbts, dts, on='run frame donor-atom')
plt.scatter(combo.d, combo.count)