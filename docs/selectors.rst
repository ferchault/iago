Selection Syntax
================

For many operations, it is crucial to know to which atoms or :term:`frames <frame>` to apply it to. This selection
is based on strings that are of relevance to the :term:`analyser` in :file:`iago-analysis.py` where all operations on
raw data are defined. For each kind of selector, a subsection describes the supported syntax. The different kinds of
selectors are not interchangeable. This means that each operation asks for either of them at appropriate places. The
documentation of each :ref:`supported calculation <supported-calculation>` highlights the required types.

.. _selection-atom:

Atom Selectors
--------------

For atom selectors, the superset of the CHARMM
`selection syntax <http://www.mdanalysis.org/mdanalysis/documentation_pages/selections.html>`_ as supported by
*MDAnalysis* is accepted.


.. _selection_plane:

Plane Selector
--------------

To select a plane where applicable, simply use *plane planename*. The following complete example creates a plane based
on three atoms, and then calculates the distance of a fourth atom to said plane.

.. code-block:: python
	:linenos:
	:emphasize-lines: 14

	import iago
	import os

	class Analyser(iago.Analyser):
		def setup(self):
			self.path = os.getcwd()

		def define_groups(self):
			self.static_group('test', 10, 11, 12)

		def calculated_columns(self):
			self.dynamic_plane('myplane', 'group test', normal=(0, 0, 1),
				comment='Plane in all frames.')
			self.dynamic_distance('mydist', 'atom 5', 'plane myplane')

	if __name__ == '__main__':
		a = Analyser()
		a.run()


.. _selection_frame:

Frame Selector
--------------

Frame selectors are regular python :func:`slice` objects. When used in any of the :class:`iago.Analyser.Analyser`
functions, the selector will be applied to all frames of the bucket, not only to a certain run. The order is defined by
the frame numbers as obtained from the raw data through the :class:`iago.Reader.Reader` subclasses.

As an example, a bucket with two runs (*run-1* and *run-2*) with 25 frames each is to be parsed. The following snipped
will fit four planes:

.. code-block:: python
	:linenos:
	:emphasize-lines: 12, 13, 14, 15, 16, 17, 18, 19

	import iago
	import os

	class Analyser(iago.Analyser):
		def setup(self):
			self.path = os.getcwd()

		def define_groups(self):
			self.static_group('test', 10, 11, 12)

		def calculated_columns(self):
			self.dynamic_plane('plane_A', 'group test', normal=(0, 0, 1),
				framesel=None, comment='Plane in all frames.')
			self.dynamic_plane('plane_B', 'group test', normal=(0, 0, 1),
				framesel=slice(2), comment='Plane in the first two frames.')
			self.dynamic_plane('plane_C', 'group test', normal=(0, 0, 1),
				framesel=slice(2, 50), comment='Plane in frames 2-50.')
			self.dynamic_plane('plane_D', 'group test', normal=(0, 0, 1),
				framesel=slice(25, None), comment='Plane after frame 25.')

	if __name__ == '__main__':
		a = Analyser()
		a.run()

The first plane, *plane_A*, has no frame selection specified via the *framesel* parameter, meaning an individual plane
will be fitted for each frame of the whole bucket, i.e. for both runs (*run-1* and *run-2*). For *plane_B*, only the
first two frames will be analysed, while *plane_C* selects the range from frame 2 to 50. This means that frames 2-25 of
*run-1* will be included, followed by all frames of the second run. If the second parameter to :func:`slice` is *None*,
as in the last example *plane_D*, then all frames until the end of the trajectory are taken into account. Note that
:func:`iago.Analyser.Analyser.dynamic_plane` fits a plane in each frame, not an averaged plane over multiple frames.