Data Structure
==============

So you have written an :term:`analyser` and ran it. You can :ref:`load the database <working-results>` but you are wondering
where to find all the properties? They are all documented here, although for some software packages, there may be no
current support for a certain quantity. Please refer to the list of :doc:`supported quantities <supported>` first.

In the following, it will be assumed that *iago* is loaded as python module in either *python* or *jupyter*. An example
of this happening could be as follows (refer to :doc:`Getting Started <gettingstarted>` for the details)

.. code-block:: python

	import iago
	lg = iago.get_location_group()
	db = lg.fetch_database('bucket name or bucket id')


Input
-----

Raw input data is found in

.. code-block:: python

	db.config

with the run name as key and a tree dictionary as value. Here are a few examples for CP2K:

.. code-block:: python

	# print time step in femtoseconds for the run 'test'
	print db.config['test']['MOTION']['MD']['TIMESTEP']

	# print the same info but with implicit navigation
	print db.config['test'].MOTION.MD.TIMESTEP

	# for each run, print all atom kind information
	for run in db.config.keys():
		print db.config[run].FORCE_EVAL.SUBSYS.KIND

While strings are converted to floats whenever applicable, no units are converted from the internal units used in the
simulation codes.

Calculated Data
---------------

Planes
++++++

Calculated plane data can be found in

.. code-block:: python

	db.planes


============= =============== ===========================================================
Column        Unit            Description
============= =============== ===========================================================
run           --              Run name. String.
frame         --              Frame number. Integer.
name          --              Plane name. String.
normal_x      --              Normal vector: x component. Float.
normal_y      --              Normal vector: y component. Float.
normal_z      --              Normal vector: z component. Float.
support_x     Angstrom        Support point: x component. Float.
support_y     Angstrom        Support point: y component. Float.
support_z     Angstrom        Support point: z component. Float.
============= =============== ===========================================================

Distances
+++++++++

Calculated distances can be found in

.. code-block:: python

	db.distances      # for atom-atom distances
	db.planedistances # for atom-plane distances

============= =============== ===========================================================
Column        Unit            Description
============= =============== ===========================================================
run           --              Run name. String.
frame         --              Frame number. Integer.
name          --              Plane name. String.
atom1         --              First atom index. Integer.
atom2         --              Second atom index. Integer.
dist          Angstrom        Distance. Float.
============= =============== ===========================================================

============= =============== ===========================================================
Column        Unit            Description
============= =============== ===========================================================
run           --              Run name. String.
frame         --              Frame number. Integer.
name          --              Plane name. String.
plane         --              Plane name. String.
atom2         --              Atom index. Integer.
dist          Angstrom        Distance. Float.
============= =============== ===========================================================
