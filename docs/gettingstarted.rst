Getting Started
===============

Installation
------------

Install the software with pip

::

	$ pip install git+https://github.com/ferchault/iago/tree/master/src

You may have to install the following requirements, depending on what you want to do:

::

	$ pip install pandas numpy pint MDAnalysis paramiko

For analysis, the use of `jupyter <http://jupyter.org/>`_ notebooks is highly recommended, although not necessary. If required, please install :command:`jupyter`

::

	$ pip install jupyter


If you want to access remote files other than SSH, you also need the :command:`rclone` `(website) <http://rclone.org/>`_ command available on your machine.

Glossary
--------
Some expressions are used throughout the document which may not be obvious at first glance. Here is an explanation of them.

.. glossary::

	Bucket
		One folder on disk. May contain multiple calculation with different software, but only for one system. This means that the number of atoms and the physical system under investigation is the same for all calculations in one bucket. The topology may change, though. Example: an ab-initio molecular dynamics trajectory of 32 waters is a bucket, even since protons may hop from water molecule to water molecule. Taking one configuration, generating a classical MD with fixed topology to compare to belongs into the same bucket. A similar simulation with 64 water molecules would be another bucket. Every bucket has an ID (currently a random MD5-like string).

	Run
		Part of a :term:`bucket`. Every :term:`bucket` has one or more runs and every run belongs to exactly one :term:`bucket`. A single run is defined as a simulation package (cp2k, namd, gromacs, Gaussian, amber, ...) started once with a single set of input files and a single set of output files. It is acceptable for runs to share input files. You could store them in an *input* folder directly in the :term:`bucket`. Restarts of continuous molecular dynamics trajectories are in fact separate runs in the same :term:`bucket`.

	Location
		Folder that may contain :term:`buckets <bucket>`. Can be on local disk or on a remote system. Iago will not distinguish by the data storage location.

	Location Group
		A set of :term:`locations <location>`. Essentially, a location group behaves as if all raw data and all databases would be available in one single local storage location.

	Analyser
		Script that defines what properties to extract from raw data. As soon as a run has completed, *iago* can calculate the derivative information like a plane fitted through a set of coordinates. The result of all this is stored in the database you can query with *jupyter*. In order to make this work, you need to tell *iago* which properties to calculate for which part of the system. This is done in the analyser script :file:`iago-analysis.py`. Per :term:`bucket`, there is exactly one analyser script.

	Frame
		Any single geometry in a :term:`bucket` is referred to as frame. For this context, it does not matter whether the :term:`run` is a single-point calculation or forms a molecular dynamics run.

This structure roughly equates to the following folder layout:

::

  bucket-name-6d78579fa1e849a2a58f794fa784c1ea
  |-- iago-analysis.py
  |-- iagodb.json
  |
  |-- input/
  |   |-- input.xyz
  |   |-- input.psf
  |
  |-- run-0/
  |   |-- run.inp
  |   |-- run.log
  |   |-- run.dcd
  |
  |-- run-1/
  |   |-- run.restart
  |   |-- run.log
  |   |-- run.dcd
  | ...

Example Bucket
--------------
If you want to start off a working example, you may copy a :term:`bucket` directory from the `tests/fixtures <https://github.com/ferchault/iago/tree/master/tests/fixtures/>`_ folder in the repository. Just configure the path as shown in the next section. The following steps, most importantly the :file:`iago-analysis.py` script is already set-up.

Configuration
-------------

Create a file :file:`.iago.conf` in your home directory. This file holds the configuration options where to look for :term:`buckets <bucket>` on both local and remote machines. One example for local file access only could look like this:

::

  [localmachine]
  url=file:///home/username/data/

Here, the path :file:`/home/{username}/data/` is the location where the :term:`bucket` folders are stored. All folders that have no MD5-hash in their name will be ignored.

If you want to set up the configuration with the example :term:`bucket`, then please make sure to specify the path to the parent folder only. If you have downloaded the sample :term:`bucket` to your :file:`~/Downloads/` directory, then add the following lines to your :file:`.iago.conf` in your home directory:

::

  [example]
  url=file:///home/username/Downloads/


Analyser
--------

Finally, *iago* needs to know what to extract from the trajectory. This is done by creating a :term:`bucket`. Since this is specific to the :term:`bucket`, the analyser script :file:`iago-analysis.py` has to be created in the top-level directory. An example of this file looks like this:

.. code-block:: python
	:linenos:

	import iago
	import os

	class Analyser(iago.Analyser):
		def setup(self):
			self.path = os.getcwd()

		def define_groups(self):
			self.static_load_groups('index.ndx')
			self.static_group('test', 1, 3, 4, 5)

		def calculated_columns(self):
			self.dynamic_plane(
				'myplane',
				'group test',
				normal=(0, 0, 1),
				framesel=slice(2),
				comment='My test plane.')

	if __name__ == '__main__':
		a = Analyser()
		a.run()

First, the *iago* module is loaded. The data to analyse is defined using the class methods as shown above, executed in that order. First, the base directory for this :term:`bucket` gets defined, followed by loading all the groups from the gromacs index file :file:`index.ndx` and defining a static group for atoms 1, 3, 4, 5. The data to calculate based on the trajectories and the meta data defined in the :class:`iago.Analyser` subclass is subsequently defined in the :func:`iago.Analyser.Analyser.calculated_columns` method. In the example, a plane with the label *myplane* is added to the database where the plane is defined by the coordinates of the atoms in group *test*. For details and a list of available methods, see the documentation of the :class:`iago.Analyser` class.

Once this :term:`analyser` file :file:`iago-analysis.py` has been created in the :term:`bucket` directory, you can run it in two ways. Locally, you can start the command line and run

::

	$ python iago-analysis.py

or from other python code (e.g., if you already have a *jupyter* notebook) with

.. code-block:: python

	import iago
	lg = iago.get_location_group()
	lg.build_database('bucket name or bucket id')
	# optionally: load it from remote
	db = lg.fetch_database('bucket name or bucket id')

The second approach also works for remote locations if they are available via SSH and have *iago* installed. For remote SSH access, passwordless and key-based login has to be set up.

Running the analysis script (ideally) gives no output on the command line and produces a database file :file:`iagodb.json`. This database file contains all calculated information from the run including the input files in a well-structured manner. You can query the results using *iago* as outlined in the next section.

Working With the Results
------------------------

So now everything is in working order. You can look into the database by loading it from *jupyter*. To do this, open a terminal and launch :command:`jupyter`

::

	$ jupyter notebook

Your browser should open. Create a new notebook and run the following code:

.. code-block:: python

	import iago
	lg = iago.get_location_group()
	db = lg.fetch_database('bucket name or bucket id')

It is always required to create a :term:`location group` first, since it caches the contents of remote repositories to speed up access. You can fetch the database by either its name or the (unique) ID. If (as in the previous example) your bucket directory is called :file:`bucket-name-6d78579fa1e849a2a58f794fa784c1ea`, then the following two lines are equivalent

.. code-block:: python

	db = lg.fetch_database('bucket-name')
	db = lg.fetch_database('6d78579fa1e849a2a58f794fa784c1ea')

Should there be two buckets of the same name, the first line will raise an error, since it is not clear which bucket the command is referring to.

The *db* object is a regular class. Its attributes are explained in detail here: :class:`iago.DatabaseProvider.DB`. E.g. if you were to inspect the configuration and then plot the z-component of the normal vector of the plane produced by the sample :file:`iago-analysis.py` above, then this could be done as follows in *jupyter*

.. code-block:: python

	import matplotlib.pyplot as plt
	%matplotlib inline

	print db.config['run-name']
	plt.plot(db.planes.frame, db.planes.normal_z)

If further data is required that currently is not part of the database, :file:`iago-analysis.py` has to be updated and re-run. Otherwise, now various data can be plotted interactively.