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

For analysis, the use of `jupyter <http://jupyter.org/>`_ notebooks is highly recommended, although not necessary. If required, please install *jupyter*

::

	$ pip install jupyter


If you want to access remote files other than SSH, you also need the `rclone <http://rclone.org/>`_ command available on your machine.

Terminology
-----------
Some expressions are used throughout the document which may not be obvious at first glance. Here is an explanation of them.

.. _whatis-bucket:

Bucket
  One folder on disk. May contain multiple calculation with different software, but only for one system. This means that the number of atoms and the physical system under investigation is the same for all calculations in one bucket. The topology may change, though. Example: an ab-initio molecular dynamics trajectory of 32 waters is a bucket, even since protons may hop from water molecule to water molecule. Taking one configuration, generating a classical MD with fixed topology to compare to belongs into the same bucket. A similar simulation with 64 water molecules would be another bucket. Every bucket has an ID (currently a random MD5-like string).

.. _whatis-run:

Run
  Part of a :ref:`bucket <whatis-bucket>`. Every :ref:`bucket <whatis-bucket>` has one or more runs and every run belongs to exactly one :ref:`bucket <whatis-bucket>`. A single run is defined as a simulation package (cp2k, namd, gromacs, Gaussian, amber, ...) started once with a single set of input files and a single set of output files. It is acceptable for runs to share input files. You could store them in an *input* folder directly in the :ref:`bucket <whatis-bucket>`. Restarts of continuous molecular dynamics trajectories are in fact separate runs in the same :ref:`bucket <whatis-bucket>`.

.. _whatis-location:

Location
  Folder that may contain :ref:`buckets <whatis-bucket>`. Can be on local disk or on a remote system. Iago will not distinguish by the data storage location.

.. _whatis-analyser:

Analyser
  Script that defines what properties to extract from raw data. As soon as a run has completed, *iago* can calculate the derivative information like a plane fitted through a set of coordinates. The result of all this is stored in the database you can query with *jupyter*. In order to make this work, you need to tell *iago* which properties to calculate for which part of the system. This is done in the analyser script *iago-analysis.py*. Per :ref:`bucket <whatis-bucket>`, there is exactly one analyser script.

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

Configuration
-------------

Create a file *.iago.conf* in your home directory. This file holds the configuration options where to look for :ref:`buckets <whatis-bucket>` on both local and remote machines. One example for local file access only could look like this:

::

  [localmachine]
  url=file:///home/username/data/

Here, the path */home/username/data/* is the location where the :ref:`bucket <whatis-bucket>` folders are stored. All folders that have no MD5-hash in their name will be ignored.

Analyser
--------

Finally, *iago* needs to know what to extract from the trajectory. This is done by creating an :ref:`bucket <whatis-analyser>`. Since this is specific to the :ref:`bucket <whatis-bucket>`, the analyser script has to be created in the top-level directory. An example of this file looks like this:

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