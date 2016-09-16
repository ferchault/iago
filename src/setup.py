#!/usr/bin/env python

from setuptools import setup
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.pardir, 'src')))
sys.path.append(os.path.abspath(os.pardir))

setup(name='iago',
      version='0.0',
      description='Implicit Archiver for Generated Observables',
      author='Guido Falk von Rudorff',
      author_email='guido@vonrudorff.de',
      url='https://github.com/ferchault/iago',
      packages=['iago', ],
      test_suite="tests",
      license='MIT',
      classifiers=['Development Status :: 3 - Alpha', ],
      )
