#!/usr/bin/env python

from setuptools import setup

setup(name = "cagecleaner",
      version = "0.0.1",
      author="Lucas De Vrieze",
      author_email="lucas.devrieze@kuleuven.be",
      license = "MIT",
      description = "Redundancy removal tool for cblaster hits",
      python_requires = ">=3.10",
      classifiers = [
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      entry_points = {"console_scripts": ['cagecleaner = cagecleaner.cagecleaner:main']},
      scripts = ['src/cagecleaner/dereplicate_assemblies.sh',
                 'src/cagecleaner/download_assemblies.sh',
                 'src/cagecleaner/get_accessions.sh']
      )
