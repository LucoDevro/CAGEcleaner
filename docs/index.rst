.. CAGEcleaner documentation master file, created by
   sphinx-quickstart on Tue Apr 14 13:00:12 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CAGEcleaner
=======================================

**Welcome to CAGECleaner's documentation!**

|bioconda|
|conda|
|manuscript|
|doi|

.. |bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: https://bioconda.github.io/recipes/cagecleaner/README.html

.. |conda| image:: https://img.shields.io/conda/dn/bioconda/CAGEcleaner.svg
   :target: https://anaconda.org/bioconda/cagecleaner/files

.. |manuscript| image:: https://img.shields.io/badge/Manuscript-Bioinformatics-darkblue?style=flat-square&maxAge=2678400
   :target: https://doi.org/10.1093/bioinformatics/btaf373

.. |doi| image:: https://zenodo.org/badge/904110273.svg
   :target: https://doi.org/10.5281/zenodo.14726119

*CAGEcleaner: A tool to reduce redundancy in gene mining hit sets.*
   
*CAGEcleaner* reduces redundancy in gene cluster hit sets, easing downstream analyses and visualisation. It features a taxonomically conservative dereplication mode that acts at the full genome level, and a more aggressive mode that acts at the level of the genomic neighbourhood of the cluster. In addition, it prevents clusters from being discarded if they show remarkable diversity based on gene cluster contents and homology scores.

In full genome mode, *CAGEcleaner* retrieves the full genome assemblies of the clusters' host genomes, performs a fast ANI-based full genome dereplication using skDER, and only keeps clusters that were part of the retained genomes.

In region mode, *CAGEcleaner* retrieves the nucleotide sequence of each cluster with an optional sequence margin on both sides, dereplicates these using MMseqs2, and only keeps clusters part of a representative region.

*CAGEcleaner* offers seamless integration with **cblaster**, as it has originally been developed as an auxiliary tool to be used in interaction with **cblaster**. Other inputs are now also possible via TSV files. Sessions filtered by *CAGEcleaner* can be plugged back in into the **cblaster** workflow.

If you find ``CAGEcleaner`` useful, please cite:

::

	De Vrieze, L., Biltjes, M., Lukashevich, S., Tsurumi, K., Masschelein, J. (2025).
	CAGEcleaner: reducing genomic redundancy in gene cluster mining. Bioinformatics
	https://doi.org/10.1093/bioinformatics/btaf373


.. toctree::
   :caption: User Guide
   :maxdepth: 2

   guide/install
   guide/usage


Comprehensive documentation for all API exposed by ``CAGEcleaner``:

.. toctree::
   :caption: API Reference
   :maxdepth: 1

   api_ref
