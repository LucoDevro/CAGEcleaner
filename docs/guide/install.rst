Installation
============

.. note::
   
   Supported Python versions are 3.12, 3.13 and 3.14.

Conda (recommended)
-----------------------

This is the recommended and most straightforward way to install CAGEcleaner. It should work on any Linux or MacOS system. Create a fresh conda environment as well to keep it from meddling with your other tools. You can install CAGEcleaner either directly from Bioconda

.. code-block:: bash

	conda create -n cagecleaner -c bioconda -c conda-forge cagecleaner

or by using the conda `yml` environment file in this repo.

.. code-block:: bash

	conda env create -f env.yml

Then start using it by activating the conda environment.

.. code-block:: bash

	conda activate cagecleaner

Docker
-------

CAGEcleaner is also available as a Docker image from DockerHub. This is one of the recommended ways to run CAGEcleaner on Windows (the other one being running it using Windows' WSL feature).

.. code-block:: bash

	docker pull lucodevro/cagecleaner

There is no entrypoint set up so running CAGEcleaner requires prepending your CAGEcleaner command with the appropriate Docker commands.

.. code-block:: bash

	docker run lucodevro/cagecleaner -v <your-cblaster-session>:session.json -v <your-output-folder>:output cagecleaner -s session.json -o output

GitHub
-------

Alternatively, it is possible to install the latest semi-stable development version by cloning this repository and running the following command at the root of your local copy of this repository.

.. code-block:: bash

	pip install .

PyPi
------

CAGEcleaner is also installable from PyPi using pip, yet we do not recommend using this approach as some dependencies are not available from PyPi (NCBI Datasets CLI, Entrez Direct, any2fasta, MMseqs2) and therefore should be installed beforehand. So either make sure you have installed these dependencies separately, or use one of the other installation options.

.. code-block:: bash

	pip install cagecleaner

.. warning::
   
   We do not recommend using this approach as some non-Python dependencies are not available from PyPi (NCBI Datasets CLI, any2fasta, MMseqs2) and therefore should be installed beforehand. Check out CAGEcleaner's dependencies in the Bioconda recipe for more details.


