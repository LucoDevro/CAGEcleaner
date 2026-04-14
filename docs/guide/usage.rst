Usage
=====

Normally, you don't have to choose between remote and local modes. CAGEcleaner automatically recognises from the cblaster session in which mode it was generated and thus should be dereplicated. However, the local dereplication mode requires you to have a local genome database, preferably the exact one you used to run your local cblaster session on.

Remote mode
------------

The remote mode does not require other file inputs than your cblaster session. A dereplication run using default settings can be started as simply as:

.. code-block:: bash

	cagecleaner -s session.json

A slightly more complex example to run CAGEcleaner at 14 cores in full-genome dereplication mode with an identity threshold of 90 % and a coverage threshold of 50 %.

.. code-block:: bash

	cagecleaner -s session.json --cores 14 -i 90 -c 50

A more complicated example to run CAGEcleaner at 14 cores in region dereplication mode with both thresholds at 95 %, adding a sequence margin of 5 kb to both sides of each cluster, discarding clusters that are at a contig edge.

.. code-block:: bash

	cagecleaner -s session.json --cores 14 -i 90 -c 90 --method regions -m 5000 --strict

A complex example that runs CAGEcleaner at 14 cores in region dereplication mode with the identity threshold at 95 %, the coverage threshold at 80 %, a sequence margin of size 10 kb, saving the results in a subfolder `results`, keeping all intermediate files, excluding a scaffold with a duplicate ID from the analysis, bypassing another one with a unique label, and disabling the hit recovery by outlier homology scores.

.. code-block:: bash

	cagecleaner -s session.json --cores 14 -i 95 -c 80 --method regions -m 10000 -o results --keep_intermediate --exs assembly1:ctg1 --bys ctg131 --no_recovery_score


Local mode
----------

To get a local mode cblaster session, you must have generated a local cblaster database from a folder of fasta or genbank files. In local mode, CAGEcleaner needs the path to that folder to be passed on via the ``-g`` flag. For each one of the examples in section above, you can just add the ``-g`` flag to run CAGEcleaner from a local search session.

A dereplication run using default settings can be started as simply as:

.. code-block:: bash

	cagecleaner -s session.json -g genomes

A slightly more complex example to run CAGEcleaner at 14 cores in full-genome dereplication mode with an identity threshold of 90 % and a coverage threshold of 50 %.

.. code-block:: bash

	cagecleaner -s session.json -g genomes --cores 14 -i 90 -c 50

A more complicated example to run CAGEcleaner at 14 cores in region dereplication mode with both thresholds at 95 %, adding a sequence margin of 5 kb to both sides of each cluster, discarding clusters that are at a contig edge.

.. code-block:: bash

	cagecleaner -s session.json -g genomes --cores 14 -i 90 -c 90 --method regions -m 5000 --strict

A complex example that runs CAGEcleaner at 14 cores in region dereplication mode with the identity threshold at 95 %, the coverage threshold at 80 %, a sequence margin of size 10 kb, saving the results in a subfolder `results`, keeping all intermediate files, excluding a scaffold with a duplicate ID from the analysis, bypassing another one with a unique label, and disabling the hit recovery by outlier homology scores.

.. code-block:: bash

	cagecleaner -s session.json -g genomes --cores 14 -i 95 -c 80 --method regions -m 10000 -o results --keep_intermediate --exs assembly1:ctg1 --bys ctg131 --no_recovery_score

Output files
------------

This tool produces at least six output files:

- ``filtered_session.json`` a filtered cblaster session file
- ``filtered_binary.txt`` a cblaster binary presence/absence table, containing only the retained hits
- ``filtered_summary.txt`` a cblaster summary file, containing only the retained hits
- ``extended_binary.txt`` a cblaster binary table extended with the following columns

	- ``Number`` Number of the cluster in the cblaster session
	- ``Strand`` Tuple representing on which strand each query gene homolog was found.
	- ``Layout_group`` Tuple representing in which layout/order the genes are ordered on the genome
	- ``representative`` the representative genome or region
	- ``dereplication_status`` status of this cluster (Was it removed/kept? Why?)
	- ``assembly_file`` fasta file of the genome or region

- ``retained_cluster_numbers.txt`` the cluster numbers of each retained hit
- ``genome_cluster_sizes.txt`` a tab-separated text file with the number of genomes in each dereplication genome cluster

In genome dereplication mode, the file ``unmapped.scaffolds.txt`` may also be generated. This file contains scaffold IDs for which no associated genome file was found. There are several reasons this might happen: ID mismatch between the session and the downloaded fasta file, a missing download...

Optionally, the downloads can be saved in an additional subfolder ``downloads`` (flag ``--keep_downloads`` or ``--keep_intermediate``), and/or the dereplication output files generated by skDER or MMseqs2 can be kept in a subfolder ``dereplication`` (flag ``--keep_dereplication`` or ``--keep_intermediate``).

In ``extended_binary.txt``, There are four possible dereplication statuses.

- ``dereplication_representative`` This hit is part of the sequence that has been selected as the representative of a sequence cluster.
- ``readded_by_content`` This hit has been kept as it has a different group layout than the dereplication representative.
- ``readded_by_score`` This hit has been kept as it has an outlier homology score.
- ``redundant`` This hit has not been retained and is therefore removed from the filtered output.


Hit recovery
------------

The diversity in genomes or genomic neighbourhoods does not always correlate well with the diversity in the gene clusters themselves. As a result, a dereplication only on the genome level may miss diversity in gene clusters. CAGEcleaner can recover cluster diversity using the two strategies outlined below.

Cluster content
~~~~~~~~~~~~~~~

This strategy assesses the number of homologs that were found for each query gene.

For example, suppose you found two groups of hits. Group 1 has one gene A and one gene B, and contains the hit that was retained due to the dereplication. On the other hand, group 2 has two genes A and one gene B. If you would stick to the hit reduction after the dereplication, group 2 would not be represented anymore.

In each content group not containing the dereplication representative hit, the hit to be retained is chosen randomly.

.. note::

   This recovery method may also recover redundant cluster hits of which a query homolog was missed by cblaster, for example because of a too high e-value.

Homology score
~~~~~~~~~~~~~~~

This strategy assesses the homology score of each hit.

For example, suppose you found a group of gene clusters with similar gene contents. However, there is a cluster X that has evolved quite rapidly, or that has swapped one of its genes for a paralog. This shift translates in a more remote homology and thus a lower cblaster score. This diverging score is identified using outlier statistics and the hit is retained. If you would stick to the hit reduction after the dereplication, cluster X would have been lost.

.. note::
   
   When recovery by cluster content is enabled, outliers are determined per content group.

.. note::

   By default, outliers scores are required to be sufficiently different from the mean score, implemented using an absolute score threshold. Otherwise, you may find an outlier having a score of 3.17 in a group of hundreds of hits having a score of 3.16.

Example hit recovery
~~~~~~~~~~~~~~~~~~~~

Let's consider a (simplified) extended binary table output table where CAGEcleaner has added the representative, and the dereplication status for each row:

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org1     scaff1   3.17  1     1     1     org1           dereplication_representative
org2     scaff2   3.16  1     1     1     org1           redundant
org3     scaff3   3.17  1     0     1     org1           redundant
org4     scaff4   3.17  1     1     1     org4           dereplication_representative
org5     scaff5   2.92  1     0     1     org4           redundant
org6     scaff6   3.12  1     0     1     org4           redundant
org7     scaff7   3.16  1     1     1     org4           redundant
org8     scaff8   3.14  1     0     1     org4           redundant
======== ======== ===== ===== ===== ===== ============== ====================

The first step is to group the rows based on their corresponding representatives. In this case, there are two groups; org1 and org4. This will result in the following two tables:

*REPRESENTATIVE = org1:*

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org1     scaff1   3.17  1     1     1     org1           dereplication_representative
org2     scaff2   3.16  1     1     1     org1           redundant
org3     scaff3   3.17  1     0     1     org1           redundant
======== ======== ===== ===== ===== ===== ============== ====================

*REPRESENTATIVE = org4:*

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org4     scaff4   3.17  1     1     1     org4           dereplication_representative
org5     scaff5   2.92  1     0     1     org4           redundant
org6     scaff6   3.12  1     0     1     org4           redundant
org7     scaff7   3.16  1     1     1     org4           redundant
org8     scaff8   3.14  1     0     1     org4           redundant
======== ======== ===== ===== ===== ===== ============== ====================

Now notice that the gene cluster composition varies within the groupings, something that might also occur in a real-life example. Some hits contain all the genes from the queries, whereas others are missing gene2. In the first grouping, org3 is missing gene2. If recovery by content is enabled, CAGEcleaner will make two groups and select a random representative for each group if it doesn't contain a dereplication representative. Hence, in this case, two hits are retained: org1 as representative for the structural group of org1 and org2, and org3 for its own structural group.

*ALL GENES PRESENT:*

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org1     scaff1   3.17  1     1     1     org1           dereplication_representative
org2     scaff2   3.16  1     1     1     org1           redundant
======== ======== ===== ===== ===== ===== ============== ====================

*GENE2 MISSING:*

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org3     scaff3   3.17  1     0     1     org1           redundant
======== ======== ===== ===== ===== ===== ============== ====================

Now let's consider the second grouping with org4 as dereplication representative and divide this table by gene cluster content:

*ALL GENES PRESENT:*

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org4     scaff4   3.17  1     1     1     org4           dereplication_representative
org7     scaff7   3.16  1     1     1     org4           redundant
======== ======== ===== ===== ===== ===== ============== ====================

*GENE2 MISSING:*

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org5     scaff5   2.92  1     0     1     org4           redundant
org6     scaff6   3.12  1     0     1     org4           redundant
org8     scaff8   3.14  1     0     1     org4           redundant
======== ======== ===== ===== ===== ===== ============== ====================

The grouping where all genes are present contains a representative, so no recovery is required. The other grouping contains hits where gene2 is missing. Notice that org5's cblaster score (2.92) significantly deviates from the scores in the rest of this grouping. If recovery by score is enabled, CAGEcleaner will retain this score outlier. Equivalent to the org1 case, a random non-outlier hit will be selected to represent this grouping if recovery by content is enabled (in this case org6).

To summarise, the (simplified) binary table for this example case will look like below:

======== ======== ===== ===== ===== ===== ============== ====================
Organism Scaffold Score Gene1 Gene2 Gene3 Representative Dereplication status
======== ======== ===== ===== ===== ===== ============== ====================
org1     scaff1   3.17  1     1     1     org1           dereplication_representative
org3     scaff3   3.17  1     0     1     org1           redundant
org4     scaff4   3.17  1     1     1     org4           dereplication_representative
org5     scaff5   2.92  1     0     1     org4           redundant
org6     scaff6   3.12  1     0     1     org4           redundant
======== ======== ===== ===== ===== ===== ============== ====================


