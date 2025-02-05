# cagecleaner

**>>> `cagecleaner` has now been integrated into [`cblaster`](https://github.com/gamcil/cblaster)! <<<**

## Outline

`cagecleaner` removes genomic redundancy from gene cluster hit sets identified by [`cblaster`](https://github.com/gamcil/cblaster). The redundancy in target databases used by `cblaster` often propagates into the result set, requiring extensive manual curation before downstream analyses and visualisation can be carried out.

Given a session file from a `cblaster` run (or from a [`CAGECAT`](https://cagecat.bioinformatics.nl/) run), `cagecleaner` retrieves all hit-associated genome assemblies, groups these into assembly clusters by ANI and identifies a representative assembly for each assembly cluster using `skDER`. In addition, `cagecleaner` can reinclude hits that are different at the gene cluster level despite the genomic redundancy, and this by different gene cluster content and/or by outlier `cblaster` scores. Finally, `cagecleaner` returns a filtered `cblaster` session file as well as a list of retained gene cluster IDs for easier downstream analysis.

![workflow](workflow.png)

## Output

This tool will produce seven final output files
    - filtered_session.json: a filtered cblaster session file
    - filtered_binary.txt: a cblaster binary presence/absence table, containing only the retained hits.
    - filtered_summary.txt: a cblaster summary file, containing only the retained hits.
    - clusters.txt: the corresponding cluster IDs from the cblaster summary file for each retained hit.
    - genome_cluster_sizes.txt: the number of genomes in a dereplication genome cluster, referred to by the dereplication representative genome.
    - genome_cluster_status.txt: a table with scaffold IDs, their representative genome assembly and their dereplication status.
    - scaffold_assembly_pairs.txt: a table with scaffold IDs and the IDs of the genome assemblies of which they are part.
    
There are four possible dereplication statuses:
    - 'dereplication_representative': this scaffold is part of the genome assembly that has been selected as the representative of a genome cluster.
    - 'readded_by_content': this scaffold has been kept as it contains a hit that is different in content from the one of the dereplication representative.
    - 'readded_by_score': this scaffold has been kept as it contains a hit that has an outlier cblaster score.
    - 'redundant': this scaffold has not been retained and is therefore removed from the final output.

## Installation

First set up a `conda` environment using the `env.yml` file in this repo, and activate the environment.

```
conda env create -f env.yml
conda activate cagecleaner
```

Then install `cagecleaner` inside this environment using `pip`. First check you have the right `pip` using `which pip`, which should point to the `pip` instance inside the `cagecleaner` environment.

```
pip install cagecleaner
```

## Dependencies

`cagecleaner` has been developed on Python 3.10. All external dependencies listed below are managed by the `conda` environment, except for the NCBI EDirect utilities, which can be installed as outlined [here](https://www.ncbi.nlm.nih.gov/books/NBK179288/).

 - NCBI EDirect utilities (>= v21.6)
 - NCBI Datasets CLI (v16.39.0)
 - skDER (v1.2.8)
 - pandas (v2.2.3)
 - scipy (v1.14.1)
 - BioPython (v1.84)
 - more-itertools (v10.5)

 ## Usage

 `cagecleaner` expects as inputs at least the `cblaster` binary and summary files containing NCBI Nucleotide accession IDs. A dereplication run using the default settings can be started as simply as:
 ```
 cagecleaner -b binary.txt -s summary.txt
 ```

 Help message:
 ```
usage: cagecleaner [-c CORES] [-h] [-v] [-s SESSION_FILE] [-o OUTPUT_DIR] [--keep_downloads] [--keep_dereplication]
                   [--keep_intermediate] [--download_batch DOWNLOAD_BATCH] [-a ANI] [--no_recovery_content]
                   [--no_recovery_score] [--min_z_score ZSCORE_OUTLIER_THRESHOLD]
                   [--min_score_diff MINIMAL_SCORE_DIFFERENCE]

    cagecleaner: A tool to remove redundancy from cblaster hits.
    
    cagecleaner reduces redundancy in cblaster hit sets by dereplicating the genomes containing the hits. 
    It can also recover hits that would have been omitted by this dereplication if they have a different gene cluster content
    or an outlier cblaster score.
    
    cagecleaner first retrieves the assembly accession IDs of each cblaster hit via NCBI Entrez-Direct utilities, 
    then downloads these assemblies using NCBI Datasets CLI, and then dereplicates these assemblies using skDER.
    By default, cblaster hits that have an alternative gene cluster content or an outlier cblaster score 
    (calculated via z-scores) are recovered as well.
                                     

General:
  -c CORES, --cores CORES
                        Number of cores to use (default: 1)
  -h, --help            Show this help message and exit
  -v, --version         show program's version number and exit

Input / Output:
  -s SESSION_FILE, --session SESSION_FILE
                        Path to cblaster session file
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output directory (default: current working directory)
  --keep_downloads      Keep downloaded genomes
  --keep_dereplication  Keep skDER output
  --keep_intermediate   Keep all intermediate data. This overrules other keep flags.

Download:
  --download_batch DOWNLOAD_BATCH
                        Number of genomes to download in one batch (default: 300)

Dereplication:
  -a ANI, --ani ANI     ANI dereplication threshold (default: 99.0)

Hit recovery:
  --no_recovery_content
                        Skip recovering hits by cluster content
  --no_recovery_score   Skip recovering hits by outlier scores
  --min_z_score ZSCORE_OUTLIER_THRESHOLD
                        z-score threshold to consider hits outliers (default: 2.0)
  --min_score_diff MINIMAL_SCORE_DIFFERENCE
                        minimum cblaster score difference between hits to be considered different. Discards outlier
                        hits with a score difference below this threshold. (default: 0.1)

    Lucas De Vrieze, 2025
    (c) Masschelein lab, VIB
 ```

## Example case

We provide two example cases in the folder `examples` in this repo. We have already provided the `cblaster` output files as well as the original query `fasta`.

In the first case, 1146 gene cluster hits from *Staphylococcus* spp. should be reduced to 22 non-redundant hits. Running the `cagecleaner` for this example is done like below

```
cd N398V589S066P61
cagecleaner -s session.json -o output -c 20
```
This should give the seven output files in a new subfolder `output`. This should take about 10' using 20 cores, depending on the download speed of your internet connection. This requires 1.2 GB of disk space and 1.7 GB of RAM.

```
$ dir -1 output
clusters.txt
filtered_binary.txt
filtered_session.txt
filtered_summary.txt
genome_cluster_sizes.txt
genome_cluster_status.txt
scaffold_assembly_pairs.txt
```

In the second case, we queried four genes from MIBiG entry BGC0000194 (actinorhodin from *Streptomyces coelicolor A3(2)*), which yielded 8934 gene cluster hits. `cagecleaner` should reduce this to 4847 hits in about 1.5 h using 20 cores. 28.5 GB of disk space and 27.6 GB of RAM are required for this example case.

```
cd actinorhodin
cagecleaner -s session.json -o output -c 20
```

## Citations

`cagecleaner` relies heavily on the `skDER` genome dereplication tool and its main dependendy `skani`, so we give these tools proper credit.
```
Salamzade, R., & Kalan, L. R. (2023). skDER: microbial genome dereplication approaches for comparative and metagenomic applications. https://doi.org/10.1101/2023.09.27.559801
Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nature Methods, 20(11), 1661–1665. https://doi.org/10.1038/s41592-023-02018-3
```

Please cite the `cagecleaner` manuscript:
```
In preparation
```

## License

`cagecleaner` is freely available under an MIT license.

Use of the third-party software, libraries or code referred to in the References section above may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.
