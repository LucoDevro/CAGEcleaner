# cagecleaner

## Outline

`cagecleaner` removes genomic redundancy from gene cluster hit sets identified by [`cblaster`](https://github.com/gamcil/cblaster). The redundancy in target databases used by `cblaster` often propagates into the result set, requiring extensive manual curation before downstream analyses and visualisation can be carried out.

Given the results files of a `cblaster` run (or a [`CAGECAT`](https://cagecat.bioinformatics.nl/) run), `cagecleaner` retrieves all hit-associated genome assemblies, groups these into assembly clusters and identifies a representative assembly for each assembly cluster using `skDER`. In addition, `cagecleaner` can reinclude hits that seem different at the gene cluster level despite the genomic redundancy, and this by different gene cluster content and/or by outlier `cblaster` scores. Finally, `cagecleaner` returns a filtered `cblaster` binary file and a list of retained gene cluster IDs for straightforward downstream analyses.

This tool has primarily been developed for `cblaster` searches against the NCBI nr database, but should work for any result set containing NCBI Nucleotide accession codes.

## Installation

First set up a `conda` environment using the `env.yml` file in this repo, and activate the environment.

```
conda env create -y -f env.yml
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
 usage: cagecleaner [-c CORES] [-h] [-v] [-o OUTPUT_DIR] [-b BINARY_FILE] [-s SUMMARY_FILE] [--validate-files]
                   [--keep-downloads] [--keep-dereplication] [--keep-intermediate]
                   [--download-batch DOWNLOAD_BATCH] [-a ANI] [--no-content-revisit] [--no-score-revisit]
                   [--min-z-score ZSCORE_OUTLIER_THRESHOLD] [--min-score-diff MINIMAL_SCORE_DIFFERENCE]

    cagecleaner: A tool to remove redundancy from cblaster hits.
    
    cagecleaner reduces redundancy in cblaster hit sets by dereplicating the genomes containing the hits. 
    It can also recover hits that would have been omitted by this dereplication if they have a different gene cluster content
    or an outlier cblaster score.
    
    cagecleaner has been designed for usage with the NCBI nr database. It first retrieves the assembly accession IDs
    of each cblaster hit via NCBI Entrez-Direct utilities, then downloads these assemblies using NCBI Datasets CLI,
    and then dereplicates these assemblies using skDER. If requested, cblaster hits that have an alternative gene cluster content
    or an outlier cblaster score (calculated via z-scores) are recovered.
                                     

General:
  -c CORES, --cores CORES
                        Number of cores to use (default: 1)
  -h, --help            Show this help message and exit
  -v, --version         show program's version number and exit

Input / Output:
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output directory (default: current working directory)
  -b BINARY_FILE, --binary BINARY_FILE
                        Path to cblaster binary file
  -s SUMMARY_FILE, --summary SUMMARY_FILE
                        Path to cblaster summary file
  --validate-files      Validate cblaster input files
  --keep-downloads      Keep downloaded genomes
  --keep-dereplication  Keep skDER output
  --keep-intermediate   Keep all intermediate data. This overrules other keep flags.

Download:
  --download-batch DOWNLOAD_BATCH
                        Number of genomes to download in one batch (default: 300)

Dereplication:
  -a ANI, --ani ANI     ANI dereplication threshold (default: 99.0)

Hit recovery:
  --no-content-revisit  Do not recover hits by cluster content
  --no-score-revisit    Do not recover hits by outlier scores
  --min-z-score ZSCORE_OUTLIER_THRESHOLD
                        z-score threshold to consider hits outliers (default: 2.0)
  --min-score-diff MINIMAL_SCORE_DIFFERENCE
                        minimum cblaster score difference between hits to be considered different. Discards outlier
                        hits with a score difference below this threshold. (default: 0.1)

    Lucas De Vrieze, 2025
    (c) Masschelein lab, VIB
 ```

## Citations

`cagecleaner` relies heavily on the `skDER` genome dereplication tool, so we give it proper credit.
```
Salamzade, R., & Kalan, L. R. (2023). skDER: microbial genome dereplication approaches for comparative and metagenomic applications. https://doi.org/10.1101/2023.09.27.559801
```

Please cite the `cagecleaner` manuscript:
```
In preparation
```

## License

`cagecleaner` is freely available under an MIT license.

Use of the third-party software, libraries or code referred to in the References section above may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.